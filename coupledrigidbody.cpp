#include "coupledrigidbody.h"
#include "util.h"
#include "array2.h"
#include "levelset.h"

//#include "cellweights.h"

#include <iostream>

using namespace std;

CoupledRigidBody::CoupledRigidBody(RigidBody& rb, int ni_, int nj_, float dx_):
ni(ni_),nj(nj_), dx(dx_),
rigidbody(rb), 
phi(ni,nj), 
u_vol(ni+1,nj), v_vol(ni,nj+1), 
u_mass(0), v_mass(0),
base_trans_x(ni,nj), base_trans_y(ni,nj), base_rot_z(ni,nj)
{
	phi.set_zero();
	u_vol.set_zero();
	v_vol.set_zero();
	update_grid_data();
}

bool CoupledRigidBody::is_inside(const Vec2f& point) {
   return interpolate_phi(point, phi, 0.5f*dx*Vec2f(1,1), dx) < 0;
}

const bool CoupledRigidBody::test_collision_and_project(const Vec2f& P, Vec2f& newPos) {
   return rigidbody.testCollisionAndProject(P, newPos);
}

void CoupledRigidBody::recompute_volumes() {

   u_vol.set_zero();
   v_vol.set_zero();
   const float full_vol = 1;

   //compute_finitevolume_area_weights(phi, full_vol, dx, u_vol, v_vol);

}

void CoupledRigidBody::recompute_phi() {
   for(int j = 0; j < nj; ++j) {
      for(int i =0; i < ni; ++i) {
         Vec2f location((i+0.5f)*dx, (j+0.5f)*dx);
         phi(i,j) = rigidbody.getSignedDist(location);
      }
   }
}

float sum_array(const Array2f& array) {
	float sum = 0;
	for(unsigned int i = 0; i < array.size(); ++i)
		sum += array.a[i];
	return sum;
}

void CoupledRigidBody::recompute_masses() {

   u_mass = rigidbody.getDensity() * sum_array(u_vol);
   v_mass = rigidbody.getDensity() * sum_array(v_vol);
}

void CoupledRigidBody::update_grid_data() {
   recompute_phi();
   recompute_volumes();
   recompute_masses();
}

void CoupledRigidBody::add_coupling(SparseMatrix<double>& matrix, Array2d& rhs) {
   build_base_vectors();
   add_coefficients(matrix);
   add_right_hand_side(rhs);
}

void CoupledRigidBody::build_base_vectors()
{
   // Translation coupling
   base_trans_x.assign(0);
   base_trans_y.assign(0);
   base_rot_z.assign(0);

   Vec2f centre_of_mass;
   rigidbody.getCOM(centre_of_mass);

   // Browsing u and v-cells together, since it's the same loop
   for (int j = 0; j != nj; ++j) {
      for (int i = 0; i != ni; ++i) {
         double u_term = ( u_vol(i,j) - u_vol(i+1,j) ) / dx;
         double v_term = ( v_vol(i,j) - v_vol(i,j+1) ) / dx;

         // Translation coupling
         base_trans_x(i, j) = -u_term / u_mass;
         base_trans_y(i, j) = -v_term / v_mass;

         // Rotation coupling
         Vec3f position((i + 0.5f) * dx, (j + 0.5f) * dx, 0); 
         Vec3f centre_3d(centre_of_mass[0], centre_of_mass[1], 0);
         Vec3f rad3d = position - centre_3d;
         Vec3f vol_terms((float)u_term, (float)v_term, 0);

         Vec3f result = -cross(rad3d,vol_terms);

         base_rot_z(i,j) = result[2];
      }
   }
}

void CoupledRigidBody::add_coefficients(SparseMatrixd& sparseMat)
{
   double val;
   Vec2f solidPosG;
   rigidbody.getCOM(solidPosG);

   const float Jinv = rigidbody.getInvInertiaModulus();

   for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {

         for (int k = 0; k < ni; ++k) {
            for (int m = 0; m < nj; ++m) {
               val = 0;

               //Translation
               val += base_trans_x(i,j) * base_trans_x(k, m) * u_mass; 
               val += base_trans_y(i,j) * base_trans_y(k, m) * v_mass;

               //Rotation
               val += base_rot_z(i,j) * base_rot_z(k,m) * Jinv;

               if(val != 0)
                  sparseMat.add_to_element(i + ni*j, k + ni*m, val);
            }
         }
      }
   }
}

void CoupledRigidBody::add_right_hand_side(Array2d& rhs) {

   Vec2f solidLinearVelocity;
   rigidbody.getLinearVelocity(solidLinearVelocity);

   float angular_velocity;
   rigidbody.getAngularVelocity(angular_velocity);

   for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
         // Translation
         rhs(i,j) += -solidLinearVelocity[0] * base_trans_x(i,j) * u_mass;
         rhs(i,j) += -solidLinearVelocity[1] * base_trans_y(i,j) * v_mass;

         //Rotation
         rhs(i,j) += -angular_velocity * base_rot_z(i,j);
      }
   }
}

void CoupledRigidBody::apply_pressure(const Array2d& pressure) {
   // Computing new velocity, from the pressure force applied on the solid

   Vec2f vGnext;
   float Lznext;

   rigidbody.getLinearVelocity(vGnext);
   rigidbody.getAngularMomentum(Lznext);

   for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {

         vGnext[0] += (float)(base_trans_x(i,j)*pressure(i,j));
         vGnext[1] += (float)(base_trans_y(i,j)*pressure(i,j));

         Lznext += (float)(base_rot_z(i,j) * pressure(i,j));
      }
   }

   rigidbody.setLinearVelocity(vGnext);
   rigidbody.setAngularMomentum(Lznext);

}
