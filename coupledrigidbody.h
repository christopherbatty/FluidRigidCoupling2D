#ifndef COUPLEDRIGIDBODY_H
#define COUPLEDRIGIDBODY_H

#include "rigidbody.h"
#include "pcgsolver/sparse_matrix.h"
#include "array2.h"

class CoupledRigidBody {
public:
   int ni, nj;
   float dx;
   
   RigidBody& rigidbody;

   Array2f phi;
   Array2f u_vol, v_vol;
   float u_mass, v_mass;
   
   CoupledRigidBody(RigidBody& rb, int nx_, int nj_, float dx_);
   void update_grid_data();
   void add_coupling(SparseMatrix<double>& matrix, Array2d& rhs);
   void apply_pressure(const Array2d& pressure);
   
   bool is_inside(const Vec2f& point);
   const bool test_collision_and_project(const Vec2f& P, Vec2f& newPos);
   
protected:
   void recompute_phi();
   void recompute_volumes();
   void recompute_masses();

   void build_base_vectors();
   void add_coefficients(SparseMatrixd& sparseMat);
   void add_right_hand_side(Array2d& rhs);

   Array2d base_trans_x, base_trans_y;    //used to compute the translation update of the solid
   Array2f base_rot_z;                    //used to compute the rotation update of the solid

};


#endif
