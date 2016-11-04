#include "rigidbody.h"
#include <iostream>
#include "mat.h"
#include "vec.h"
#include "gluvi.h"

using namespace std;

RigidBody::RigidBody(const float d, RigidGeometry& geom) : 
   density(d), geometry(geom),
   velocity(0,0), angularMomentum(0), angularVelocity(0)
{
  
  localFrame.set_position(Vec2f(0,0));
  localFrame.set_orientation(0);

  // Initializing inertial matrix
  imatrix = 0;
  geometry.get_inertia_tensor(density, imatrix);
  mass = geometry.get_mass(density);
  // Computing the inverse
  imatrixInv = 1.0f/imatrix;
  
}

void RigidBody::getWorldInvInertiaMatrix(float& worldInvImatrix) const
{
  worldInvImatrix = imatrixInv;
}

void RigidBody::get2DVertices(std::vector<Vec2f>&verts) const
{
 
   geometry.get_vertices(verts);
   for (unsigned int i = 0; i < verts.size(); ++i)
      verts[i] = localFrame.local_to_global(verts[i]);

}

const float RigidBody::getMass() const
{
  return mass;
}

CoordFrame2D RigidBody::getCoordFrame() const{
   return localFrame;
}

void RigidBody::setCoordFrame(CoordFrame2D new_frame) {
   localFrame = new_frame;
} 

void RigidBody::getCOM(Vec2f& pos) const
{
   pos = localFrame.position;
}

void RigidBody::setCOM(const Vec2f& pos) 
{
  localFrame.set_position(pos);
}

void RigidBody::setAngle(float angle) {
   localFrame.set_orientation(angle);
}

float RigidBody::getAngle() {
   return localFrame.orientation;
}

const float RigidBody::getDensity() const
{
  return density;
}

const float RigidBody::getInertiaModulus() const
{
  return imatrix;
}

const float RigidBody::getInvInertiaModulus() const
{
  return imatrixInv;
}

void RigidBody::getLinearVelocity(Vec2f& vel) const
{
  for (int i = 0; i != 2; ++i)
    vel[i] = velocity[i];
}

void RigidBody::getAngularMomentum(float& Lz) const
{
  Lz = angularMomentum;
}

// Setting angular momentum
void RigidBody::setAngularMomentum(const float& Lz)
{
  angularMomentum = Lz;
  recomputeAngularVelocity();
}

void RigidBody::setLinearVelocity(const Vec2f& vel)
{
  for (int i = 0; i != 2; ++i)
    velocity[i] = vel[i];
}

Vec2f RigidBody::getPointVelocity(const Vec2f& point) {
   Vec2f linear;
   float angular;
   getLinearVelocity(linear);
   getAngularVelocity(angular);
   Vec2f centre;
   getCOM(centre);
   Vec3f radius(point[0]-centre[0], point[1]-centre[0], 0);
   Vec3f point_vel = Vec3f(linear[0], linear[1], 0) + cross(Vec3f(0,0,angular), radius);
   
   return Vec2f(point_vel[0], point_vel[1]);
}

const bool RigidBody::isInside(const Vec2f& P) const
{
  // Convert P to the local frame
  Vec2f Plocal = localFrame.global_to_local(P);

  // Return true if Plocal belongs to the geometry
  return geometry.is_inside(Plocal);
}

const bool RigidBody::testCollisionAndProject(const Vec2f& P, Vec2f& newPos)
{

   // Convert P to the local frame
   Vec2f Plocal = localFrame.global_to_local(P);
   geometry.project_out(Plocal);
   newPos = localFrame.local_to_global(Plocal);
   return true;
}


const float RigidBody::getSignedDist(const Vec2f& P) const
{
   // Convert P to the local frame
   Vec2f Plocal = localFrame.global_to_local(P);

   return geometry.get_signed_distance(Plocal);
}

void RigidBody::advance(const float dt)
{
  
  // Updating position/orientation
  updatePosition(dt);
  updateOrientation(dt);

  //Add force and torque
  applyForce(Vec2f(0, -0.1f), dt);
  applyTorque(0.0, dt);
  
}
void RigidBody::getAngularVelocity(float&wz) 
{
   recomputeAngularVelocity(); //just to be safe
   wz = angularVelocity;
}
void RigidBody::applyForce(const Vec2f& force, const float dt)
{
  //Updating velocity with external force
  velocity[0] += force[0] * dt;
  velocity[1] += force[1] * dt;
}

void RigidBody::applyTorque(float torque, const float dt)
{
  // Updating angular momentum
  angularMomentum += torque * dt;
  recomputeAngularVelocity();
}


void RigidBody::updatePosition(const float dt)
{
  // Updating center of mass
  localFrame.position += velocity*dt;
}


void RigidBody::recomputeAngularVelocity()
{
   // Updating angular velocity omega
   // NB: omega = (I^-1)L
   angularVelocity = imatrixInv * angularMomentum; 
   
}

void RigidBody::updateOrientation(const float dt)
{
  float rotation = angularVelocity*dt;
  localFrame.rotate(rotation);
  recomputeAngularVelocity();
}
