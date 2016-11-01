#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include "array2.h"
#include "vec.h"
#include "coordframe2d.h"
#include "rigidgeometry.h"

class RigidBody
{
protected:
   float density;
   float mass;

   RigidGeometry& geometry;

   // Local frame of the solid (tracks orientation and position)
   CoordFrame2D localFrame;

   //Linear velocity
   Vec2f velocity;

   // Angular momentum (= moment cinetique)
   float angularMomentum;

   //Computed from momentum: it's the angular velocity!
   float angularVelocity;

   // Inertia matrix (and its inverse), expressed in local frame, computed from mass and geometric properties
   float imatrix;
   float imatrixInv;
  
public:
   // Constructor
   RigidBody(const float d, RigidGeometry& geom);

   // Compute inverse inertia matrix at point G in the world frame
   void getWorldInvInertiaMatrix(float &imatrix) const;

   // Getting 2D vertices of the box, COUNTER-WISE CLOCK ORDER
   void get2DVertices(std::vector<Vec2f>&verts) const;

   CoordFrame2D getCoordFrame() const;
   void setCoordFrame(CoordFrame2D new_frame);

   const float getMass() const;
   const float getDensity() const;

   // Access to the inertia modulus Jz and inverse
   const float getInertiaModulus() const;
   const float getInvInertiaModulus() const;

   // Access to position of the center of mass
   void getCOM(Vec2f& pos) const;
   void setCOM(const Vec2f& pos);
  
   //Get/set orientation
   float getAngle();
   void setAngle(float angle);

   //Access to linear velocity
   void getLinearVelocity(Vec2f& vel) const;
   void setLinearVelocity(const Vec2f& vel);

   // Access to angular momentum of the solid (in 2d, just the z component)
   void getAngularMomentum(float& Lz) const;
   void setAngularMomentum(const float& Lz);

   Vec2f getPointVelocity(const Vec2f& point);

   // Indicates if point P (in world frame) belongs to the solid
   const bool isInside(const Vec2f& P) const;
   const bool testCollisionAndProject(const Vec2f& P, Vec2f& newPos);
   const float getSignedDist(const Vec2f& P) const;

   // Computing one time step
   void advance(const float dt);

   //integration of position/angle
   void updatePosition(const float dt);
   void updateOrientation(const float dt);
   
   //accelerations
   void applyForce(const Vec2f& force, const float dt);
   void applyTorque(float torque, const float dt);

   void getAngularVelocity(float & wz);



protected:
   //just for internal bookkeeping
   void recomputeAngularVelocity();
  
};

#endif
