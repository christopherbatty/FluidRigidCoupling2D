#ifndef RIGIDGEOMETRY_H
#define RIGIDGEOMETRY_H

#include "vec.h"
#include "mat.h"
#include <vector>

class RigidGeometry {
   //defines the geometry of an object
   //and various operations on that geometry

public:

   // Virtual destructor
   virtual ~RigidGeometry() {};
  
   //check if a point in the object's frame of reference is inside or not
   virtual bool is_inside(const Vec2f& point) = 0;
   
   //return an estimate of the signed distance from this point in object space
   //to the object surface (zero isocontour)
   virtual float get_signed_distance(const Vec2f& point) = 0;
   
   //project the point to the surface if it's inside
   virtual void project_out(Vec2f& point) = 0;

   //return sample points on the surface, typically for collision processing
   virtual void get_sample_points(std::vector<Vec2f>& points) = 0;

   //compute the inertia tensor matrix for the object
   virtual float get_mass(float density) = 0;
   virtual void get_inertia_tensor(float density, float& imatrix) = 0;

   //for rendering
   virtual void get_segments() = 0;
   virtual void get_vertices(std::vector<Vec2f>& verts) = 0;
};


#endif
