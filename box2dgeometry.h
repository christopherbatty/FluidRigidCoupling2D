#ifndef BOX2DGEOMETRY_H
#define BOX2DGEOMETRY_H

#include "rigidgeometry.h"

struct Box2DGeometry : public RigidGeometry {
   
   float width;
   float height;
   //(0,0) = origin is at point 0.5*(width,height)
   
   Box2DGeometry(float w, float h):width(w), height(h) {}
   
   bool is_inside(const Vec2f& point);
   float get_signed_distance(const Vec2f& point);
   void project_out(Vec2f& point);
   void get_sample_points(std::vector<Vec2f>& points);
   
   float get_mass(float density);
   void get_inertia_tensor(float density, float& imatrix);
   
   void get_vertices(std::vector<Vec2f>& verts);
   void get_segments(){}
   void get_triangles() {}
};
#endif
