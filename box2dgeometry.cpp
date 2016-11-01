#include "box2dgeometry.h"
#include <vector>
#include "util.h"

//check if a point in the object's frame of reference is inside or not
bool Box2DGeometry::is_inside(const Vec2f& point) {
   return point[0] >= -width/2 && point[0] <= width/2 && point[1] >= -height/2 && point[1] <= height/2;
}

//determine the signed distance value (ie. distance to surface, negative if inside)
float Box2DGeometry::get_signed_distance(const Vec2f& point) {
   static std::vector<float> positives;
   positives.clear();

   //Determine distance from all walls
   float x_left = -width/2.0f - point[0];
   float x_right = point[0] - width/2.0f;
   float y_bottom = -height/2.0f - point[1];
   float y_top = point[1] - height/2.0f;

   //Figure out how many are positive
   if(x_left >= 0)
      positives.push_back(x_left);
   if(x_right >= 0)
      positives.push_back(x_right);
   if(y_bottom >= 0)
      positives.push_back(y_bottom);
   if(y_top >= 0)
      positives.push_back(y_top);

   if(positives.size() == 0)
      return max(max(x_left,x_right),max(y_bottom, y_top));
   else if(positives.size() == 1)
      return positives[0];
   else if(positives.size() == 2)
      return sqrt(sqr(positives[0])+ sqr(positives[1]));
   else {
      printf("I screwed this up.");
      return 0;
   }
}

void Box2DGeometry::project_out(Vec2f& point) {
   //Project point to the closest boundary
   
   if(is_inside(point)) {
      const float leftDist = point[0] + width/2;
      const float rightDist = width/2 - point[0];
      const float bottomDist = point[1] + height/2;
      const float topDist = height/2 - point[1];

      if ((bottomDist<=rightDist) && (bottomDist<=leftDist) && (bottomDist<=topDist))
         point[1] = -height/2;
      else if ((rightDist<=leftDist) && (rightDist<=bottomDist) && (rightDist<=topDist))
         point[0] = width/2;
      else if ((leftDist<=rightDist) && (leftDist<=bottomDist) && (leftDist<=topDist))
         point[0] = -width/2;
      else if ((topDist<=leftDist) && (topDist<=bottomDist) && (topDist<=rightDist))
         point[1] = height/2;
   }
   
}

void Box2DGeometry::get_sample_points(std::vector<Vec2f>& verts) {
   float divisions = 3;
   float dw = width / divisions;
   float dh = height / divisions;
   verts.clear();
   for(int i = 0; i <= divisions; ++i) {
      Vec2f point1(-width/2 + i*dw,-height/2);
      Vec2f point2(-width/2 + i*dw, height/2);
      verts.push_back( point1 );
      verts.push_back( point2 );

      if(i != 0 && i != divisions) {
         Vec2f point3(-width/2,-height/2 + i*dh);
         Vec2f point4(width/2,-height/2 + i*dh);
         verts.push_back( point3 );
         verts.push_back( point4 );
      }
   }
}

float Box2DGeometry::get_mass(float density) {
   return density*width*height;
}

void Box2DGeometry::get_inertia_tensor(float density, float& imatrix){
   //zero(imatrix);
   const float mass = get_mass(density);

   const float widthSq = width*width;
   const float heightSq = height*height;
   imatrix = mass / 12.0f * (heightSq + widthSq);
}

void Box2DGeometry::get_vertices(std::vector<Vec2f>& verts) {
   verts.push_back(Vec2f(-width/2, height/2));
   verts.push_back(Vec2f(width/2, height/2));
   verts.push_back(Vec2f(width/2, -height/2));
   verts.push_back(Vec2f(-width/2, -height/2));
}
