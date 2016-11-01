#ifndef COORDFRAME_H
#define COORDFRAME_H

#include "quaternion.h"
#include "vec.h"
#include "mat.h"
#include "gluvi.h"

//Personally, I'd rather not have gluvi.h included here, 
//so we avoid having opengl-specific stuff in our general graphics data structure routines
//-CB

struct CoordFrame {

  Vec3f position;
  Quaternion<float> orientation;

  //It's a little silly to have these accessors when the data is public
  void set_position(const Vec3f& pos)
  {
    position = pos;
  }

  void set_orientation(const Quaternion<float>& q)
  {
    orientation = q;
  }
  
  //TODO What is the difference between these two functions?
  void rotate(const Quaternion<float>& rotation)
  {
    Quaternion<float> qbis = rotation;
    rotate(qbis);
    //orientation *= rotation;
  }

  void rotate(Quaternion<float>& q)
  {
    orientation *= q;
    orientation.normalize(); // Prevents numerical drift
  }

  Vec3f global_to_local(Vec3f point) const { 
    return orientation.inverseRotate(point - position);
  }

  Vec3f local_to_global(Vec3f point) const { 
    return orientation.rotate(point) + position;
  }

  Mat44f mat44f() const {

     Mat44f translation;
     zero(translation);
     translation(0,0) = 1;
     translation(1,1) = 1;
     translation(2,2) = 1;
     translation(3,3) = 1;
     translation(0,3) = position[0];
     translation(1,3) = position[1];
     translation(2,3) = position[2];

     Mat33f r3x3;
     orientation.get_rotation_matrix(r3x3);
     Mat44f rotation;
     zero(rotation);
     for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
           rotation(i,j) = r3x3(i,j);
        }
     }
     rotation(3,3) = 1;

     return translation * rotation;
  }

  const GLdouble* glmatrix() const
  {
     static GLdouble m[4][4];
     get_glmatrix(m);
     return (const GLdouble*)(m);
  }

  // GLdouble[4][4] version of matrix(). See also getWorldMatrix() and matrix(). */
  void get_glmatrix(GLdouble m[4][4]) const
  {
     orientation.get_opengl_rotation_matrix(m);

     m[3][0] = position[0];
     m[3][1] = position[1];
     m[3][2] = position[2];
  }
  
};

#endif
