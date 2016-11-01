#ifndef COORDFRAME2D_H
#define COORDFRAME2D_H

#include "vec.h"

struct CoordFrame2D {

  Vec2f position;
  float orientation;

  CoordFrame2D():position(0,0),orientation(0){}
  CoordFrame2D(Vec2f pos, float ang):position(pos),orientation(ang){}

  void set_position(const Vec2f& _position) {
     position = _position;
  }
  void set_orientation(float f) {
     orientation = f;
  }

  void translate(const Vec2f& vector)
  {
      position += vector;
  }

  void rotate(float angle)
  {
    orientation += angle;
  }
  
  Vec2f global_to_local(const Vec2f& point) const { 
     return ::rotate(point - position, -orientation);
  }

  Vec2f local_to_global(const Vec2f& point) const { 
     return ::rotate(point, orientation) + position;
  }

};

#endif
