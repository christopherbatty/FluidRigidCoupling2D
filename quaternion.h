#ifndef QUATERNION_H
#define QUATERNION_H

#include "vec.h"
#include "mat.h"

template<class T>
struct Quaternion {
   //See Matrix/Quaternion FAQ for some details
   
   //Quaternion data in [x,y,z,w] order
   Vec<4,T> data;
   
   Quaternion() 
   {
      set_axis_angle(Vec<3,T>(1,0,0), 0);
   }

   Quaternion<T>(Vec<4,T> data_) 
   {
      data = data_;
   }

   Quaternion<T>(const float q0, const float q1, const float q2, const float q3)
   {
         data[0] = q0;
         data[1] = q1;
         data[2] = q2;
         data[3] = q3;
   }
  
   Quaternion(Vec<3,T> axis, T angle) {
      set_axis_angle(axis,angle);
   }

   explicit Quaternion(Vec3f vec) {
      data[0] = vec[0];
      data[1] = vec[1];
      data[2] = vec[2];
      data[3] = 0;
   }

   void set_axis_angle(Vec<3,T> axis, T angle) {
      
      const double norm = mag(axis);
      
      if (norm < 1E-8) {
	      data[0] = 0.0;
         data[1] = 0.0;
         data[2] = 0.0;
         data[3] = 1.0;
      }
      else
      {
         T half_angle = angle / (T)2.0;
	      const double sin_half_angle_over_norm = sin(half_angle)/norm;
	      data[0] = (T)sin_half_angle_over_norm*axis[0];
	      data[1] = (T)sin_half_angle_over_norm*axis[1];
	      data[2] = (T)sin_half_angle_over_norm*axis[2];
	      data[3] = (T)cos(half_angle);
      }
   }

   void get_axis_angle(Vec<3,T>& axis, T& angle) const {
      T cos_a = data[3];
      angle = acos(cos_a) * 2;
      if(angle < 1E-8) {
         axis = Vec<3,T>(1,0,0);
      }
      else {
         T sin_a = sqrt(1-sqr(cos_a));
         if(fabs(sin_a)< 0.0005) sin_a = 1;
         axis[0] = data[0] / sin_a;
         axis[1] = data[1] / sin_a;
         axis[2] = data[2] / sin_a;
      }
   }

   // American format (transpose of European format)
   void get_rotation_matrix(Mat33f& matrix) const {
      T xx = 2 * sqr(data[0]);
      T yy = 2 * sqr(data[1]);
      T zz = 2 * sqr(data[2]);
      T xy = 2 * data[0]*data[1];
      T xz = 2 * data[0]*data[2];
      T yz = 2 * data[1]*data[2];
      T xw = 2 * data[0]*data[3];
      T yw = 2 * data[1]*data[3];
      T zw = 2 * data[2]*data[3];

      matrix(0,0) = 1 - (yy + zz);
      matrix(0,1) = (xy + zw);
      matrix(0,2) = (xz - yw);

      matrix(1,0) = (xy - zw);
      matrix(1,1) = 1 - (xx + zz);
      matrix(1,2) = (yz + xw);

      matrix(2,0) = (xz + yw);
      matrix(2,1) = (yz - xw);
      matrix(2,2) = 1 - (xx + yy);

   }

   // Same get_rotation_matrix: American format
   void get_opengl_rotation_matrix(double matrix[4][4]) const {
      T xx = 2 * sqr(data[0]);
      T yy = 2 * sqr(data[1]);
      T zz = 2 * sqr(data[2]);
      T xy = 2 * data[0]*data[1];
      T xz = 2 * data[0]*data[2];
      T yz = 2 * data[1]*data[2];
      T xw = 2 * data[0]*data[3];
      T yw = 2 * data[1]*data[3];
      T zw = 2 * data[2]*data[3];

      matrix[0][0] = 1 - (yy + zz);
      matrix[0][1] = (xy + zw);
      matrix[0][2] = (xz - yw);

      matrix[1][0] = (xy - zw);
      matrix[1][1] = 1 - (xx + zz);
      matrix[1][2] = (yz + xw);

      matrix[2][0] = (xz + yw);
      matrix[2][1] = (yz - xw);
      matrix[2][2] = 1 - (xx + yy);

      matrix[3][0] = 0.00;
      matrix[3][1] = 0.00;
      matrix[3][2] = 0.00;

      matrix[0][3] = 0;
      matrix[1][3] = 0;
      matrix[2][3] = 0;
      matrix[3][3] = 1.00;
   }

   // European format (transpose of American format, columns are vectors of the eigen basis)
   void get_EU_rotation_matrix(Mat33f& matrix) const {
      get_inv_rotation_matrix(matrix);
   }
  
   void get_inv_rotation_matrix(Mat33f& matrix) const {
      get_rotation_matrix(matrix);
      matrix = matrix.transpose();
   }

   
   Quaternion<T> operator+=(Quaternion<T> &w)
   {
      data+=w.data;
      return *this;
   }
   
   Quaternion<T> operator+(Quaternion<T> a) const
   {
      Quaternion<T> result(*this);
      result += a;
      return result;
   }
   Quaternion<T> operator-=(Quaternion<T> &w)
   {
      data-=w.data;
      return *this;
   }
   
   Quaternion<T> operator-(Quaternion<T> a) const
   {
      Quaternion<T> result(*this);
      result -= a;
      return result;
   }

   Quaternion<T> operator*=(Quaternion<T> a)
   {
      *this = (*this) * a;
      return *this;
   }

   Quaternion<T> operator*(Quaternion<T> a) const
   {
      //TODO Optimize by stripping out intermediate data structures?
      Quaternion<T> result;
      Vec3f v1(data[0],data[1],data[2]);
      Vec3f v2(a.data[0],a.data[1],a.data[2]);
      
      Vec3f c = data[3] * v2 + a.data[3] * v1 + cross(v1, v2);
      result.data[0] = c[0];
      result.data[1] = c[1];
      result.data[2] = c[2];
      result.data[3] = data[3]*a.data[3] - dot(v1,v2);
      
      return result;
   }

   float magnitude() {
      return sqrt(sqr(data[0]) + sqr(data[1])+sqr(data[2]) + sqr(data[3]));
   }
   
   inline void normalize() {
      float mag = magnitude();
      data[0]/=mag;
      data[1]/=mag;
      data[2]/=mag;
      data[3]/=mag;
   }

   Quaternion<T> inverse() const { return Quaternion<T>(-data[0], -data[1], -data[2], data[3]); }

  
   Vec3f rotate(const Vec3f& v) const {
     const double q00 = 2.0l * data[0] * data[0];
     const double q11 = 2.0l * data[1] * data[1];
     const double q22 = 2.0l * data[2] * data[2];

     const double q01 = 2.0l * data[0] * data[1];
     const double q02 = 2.0l * data[0] * data[2];
     const double q03 = 2.0l * data[0] * data[3];

     const double q12 = 2.0l * data[1] * data[2];
     const double q13 = 2.0l * data[1] * data[3];

     const double q23 = 2.0l * data[2] * data[3];

     return Vec3f(
        (float)(1.0 - q11 - q22)*v[0] + (float)(      q01 - q23)*v[1] + (float)(      q02 + q13)*v[2],
	     (float)(      q01 + q23)*v[0] + (float)(1.0 - q22 - q00)*v[1] + (float)(      q12 - q03)*v[2],
	     (float)(      q02 - q13)*v[0] + (float)(      q12 + q03)*v[1] + (float)(1.0 - q11 - q00)*v[2] );
   }

   Vec3f inverseRotate(const Vec3f& v) const
   {
    return inverse().rotate(v);
   }

};

typedef Quaternion<float> Quaternionf;
typedef Quaternion<double> Quaterniond;

template<class T, class U>
inline Quaternion<T> slerp(const Quaternion<T> &value0, const Quaternion<T> &value1, U f) {
   Quaternion<T> result;
   U w1, w2;

   U cos_theta = dot(value0.data, value1.data);
   U theta = acos(cos_theta);
   U sin_theta = sin(theta);
   if(sin_theta > 0.001f) {
      w1 = (U)sin((1.0 - f)*theta) / sin_theta;
      w2 = (U)sin(f*theta) / sin_theta;
   }
   else {
      w1 = (U)(1.0 - f);
      w2 = f;
   }
   result.data = w1*value0.data + w2*value1.data;
   return result;
}

#endif
