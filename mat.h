#ifndef MAT_H
#define MAT_H
 
#include "vec.h"

template<unsigned int M, unsigned int N, class T>
struct Mat
{
   T a[M*N]; // entries stored column by column (FORTRAN layout)

   Mat<M,N,T>(void)
   {}

   template<class S>
   Mat<M,N,T>(const S *source)
   {
      for(unsigned int i=0; i<M*N; ++i) a[i]=source[i];
   }

   Mat<M,N,T>(T a0, T a1, T a2, T a3)
   {
      assert(M*N==4);
      a[0]=a0; a[1]=a1; a[2]=a2; a[3]=a3;
   }

   Mat<M,N,T>(T a0, T a1, T a2, T a3, T a4, T a5)
   {
      assert(M*N==6);
      a[0]=a0; a[1]=a1; a[2]=a2; a[3]=a3; a[4]=a4; a[5]=a5;
   }

   Mat<M,N,T>(T a0, T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8)
   {
      assert(M*N==9);
      a[0]=a0; a[1]=a1; a[2]=a2; a[3]=a3; a[4]=a4; a[5]=a5; a[6]=a6; a[7]=a7; a[8]=a8;
   }

   Mat<M,N,T>(const Vec<M,T> &col0, const Vec<M,T> &col1)
   {
      assert(N==2);
      setcol(0,col0); setcol(1,col1);
   }

   Mat<M,N,T>(const Vec<M,T> &col0, const Vec<M,T> &col1, const Vec<M,T> &col2)
   {
      assert(N==3);
      setcol(0,col0); setcol(1,col1); setcol(2,col2);
   }

   Mat<M,N,T>(const Vec<M,T> &col0, const Vec<M,T> &col1, const Vec<M,T> &col2, const Vec<M,T> &col3)
   {
      assert(N==4);
      setcol(0,col0); setcol(1,col1); setcol(2,col2); setcol(3,col3);
   }

   T &operator()(int i, int j)
   {
      assert(0<=i && (unsigned int)i<M && 0<=j && (unsigned int)j<N);
      return a[i+M*j];
   }

   const T &operator()(int i, int j) const
   {
      assert(0<=i && (unsigned int)i<M && 0<=j && (unsigned int)j<N);
      return a[i+M*j];
   }

   Vec<M,T> col(int j) const
   {
      assert(0<=j && (unsigned int)j<N);
      return Vec<M,T>(a+j*M);
   }

   Vec<N,T> row(int i) const
   {
      assert(0<=i && i<M);
      Vec<N,T> v;
      for(unsigned int j=0; j<N; ++j) v[j]=a(i,j);
      return v;
   }

   Mat<M,N,T> operator+=(const Mat<M,N,T> &b)
   {
      for(unsigned int i=0; i<M*N; ++i) a[i]+=b.a[i];
      return *this;
   }

   Mat<M,N,T> operator+(const Mat<M,N,T> &b) const
   {
      Mat<M,N,T> sum(*this);
      sum+=b;
      return sum;
   }

   Mat<M,N,T> operator-=(const Mat<M,N,T> &b)
   {
      for(unsigned int i=0; i<M*N; ++i) a[i]-=b.a[i];
      return *this;
   }

   Mat<M,N,T> operator-(const Mat<M,N,T> &b) const
   {
      Mat<M,N,T> diff(*this);
      diff-=b;
      return diff;
   }

   Mat<M,N,T> operator*=(T scalar)
   {
      for(unsigned int i=0; i<M*N; ++i) a[i]*=scalar;
      return *this;
   }

   Mat<M,N,T> operator*(T scalar) const
   {
      Mat<M,N,T> b(*this);
      b*=scalar;
      return b;
   }

   Vec<M,T> operator*(const Vec<N,T> v) const
   {
      Vec<M,T> r;
      unsigned int i, j;
      const T *pa, *pv;
      T s, *pr=r.v;
      for(i=0; i<M; ++i, ++pr){
         pa=a+i;
         pv=v.v;
         s=0;
         for(j=0; j<N; ++j, pa+=M, ++pv)
            s+=*pa*(*pv);
         *pr=s;
      }
	  return r;
   }

   template<unsigned int P>
   Mat<M,P,T> operator*(const Mat<N,P,T> b) const
   {
      Mat<M,P,T> c;
      unsigned int i, j, k;
      const T *pa, *pb;
      T s, *pc=c.a;
      for(k=0; k<P; ++k){
         for(i=0; i<M; ++i, ++pc){
            pa=a+i;
            pb=b.a+N*k;
            s=0;
            for(j=0; j<N; ++j, pa+=M, ++pb)
               s+=*pa*(*pb);
            *pc=s;
         }
      }
	  return c;
   }

   Mat<M,N,T> operator/=(T scalar)
   {
      for(unsigned int i=0; i<M*N; ++i) a[i]/=scalar;
      return *this;
   }

   Mat<M,N,T> operator/(T scalar) const
   {
      Mat<M,N,T> b(*this);
      b/=scalar;
      return b;
   }

   Mat<N,M,T> transpose() const {
      Mat<N,M,T> result;
      
      for(unsigned int i = 0; i < M; ++i) {
         for(unsigned int j = 0; j < N; ++j) {
            result(j,i) = (*this)(i,j);
         }
      }
      return result;
   }
};

typedef Mat<2,2,double> Mat22d;
typedef Mat<2,2,float>  Mat22f;
typedef Mat<2,2,int>    Mat22i;
typedef Mat<3,2,double> Mat32d;
typedef Mat<3,2,float>  Mat32f;
typedef Mat<3,2,int>    Mat32i;
typedef Mat<2,3,double> Mat23d;
typedef Mat<2,3,float>  Mat23f;
typedef Mat<2,3,int>    Mat23i;
typedef Mat<3,3,double> Mat33d;
typedef Mat<3,3,float>  Mat33f;
typedef Mat<3,3,int>    Mat33i;
typedef Mat<4,4,double> Mat44d;
typedef Mat<4,4,float>  Mat44f;
typedef Mat<4,4,int>    Mat44i;

// more for human eyes than a good machine-readable format
template<unsigned int M, unsigned int N, class T>
std::ostream &operator<<(std::ostream &out, const Mat<M,N,T> &a)
{
   for(unsigned int i=0; i<M; ++i){
      out<<(i==0 ? '[' : ' ');
      for(unsigned int j=0; j<N-1; ++j)
         out<<a(i,j)<<',';
      out<<a(i,N-1);
      if(i<M-1) out<<';'<<std::endl;
      else out<<']';
   }
   return out;
}

template<unsigned int M, unsigned int N, class T>
inline Mat<M,N,T> operator*(T scalar, const Mat<M,N,T> &a)
{
   Mat<M,N,T> b(a);
   b*=scalar;
   return b;
}

template<unsigned int M, unsigned int N, class T>
inline Mat<M,N,T> outer(const Vec<M,T> &x, const Vec<N,T> &y)
{
   Mat<M,N,T> r;
   T *pr=r.a;
   for(unsigned int j=0; j<N; ++j)
      for(unsigned int i=0; i<M; ++i, ++pr)
         *pr=x[i]*y[j];
   return r;
}

template<unsigned int M, unsigned int N, class T>
inline void zero(Mat<M,N,T>& mat) {
   std::memset(mat.a, 0, N*M*sizeof(T));
}

template<unsigned int N, class T>
inline T trace(Mat<N,N,T>& mat)
{
   T t=0;
   for(unsigned int i=0; i<N; ++i)
      t+=mat.a[(N+1)*i];
   return t;
}

template<class T>
inline Mat<3,3,T> star_matrix(const Vec<3,T> &w)
{
   return Mat<3,3,T>(0, -w.v[2], w.v[1],
                     w.v[2], 0, -w.v[0],
                     -w.v[1], w.v[0], 0);
}

#endif
