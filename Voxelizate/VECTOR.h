#ifndef	__VECTORDEF__

#define	__VECTORDEF__

#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

#define PI 3.141592653589793

struct	VECTOR {
	double	x,y,z;

	VECTOR(double c=0) : x(c),y(c),z(c) {}
	VECTOR(double a,double b,double c) : x(a), y(b), z(c) {}
	VECTOR  operator + (const VECTOR& b);
	VECTOR  operator - (const VECTOR& b);
	VECTOR	operator - ();
	VECTOR  operator * (const double c);
	VECTOR	operator * (const VECTOR& b);
	VECTOR  operator / (const double c);
	VECTOR&	operator = (const double c);
	VECTOR&	operator += (const VECTOR& b);
	VECTOR&	operator -= (const VECTOR& b);
	VECTOR&	operator /=(const double c);
	friend VECTOR	operator * (const double c,const VECTOR& b);
	friend VECTOR	operator * (const double R[3][3],const VECTOR& a);
};


inline	VECTOR	VECTOR::operator + (const VECTOR& b) {
	VECTOR	c;
	c.x=x+b.x;
	c.y=y+b.y;
	c.z=z+b.z;
	return c;
}

inline	VECTOR	VECTOR::operator - (const VECTOR& b) {
	VECTOR	c;
	c.x=x-b.x;
	c.y=y-b.y;
	c.z=z-b.z;
	return c;
}

inline  VECTOR  VECTOR::operator - () {
	VECTOR	c(*this);
	c.x=-c.x;
	c.y=-c.y;
	c.z=-c.z;
	return c;
}

inline	VECTOR	VECTOR::operator * (const double c) {
	VECTOR	b;
	b.x=x*c;
	b.y=y*c;
	b.z=z*c;
	return b;
}

inline	VECTOR	VECTOR::operator * (const VECTOR& b) {
	VECTOR	c;
	c.x=y*b.z-z*b.y;
	c.y=z*b.x-x*b.z;
	c.z=x*b.y-y*b.x;
	return	c;
}

inline  VECTOR  VECTOR::operator / (const double c) {
	VECTOR	b;
	b.x=x/c;
	b.y=y/c;
	b.z=z/c;
	return b;
}

inline	VECTOR&	VECTOR::operator = (const double c) {
	x=c;
	y=c;
	z=c;
	return *this;
}

inline	VECTOR&	VECTOR::operator += (const VECTOR& b) {
	x+=b.x;
	y+=b.y;
	z+=b.z;
	return *this;
}

inline	VECTOR&	VECTOR::operator -= (const VECTOR& b) {
	x-=b.x;
	y-=b.y;
	z-=b.z;
	return *this;
}

inline	VECTOR&	VECTOR::operator /=(const double c) {
	x/=c;
	y/=c;
	z/=c;
	return *this;
}

inline	VECTOR	operator * (const double c,const VECTOR& a) {
	VECTOR	b;
	b.x=a.x*c;
	b.y=a.y*c;
	b.z=a.z*c;
	return b;
}

inline	VECTOR	operator * (const double R[3][3],const VECTOR& a) {
	VECTOR	b;
	b.x=R[0][0]*a.x+R[0][1]*a.y+R[0][2]*a.z;
	b.y=R[1][0]*a.x+R[1][1]*a.y+R[1][2]*a.z;
	b.z=R[2][0]*a.x+R[2][1]*a.y+R[2][2]*a.z;
	return	b;
}


inline	double	dot(const VECTOR& a,const VECTOR& b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}


inline  double  getDistance(const VECTOR& a,const VECTOR& b) {
	double	x=a.x-b.x, y=a.y-b.y, z=a.z-b.z;
	return sqrt(x*x+y*y+z*z);
}


inline  double  getDistanceSquare(const VECTOR& a,const VECTOR& b) {
	double	x=a.x-b.x, y=a.y-b.y, z=a.z-b.z;
	return x*x+y*y+z*z;
}


inline	double	lengthof(const VECTOR& a) {
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

double	getAngle(VECTOR& a,VECTOR& b,VECTOR& c);
double	getDihedralAngle(VECTOR& a,VECTOR& b,VECTOR& c,VECTOR& d);
void	transpose(double a[3][3]);
void	matrix_multiply(const double a[3][3],const double b[3][3],double R[3][3]);
ostream& operator << (ostream& output, const VECTOR& a);
ostream& operator << (ostream& output, double a[3][3]);

inline double myabs(double x) {
	if(x>0) return x;
	return -x;
}

#endif

