#include <iostream>
#include <iomanip>
#include "VECTOR.h"

double	getAngle(VECTOR& a,VECTOR& b,VECTOR& c) {
	double	d=dot(a-b,c-b);
	double	A=getDistance(b,c);
	double	C=getDistance(b,a);
	return	acos(d/A/C)*180.0/PI;
}


double	getDihedralAngle(VECTOR& a,VECTOR& b,VECTOR& c,VECTOR& d) {
	VECTOR	v12=b-a, v23=c-b, v34=d-c;
	VECTOR	n1=v12*v23, n2=v23*v34;
	double	d1=lengthof(n1), d2=lengthof(n2);
	double	dd=dot(n1,n2);
	double	fai=acos( dd/(d1*d2) )*180./PI;
	double	chirality=dot(n1,v34);
	if(chirality<0)
		//fai=-fai;
		fai=360.0-fai;
	return	fai;
}


void	transpose(double a[3][3]) {
	double	t;
	t=a[0][1];
	a[0][1]=a[1][0];
	a[1][0]=t;
	t=a[0][2];
	a[0][2]=a[2][0];
	a[2][0]=t;
	t=a[1][2];
	a[1][2]=a[2][1];
	a[2][1]=t;
}


void	matrix_multiply(const double a[3][3],const double b[3][3],double R[3][3]) {
	double	t[3][3];	//just in case when one of input matrix is the same with output matrix R.
	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j)
			t[i][j]=a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];

	for(int i=0; i<3; ++i)
		for(int j=0; j<3; ++j)
			R[i][j]=t[i][j];
}


ostream& operator << (ostream& output, const VECTOR& a) {
	output.precision(5);
	output<<a.x<<","<<a.y<<","<<a.z;

	return output;
}


ostream& operator << (ostream& output, double a[3][3]) {
	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) cout<<setw(8)<<setprecision(3)<<a[i][j]<<" ";
		cout<<endl;
	}
	return output;
}



