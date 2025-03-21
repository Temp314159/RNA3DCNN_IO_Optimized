#include <iostream>
#include <iomanip>
#include <math.h>

#include "calc_crmsd_VECTOR.h"

//these two head files must be included AFTER the others, otherwise it will cause wierd errors.
#include "f2c.h"
#include "clapack.h"

#define NTATOM 15000

//the 1st dimension is from 0-4, but the 0 position is useless, this is to be compatible with fortran.
static	double	vm[4][NTATOM],vp[4][NTATOM];
static	double	A[5][5];	//the first position 0 is useless, in order to be compatible with FORTRAN.
static	double	B[4*4];

double	calc_crmsd_VECTOR(vector<VECTOR>& s1, vector<VECTOR>& s2,const int rotation){
	if(s1.size()!=s2.size()){
		cout<<"not equal atoms."<<endl;
		exit(0);
	}
	
	int natom=s1.size();

	if(natom>NTATOM) {cout<<"Error too small NTATOM."<<endl; exit(0);}

	int i,j;

	VECTOR c1,c2;
	c1.x=0; c1.y=0; c1.z=0;
	c2.x=0; c2.y=0; c2.z=0;
	for(j=0;j<natom;++j) {
		c1.x += s1[j].x;
		c1.y += s1[j].y;
		c1.z += s1[j].z;
		c2.x += s2[j].x;
		c2.y += s2[j].y;
		c2.z += s2[j].z;
	}

	c1.x/=natom; c1.y/=natom; c1.z/=natom;
	c2.x/=natom; c2.y/=natom; c2.z/=natom;
	for(j=0;j<natom;++j) {
		s1[j].x -= c1.x;
		s1[j].y -= c1.y;
		s1[j].z -= c1.z;
		s2[j].x -= c2.x;
		s2[j].y -= c2.y;
		s2[j].z -= c2.z;
	}

	for(i=0;i<natom;++i) {
		vm[1][i]=s1[i].x-s2[i].x;
		vp[1][i]=s1[i].x+s2[i].x;
		vm[2][i]=s1[i].y-s2[i].y;
		vp[2][i]=s1[i].y+s2[i].y;
		vm[3][i]=s1[i].z-s2[i].z;
		vp[3][i]=s1[i].z+s2[i].z;
	}

	double tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vm[1][i]*vm[1][i]+vm[2][i]*vm[2][i]+vm[3][i]*vm[3][i];
	A[1][1]=tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vp[2][i]*vm[3][i]-vm[2][i]*vp[3][i];
	A[2][1]=tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vm[1][i]*vp[3][i]-vp[1][i]*vm[3][i];
	A[3][1]=tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vp[1][i]*vm[2][i]-vm[1][i]*vp[2][i];
	A[4][1]=tmp;

	A[1][2]=A[2][1];

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vp[2][i]*vp[2][i]+vp[3][i]*vp[3][i]+vm[1][i]*vm[1][i];
	A[2][2]=tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vm[1][i]*vm[2][i]-vp[1][i]*vp[2][i];
	A[3][2]=tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vm[1][i]*vm[3][i]-vp[1][i]*vp[3][i];
	A[4][2]=tmp;

	A[1][3]=A[3][1];

	A[2][3]=A[3][2];

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vp[1][i]*vp[1][i]+vp[3][i]*vp[3][i]+vm[2][i]*vm[2][i];
	A[3][3]=tmp;

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vm[2][i]*vm[3][i]-vp[2][i]*vp[3][i];
	A[4][3]=tmp;

	A[1][4]=A[4][1];

	A[2][4]=A[4][2];

	A[3][4]=A[4][3];

	tmp=0.0;
	for(i=0;i<natom;++i)
		tmp+= vp[1][i]*vp[1][i]+vp[2][i]*vp[2][i]+vm[3][i]*vm[3][i];
	A[4][4]=tmp;

	int	ip=0;
	for(i=1;i<=4;++i)
		for(j=1;j<=4;++j) {
			B[ip]=A[j][i];
			ip++;
		}

	char	JOBZ,UPLO='L';
	if(rotation==1) JOBZ='V';
	else JOBZ='N';
	integer	N=4,LWORK=1000,outinfo;
	double	W[4],WORK[1000];

	dsyev_(&JOBZ,&UPLO,&N,B,&N,W,WORK,&LWORK,&outinfo);
	if(outinfo!=0) {cout<<"calling lapack failed!"; exit(0);}

/*	
	   cout<<"Matrix:"<<endl;
	   for(i=1;i<=4;++i) {
	   for(j=1;j<=4;++j) cout<<A[j][i]<<" ";
	   cout<<endl;
	   }
	   cout<<"Eigenvalues:"<<endl;
	   for(i=0;i<4;++i) cout<<W[i]<<" "; cout<<endl;*/
	   

	double	rmsd=1e10;
	int	imin;
	for(i=0;i<N;++i) {
		if(rmsd>abs(W[i])) {
			rmsd=abs(W[i]);
			imin=i;
		}
	}
	rmsd=sqrt(rmsd/natom);
	//cout<<"rmsd, imin "<<rmsd<<" "<<imin<<endl;

	if(rotation==1) {
		//convert the matrix from 1D to 2D.
		ip=0;
		for(i=1;i<=4;++i)
			for(j=1;j<=4;++j) {
				A[j][i]=B[ip];
				ip++;
			}	
		//get the eigenvector
		double	Q[5];	//the zero position is useless
		for(i=1;i<=4;++i) Q[i]=A[i][imin+1]; //imin+1 because fortran matrix starts from 1

		double	R[4][4];
		R[1][1]=Q[1]*Q[1]+Q[2]*Q[2]-Q[3]*Q[3]-Q[4]*Q[4];
		R[2][2]=Q[1]*Q[1]-Q[2]*Q[2]+Q[3]*Q[3]-Q[4]*Q[4];
		R[3][3]=Q[1]*Q[1]-Q[2]*Q[2]-Q[3]*Q[3]+Q[4]*Q[4];

		R[1][2]=2.0*( Q[2]*Q[3]+Q[1]*Q[4] );
		R[2][1]=2.0*( Q[2]*Q[3]-Q[1]*Q[4] );

		R[1][3]=2.0*( Q[2]*Q[4]-Q[1]*Q[3] );
		R[3][1]=2.0*( Q[2]*Q[4]+Q[1]*Q[3] );

		R[2][3]=2.0*( Q[3]*Q[4]+Q[1]*Q[2] );
		R[3][2]=2.0*( Q[3]*Q[4]-Q[1]*Q[2] );

		for(i=0;i<natom;++i) {
			double x=s2[i].x, y=s2[i].y, z=s2[i].z;
			s2[i].x=R[1][1]*x+R[1][2]*y+R[1][3]*z;
			s2[i].y=R[2][1]*x+R[2][2]*y+R[2][3]*z;
			s2[i].z=R[3][1]*x+R[3][2]*y+R[3][3]*z;
			//shift the center of the second structure to that of the first structure.
			s2[i].x+=c1.x;
			s2[i].y+=c1.y;
			s2[i].z+=c1.z;
		}
	}

	return rmsd;
}


