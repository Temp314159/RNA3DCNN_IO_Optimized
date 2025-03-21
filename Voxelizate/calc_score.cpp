#include <iostream>
#include <iomanip>
#include <math.h>

#include "calc_score.h"

//these two head files must be included AFTER the others, otherwise it will cause wierd errors.
#include "f2c.h"
#include "clapack.h"

#define NTATOM 15000

//the 1st dimension is from 0-4, but the 0 position is useless, this is to be compatible with fortran.
static	double	vm[4][NTATOM],vp[4][NTATOM];
static	double	A[5][5];	//the first position 0 is useless, in order to be compatible with FORTRAN.
static	double	B[4*4];

void calc_score(vector<VECTOR>& s1,vector<VECTOR>& s2,vector<VECTOR>& t1,vector<VECTOR>& t2,double& rmsd1,double& rmsd2){
	if(s1.size()!=s2.size() || t1.size()!=t2.size()){cout<<"not equal atom number."<<endl; exit(0);}

	int natom1=s1.size(),natom2=t1.size();
	if(natom1>NTATOM) {cout<<"Error too small NTATOM."<<endl; exit(0);}

	int i,j;

	VECTOR c1(0.0),c2(0.0);
	for(j=0;j<natom1;++j) {
		c1+=s1[j];
		c2+=s2[j];
	}
	c1/=(double)natom1;
	c2/=(double)natom1;
	for(j=0;j<natom1;++j) {
		s1[j]-=c1;
		s2[j]-=c2;
	}
	for(j=0;j<natom2;j++){
		t1[j]-=c1;
		t2[j]-=c2;
	}		

	for(i=0;i<natom1;++i) {
		vm[1][i]=s1[i].x-s2[i].x;
		vp[1][i]=s1[i].x+s2[i].x;
		vm[2][i]=s1[i].y-s2[i].y;
		vp[2][i]=s1[i].y+s2[i].y;
		vm[3][i]=s1[i].z-s2[i].z;
		vp[3][i]=s1[i].z+s2[i].z;
	}

	double tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vm[1][i]*vm[1][i]+vm[2][i]*vm[2][i]+vm[3][i]*vm[3][i];
	A[1][1]=tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vp[2][i]*vm[3][i]-vm[2][i]*vp[3][i];
	A[2][1]=tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vm[1][i]*vp[3][i]-vp[1][i]*vm[3][i];
	A[3][1]=tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vp[1][i]*vm[2][i]-vm[1][i]*vp[2][i];
	A[4][1]=tmp;

	A[1][2]=A[2][1];

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vp[2][i]*vp[2][i]+vp[3][i]*vp[3][i]+vm[1][i]*vm[1][i];
	A[2][2]=tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vm[1][i]*vm[2][i]-vp[1][i]*vp[2][i];
	A[3][2]=tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vm[1][i]*vm[3][i]-vp[1][i]*vp[3][i];
	A[4][2]=tmp;

	A[1][3]=A[3][1];

	A[2][3]=A[3][2];

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vp[1][i]*vp[1][i]+vp[3][i]*vp[3][i]+vm[2][i]*vm[2][i];
	A[3][3]=tmp;

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vm[2][i]*vm[3][i]-vp[2][i]*vp[3][i];
	A[4][3]=tmp;

	A[1][4]=A[4][1];

	A[2][4]=A[4][2];

	A[3][4]=A[4][3];

	tmp=0.0;
	for(i=0;i<natom1;++i)
		tmp+= vp[1][i]*vp[1][i]+vp[2][i]*vp[2][i]+vm[3][i]*vm[3][i];
	A[4][4]=tmp;

	int	ip=0;
	for(i=1;i<=4;++i)
		for(j=1;j<=4;++j) {
			B[ip]=A[j][i];
			ip++;
		}

	char	JOBZ='V',UPLO='L';
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

	rmsd1=1e10;
	int	imin;
	for(i=0;i<N;++i) {
		if(rmsd1>abs(W[i])) {
			rmsd1=abs(W[i]);
			imin=i;
		}
	}
	rmsd1=sqrt(rmsd1/natom1);

	{
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

		for(i=0;i<natom2;++i) {
			double x=t2[i].x, y=t2[i].y, z=t2[i].z;
			t2[i].x=R[1][1]*x+R[1][2]*y+R[1][3]*z;
			t2[i].y=R[2][1]*x+R[2][2]*y+R[2][3]*z;
			t2[i].z=R[3][1]*x+R[3][2]*y+R[3][3]*z;
		}
	
		/*
		for(i=1;i<=3;i++){
			for(j=1;j<=3;j++){
				cout<<R[i][j]<<" ";
			}
			cout<<endl;
		}*/
	}

	rmsd2=0.0;

	for(i=0;i<natom2;++i){
		rmsd2+=getDistanceSquare(t1[i],t2[i]);
	}
	rmsd2=sqrt(rmsd2/natom2);
		
}


