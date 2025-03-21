#include "reshapeResidue.h"
#include "pdbWrite.h"

double mP=30.9738,
	   mO=15.9994,
	   mC=12.0107,
	   mN=14.0067;

double atommass[27]={mP,mO,mO,mO,mC,mC,mC,mO,mO,mC,mC,mO,mN,mC,mN,mC,mC,mO,mN,mC,mN,mN,mC,mN,mO,mO,mN};

double qP=1.1662,
       qOP1=-0.776,
       qOP2=-0.776,
       qO5s_1=-0.4989,qO5s_2=-0.6223,
       qC5s=0.0558,
       qC4s=0.1065,
       qC3s=0.2022,
       qO3s_1=-0.5246,qO3s_2=-0.6541,
       qO4s=-0.3548,
       qC1s_A=0.0394,qC1s_G=0.0191,qC1s_U=0.0674,qC1s_C=0.0066,
       qC2s=0.067,
       qO2s=-0.6139,
       qN9_A=-0.0251,qN9_G=0.0492,
       qC8_A=0.2006,qC8_G=0.1374,
       qN7_A=-0.6073,qN7_G=-0.5709,
       qC5_A=0.0515,qC5_G=0.1744,qC5_U=-0.3635,qC5_C=-0.5215,
       qC6_A=0.7009,qC6_G=0.477,qC6_U=-0.1126,qC6_C=0.0053,
       qO6=-0.5597,
       qN1_A=-0.7615,qN1_G=-0.4787,qN1_U=0.0418,qN1_C=-0.0484,
       qC2_A=0.5875,qC2_G=0.7657,qC2_U=0.4687,qC2_C=0.7538,
       qN2=-0.9672,
       qN3_A=-0.6997,qN3_G=-0.6323,qN3_U=-0.3549,qN3_C=-0.7584,
       qC4_A=0.3053,qC4_G=0.1222,qC4_U=0.5952,qC4_C=0.8185,
       qN6=-0.9019,
       qO4=-0.5761,
       qO2_U=-0.5477,qO2_C=-0.6252,
       qN4=-0.953;

void writeMassChargeInfo(vector<RESIDUE>& rna){
	for(int i=0;i<rna.size();i++){
		for(int j=0;j<27;j++){
			rna[i].atoms[j].mass=atommass[j];
		}

		RESIDUE& res=rna[i];
		res.atoms[p].charge=qP;
		res.atoms[op1].charge=qOP1;
		res.atoms[op2].charge=qOP2;
		
		bool isHead=false;
		bool isTail=false;

		if(res.atoms[p].name[0]=='\0'){
			if(i==0)
				isHead=true;
			else{
				if(res.chain!=rna[i-1].chain)
					isHead=true;
			}
		}

		if(i==rna.size()-1)
			isTail=true;
		else{
			if(res.chain!=rna[i+1].chain)
				isTail=true;
		}

		if(!isHead)
			res.atoms[o5s].charge=qO5s_1;
		else
			res.atoms[o5s].charge=qO5s_2;

		if(!isTail)
			res.atoms[o3s].charge=qO3s_1;
		else
			res.atoms[o3s].charge=qO3s_2;

		res.atoms[c5s].charge=qC5s;
		res.atoms[c4s].charge=qC4s;
		res.atoms[c3s].charge=qC3s;
		res.atoms[o4s].charge=qO4s;
		res.atoms[c2s].charge=qC2s;
		res.atoms[o2s].charge=qO2s;
		res.atoms[o6].charge=qO6;
		res.atoms[n2].charge=qN2;
		res.atoms[n6].charge=qN6;
		res.atoms[o4].charge=qO4;
		res.atoms[n4].charge=qN4;

		if(res.name[2]=='A'){
			res.atoms[c1s].charge=qC1s_A;
			res.atoms[n9].charge=qN9_A;
			res.atoms[c8].charge=qC8_A;
			res.atoms[n7].charge=qN7_A;
			res.atoms[c5].charge=qC5_A;
			res.atoms[c6].charge=qC6_A;
			res.atoms[n1].charge=qN1_A;
			res.atoms[c2].charge=qC2_A;
			res.atoms[n3].charge=qN3_A;
			res.atoms[c4].charge=qC4_A;
		}
		else if(res.name[2]=='G'){
			res.atoms[c1s].charge=qC1s_G;
			res.atoms[n9].charge=qN9_G;
			res.atoms[c8].charge=qC8_G;
			res.atoms[n7].charge=qN7_G;
			res.atoms[c5].charge=qC5_G;
			res.atoms[c6].charge=qC6_G;
			res.atoms[n1].charge=qN1_G;
			res.atoms[c2].charge=qC2_G;
			res.atoms[n3].charge=qN3_G;
			res.atoms[c4].charge=qC4_G;
		}
		else if(res.name[2]=='U'){
			res.atoms[c1s].charge=qC1s_U;
			res.atoms[c5].charge=qC5_U;
			res.atoms[c6].charge=qC6_U;
			res.atoms[n1].charge=qN1_U;
			res.atoms[c2].charge=qC2_U;
			res.atoms[n3].charge=qN3_U;
			res.atoms[c4].charge=qC4_U;
			res.atoms[o2].charge=qO2_U;
		}
		else if(res.name[2]=='C'){
			res.atoms[c1s].charge=qC1s_C;
			res.atoms[c5].charge=qC5_C;
			res.atoms[c6].charge=qC6_C;
			res.atoms[n1].charge=qN1_C;
			res.atoms[c2].charge=qC2_C;
			res.atoms[n3].charge=qN3_C;
			res.atoms[c4].charge=qC4_C;
			res.atoms[o2].charge=qO2_C;
		}
		else{
			cout<<"wrong residue name "<<res.name[2]<<endl;
			exit(0);	
		}
	}
}
void calc_local_reference(RESIDUE& res,VECTOR& O,VECTOR& vx,VECTOR& vy,VECTOR& vz){
	O=res.atoms[c1s].r;
	if(res.name[2]=='A' || res.name[2]=='G')
		vx=res.atoms[n9].r-O;
	else
		vx=res.atoms[n1].r-O;
	vy=(res.atoms[o5s].r+res.atoms[c5s].r)/2.0-O;
	vz=vx*vy;
	vz/=lengthof(vz);
	vx/=lengthof(vx);
	vy=vz*vx;
	vy/=lengthof(vy);
}

void ExtractAtomsInResidueXCenteredBox(vector<RESIDUE>& native,vector<RESIDUE>& rna,int res_index,vector<ATOM>& atomsinbox,double& score){
	if(res_index<0 || res_index>=rna.size()){
		cout<<"illegal res_index "<<res_index<<endl;
		exit(0);
	}

	RESIDUE& res1=rna[res_index];

	VECTOR O,vx,vy,vz;
	calc_local_reference(res1,O,vx,vy,vz);

	VECTOR r;
	double dx,dy,dz;

//	vector<RESIDUE> newrna=rna;

	vector<VECTOR> s1,s2;
	vector<VECTOR> t1,t2;

	for(int i=0;i<rna.size();i++){
		RESIDUE& res2=rna[i];

		r=res2.atoms[c1s].r-O;
		dx=fabs(dot(r,vx));
		if(dx>TrialHalfBoxSideLength)
			continue;
		dy=fabs(dot(r,vy));
		if(dy>TrialHalfBoxSideLength)
			continue;
		dz=fabs(dot(r,vz));
		if(dz>TrialHalfBoxSideLength)
			continue;
		
		for(int j=0;j<27;j++){
			ATOM& a=res2.atoms[j];
			if(a.name[0]=='\0')
				continue;
			
			r=a.r-O;
			dx=dot(r,vx);
			if(fabs(dx)>=HalfBoxSideLength)
				continue;
			dy=dot(r,vy);
			if(fabs(dy)>=HalfBoxSideLength)
				continue;
			dz=dot(r,vz);
			if(fabs(dz)>=HalfBoxSideLength)
				continue;

			ATOM b=a;
			b.r.x=dx+HalfBoxSideLength;
			b.r.y=dy+HalfBoxSideLength;
			b.r.z=dz+HalfBoxSideLength;

			atomsinbox.push_back(b);

			if(i==res_index){
				s1.push_back(native[i].atoms[j].r);
				s2.push_back(a.r);
			}
			else{
				t1.push_back(native[i].atoms[j].r);
				t2.push_back(a.r);
			}
			
			/*
			newrna[i].atoms[j].r.x=dx;
			newrna[i].atoms[j].r.y=dy;
			newrna[i].atoms[j].r.z=dz;*/
            //cout<<a.name<<"  "<<a.mass<<"  "<<a.charge<<" "<<res2.ridx<<endl;
		}
	}

	double Score[2];
	calc_score(s1,s2,t1,t2,Score[0],Score[1]); 
	score=Score[0]+Score[1];

	/*
	PDBW pdbw("newpdb.pdb");
	for(int i=0;i<newrna.size();i++){
		pdbw.writeresidue(newrna[i],newrna[i].chain,newrna[i].ridx);
	}
	pdbw.outTerm();	*/
}

void lattice1Dpoint(double& r,int x[2],double p[2]){
	int x1,x2;
	double p1,p2;
	x1=int((r+HalfBinWidth)/BinWidth)-1;
	p1=((x1+1)*BinWidth+HalfBinWidth-r)/BinWidth;
	x2=x1+1;
	p2=1.0-p1;

	x[0]=x1;
	x[1]=x2;
	p[0]=p1;
	p[1]=p2;

	if(x1<-1 || x2>nBins){
		cout<<"Atoms overflow the box."<<endl;
		exit(0);
	}
}

void lattice3Dpoint(VECTOR& r,int x[2],int y[2],int z[2],double px[2],double py[2],double pz[2]){
	lattice1Dpoint(r.x,x,px);
	lattice1Dpoint(r.y,y,py);
	lattice1Dpoint(r.z,z,pz);
}	

void pixelateAtomsInBox(vector<ATOM>& atomsinbox, double pixels[3][nBins][nBins][nBins]){
	for(int i=0;i<3;i++){
		for(int j=0;j<nBins;j++){
			for(int k=0;k<nBins;k++){
				for(int m=0;m<nBins;m++){
					pixels[i][j][k][m]=0.0;
				}
			}
		}
	}

	int x[2],y[2],z[2];
	double px[2],py[2],pz[2];

	for(int i=0;i<atomsinbox.size();i++){
		ATOM& a=atomsinbox[i];

		lattice3Dpoint(a.r,x,y,z,px,py,pz);

        //cout<<a.name<<"  "<<x[0]<<"  "<<x[1]<<"  "<<y[0]<<"  "<<y[1]<<"  "<<z[0]<<"  "<<z[1]<<"  "<<px[0]<<"  "<<px[1]<<"  "<<py[0]<<"  "<<py[1]<<"  "<<pz[0]<<"  "<<pz[1]<<" "<<a.mass<<" "<<a.charge<<endl;

		for(int j=0;j<2;j++){
			if(x[j]<0 || x[j]>=nBins)
				continue;
			for(int k=0;k<2;k++){
				if(y[k]<0 || y[k]>=nBins)
					continue;
				for(int m=0;m<2;m++){
						if(z[m]<0 || z[m]>=nBins)
							continue;

						pixels[0][ x[j] ][ y[k] ][ z[m] ]+=px[j]*py[k]*pz[m];
						pixels[1][ x[j] ][ y[k] ][ z[m] ]+=px[j]*py[k]*pz[m]*a.mass;
						pixels[2][ x[j] ][ y[k] ][ z[m] ]+=px[j]*py[k]*pz[m]*a.charge;	

                        /*
                        if(x[j]==23 && y[k]==6 && z[m]==3){
                            cout<<a.name<<" "<<a.mass<<" "<<a.charge<<" "<<px[j]<<" "<<py[k]<<" "<<pz[m]<<endl;
                        }*/
				}
			}
		}
	}
}

void pixelateResidue(vector<RESIDUE>& native,vector<RESIDUE>& rna,int res_index,double pixels[3][nBins][nBins][nBins],double& score){
	vector<ATOM> atomsinbox;
	ExtractAtomsInResidueXCenteredBox(native,rna,res_index,atomsinbox,score);
	pixelateAtomsInBox(atomsinbox,pixels);
}
		
