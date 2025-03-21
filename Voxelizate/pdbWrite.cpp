#include "pdbWrite.h"
#include <iostream>
#include <iomanip>

PDBW::PDBW(char* filename) {
	out.open(filename);
	if(!out) {
		cout<<"Error open file to write pdb."<<endl;
		exit(0);
	}
	iatom=1;
	iresinit=1;
	ires=iresinit;
	imodel=1;
}


PDBW::PDBW(char* filename,int ir,int im) {
	out.open(filename);
	if(!out) {
		cout<<"Error open file to write pdb."<<endl;
		exit(0);
	}
	iatom=1;
	iresinit=ir;
	ires=ir;
	imodel=im;
}


PDBW::~PDBW() {
	out.close();
}


void	PDBW::reset() {
	iatom=1;
	ires=iresinit;
	//not reset imodel
}

void PDBW::writeatom(ATOM& a,char* rname,char chainname,int ridx) {
	if(a.name[0]=='\0') return;
	if(a.name[3]=='*') a.name[3]='\'';
	out.setf(ios::fixed);
	out<<"ATOM  "<<setw(5)<<iatom<<" "<<a.name<<" "<<rname<<setw(2)<<chainname<<setw(4)<<ridx<<"    "
		<<setw(8)<<setprecision(3)<<a.r.x
		<<setw(8)<<setprecision(3)<<a.r.y
		<<setw(8)<<setprecision(3)<<a.r.z<<endl;
	iatom++;
}

void	PDBW::writeresidue(RESIDUE& res) {
	for(int j=0; j<enumAtomTotNumconst; ++j) writeatom(res.atoms[j],res.name,res.chain,res.ridx);
	ires++;
}

void    PDBW::writeresidue(RESIDUE& res,char chainname,int ridx) {
	for(int j=0; j<enumAtomTotNumconst; ++j) writeatom(res.atoms[j],res.name,chainname,ridx);
	ires++;
}

void    PDBW::writeRNA(vector<RESIDUE>& rna,vector<char>& chainname,vector<int>& ridx){
	outModel();
	for(int j=0;j<rna.size();++j){
		writeresidue(rna[j],chainname[j],ridx[j]);
	}
	outTerm();
}

void    PDBW::writeRNA(vector<RESIDUE>& rna){
	outModel();
	for(int j=0;j<rna.size();++j){
		writeresidue(rna[j]);
	}
	outTerm();
}

void	PDBW::outModel() {
	out.setf(ios::fixed);
	out<<"MODEL     "<<setw(4)<<imodel<<endl;
	iatom=1;
}

void	PDBW::outTerm() {
	out<<"TER "<<endl<<"ENDMDL"<<endl;
	imodel++;
}

