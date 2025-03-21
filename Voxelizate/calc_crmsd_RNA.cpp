#include "calc_crmsd_RNA.h"
#include "calc_crmsd_VECTOR.h"

double  calc_crmsd(RESIDUE* native,RESIDUE* conf,int nres){
	vector<VECTOR> s1,s2;

	int i,j,natom=0;

	for(i=0;i<nres;++i){
		for(j=0;j<enumAtomTotNumconst;++j) {
			if( conf[i].atoms[j].name[0]!='\0' && native[i].atoms[j].name[0]!='\0' ) {
				s1.push_back(native[i].atoms[j].r);
				s2.push_back(conf[i].atoms[j].r);
				natom++;
			}
		}
	}

	if(natom==0){
		cout<<"Error in calc_rmsd(), natom=0."<<endl;
		exit(0);
	}

	double rmsd = calc_crmsd_VECTOR(s1,s2,0);

	return rmsd;
}

