#include <iostream>
#include <iomanip>
#include <string.h>
#include "pdb.h"

enumATOMname AnucleotideAtoms[AnucleotideAtomNum]= {P,OP1,OP2,O5s,C5s,C4s,C3s,O3s,O4s,C1s,C2s,O2s,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4};
enumATOMname GnucleotideAtoms[GnucleotideAtomNum]= {P,OP1,OP2,O5s,C5s,C4s,C3s,O3s,O4s,C1s,C2s,O2s,N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4};
enumATOMname UnucleotideAtoms[UnucleotideAtomNum]= {P,OP1,OP2,O5s,C5s,C4s,C3s,O3s,O4s,C1s,C2s,O2s,N1,C6,C5,C4,O4,N3,C2,O2};
enumATOMname CnucleotideAtoms[CnucleotideAtomNum]= {P,OP1,OP2,O5s,C5s,C4s,C3s,O3s,O4s,C1s,C2s,O2s,N1,C6,C5,C4,N4,N3,C2,O2};
enumATOMname BkboneSugarAtoms[BkboneSugarAtomNum]= {P,OP1,OP2,O5s,C5s,C4s,C3s,O3s,O4s,C1s,C2s,O2s};
enumATOMname RingAtoms[RingAtomNum]= {N1,C2,N3,C4,C5,C6};
enumATOMname ImidazoleRingAtoms[ImidazoleRingAtomNum]= {C4,N9,C8,N7,C5};

enumATOMname AbaseAtoms[AbaseAtomNum]= {N9,C8,N7,C5,C6,N6,N1,C2,N3,C4};
enumATOMname GbaseAtoms[GbaseAtomNum]= {N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4};
enumATOMname UbaseAtoms[UbaseAtomNum]= {N1,C6,C5,C4,O4,N3,C2,O2};
enumATOMname CbaseAtoms[CbaseAtomNum]= {N1,C6,C5,C4,N4,N3,C2,O2};

enumATOMname p=P,op1=OP1,op2=OP2,o5s=O5s,c5s=C5s,c4s=C4s,c3s=C3s,o3s=O3s,o4s=O4s,c1s=C1s,c2s=C2s,o2s=O2s,n9=N9,c8=C8,n7=N7,c5=C5,c6=C6,o6=O6,n1=N1,c2=C2,n2=N2,n3=N3,c4=C4,n6=N6,o4=O4,o2=O2,n4=N4;

static int residueType(RESIDUE& res);
static int residueType(RESIDUEpdb& res);

PDB::PDB() {
	if(enumAtomTotNumconst!=enumAtomEND+1) {
		cout<<"error init PDB()"<<endl;
		cout<<enumAtomTotNumconst<<" "<<enumAtomEND+1<<endl;
		exit(0);
	}
	//check consistency
	if(nucA!=0 || nucG!=1 || nucU!=2 || nucC!=3) {
		cout<<"program in danger, unexpected default behaviors of the current compiler"<<endl;
		exit(0);
	}
}

int     PDB::readpdb(const char *pdbname,vector<RESIDUE>& rlist,vector<int>& chain) {
	/*
	   Read atoms, store these in the same residue together, and check how many chains are there in the pdb. The chain[i] stores of the position of the first residue in the i-th chain.
	   the chain[0] is always set to 0, and the chain[nres] is always set to nres, where nres is the number of residues. Therefore, the residues belonging to the i-th chain will be
	   indexed between chain[i] to chain[i+1]-1 in the vector<RESIDUEpdb> residues.
	   */

	vector<RESIDUEpdb>  residues;
	residues.clear();
	chain.clear();

	readpdbresidues(pdbname,residues);

	chain.push_back(0);
	int     nres=residues.size(),ip=0;
	for(int i=1; i<nres; ++i) {
		if(residues[i].atoms[0].chain!=residues[ip].atoms[0].chain) {
			chain.push_back(i);
			ip=i;
		}
	}
	chain.push_back(nres);

	//rearrange atoms, so that we can access them directly according to their names
	for(int i=0; i<nres; ++i) {
		RESIDUEpdb&   rold=residues[i];
		RESIDUE       r;
		strcpy(r.name,rold.name);
		for(int j=0; j<enumAtomTotNumconst; ++j) {
			strcpy(r.atoms[j].name,"\0");
			r.atoms[j].r=0;
		}
		if(1) {
			strcpy(r.pn,pdbname);
			r.ridx=rold.atoms[0].ridx;
			int natom=rold.atoms.size();
			for(int j=0; j<natom; ++j) {
				enumATOMname ename=name2enum(rold.atoms[j].name);
				if(ename<enumAtomTotNumconst-1) {
					strcpy(r.atoms[ename].name,rold.atoms[j].name);
					r.atoms[ename].r=rold.atoms[j].r;
					r.atoms[ename].chain=rold.atoms[j].chain;
					r.chain=rold.atoms[j].chain;
				}
			}
		}
		rlist.push_back(r);
	}

	for(int i=0;i<rlist.size();++i){
		if(strcmp(rlist[i].name,"  A")==0 || strcmp(rlist[i].name,"  G")==0)
			rlist[i].isPurine=true;
		else if(strcmp(rlist[i].name,"  U")==0 || strcmp(rlist[i].name,"  C")==0)
			rlist[i].isPurine=false;
		else{
			cout<<"Wrong name of residue "<<rlist[i].ridx<<" in PDB "<<rlist[i].pn<<endl;
			exit(0);
		}
	}

	//check if residue index is continuous
	
	for(int i=0;i<rlist.size()-1;++i){
		RESIDUE& r1=rlist[i];
		RESIDUE& r2=rlist[i+1];
		if(r1.chain==r2.chain){
			if(r2.ridx-r1.ridx!=1){
				//cout<<"residue "<<r1.ridx<<" and residue "<<r2.ridx<<" in chain "<<r1.chain<<" are not continuous in PDB "<<r1.pn<<endl;
				//exit(0);
			}
		}
	}


	return 0;
}

int	PDB::readpdbresidues(const char *pdbname,vector<RESIDUEpdb>& residues) {
	/*
	   read atoms in the pdb, and put atoms in the same residue together.
	   do not care if these residues in the same chain.
	   */
	residues.clear();

	vector<ATOMpdb>	atoms;
	readpdbatoms(pdbname,atoms);

	int	natom=atoms.size();
	if(natom==0) {
		cout<<"Number of atoms is zero."<<endl;
		return 0;
	}

	int	ip=0;
	RESIDUEpdb	res;
	for(int i=1; i<natom; ++i) {
        
		if(atoms[i].ridx!=atoms[ip].ridx || atoms[i].chain != atoms[ip].chain) {
			strcpy(res.name,atoms[ip].rname);
			res.atoms.clear();
			for(int j=ip; j<i; ++j)
				res.atoms.push_back(atoms[j]);
			residues.push_back(res);
			ip=i;
		}
	}
	//last nucleotide
	{
		strcpy(res.name,atoms[ip].rname);
		res.atoms.clear();
		for(int j=ip; j<natom; ++j)
			res.atoms.push_back(atoms[j]);
		residues.push_back(res);
	}

	return 0;
}


int	PDB::readpdbatoms(const char* pdbname,vector<ATOMpdb>& atoms) {
	/* Read all atoms in the pdb and store them in vector<> atoms.
	   Do not check if these atoms are in the same chain or not.
	   */
	atoms.clear();

	ifstream	in(pdbname);
	if(!in) {
		cout<<"Error open file!!! "<<pdbname<<endl;
		exit(0);
	}
	int	tag;
	ATOMpdb	a;
	while( (tag=getatom(in,a))!=-1) {
		//cout<<"reading "<<a.rname<<" "<<a.ridx<<" "<<tag<<" || "<<a<<endl;
		if(tag==1) atoms.push_back(a);
		if(tag==3) break;
	}
	in.close();

	return 0;
}


void	PDB::stripNonRNAchains(vector<RESIDUE>& residues,vector<int>& chain) {
	/*
	   if a chain contains a single non-RNA residues, this chain will be deleted
	   only pure RNA chains left.
	   */
	vector<RESIDUE>	residuesnew;
	vector<int>	chainnew;

	chainnew.push_back(0);
	int	nchain=chain.size()-1;
	for(int i=0; i<nchain; ++i) {
		int j1=chain[i],j2=chain[i+1],j,nc=0;
		for(j=j1; j<j2; ++j) if( residueType(residues[j])!=2 ) nc++;
		if( nc!=0) { //too many non-RNA residues, skip this chain
			//cout<<"percentage of non-RNA residues "<<(double)nc/(j2-j1)*100<<endl;
			continue;
		} else {
			for(int j=j1; j<j2; ++j) residuesnew.push_back(residues[j]);
			chainnew.push_back(residuesnew.size());
		}
	}
	residues=residuesnew;
	chain=chainnew;
}


void	PDB::checkIntegrity(vector<RESIDUE>& residues,vector<int>& chain) {
	/*
	   The input vector<> residues has to be rearranged first by calling rearrangeatoms().
	   If a residue is deemed to be incomplete, its name residues[i].name will be set to a NULL string.
	   */
	cout<<"Checking Integrity: "<<endl;
	if(residues.size()==0) return;

	int	nc=0;
	int	nres=residues.size();
	for(int j=0; j<nres; ++j) {
		ATOM*	list=residues[j].atoms;
		int	tag=0;
		if(strcmp(residues[j].name,"  A")==0) {
			for(int k=0; k<AnucleotideAtomNum; ++k) {
				if( list[ AnucleotideAtoms[k] ].name[0] == '\0' ) {
					cout<<"Residue A "<<j<<" incomplete, missing atoms: "<<enum2name(AnucleotideAtoms[k])<<endl;
					tag=1;
				}
			}
		} else if (strcmp(residues[j].name,"  G")==0) {
			for(int k=0; k<GnucleotideAtomNum; ++k) {
				if( list[ GnucleotideAtoms[k] ].name[0] == '\0' ) {
					cout<<"Residue G "<<j<<" incomplete, missing atoms: "<<enum2name(GnucleotideAtoms[k])<<endl;
					tag=1;
				}
			}
		} else if (strcmp(residues[j].name,"  U")==0) {
			for(int k=0; k<UnucleotideAtomNum; ++k) {
				if( list[ UnucleotideAtoms[k] ].name[0] == '\0' ) {
					cout<<"Residue U "<<j<<" incomplete, missing atoms: "<<enum2name(UnucleotideAtoms[k])<<endl;
					tag=1;
				}
			}
		} else if (strcmp(residues[j].name,"  C")==0) {
			for(int k=0; k<CnucleotideAtomNum; ++k) {
				if( list[ CnucleotideAtoms[k] ].name[0] == '\0' ) {
					cout<<"Residue C "<<j<<" incomplete, missing atoms: "<<enum2name(CnucleotideAtoms[k])<<endl;
					tag=1;
				}
			}
		} else {
			cout<<"Residue "<<residues[j].name<<" "<<j<<", Unknow type"<<endl;
			tag=2;
		}
		if(tag==1) {
			//strcpy(residues[j].name,"");
			residues[j].name[0]='z'; //tag it with a 'z'
			residues[j].name[1]='z';
			nc++;
		} else if(tag==2) {
			strcpy(residues[j].name,"XXX");
			nc++;
		}
	}
	printf("Incomplete residues: %5.1f\n\n",(float)nc/residues.size()*100);
}

const char* RNAname1[4]= {"  A","  G","  U","  C"};

const int noncan=88; //number of noncanonical nucleotide
const char* RNAname2[noncan]= {" RA"," RG"," RU"," RC",
	"2MG","H2U","M2G","OMG"," YG","PSU","5MC","7MG","5MU","1MA",
	"OMC","  I","1MG","GDP","A23","5BU","5IC","CB2","GTP","DHU",
	"AET","G7M","4SU","CCC","S4U","  T","FHU","AVC","OMU","UR3",
	"T6A","RIA","PGP","BRU","U34","YYG","CBR","A2M","BGM","UMS",
	"CSL"," IU","UD5","S4C","FMU","5FU"," DU","XUG","TM2","G46",
	"1SC","CFL","UFT","SUR","MTU","6FC"," CH","U8U","RUS"," IG",
	" IC","6MZ","CM0","MIA"," 0C"," 0U"," 0G"," DG","AP7","LCA",
    "10C","SSU","CBV","RA5","RG5","RU5","RC5","RA3","RG3","RU3",
    "RC3","PPU","6IA","  N"
};

const char* RNAname3[noncan]= {"  A","  G","  U","  C",
	"  G","  U","  G","  G","  G","  U","  C","  G","  U","  A",
	"  C","  G","  G","  G","  A","  U","  C","  C","  G","  U",
	"  A","  G","  U","  C","  U","  U","  U","  A","  U","  U",
	"  A","  A","  G","  U","  U","  G","  C","  A","  G","  U",
	"  C","  U","  U","  C","  U","  U","  U","  G","  U","  G",
	"  C","  C","  U","  U","  G","  C","  C","  U","  U","  G",
	"  C","  A","  U","  A","  C","  U","  G","  G","  A","  A",
    "  C","  U","  C","  A","  G","  U","  C","  A","  G","  U",
    "  C","  A","  A","  C"
};


int PDB::getatom(ifstream &in,ATOMpdb &a) {
	/*return:
	  -1:	end of file
0	:	line skipped
1	:	RNA ATOM read
2	:	TER
3   :       ENDMDL
*/
	char	line[200],str[100];
	int	tag=-1;
	if(in.getline(line,100)) {
		if( strncmp(line,"ATOM",4)==0 || strncmp(line,"HETATM",6)==0 ) {
			strmncpy(str,line,17,19);
			int	nlen=strlen(str);
			for(int i=0; i<nlen; ++i) if('a'<=str[i]&&str[i]<='z') str[i]+='A'-'a'; //capitalization

			bool ifRNA=false;
			for(int i=0; i<4; ++i) {
				if(strcmp(str,RNAname1[i])==0) {
					ifRNA=true;
					strcpy(a.rname,RNAname1[i]);
					break;
				}
			}

			if(ifRNA==false) {
				for(int i=0; i<noncan; ++i) {
					if(strcmp(str,RNAname2[i])==0) {
						ifRNA=true;
						strcpy(a.rname,RNAname3[i]);
						break;
					}
				}
			}

			if(ifRNA==false)
				return 0;

			strmncpy(str,line,6,10);
			a.idx=atoi(str);
			strmncpy(a.name,line,12,15);
			if(a.name[3]=='\'') a.name[3]='*';

			if(strcmp(a.name," O1P")==0) strcpy(a.name," OP1");
			if(strcmp(a.name," O2P")==0) strcpy(a.name," OP2");

			a.chain=line[21];

			strmncpy(str,line,22,25);
			a.ridx=atoi(str);	//residue index
			strmncpy(str,line,30,37);
			a.r.x=(double)atof(str); //x coordinate
			strmncpy(str,line,38,45);
			a.r.y=(double)atof(str); //y coordinate
			strmncpy(str,line,46,53);
			a.r.z=(double)atof(str); //z coordinate
			tag=1;
		} else if(strncmp(line,"TER ",4)==0) {
			tag=2;
		} else if(strncmp(line,"ENDMDL",6)==0) {
			tag=3;
		} else
			tag=0;
	}
	return tag;
}


void	PDB::strmncpy(char* str,char* line,int m, int n) {
	int	i;
	for(i=0; i<n-m+1; ++i) str[i]=line[i+m];
	str[i]='\0';
}

void	PDB::outputlevel0(vector<RESIDUE>& residues,vector<int>& chains) {
	cout<<"Number of chains: "<<chains.size()-1<<endl;
	cout<<"Number of residues: "<<residues.size()<<endl;
	for(int i=0; i<chains.size()-1; ++i) {
		cout<<"Chain: "<<i<<", number of residues: "<<chains[i+1]-chains[i]<<endl;
		for(int ic=1,j=chains[i]; j<chains[i+1]; ++j) {
			cout<<residues[j].name<<" ";
			if(ic%20==0) cout<<endl;
			ic++;
		}
		cout<<endl;
	}
}


void	PDB::outputlevel1(vector<RESIDUE>& residues,vector<int>& chains) {
	cout<<"Number of chains: "<<chains.size()-1<<endl;
	cout<<"Number of residues: "<<residues.size()<<endl;
	for(int i=0; i<chains.size()-1; ++i) {
		cout<<"--------------------> Chain: "<<i<<", number of residues: "<<chains[i+1]-chains[i]<<endl;
		for(int j=chains[i]; j<chains[i+1]; ++j) {
			cout<<"    "<<setw(4)<<j<<" Residue: "<<residues[j].name<<endl;
			for(int k=0; k<enumAtomTotNumconst; ++k)
				cout<<"        "<<setw(4)<<k<<" "<<residues[j].atoms[k].name<<endl;
		}
	}
}


enumATOMname	name2enum(char* name) {
	enumATOMname	ename;
	//write the code explicitly, for the sake of safety, though a little bit ugly.
	if( strcmp(name," P  ")==0 )	ename=P;
	else if ( strcmp(name," OP1")==0 ) ename=OP1;
	else if ( strcmp(name," OP2")==0 ) ename=OP2;
	else if ( strcmp(name," O5*")==0 ) ename=O5s;
	else if ( strcmp(name," C5*")==0 ) ename=C5s;
	else if ( strcmp(name," C4*")==0 ) ename=C4s;
	else if ( strcmp(name," C3*")==0 ) ename=C3s;
	else if ( strcmp(name," O3*")==0 ) ename=O3s;
	else if ( strcmp(name," O4*")==0 ) ename=O4s;
	else if ( strcmp(name," C1*")==0 ) ename=C1s;
	else if ( strcmp(name," C2*")==0 ) ename=C2s;
	else if ( strcmp(name," O2*")==0 ) ename=O2s;
	else if ( strcmp(name," N9 ")==0 ) ename=N9;
	else if ( strcmp(name," C8 ")==0 ) ename=C8;
	else if ( strcmp(name," N7 ")==0 ) ename=N7;
	else if ( strcmp(name," C5 ")==0 ) ename=C5;
	else if ( strcmp(name," C6 ")==0 ) ename=C6;
	else if ( strcmp(name," O6 ")==0 ) ename=O6;
	else if ( strcmp(name," N1 ")==0 ) ename=N1;
	else if ( strcmp(name," C2 ")==0 ) ename=C2;
	else if ( strcmp(name," N2 ")==0 ) ename=N2;
	else if ( strcmp(name," N3 ")==0 ) ename=N3;
	else if ( strcmp(name," C4 ")==0 ) ename=C4;
	else if ( strcmp(name," N6 ")==0 ) ename=N6;
	else if ( strcmp(name," O4 ")==0 ) ename=O4;
	else if ( strcmp(name," O2 ")==0 ) ename=O2;
	else if ( strcmp(name," N4 ")==0 ) ename=N4;
	else {
		ename=enumAtomEND;
	}
	return ename;
}


char*	enum2name(enumATOMname	ename) {
	static	char	name[5];
	//write the code explicitly, for the sake of safety, though a little bit ugly.
	if( ename==P )	strcpy(name," P  ");
	else if ( ename==OP1 ) strcpy(name," OP1");
	else if ( ename==OP2 ) strcpy(name," OP2");
	else if ( ename==O5s ) strcpy(name," O5*");
	else if ( ename==C5s ) strcpy(name," C5*");
	else if ( ename==C4s ) strcpy(name," C4*");
	else if ( ename==C3s ) strcpy(name," C3*");
	else if ( ename==O3s ) strcpy(name," O3*");
	else if ( ename==O4s ) strcpy(name," O4*");
	else if ( ename==C1s ) strcpy(name," C1*");
	else if ( ename==C2s ) strcpy(name," C2*");
	else if ( ename==O2s ) strcpy(name," O2*");
	else if ( ename==N9 ) strcpy(name," N9 ");
	else if ( ename==C8 ) strcpy(name," C8 ");
	else if ( ename==N7 ) strcpy(name," N7 ");
	else if ( ename==C5 ) strcpy(name," C5 ");
	else if ( ename==C6 ) strcpy(name," C6 ");
	else if ( ename==O6 ) strcpy(name," O6 ");
	else if ( ename==N1 ) strcpy(name," N1 ");
	else if ( ename==C2 ) strcpy(name," C2 ");
	else if ( ename==N2 ) strcpy(name," N2 ");
	else if ( ename==N3 ) strcpy(name," N3 ");
	else if ( ename==C4 ) strcpy(name," C4 ");
	else if ( ename==N6 ) strcpy(name," N6 ");
	else if ( ename==O4 ) strcpy(name," O4 ");
	else if ( ename==O2 ) strcpy(name," O2 ");
	else if ( ename==N4 ) strcpy(name," N4 ");
	else {
		cout<<"Unknow atom enumType "<<ename<<endl;
		exit(0);
	}
	return name;
}


char*    enum2name(enumNucType nuc) {
	static	char	name[4];
	if(nuc==nucA) strcpy(name,"  A");
	else if (nuc==nucG) strcpy(name,"  G");
	else if (nuc==nucU) strcpy(name,"  U");
	else if (nuc==nucC) strcpy(name,"  C");
	else {
		cout<<"Unknow nuc Type "<<nuc<<endl;
		exit(0);
	}
	return name;
}


char    enum2char(enumNucType nuc) {
	char	name;
	if(nuc==nucA) name='A';
	else if (nuc==nucG) name='G';
	else if (nuc==nucU) name='U';
	else if (nuc==nucC) name='C';
	else {
		cout<<"Unknow nuc Type "<<nuc<<endl;
		exit(0);
	}
	return name;
}

char name2char(RESIDUE& res) {
	if(strcmp(res.name,"  A")==0)
		return 'A';
	else if(strcmp(res.name,"  G")==0)
		return 'G';
	else if(strcmp(res.name,"  U")==0)
		return 'U';
	else if(strcmp(res.name,"  C")==0)
		return 'C';
	else {
		cout<<"wrong residue name. "<<res.name<<endl;
		cout<<res.pn<<" "<<res.ridx<<endl;
		exit(0);
	}
}

int	residueType(RESIDUE& res) {
	/* return:
0: unknown
1: protein
2: RNA
*/
	static const char*	proteins[]= {"GLY","ALA","VAL","LEU","ILE","PHE","TYR","TRP","SER","THR",
		"CYS","MET","LYS","ARG","HIS","ASP","GLU","PRO","ASN","GLN"
	};
	static const char*	RNAs[]= {"  A","  U","  G","  C","RA ","RU ","RG ","RC ","PSA","PSU","PSG","PSC","H2U"};

	int	i;

	for(i=0; i<13; ++i) if(strcmp(res.name,RNAs[i])==0) break;
	if(i<13) return 2;

	for(i=0; i<20; ++i) if(strcmp(res.name,proteins[i])==0) break;
	if(i<20) return 1;

	//test the atom names inside this residue
	int	tag=0;
	for(int i=0; i<enumAtomTotNumconst; ++i) {
		char *name=res.atoms[i].name;
		if(strcmp(name," P  ")==0||strcmp(name," O5*")==0||strcmp(name," C5*")==0||
				strcmp(name," C4*")==0||strcmp(name," O3*")==0) tag++;
		if(tag==5) break;
	}
	if(tag==5) return 2;

	return 0;
}


int	residueType(RESIDUEpdb& res) {
	/* return:
0: unknown
1: protein
2: RNA
*/
	static const char*	proteins[]= {"GLY","ALA","VAL","LEU","ILE","PHE","TYR","TRP","SER","THR",
		"CYS","MET","LYS","ARG","HIS","ASP","GLU","PRO","ASN","GLN"
	};
	static const char*	RNAs[]= {"  A","  U","  G","  C","RA ","RU ","RG ","RC ","PSA","PSU","PSG","PSC","H2U","1MG","5MC","5MU"};

	int	i;

	for(i=0; i<16; ++i) if(strcmp(res.name,RNAs[i])==0) break;
	if(i<16) return 2;

	for(i=0; i<20; ++i) if(strcmp(res.name,proteins[i])==0) break;
	if(i<20) return 1;

	//test the atom names inside this residue
	/*int	tag=0,natom=res.atoms.size();
	  for(int i=0;i<natom;++i) {
	  char *name=res.atoms[i].name;
	  if(strcmp(name," P  ")==0||strcmp(name," O5*")==0||strcmp(name," C5*")==0||
	  strcmp(name," C4*")==0||strcmp(name," O3*")==0) tag++;
	  if(tag==5) break;
	  }
	  if(tag==5) return 2;*/

	return 0;
}


int	isRNAbackboneAtom(ATOM& a) {
	static const  char* backbone[]= {" P  "," O5*"," C5*"," C4*"," C4*"," O3*"};
	int	i;
	for(i=0; i<6; ++i) {
		if(strcmp(a.name,backbone[i])==0) break;
	}
	if(i<6) return 1;
	else return 0;
}

int	isPurine(RESIDUE& res) {
	//return 1 if purine, return 0 if pyrimidine, otherwise, return -1
	if(strcmp(res.name,"  A")==0||strcmp(res.name,"  G")==0) return 1;
	else if(strcmp(res.name,"  U")==0||strcmp(res.name,"  C")==0) return 0;
	else return -1;
}


void	convertNucname2Num(vector<RESIDUE>& residues,int* sequenceidx) {
	static int tagCheck=0;
	if(tagCheck==0) {
		if(nucA!=0 || nucG!=1 || nucU!=2 || nucC!=3) {
			cout<<"program in danger, unexpected default behaviors of the current compiler"<<endl;
			exit(0);
		}
		tagCheck=1;
	}
	int nres=residues.size();
	for(int i=0; i<nres; ++i) {
		if( strcmp(residues[i].name,"  A" )==0 ) sequenceidx[i]=0;
		else if( strcmp(residues[i].name,"  G" )==0 ) sequenceidx[i]=1;
		else if( strcmp(residues[i].name,"  U" )==0 ) sequenceidx[i]=2;
		else if( strcmp(residues[i].name,"  C" )==0 ) sequenceidx[i]=3;
		else {
			//cout<<"error sequence information: "<<i<<" "<<residues[i].name<<endl;
			sequenceidx[i]=-1;
		}
	}
}


void	getSequenceInfo(vector<RESIDUE>& list,enumNucType* seq,int nres) {
	int n=list.size();
	if(n>nres) n=nres;
	/*
	   cout<<"Warning: I am changing the name of residues..."<<endl;
	   for(int i=0;i<n;++i) {
	   cout<<i<<"("<<list[i].name<<")"<<endl;
	   list[i].name[0]=' ';
	   list[i].name[1]=' ';
	   }*/
	for(int i=0; i<n; ++i) {
		if( strcmp(list[i].name,"  A" )==0 ) seq[i]=nucA;
		else if( strcmp(list[i].name,"  G" )==0 ) seq[i]=nucG;
		else if( strcmp(list[i].name,"  U" )==0 ) seq[i]=nucU;
		else if( strcmp(list[i].name,"  C" )==0 ) seq[i]=nucC;
		else {
			cout<<"error: not an RNA nucleotide -->"<<i<<"("<<list[i].name<<")"<<endl;
			exit(0);
		}
	}
}


enumNucType	res2enum(RESIDUE& res) {
	enumNucType	ename;
	if( strcmp(res.name,"  A" )==0 ) ename=nucA;
	else if( strcmp(res.name,"  G" )==0 ) ename=nucG;
	else if( strcmp(res.name,"  U" )==0 ) ename=nucU;
	else if( strcmp(res.name,"  C" )==0 ) ename=nucC;
	else {
		cout<<"Error converting residue: "<<res<<endl;
		exit(0);
	}
	return ename;
}

void read_rna_from_pdb(const char* dataset,const char* pdbname,vector<RESIDUE>& rna) {
	PDB pdb;
	vector<int> ichains;
	char pn[300];
	strcpy(pn,dataset);
	strcat(pn,pdbname);

	pdb.readpdb(pn,rna,ichains);
}

void read_rna_from_pdb(const char* dataset,const char* pdbname,vector<RESIDUE>& rna,vector<int>& ichains) {
	PDB pdb;
	char pn[300];
	strcpy(pn,dataset);
	strcat(pn,pdbname);

	pdb.readpdb(pn,rna,ichains);
}

bool isSameSequence(vector<RESIDUE>& native,vector<RESIDUE>& decoy){
	if(native.size()!=decoy.size())
		return false;

	for(int i=0;i<native.size();++i){
		if(strcmp(native[i].name,decoy[i].name)!=0)
			return false;
	}

	return true;
}

ostream& operator << (ostream& output, const ATOMpdb& a) {
	output.setf(ios::fixed);
	output<<a.name<<" "<<a.chain<<" "<< setw(8)<<setprecision(3)<<a.r.x
		<<setw(8)<<setprecision(3)<<a.r.y<<setw(8)<<setprecision(3)<<a.r.z;

	return output;
}


ostream& operator << (ostream& output, const RESIDUE& res) {
	cout<<" Residue: "<<res.name<<endl;
	for(int k=0; k<enumAtomTotNumconst; ++k) {
		if(res.atoms[k].name[0]=='\0') continue;
		cout<<setw(2)<<k<<" "<<res.atoms[k].name<<" "<<res.atoms[k].r.x<<" "<<res.atoms[k].r.y<<" "<<res.atoms[k].r.z<<endl;
	}
	return output;
}
