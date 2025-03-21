#ifndef	__PDBALLATOM_H__

#define	__PDBALLATOM_H__

#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include "VECTOR.h"

using namespace std;

//------------------------------------------------------------------------------------------------------------------------
enum enumNucType {nucA,nucG,nucU,nucC,otherNucleotideType};

enum enumATOMname {P,OP1,OP2,O5s,C5s,C4s,C3s,O3s,O4s,C1s,C2s,O2s,N9,C8,N7,C5,C6,O6,N1,C2,N2,N3,C4,N6,O4,O2,N4,enumAtomEND};

extern enumATOMname p,op1,op2,o5s,c5s,c4s,c3s,o3s,o4s,c1s,c2s,o2s,n9,c8,n7,c5,c6,o6,n1,c2,n2,n3,c4,n6,o4,o2,n4;
//------------------------------------------------------------------------------------------------------------------------
#define	enumAtomTotNumconst	28

#define	AnucleotideAtomNum 22
#define	GnucleotideAtomNum 23
#define	UnucleotideAtomNum 20
#define	CnucleotideAtomNum 20
#define AbaseAtomNum 10
#define GbaseAtomNum 11
#define UbaseAtomNum 8
#define CbaseAtomNum 8

#define BkboneSugarAtomNum 12
#define	RingAtomNum	6
#define ImidazoleRingAtomNum 5

extern enumATOMname AbaseAtoms[AbaseAtomNum];
extern enumATOMname GbaseAtoms[GbaseAtomNum];
extern enumATOMname UbaseAtoms[UbaseAtomNum];
extern enumATOMname CbaseAtoms[CbaseAtomNum];

extern	enumATOMname AnucleotideAtoms[AnucleotideAtomNum], GnucleotideAtoms[GnucleotideAtomNum];
extern	enumATOMname UnucleotideAtoms[UnucleotideAtomNum], CnucleotideAtoms[CnucleotideAtomNum];
extern	enumATOMname BkboneSugarAtoms[BkboneSugarAtomNum];
extern	enumATOMname RingAtoms[RingAtomNum];

extern  enumATOMname ImidazoleRingAtoms[ImidazoleRingAtomNum];

//------------------------------------------------------------------------------------------------------------------------
struct ORIENTATION {
	VECTOR O,vx,vy,vz;
};

struct	ATOM {
	char name[5];
	char chain;
	VECTOR r;
	double mass;
	double charge;
};

struct	RESIDUE {
	char pn[200];
	int ridx;
	char chain;
	char name[4];
	ATOM atoms[enumAtomTotNumconst];
	double charge[enumAtomTotNumconst];
	int atoms_type[enumAtomTotNumconst];
	ORIENTATION ori;
	bool isPurine;
};

struct ATOMpdb {
	int	idx,ridx;
	char name[5];
	char rname[4],chain;
	VECTOR r;
};

struct RESIDUEpdb {
	char name[4];
	vector<ATOMpdb>	atoms;
};

class PDB {

	int	readpdbresidues(const char *pdbname,vector<RESIDUEpdb>& residues);
	int	readpdbatoms(const char* pdbname,vector<ATOMpdb>& atoms);
	int	getatom(ifstream &in,ATOMpdb &a);
	void strmncpy(char* str,char* line,int m, int n);

	public:
	PDB();
	~PDB() {};
	int	readpdb(const char *pdbname,vector<RESIDUE>& residues,vector<int>& chain);
	void modifyResname(vector<RESIDUE>& residues,vector<int>& chains);
	void stripNonRNAchains(vector<RESIDUE>& residues,vector<int>& chain);
	void outputlevel0(vector<RESIDUE>& residues,vector<int>& chains);
	void outputlevel1(vector<RESIDUE>& residues,vector<int>& chains);
	void extractRNABackbone(vector<RESIDUE>& residues,vector<int>& chain,vector<ATOMpdb>& atoms,vector<int>& atomschain);
	void checkIntegrity(vector<RESIDUE>& residues,vector<int>& chain);
};
//------------------------------------------------------------------------------------------------------------------------
int	incompleteResidue(RESIDUE& res);
int	isPurine(RESIDUE& res);
int	isRNAbackboneAtom(ATOM& a);

void getSequenceInfo(vector<RESIDUE>& list,enumNucType* seq,int nres);
void convertNucname2Num(vector<RESIDUE>& residues,int* sequenceidx);

char* enum2name(enumATOMname	ename);
char* enum2name(enumNucType nuc);
char enum2char(enumNucType nuc);
enumATOMname name2enum(char* name);
enumNucType res2enum(RESIDUE& res);
char name2char(RESIDUE& res);

void read_rna_from_pdb(const char* dataset,const char* pdbname,vector<RESIDUE>& rna);
void read_rna_from_pdb(const char* dataset,const char* pdbname,vector<RESIDUE>& rna,vector<int>& ichains);

bool isSameSequence(vector<RESIDUE>& native,vector<RESIDUE>& decoy);

ostream& operator << (ostream& output, const ATOMpdb& a);
ostream& operator << (ostream& output, const RESIDUE& r);

#endif

