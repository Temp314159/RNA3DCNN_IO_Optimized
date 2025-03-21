#ifndef	__PDBALLATOMWRITE__

#define	__PDBALLATOMWRITE__

#include <fstream>
#include <vector>
#include "pdb.h"

using namespace std;


class PDBW {
	private:
		int	iresinit;
		int	iatom,ires,imodel;
		ofstream out;

	public:
		explicit	PDBW(char* filename);
		PDBW(char* filename,int ir,int im);
		~PDBW();
		void	reset();
		void 	writeatom(ATOM& a,char* rname,char chainname,int ridx);
		void 	writeresidue(RESIDUE& res,char chainname,int ridx);
		void	writeresidue(RESIDUE& res);
		void	writeRNA(vector<RESIDUE>& rna);
		void 	writeRNA(vector<RESIDUE>& rna,vector<char>& chainname,vector<int>& ridx);
		void	outModel();
		void	outTerm();
};

#endif

