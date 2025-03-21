#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>
#include<ctime>
#include "pdb.h"
#include "pdbWrite.h"
#include "calc_crmsd_RNA.h"
#include "reshapeResidue.h"

using namespace std;

int main(int argc, char** argv) {
	if(argc!=5){
		cout<<"Error in input. Please input as: ./ExtractList pdblist datapath outputpath t_or_v_for_Train_or_Val"<<endl;
		//for example: /home/ExtractList /home/pdblist /home/data/ /home/datalist/ t
		//it means to calculate all nucleotides in pdbs in pdblist, and extract tdatalist or vdatalist into datalist folder
		return 1;
	}

	char* pdblist = argv[1];
	char* dpath = argv[2];
	char* opath = argv[3];
	char* tvtype = argv[4];
	char tvt[2] = "t";
	char tvv[2] = "v";

	double train_prob[30] = { 1,0.8343,0.14964,0.091922,0.092957,0.11447,0.14888,0.1913,0.23598,0.28425,0.33259,0.38958,0.44591,0.50788,0.57425,0.64178,0.71234,0.7882,0.84625,0.92125,0.98804,1,1,1,1,1,1,1,1,1 };
	double val_prob[30] = { 1,0.64186,0.12987,0.076907,0.072622,0.087519,0.11786,0.15245,0.18515,0.22447,0.26815,0.31805,0.37072,0.41574,0.4674,0.52372,0.58824,0.65021,0.7088,0.78359,0.8385,0.90323,1,1,1,1,1,1,1,1 };
	double tv_prob[30] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	if (strcmp(tvt, tvtype) == 0) {
		for (int tvn = 0; tvn < 30; tvn++) {
			tv_prob[tvn] = train_prob[tvn];
		}
	}
	else if (strcmp(tvv, tvtype) == 0) {
		for (int tvn = 0; tvn < 30; tvn++) {
			tv_prob[tvn] = val_prob[tvn];
		}
	}
	else {
		cout << "Error in setting Train_or_Val" << endl;
		return 2;
	}

	ifstream fin(pdblist);
	if (!fin) {
		cout << "Error in opening " << pdblist << endl;
		return 3;
	}
	
	vector<string> pdbname;
	string word;
	while (fin >> word) {
		pdbname.push_back(word);
	}
	fin.close();

	char listname[100] = "";
	sprintf(listname, "%s%sdatalist", opath, tvtype);
	ofstream outfile(listname);

	long int seed = time(NULL);
	srand48(seed);
	double prob;
	double pixels[3][nBins][nBins][nBins];
	double score;
	int count = 1;
	
	char nativeDB[100] = "";
	sprintf(nativeDB, "%sFilteredDataset/", dpath);

	for (int npdb = 0; npdb < pdbname.size(); npdb++) {
		char nativePDB[20];
		sprintf(nativePDB, "%s.pdb", pdbname[npdb].c_str());

		for (int ndecoy = 0; ndecoy < 300; ndecoy++) {
			char decoyDB[100] = "";
			sprintf(decoyDB, "%sMDPDB/%s/SelectedPDB/", dpath, pdbname[npdb].c_str());
			char decoyPDB[20];
			sprintf(decoyPDB, "%d.pdb", ndecoy);

			vector<RESIDUE> native, decoy;
			vector<int> ichains;
			read_rna_from_pdb(nativeDB, nativePDB, native, ichains);
			read_rna_from_pdb(decoyDB, decoyPDB, decoy);
			if (!isSameSequence(native, decoy)) {
				cout << "not the same sequence" << endl;
				cout << nativeDB << nativePDB << " " << decoyDB << decoyPDB << endl;
				return 4;
			}

			writeMassChargeInfo(decoy);

			for (int res_index = 0; res_index < decoy.size(); res_index++) {
				pixelateResidue(native, decoy, res_index, pixels, score);
				if (score >= 30.0)
					continue;

				int index = (int)(score / 1.0);
				double tp = tv_prob[index];
				prob = drand48();
				if (prob > tp)
					continue;

				outfile.setf(ios::fixed);
				outfile.precision(3);
				outfile << setw(10) << count << setw(10) << pdbname[npdb] << setw(10) << ndecoy << setw(10) << res_index << setw(10) << score << endl;

				count++;
			}
		}
	}
	outfile.close();
	return 0;
}
