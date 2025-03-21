#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>
#include<ctime>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
#include<cassert>

#include "pdb.h"
#include "pdbWrite.h"
#include "calc_crmsd_RNA.h"
#include "reshapeResidue.h"

using namespace std;

int main(int argc, char** argv) {

	//When testing, you can uncomment error report lines to turn on error feedback. (Meanwhile you should uncomment corresponding lines in python program.)

	/*
	if (argc != 5) {
		cout << "Error in input. Please input as :./voxelate.out data_path datalist start_number batch_size" << endl;
		//for example: /home/Voxelizate /home/data/ /home/datalist/tdatalist 129 128
		//it means voxelizate No.129 to No.256 nucleotides in tdatalist, pdb in 
		return 1;
	}
	*/

	char* dpath = argv[1];
	char* datalist = argv[2];
	int ndata = atoi(argv[3]);
	int batchsize = atoi(argv[4]);
	ifstream fin(datalist);

	/*
	if (!fin) {
		cout << "Error in opening " << datalist << endl;
		return 2;
	}
	*/

	string word;
	for (int ndt = 0; ndt < ndata; ndt++) {
		getline(fin, word);
	}

	char nativePDB[20] = "";
	char decoyPDB[20] = "";
	char resindexstr[20] = "";
	char tempw[20] = "";
	char nativeDB[100] = "";
	char decoyDB[100] = "";
	sprintf(nativeDB, "%sFilteredDataset/", dpath);

	char tword[60];
	char spc = ' ';
	char blank[2] = "";

	double pixels[3][nBins][nBins][nBins];
	double score;
	int res_index;

	for (int bn = 0; bn < batchsize; bn++) {
		strcpy(nativePDB, blank);
		strcpy(decoyPDB, blank);
		strcpy(decoyDB, blank);
		strcpy(resindexstr, blank);
		strcpy(tword, word.c_str());
		res_index = atoi(resindexstr);

		for (int nt = 0; nt < 10; nt++) {
			if (tword[nt + 10] != spc) {
				sprintf(tempw, "%c", tword[nt + 10]);
				strcat(nativePDB, tempw);
			}
			if (tword[nt + 20] != spc) {
				sprintf(tempw, "%c", tword[nt + 20]);
				strcat(decoyPDB, tempw);
			}
			if (tword[nt + 30] != spc) {
				sprintf(tempw, "%c", tword[nt + 30]);
				strcat(resindexstr, tempw);
			}
		}//This is a bad way to split the string by space. Try strtok() to rewrite it for higher efficiency.
		//It seems sprintf() may waste a lot of time. Use it as little as possible.

		sprintf(decoyDB, "%sMDPDB/%s/SelectedPDB/", dpath, nativePDB);
		strcat(decoyPDB, ".pdb");
		strcat(nativePDB, ".pdb");

		vector<RESIDUE> native, decoy;
		vector<int> ichains;

		read_rna_from_pdb(nativeDB, nativePDB, native, ichains);
		read_rna_from_pdb(decoyDB, decoyPDB, decoy);

		/*
		if (!isSameSequence(native, decoy)) {
			cout << "not the same sequence" << endl;
			cout << nativeDB << nativePDB << " " << decoyDB << decoyPDB << endl;
			return 3;
		}
		*/

		//calculate features
		writeMassChargeInfo(decoy);
		pixelateResidue(native, decoy, res_index, pixels, score);

		cout.setf(ios::fixed);
		cout.precision(3);

		for (int m = 0; m < 3; m++) {
			for (int i = 0; i < nBins; i++) {
				for (int j = 0; j < nBins; j++) {
					for (int k = 0; k < nBins; k++) {
						if (fabs(pixels[m][i][j][k]) >= 0.001) {
							cout << setw(4) << bn << setw(3) << m << setw(3) << i << setw(3) << j << setw(3) << k << setw(10) << pixels[m][i][j][k] << endl;
						}
					}
				}
			}
		}

		getline(fin, word);
	}

	fin.close();
	return 0;
}

