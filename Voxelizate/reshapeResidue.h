#ifndef __RESHAPERESIDUE_H__
#define __RESHAPERESIDUE_H__

#include <vector>
#include <iomanip>
#include "pdb.h"
#include "VECTOR.h"
#include "calc_crmsd_VECTOR.h"
#include "calc_score.h"

static const double BoxSideLength=32.0;
static const double HalfBoxSideLength=BoxSideLength/2.0;
static const double TrialHalfBoxSideLength=HalfBoxSideLength+5.0;
static const double BinWidth=1.0;
static const double HalfBinWidth=BinWidth/2.0;
static const int nBins=(int)32.0;

void calc_local_reference(RESIDUE& res,VECTOR& O,VECTOR& vx,VECTOR& vy,VECTOR& vz);

void writeMassChargeInfo(vector<RESIDUE>& rna);

void pixelateResidue(vector<RESIDUE>& native,vector<RESIDUE>& rna,int res_index,double pixels[3][nBins][nBins][nBins],double& score);

#endif
