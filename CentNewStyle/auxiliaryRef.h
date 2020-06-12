#ifndef AUXILIARYREF
#define AUXILIARYREF


// fs/fu
// PDG
double BsFrag = 10.3;
double BsErr  = 0.5;
double BPFrag = 40.5;
double BPErr  = 0.6;
const int BandBin = 1;
//double BandY[BandBin] = {BsFrag/BPFrag};
//double BandY[BandBin] = {0.256};

double BandY[BandBin] = {0.244};

//double bin_min = 5;
//double bin_max = 50;

double bin_min = 4;
double bin_max = 25;




double BandX[BandBin] = {(bin_max+bin_min)/2};
double BandXErr[BandBin] = {(bin_max-bin_min)/2};


//forCent
double binCent_min = 0;
double binCent_max = 450;
//double binCent_min = 1.5;
//double binCent_max = 2.5;
double BandXCent[BandBin] = {(binCent_max+binCent_min)/2};
double BandXErrCent[BandBin] = {(binCent_max-binCent_min)/2};

//double BandYErr[BandBin] = { BsFrag/BPFrag*TMath::Sqrt((BsErr/BsFrag)*(BsErr/BsFrag)+(BPErr/BPFrag)*(BPErr/BPFrag)) };
//double BandYErr[BandBin] = { 0.020 };

double BandYErr[BandBin] = { 0.012 };

#endif
