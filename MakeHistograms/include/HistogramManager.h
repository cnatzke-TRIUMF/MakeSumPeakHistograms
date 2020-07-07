#ifndef HISTOGRAM_MANAGER_H
#define HISTOGRAM_MANAGER_H

#include "TH1.h"
#include "TH2.h"
#include "TGriffin.h"
#include "TChain.h"
#include "TVector3.h"

class HistogramManager
{
public:
void MakeHistogramFile(TChain *inputChain, const char * linearParamFile);
void InitializeHistograms(int verbose = 0);
void FillHistograms(TChain *gChain);

private:
void PreProcessData();
void WriteHistogramsToFile();
bool EnergyGate(float gate, float hit_1_energy, float hit_2_energy, float threshold);

TGriffin *fGrif = NULL;
TGriffinBgo *fGriffinBgo = NULL;

std::vector<double> fAngleCombinations = {15.442, 21.9054, 29.1432, 33.1433, 38.382, 44.57, 47.4453, 48.7411, 51.4734, 55.1704, 59.9782, 60.1024, 62.3396, 62.4924, 63.4231, 68.9567, 71.4314, 73.3582, 73.6291, 75.7736, 80.9423, 81.5464, 83.8936, 86.868, 88.9658, 91.0342, 93.132, 96.1064, 98.4536, 99.0577, 104.226, 106.371, 106.642, 108.569, 111.043, 116.577, 117.508, 117.66, 119.898, 120.022, 124.83, 128.527, 131.259, 132.555, 135.43, 141.618, 146.857, 150.857, 158.095, 164.558, 180.0};

int num_crystals = 64;

double offsets[64];
double gains[64];


std::vector<TH1F*> hist_vec_1D;
std::vector<TH2F*> hist_vec_2D;

std::vector<float> energy_vec; // vector which contains the energy values
std::vector<long> time_vec; // vector which contains the time values
std::vector<TVector3> pos_vec; // vector which contains the position values
//std::vector<int> detector_vec; // vector which contains the detector values

int checkMix, lgsize, event_mixing_depth = 11;
std::vector<std::vector<float>> lastgrifEnergy;
std::vector<std::vector<TVector3>> lastgrifPosition;

};

#endif
