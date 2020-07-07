#ifndef SUM_HISTOGRAM_MANAGER_H
#define SUM_HISTOGRAM_MANAGER_H

#include "TH1.h"
#include "TH2.h"
#include "TGriffin.h"
#include "TGriffinBgo.h"
#include "TChain.h"
#include "TVector3.h"

class HistogramManager
{
public:
    void MakeHistogramFile(TChain *inputChain, std::string linearParamFile);
    void InitializeHistograms(int verbose = 0);
    void FillHistograms(TChain *gChain);

private:
    void PreProcessData();
    void WriteHistogramsToFile();
    bool EnergyGate(float gate, float hit_1_energy, float hit_2_energy, float threshold);
    int GetAngleIndex(double angle, std::vector<double> vec);
    int GetClosest(int val_1, int val_2, std::vector<double> vec, double target);
    void DisplayLoadingMessage();

    TGriffin *fGrif = NULL;
    TGriffinBgo *fGriffinBgo = NULL;

    int num_crystals = 64;

    double offsets[64];
    double gains[64];

    std::vector<TH1D*> hist_vec_1D;
    std::vector<TH2D*> hist_vec_2D;

    std::vector<float> energy_vec; // vector which contains the energy values
    std::vector<long> time_vec; // vector which contains the time values
    std::vector<TVector3> pos_vec; // vector which contains the position values
    //std::vector<int> detector_vec; // vector which contains the detector values

    int checkMix, lgsize, event_mixing_depth = 11;
    std::vector<std::vector<float>> lastgrifEnergy;
    std::vector<std::vector<TVector3>> lastgrifPosition;
};

#endif
