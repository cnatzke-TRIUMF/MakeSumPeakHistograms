//////////////////////////////////////////////////////////////////////////////////
// Creates and fills histograms
//
// Author:          Connor Natzke (cnatzke@triumf.ca)
// Creation Date:   Friday July 3, 2020	T10:56:43-07:00
// Last Update:     Friday July 3, 2020	T10:56:43-07:00
// Usage:
//
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "HistogramManager.h"
#include "progress_bar.h"

/************************************************************//**
 * Creates and Fills histograms
 *
 * @param verbose Verbosity level
 ***************************************************************/
void HistogramManager::GenerateHistogramFile(TChain *inputChain, const char * linearParamFile)
{
	int verbose = 0;

	std::cout << "Reading secondary calibration parameters ... " << std::endl;

	std::ifstream lin_coeff_file(linearParamFile);
	for (int i = 0; i < num_crystals; i++) {
		lin_coeff_file >> gains[i];
		lin_coeff_file >> offsets[i];
	}
	lin_coeff_file.close();
	std::cout << "Reading secondary calibration parameters ... [DONE]" << std::endl;

	if (verbose > 0) {
		for (int i = 0; i < (sizeof(gains)/sizeof(*gains)); i++) {
			std::cout << "Channel " << i << " Gain " << gains[i] << " Offset: "<< offsets[i] << std::endl;
		}
	}

	std::cout << "Generating Post-Calibrated Histograms ..." << std::endl;

	InitializeHistograms();
	FillHistograms(inputChain);
	WriteHistogramsToFile();

	std::cout << "Generating Post-Calibrated Histograms ... [DONE]" << std::endl;
} // GenerateHistogramFile
/************************************************************//**
 * Initializes histograms to be filled
 *
 * @param verbose Verbosity level
 ***************************************************************/
void HistogramManager::InitializeHistograms(int verbose)
{
	int g_bins = 4000;
	int g_min = 0;
	int g_max = g_bins;
	int angle_bins = 52;

	// 1D Histograms
	hist_vec_1D.push_back(new TH1D("sum_energy", "Energy Sum", g_bins, g_min, g_max));
	//for (int i = 0; i < num_crystals; ++i) {
	//	hist_vec_1D.push_back(new TH1D(Form("pos_%i", i), Form("#gamma singles crystal %1.1i", i), g_bins, g_min, g_max));
	//}

	// 2D Histograms
	hist_vec_2D.push_back(new TH2D("AB_AngDiff_Supp_1022", "#gamma_1 Energy Angle;Angle", angle_bins, 0, angle_bins, g_bins, g_min, g_max));
	hist_vec_2D.push_back(new TH2D("AB_AngDiff_Supp_1022_mixed", "#gamma_1 Energy Angle;Angle", angle_bins, 0, angle_bins, g_bins, g_min, g_max));

	if (verbose > 0 ) std::cout << "Histograms created" << std::endl;

} // InitializeHistograms

/************************************************************//**
 * Fills histograms
 *
 * @param gChain Data chain
 ***************************************************************/
void HistogramManager::FillHistograms(TChain *gChain)
{

	float ggPrompt = 30.; // max time difference for gamma gamma; 30 ns
	float bgLow = 500.; // min time difference for gamma gamma time random
	float bgHigh = 2000.; // max time diff for gamma gamma time random

	if (gChain->FindBranch("TGriffin")) {
		gChain->SetBranchAddress("TGriffin", &fGrif);
		if (fGrif != NULL) {
			std::cout << "Succesfully found TGriffin branch" << std::endl;
		} else {
			std::cout << "Could not find TGriffin branch ... exiting" << std::endl;
		}
	}

	DisplayLoadingMessage();

	long analysis_entries = gChain->GetEntries();

	/* Creates a progress bar that has a width of 70,
	 * shows '=' to indicate completion, and blank
	 * space for incomplete
	 */
	ProgressBar progress_bar(analysis_entries, 70, '=', ' ');
	for (auto i = 0; i < analysis_entries; i++) {
		gChain->GetEntry(i);

		// Applies secondary energy calculation
		PreProcessData();

		// Filling histograms
		if (energy_vec.size() > 0) {
			for (unsigned int g1 = 0; g1 < energy_vec.size(); ++g1) {
				// gamma-gamma matrices
				for(unsigned int g2 = 0; g2 < energy_vec.size(); ++g2) {
					if (g1 == g2) continue;

					double angle = pos_vec.at(g1).Angle(pos_vec.at(g2)) * 180. / TMath::Pi();
					if (angle < 0.0001) continue;

					int angleIndex = GetAngleIndex(angle, fAngleCombinations);
					double ggTime = TMath::Abs(time_vec.at(g1) - time_vec.at(g2));

					// check for bad angles
					if (angleIndex == -1) {
						std::cout << "Bad Angle" << std::endl;
						continue;
					}

					// Prompt coincidences
					if (ggTime < ggHigh) {
						hist_vec_1D.at(0)->Fill(energy_vec.at(g1) + energy_vec.at(g2));
						if (EnergyGate(1022, energy_vec.at(g1), energy_vec.at(g2), threshold)) {
							hist_vec_2D.at(0)->Fill(angle, energy_vec.at(g1));
							hist_vec_2D.at(0)->Fill(angle, energy_vec.at(g2));
						}
					}
					// Background subtraction
					else if (bgLow < ggTime && ggTime < bgHigh) {
						hist_vec_1D.at(0)->Fill(energy_vec.at(g1) + energy_vec.at(g2), -ggHigh/(bgHigh-bgLow));
						if (EnergyGate(1022, energy_vec.at(g1), energy_vec.at(g2), threshold)) {
							hist_vec_2D.at(0)->Fill(angle, energy_vec.at(g1), -ggHigh/(bgHigh-bgLow));
							hist_vec_2D.at(0)->Fill(angle, energy_vec.at(g2), -ggHigh/(bgHigh-bgLow));
						}
					}
				} // grif2

				// EVENT MIXED MATRICES
				// event mixing, we use the last event as second griffin
				checkMix = (int)lastgrifEnergy.size();
				if(checkMix < event_mixing_depth) continue;
				for(auto lg = 0; lg < (checkMix - 1); ++lg) {
					unsigned int multLG = lastgrifEnergy.at(lg).size();

					for(unsigned int g3 = 0; g3 < multLG; ++g3) {
						double angle = pos_vec.at(g1).Angle(lastgrifPosition.at(lg).at(g3)) * 180. / TMath::Pi();
						if (angle < 0.0001) continue;
						int angleIndex = GetAngleIndex(angle, fAngleCombinations);

						// Filling histogram
						hist_vec_2D.at(1)->Fill(energy_vec.at(g1), lastgrifEnergy.at(lg).at(g3));
					} // end g3
				} // end LG
			} // grif1
		} // suppr_en.size() filter

		if (i % 10000 == 0) {
			progress_bar.display();
		}
		++progress_bar; // iterates progress_bar

		// update "last" event for event mixing
		lastgrifEnergy.push_back(energy_vec);
		lastgrifPosition.push_back(pos_vec);
		lgsize = static_cast<int>(lastgrifEnergy.size());
		if (lgsize > event_mixing_depth) {
			lastgrifEnergy.erase(lastgrifEnergy.begin());
			lastgrifPosition.erase(lastgrifPosition.begin());
		}

		// Cleaning up for next event
		energy_vec.clear();
		pos_vec.clear();
		time_vec.clear();
		//detector_vec.clear();
	} // end TChain loop

	progress_bar.done();
} // FillHistograms

/***************************************************************
* Applies linear calibration to data points
*
***************************************************************/
void HistogramManager::PreProcessData()
{
	float energy_temp = 0;
	int det_id = -1;
	int multiplicity_limit = 2;

	if (fGrif->GetSuppressedMultiplicity(fGriffinBgo) == 2) { // multiplicity filter
		for (auto j = 0; j < fGrif->GetSuppressedMultiplicity(fGriffinBgo); ++j) {
			det_id = fGrif->GetSuppressedHit(j)->GetArrayNumber();
			if (det_id == -1) {
				std::cout << "BAD DETECTOR" << std::endl;
				continue;
			}
			if(fGrif->GetSuppressedHit(j)->GetKValue()!=700) {continue;} // removes GRIFFIN hits pileup events

			// Applying secondary linear energy calibration
			energy_temp = offsets[det_id-1] + gains[det_id-1]*fGrif->GetSuppressedHit(j)->GetEnergy();
			//energyTmp += ((double) rand() / RAND_MAX - 0.5);

			// No secondary calibration
			//energyTmp = fGrif->GetSuppressedHit(j)->GetEnergy();

			energy_vec.push_back(energyTmp);
			pos_vec.push_back(fGrif->GetSuppressedHit(j)->GetPosition(145.0));
			time_vec.push_back(fGrif->GetSuppressedHit(j)->GetTime());
//		detector_vec.push_back(det);
		}
	} // multiplicity filter

} // PreProcessData

/***************************************************************
* Filter for energy gating
*
* @param gate Energy of interest (keV)
* @param hit_1_energy Energy of gamma 1 (keV)
* @param hit_2_energy Energy of gamma 2 (keV)
* @param threshold Accepted energy difference from gate (keV)
***************************************************************/
bool HistogramManager::EnergyGate(float gate, float hit_1_energy, float hit_2_energy, float threshold)
{
	double gate_upper = gate + threshold;
	double gate_lower = gate - threshold;
	return (hit_1_energy + hit_2_energy > gate_lower) && (hit_1_energy + hit_2_energy < gate_upper);
} // EnergyGate()

/***************************************************************
* Writes output ROOT file
*
***************************************************************/
void HistogramManager::WriteHistogramsToFile()
{
	TFile *out_file = new TFile("secondary_calibrated_histograms.root", "RECREATE");
	std::cout << "Writing output file: " << out_file->GetName() << std::endl;

	out_file->cd();
	for (auto const &h : hist_vec_1D) {
		h->Write();
	}
	for (auto const &h : hist_vec_2D) {
		h->Write();
	}
	out_file->Close();
	delete out_file;

} // WriteHistogramsToFile
