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
#include "globals.h"

std::vector<float> energy_gates_vec = {700, 1022, 1460, 1760, 2614};

/************************************************************//**
 * Creates and Fills histograms
 *
 * @param verbose Verbosity level
 ***************************************************************/
void HistogramManager::MakeHistogramFile(TChain *inputChain, std::string linearParamFile)
{
	int verbose = 0;

	if (verbose > 0) {std::cout << "Reading secondary calibration parameters ... " << std::endl;}

	std::ifstream lin_coeff_file(linearParamFile);
	for (int i = 0; i < num_crystals; i++) {
		lin_coeff_file >> gains[i];
		lin_coeff_file >> offsets[i];
	}
	lin_coeff_file.close();
	if (verbose > 0) {std::cout << "Reading secondary calibration parameters ... [DONE]" << std::endl;}

	if (verbose > 0) {
		for (unsigned int i = 0; i < (sizeof(gains)/sizeof(*gains)); i++) {
			std::cout << "Channel " << i << " Gain " << gains[i] << " Offset: "<< offsets[i] << std::endl;
		}
	}

	if (verbose > 0) {std::cout << "Generating Post-Calibrated Histograms ..." << std::endl;}

	InitializeHistograms();
	FillHistograms(inputChain);
	WriteHistogramsToFile();

	if (verbose > 0) {std::cout << "Generating Post-Calibrated Histograms ... [DONE]" << std::endl;}
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
	for (auto const &energy_gate : energy_gates_vec) {
    	hist_vec_2D.push_back(new TH2D(Form("energy_angle_%i", (int)energy_gate), Form("#gamma_1 Energy Angle %i keV;Angle; #gamma Energy (keV)", (int)energy_gate), angle_bins, 0, angle_bins, g_bins, g_min, g_max));
    	hist_vec_2D_mixed.push_back(new TH2D(Form("energy_angle_%i_mixed", (int)energy_gate), Form("#gamma_1 Energy Angle %i keV Mixed;Angle; #gamma Energy (keV)", (int)energy_gate), angle_bins, 0, angle_bins, g_bins, g_min, g_max));
	}
	//hist_vec_2D.push_back(new TH2D("AB_AngDiff_Supp_1022", "#gamma_1 Energy Angle;Angle; #gamma Energy (keV)", angle_bins, 0, angle_bins, g_bins, g_min, g_max));
	//hist_vec_2D.push_back(new TH2D("AB_AngDiff_Supp_1022_mixed", "#gamma_1 Energy Angle (Mixed);Angle; #gamma Energy (keV)", angle_bins, 0, angle_bins, g_bins, g_min, g_max));

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
	float gate_threshold = 3; // bounds of energy gate (-+ 3)

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

					int angleIndex = GetAngleIndex(angle, angle_combinations_vec);
					double ggTime = TMath::Abs(time_vec.at(g1) - time_vec.at(g2));

					// check for bad angles
					if (angleIndex == -1) {
						std::cout << "Bad Angle" << std::endl;
						continue;
					}

					// Prompt coincidences
					if (ggTime < ggPrompt) {
                        // 1D
//                        for (unsigned int gate_iter = 0; gate_iter < energy_gates_vec.size(); gate_iter++){
//    						if (EnergyGate(energy_gates_vec.at(gate_iter), energy_vec.at(g1), energy_vec.at(g2), gate_threshold)) {
//        						hist_vec_1D.at(gate_iter)->Fill(energy_vec.at(g1) + energy_vec.at(g2));
//                            }
//                        }
						hist_vec_1D.back()->Fill(energy_vec.at(g1) + energy_vec.at(g2));

                        // 2D
                        for (unsigned int gate_iter = 0; gate_iter < energy_gates_vec.size(); gate_iter++){
    						if (EnergyGate(energy_gates_vec.at(gate_iter), energy_vec.at(g1), energy_vec.at(g2), gate_threshold)) {
        						hist_vec_2D.at(gate_iter)->Fill(angleIndex, energy_vec.at(g1));
                            }
                        }
						//if (EnergyGate(1022, energy_vec.at(g1), energy_vec.at(g2), gate_threshold)) {
							//hist_vec_2D.at(0)->Fill(angleIndex, energy_vec.at(g1));
							//hist_vec_2D.at(0)->Fill(angleIndex, energy_vec.at(g2));
						//}
					}
					// Background subtraction
					else if (bgLow < ggTime && ggTime < bgHigh) {
                        // 1D
						hist_vec_1D.back()->Fill(energy_vec.at(g1) + energy_vec.at(g2), -ggPrompt/(bgHigh-bgLow));

                        // 2D
                        for (unsigned int gate_iter = 0; gate_iter < energy_gates_vec.size(); gate_iter++){
    						if (EnergyGate(energy_gates_vec.at(gate_iter), energy_vec.at(g1), energy_vec.at(g2), gate_threshold)) {
        						hist_vec_2D.at(gate_iter)->Fill(angleIndex, energy_vec.at(g1), -ggPrompt/(bgHigh-bgLow));
                            }
                        } // energy_gates_vec
						//if (EnergyGate(1022, energy_vec.at(g1), energy_vec.at(g2), gate_threshold)) {
							//hist_vec_2D.at(0)->Fill(angleIndex, energy_vec.at(g1), -ggPrompt/(bgHigh-bgLow));
							//hist_vec_2D.at(0)->Fill(angleIndex, energy_vec.at(g2), -ggPrompt/(bgHigh-bgLow));
						//}
					}
				} // grif2

				// EVENT MIXED MATRICES
				// event mixing, we use the last event as second griffin
				checkMix = (int)last_grif_energy.size();
				if(checkMix < event_mixing_depth) continue;
				for(auto lg = 0; lg < (checkMix - 1); ++lg) {
					unsigned int multLG = last_grif_energy.at(lg).size();

					for(unsigned int g3 = 0; g3 < multLG; ++g3) {
						double angle = pos_vec.at(g1).Angle(last_grif_pos.at(lg).at(g3)) * 180. / TMath::Pi();
						if (angle < 0.0001) continue;
						int angleIndex = GetAngleIndex(angle, angle_combinations_vec);

						// Filling histogram
                        for (unsigned int gate_iter = 0; gate_iter < energy_gates_vec.size(); gate_iter++){
    						if (EnergyGate(energy_gates_vec.at(gate_iter), energy_vec.at(g1), last_grif_energy.at(lg).at(g3), gate_threshold)) {
        						hist_vec_2D_mixed.at(gate_iter)->Fill(angleIndex, energy_vec.at(g1));
                            }
                        } // energy_gates_vec
					} // end g3
				} // end LG
			} // grif1
		} // suppr_en.size() filter

		if (i % 10000 == 0) {
			progress_bar.display();
		}
		++progress_bar; // iterates progress_bar

		// update "last" event for event mixing
		last_grif_energy.push_back(energy_vec);
		last_grif_pos.push_back(pos_vec);
		lgsize = static_cast<int>(last_grif_energy.size());
		if (lgsize > event_mixing_depth) {
			last_grif_energy.erase(last_grif_energy.begin());
			last_grif_pos.erase(last_grif_pos.begin());
		}

		// Cleaning up for next event
		energy_vec.clear();
		pos_vec.clear();
		time_vec.clear();
		//detector_vec.clear();
	} // end TChain loop

	progress_bar.done();
} // FillHistograms

/************************************************************//**
 * Applies linear calibration to data points
 *
 ***************************************************************/
void HistogramManager::PreProcessData()
{
	float energy_temp = 0;
	int det_id = -1;
	int multiplicity_limit = 2;

	if (fGrif->GetSuppressedMultiplicity(fGriffinBgo) == multiplicity_limit) { // multiplicity filter
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

			energy_vec.push_back(energy_temp);
			pos_vec.push_back(fGrif->GetSuppressedHit(j)->GetPosition(145.0));
			time_vec.push_back(fGrif->GetSuppressedHit(j)->GetTime());
//		detector_vec.push_back(det);
		}
	} // multiplicity filter

} // PreProcessData

/************************************************************//**
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


/************************************************************//**
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
	for (auto const &h : hist_vec_2D_mixed) {
		h->Write();
	}
	out_file->Close();
	delete out_file;

} // WriteHistogramsToFile


/************************************************************//**
 * Returns the angular index
 *
 * @param angle The angle between two gammas
 * @param vec Vector of angles
 *****************************************************************************/
int HistogramManager::GetAngleIndex(double angle, std::vector<double> vec)
{

	// corner cases
	if (angle <= vec.front()) { return 0;}
	if (angle >= vec.back() - 1.) { return vec.size() - 1;}

	// binary search
	unsigned int i = 0, j = vec.size(), mid = 0;
	while ( i < j ) {
		mid = (i + j) / 2;

		if (vec[mid] == angle) return vec[mid];

		// searching left half
		if (angle < vec[mid]) {
			// if angle is greater than previous to mid, return closest of two
			if (mid > 0 && angle > vec[mid - 1]) {
				return GetClosest(mid - 1, mid, angle_combinations_vec, angle);
			}

			// repeat for left half
			j = mid;
		}
		// if angle is greater than mid
		else{
			if (mid < vec.size() - 1 && angle < vec[mid + 1]) {
				return GetClosest(mid, mid + 1, angle_combinations_vec, angle);
			}

			// update i
			i = mid + 1;
		}
	}
	// Only single element left after search
	return mid;
} // GetAngleIndex

/************************************************************//**
 * Returns the value closest to the target
 * Assumes val2 is greater than val1 and target lies inbetween the two
 *
 * @param val1 First value to compare
 * @param val2 Second value to compare
 * @param vec Vector of values
 * @param target Target value
 *****************************************************************************/
int HistogramManager::GetClosest(int val1, int val2, std::vector<double> vec, double target)
{
	if ((target - vec[val1]) >= (vec[val2] - target))
		return val2;
	else
		return val1;
} // GetClosest

/************************************************************//**
 * Displays humourous loading message
 *
 ***************************************************************/
void HistogramManager::DisplayLoadingMessage() {

	std::string line;
	int random = 0;
	int numOfLines = 0;
	std::ifstream File("loadingQuotes.txt");

	random = rand() % 131;

	while(getline(File, line))
	{
		++numOfLines;

		if(numOfLines == random)
		{
			std::cout << line << std::endl;
		}
	}

} // DisplayLoadingMessage
