//////////////////////////////////////////////////////////////////////////////////
// Extracts the sum peak hits and creates angular histograms
//
// Author:        Connor Natzke (cnatzke@triumf.ca)
// Creation Date: Friday July 3, 2020	T10:19:56-07:00
// Last Update:   Friday July 3, 2020	T10:19:56-07:00
// Usage:
//
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TFile.h" // needed for GetRunNumber
#include "TGRSIUtilities.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "HistogramManager.h"
#include "SumPeakHistograms.h"
#include "Notifier.h"

Notifier *notifier = new Notifier;
/************************************************************//**
 * Lists all of the files in a directory matching the run number
 * (i.e. extracts all subruns for processing)
 *
 * @param dir_name   Directory name
 * @param ext        Extension of the files (default = ".root")
 ***************************************************************/
int main(int argc, char **argv) {

	if (argc == 1) { // no inputs given
		PrintUsage(argv);
		return 0;
	}

	// makes time retrival happy and loads GRSIEnv
	std::string grsi_path = getenv("GRSISYS");
	if(grsi_path.length() > 0) {
		grsi_path += "/";
	}
	grsi_path += ".grsirc";
	gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

	TParserLibrary::Get()->Load();

	for (auto i = 1; i < argc; i++) AutoFileDetect(argv[i]);

	if (!gChain) std::cout << "No gChain found" << std::endl;
	if (!gChain->GetEntries()) std::cout << "Found gChain, but no entries retrieved" << std::endl;

	if (!gChain || !gChain->GetEntries()) {
		std::cerr << "Failed to find anything. Exiting" << std::endl;
		return 1;
	}

	int process_check;
	process_check = ProcessData();

	if (process_check != 0) {
		std::cerr << "Data did not process correctly ... exiting" << std::endl;
		return 1;
	}
	return 0;

} // main

/******************************************************************************
 * Process events from gamma-gamma tree
 *
 *****************************************************************************/
int ProcessData(){
	std::string fName = gChain->GetCurrentFile()->GetName();
	int run_number = GetRunNumber(fName.c_str());

	std::cout << "Processing run " << run_number << " with " << gChain->GetNtrees() << " file(s)" << std::endl;

    HistogramManager histo_man;
    histo_man.MakeHistogramFile(gChain, lin_coeff_file);

	return 0;
} // ProcessData

/******************************************************************************
 * Detects type of input file
 *
 *****************************************************************************/
void AutoFileDetect(std::string fileName){
	size_t dot_pos = fileName.find_last_of('.');
	std::string ext = fileName.substr(dot_pos + 1);

	if (ext == "root") {
		OpenRootFile(fileName);
	}
	else if (ext == "cal") {
		notifier->AddCalFile(fileName);
    }
	else if (ext == "txt") {
		lin_coeff_file = fileName;
		std::cout << "Found secondary calibration file: " << fileName.c_str() << std::endl;
	} else {
		std::cerr << "Discarding unknown file: " << fileName.c_str() << std::endl;
	}
} // End AutoFileDetect

/************************************************************//**
 * Opens Root files
 *
 * @param dir_name   Directory name
 * @param ext        Extension of the files (default = ".root")
 ***************************************************************/
void OpenRootFile(std::string fileName){
	TFile f(fileName.c_str());
	if (f.Get("AnalysisTree")) {
		if (!gChain) {
			gChain = new TChain("AnalysisTree");
			notifier->AddChain(gChain);
			gChain->SetNotify(notifier);
		}
		gChain->Add(fileName.c_str());
		std::cout << "Added: " << fileName << std::endl;
	}
} // end OpenRootFile

/******************************************************************************
 * Prints usage message and version
 *****************************************************************************/
void PrintUsage(char* argv[]){
	std::cerr << argv[0] << " Version: " << SumPeakHistograms_VERSION_MAJOR
	          << "." << SumPeakHistograms_VERSION_MINOR << "\n"
	          << "usage: " << argv[0] << " calibration_file analysis_tree [analysis_tree_2 ... ] linear_parameter_file\n"
	          << " calibration_file:       calibration file (must end with .cal)\n"
	          << " analysis_tree:          analysis tree to process (must end with .root)\n"
	          << " linear_parameter_file:  contains secondary linear parameters (must end with .txt)"
	          << std::endl;
} // end PrintUsage
