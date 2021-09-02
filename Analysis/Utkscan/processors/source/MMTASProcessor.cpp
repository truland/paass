/** @file MMTASProcessor.cpp
 *  @brief Basic ROOT output. Fills root tree to mimic CAEN's CoMPASS output. It has NO damm output
 *  @authors T. Ruland
 *  @date 07/27/2021
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Globals.hpp"
#include "DetectorDriver.hpp"
#include "DetectorLibrary.hpp"
#include "HelperFunctions.hpp"
#include "RawEvent.hpp"
#include "MMTASProcessor.hpp"

using namespace std;

MMTASProcessor::MMTASProcessor(bool recordtrace) : EventProcessor() {
	associatedTypes.insert("MMTAS");
	
	RecordTraces = recordtrace;

	Rev = Globals::get()->GetPixieRevision();
	std::string name = Globals::get()->GetOutputPath() + Globals::get()->GetOutputFileName() + "_MMTAS.root";
	outputfile = new TFile(name.c_str(),"RECREATE");
	outputtree = new TTree("Event","Output from MMTASProcess in utkscan. For ORNL Pixie setup");
	
	for( auto ii = 0; ii < NumNaILeft; ++ii ){
		string currdet = "NaI"_+"L_"+to_string(ii);
		HitsMap[currdet] = MMTASSingleDetector();
		outputtree->Branch((currdet+"_Energy").c_str(),&HitsMap[currdet].energy);
		outputtree->Branch((currdet+"_Energy_RAW").c_str(),&HitsMap[currdet].rawenergy);
		outputtree->Branch((currdet+"_Timestamp").c_str(),&HitsMap[currdet].timestamp);
		outputtree->Branch((currdet+"_FineTimestamp").c_str(),&HitsMap[currdet].finetimestamp);
		outputtree->Branch((currdet+"_PSD").c_str(),&HitsMap[currdet].psd);
		if( RecordTraces)
			outputtree->Branch((currdet+"_Trace").c_str(),&HitsMap[currdet].trace);
	}
	for( auto ii = 0; ii < NumNaIRight; ++ii ){
		string currdet = "NaI"_+"R_"+to_string(ii);
		HitsMap[currdet] = MMTASSingleDetector();
		outputtree->Branch((currdet+"_Energy").c_str(),&HitsMap[currdet].energy);
		outputtree->Branch((currdet+"_Energy_RAW").c_str(),&HitsMap[currdet].rawenergy);
		outputtree->Branch((currdet+"_Timestamp").c_str(),&HitsMap[currdet].timestamp);
		outputtree->Branch((currdet+"_FineTimestamp").c_str(),&HitsMap[currdet].finetimestamp);
		outputtree->Branch((currdet+"_PSD").c_str(),&HitsMap[currdet].psd);
		if( RecordTraces)
			outputtree->Branch((currdet+"_Trace").c_str(),&HitsMap[currdet].trace);
	}
	for( auto ii = 0; ii < NumBeesLeft; ++ii ){
		string currdet = "BeES_"_+"L_"+to_string(ii);
		HitsMap[currdet] = MMTASSingleDetector();
		outputtree->Branch((currdet+"_Energy").c_str(),&HitsMap[currdet].energy);
		outputtree->Branch((currdet+"_Energy_RAW").c_str(),&HitsMap[currdet].rawenergy);
		outputtree->Branch((currdet+"_Timestamp").c_str(),&HitsMap[currdet].timestamp);
		outputtree->Branch((currdet+"_FineTimestamp").c_str(),&HitsMap[currdet].finetimestamp);
		outputtree->Branch((currdet+"_PSD").c_str(),&HitsMap[currdet].psd);
		if( RecordTraces)
			outputtree->Branch((currdet+"_Trace").c_str(),&HitsMap[currdet].trace);
	}
	for( auto ii = 0; ii < NumBeesRight; ++ii ){
		string currdet = "BeES"_+"R_"+to_string(ii);
		HitsMap[currdet] = MMTASSingleDetector();
		outputtree->Branch((currdet+"_Energy").c_str(),&HitsMap[currdet].energy);
		outputtree->Branch((currdet+"_Energy_RAW").c_str(),&HitsMap[currdet].rawenergy);
		outputtree->Branch((currdet+"_Timestamp").c_str(),&HitsMap[currdet].timestamp);
		outputtree->Branch((currdet+"_FineTimestamp").c_str(),&HitsMap[currdet].finetimestamp);
		outputtree->Branch((currdet+"_PSD").c_str(),&HitsMap[currdet].psd);
		if( RecordTraces)
			outputtree->Branch((currdet+"_Trace").c_str(),&HitsMap[currdet].trace);
	}
}

bool MMTASProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::PreProcess(event))
		return false;

	return true;
}

bool MMTASProcessor::Process(RawEvent &event) {
	if (!EventProcessor::Process(event))
		return false;

	static const auto &Events = event.GetSummary("MMTAS", true)->GetList();
	string currdet;
	for (auto it = Events.begin(); it != Events.end(); it++) {
		if (Rev == "F") {
			timestamp =  (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds((*it)->GetChanID().GetModFreq()) * 1e12;
		} else {
			timestamp =  (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds() * 1e12;
		}

		energy =  (*it)->GetEnergy();
		channel = (*it)->GetChannelNumber();
		board = (*it)->GetModuleNumber();

		if( RecordTraces ){
			if ((*it)->GetTrace().size() > 0) {
				unsigned int len_arr = (*it)->GetTrace().size();
				samples->Set(len_arr);
				auto trace = (*it)->GetTrace();
				for( unsigned int ii = 0; ii < len_arr; ++ii )
					samples->SetAt(trace.at(ii),ii);
			}
		}
		//if (!(*it)->GetQdc().empty()) {
			//set energyshort in here
			//currently is a no-op
		//}
	}
	outputtree->Fill();
	EndProcess();
	return true;
}
