/** @file CompassProcessor.cpp
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
#include "CompassProcessor.hpp"

using namespace std;

CompassProcessor::CompassProcessor() : EventProcessor() {
	associatedTypes.insert("Compass");
	Rev = Globals::get()->GetPixieRevision();
	std::string name = Globals::get()->GetOutputPath() + Globals::get()->GetOutputFileName() + "_Compass.root";
	outputfile = new TFile(name.c_str(),"RECREATE");
	outputtree = new TTree("Data","Output from CompassProcess in utkscan to mimic CAEN CoMPASS");
	outputtree->Branch("Channel",&channel);
	outputtree->Branch("Board",&board);
	outputtree->Branch("Energy",&energy);
	outputtree->Branch("EnergyShort",&energyshort);
	outputtree->Branch("Flags",&flags);
	outputtree->Branch("Timestamp",&timestamp);
	samples = new TArrayS(1024);
	outputtree->Branch("Samples",&samples);
}

bool CompassProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::PreProcess(event))
		return false;

	return true;
}

bool CompassProcessor::Process(RawEvent &event) {
	if (!EventProcessor::Process(event))
		return false;

	static const auto &Events = event.GetSummary("Compass", true)->GetList();
	for (auto it = Events.begin(); it != Events.end(); it++) {
		if (Rev == "F") {
			timestamp =  (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds((*it)->GetChanID().GetModFreq()) * 1e12;
		} else {
			timestamp =  (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds() * 1e12;
		}
		energy =  (*it)->GetEnergy();
		energyshort = 0; //temporary while the qdcSums are off
		channel = (*it)->GetChannelNumber();
		board = (*it)->GetModuleNumber();
		flags = 0;
		if( (*it)->IsPileup() )
			flags |= 0x8000;
		if( (*it)->IsSaturated() )
			flags |= 0x80;

		if ((*it)->GetTrace().size() > 0) {
			unsigned int len_arr = (*it)->GetTrace().size();
			samples->Set(len_arr);
			for( unsigned int ii = 0; ii < len_arr; ++ii )
				samples->SetAt((*it)->GetTrace().at(ii),ii);
		}
		if (!(*it)->GetQdc().empty()) {
			//set energyshort in here
			//currently is a no-op
		}
		outputtree->Fill();
	}
	//outputtree->Write();
        //PixieFile = PTree->GetCurrentFile();
        //PixieFile->Write();
        //PixieFile->Close();
	EndProcess();
	return true;
}
