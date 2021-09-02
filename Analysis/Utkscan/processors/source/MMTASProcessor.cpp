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

MMTASProcessor::MMTASProcessor() : EventProcessor() {
	associatedTypes.insert("MMTAS");
	
	Rev = Globals::get()->GetPixieRevision();
	std::string name = Globals::get()->GetOutputPath() + Globals::get()->GetOutputFileName() + "_MMTAS.root";
	outputfile = new TFile(name.c_str(),"RECREATE");
	outputtree = new TTree("Event","Output from MMTASProcess in utkscan. For ORNL Pixie setup");
	outputtree->Branch("mmtas_vec_",&mmtas_vec_);
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
		MMstruct.energy = (*it)->GetCalibratedEnergy();
		MMstruct.rawEnergy = (*it)->GetEnergy();
		if (Rev == "F") {
			MMstruct.timeSansCfd = (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds((*it)->GetChanID().GetModFreq()) * 1e9;
			MMstruct.time = (*it)->GetTime() * Globals::get()->GetAdcClockInSeconds((*it)->GetChanID().GetModFreq()) * 1e9;
		} else {
			MMstruct.timeSansCfd = (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds() * 1e9;
			MMstruct.time = (*it)->GetTime() * Globals::get()->GetAdcClockInSeconds() * 1e9;
		}
		MMstruct.detNum = (*it)->GetChanID().GetLocation();
		MMstruct.modNum = (*it)->GetModuleNumber();
		MMstruct.chanNum = (*it)->GetChannelNumber();
		MMstruct.subtype = (*it)->GetChanID().GetSubtype();
		MMstruct.group = (*it)->GetChanID().GetGroup();
		MMstruct.pileup = (*it)->IsPileup();
		MMstruct.saturation = (*it)->IsSaturated();

		if ((*it)->GetTrace().size() > 0) {
			auto trace = (*it)->GetTrace();
			MMstruct.hasValidTimingAnalysis = trace.HasValidTimingAnalysis();
			MMstruct.hasValidWaveformAnalysis = trace.HasValidWaveformAnalysis();
			MMstruct.baseline = trace.GetBaselineInfo().first;
			MMstruct.stdBaseline = trace.GetBaselineInfo().second;
			MMstruct.trace = trace;
			MMstruct.maxPos = trace.GetMaxInfo().first;
			MMstruct.maxVal = trace.GetMaxInfo().second;
			MMstruct.extMaxVal = trace.GetExtrapolatedMaxInfo().second;
			MMstruct.tqdc = trace.GetQdc();
			MMstruct.highResTime = (*it)->GetHighResTimeInNs();
			if (Rev == "F") {
				MMstruct.phase = trace.GetPhase() * Globals::get()->GetAdcClockInSeconds((*it)->GetChanID().GetModFreq()) * 1e9;
			} else {
				MMstruct.phase = trace.GetPhase() * Globals::get()->GetAdcClockInSeconds() * 1e9;
			}
		}
		if (!(*it)->GetQdc().empty()) {
			MMstruct.qdcSums = (*it)->GetQdc();
		}
		mmtas_vec_.push_back(MMstruct);
		//pixie_tree_event_->root_dev_vec_.emplace_back(MMstruct);
		MMstruct = MMTAS_DEFAULT_STRUCT;

	}
	outputtree->Fill();
	mmtas_vec_.clear();
	EndProcess();
	return true;
}
