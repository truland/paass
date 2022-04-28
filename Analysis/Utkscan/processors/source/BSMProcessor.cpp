/** @file BSMProcessor.cpp
 * @brief  Basic BSMProcessor
 * @author T. Ruland
 * @date 04/25/2022
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include "DetectorDriver.hpp"
#include "DetectorLibrary.hpp"
#include "Globals.hpp"
#include "HelperFunctions.hpp"
#include "StringManipulationFunctions.hpp"
#include "RawEvent.hpp"
#include "BSMProcessor.hpp"
#include "DammPlotIds.hpp"

namespace dammIds {
	namespace bsm {
		//6500 is beginning, 500 range
		const unsigned TOTAL_OFFSET = 0;
		const unsigned DD_OFFSET = 250;

		const unsigned D_TIMING_OFFSET = 300;

		const unsigned D_BSM_TOTAL = TOTAL_OFFSET; 
		const unsigned D_BSM_MTAS_SUM = TOTAL_OFFSET + 1;
		const unsigned D_BSM_ZERO_MTAS = TOTAL_OFFSET + 2;

		const unsigned DD_BSM_MTAS_TOTAL = DD_OFFSET; 

		const unsigned D_TDIFF_EVENTS = D_TIMING_OFFSET;
	}
}
using namespace std;
using namespace dammIds::bsm;

void BSMProcessor::DeclarePlots(void){
	DeclareHistogram1D(D_BSM_TOTAL,SE,"BSM Total");
	DeclareHistogram1D(D_BSM_MTAS_SUM,SE,"BSM Total + MTAS Total");
	DeclareHistogram1D(D_BSM_ZERO_MTAS,SE,"BSM Total No MTAS");
	DeclareHistogram2D(DD_BSM_MTAS_TOTAL,SD,SD,"BSM Total vs MTAS Total");
	DeclareHistogram1D(D_TDIFF_EVENTS,SE,"Tdiff between bsm events");
}


BSMProcessor::BSMProcessor(int numsegments,bool zerosuppress,bool alone) : EventProcessor(OFFSET, RANGE, "BSMProcessor") {
	associatedTypes.insert("bsm");
	PixieRev = Globals::get()->GetPixieRevision();
	NumSegments = numsegments;
	HasZeroSuppression = zerosuppress;
	StandAlone = alone;
	FoundFirst = false;
}

bool BSMProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::PreProcess(event))
		return false;

	static const auto &chanEvents = event.GetSummary("bsm", true)->GetList();
	vector<BSMSegment> BSMSegVec(NumSegments,BSMSegment(HasZeroSuppression));
	vector<short> BSMSegMulti(2*NumSegments,0); // MTAS segment multiplicity "map"

	double EarliestTime = 1.0e99;
	for (auto chanEvtIter = chanEvents.begin(); chanEvtIter != chanEvents.end(); ++chanEvtIter){
		int segmentNum = stoi((*chanEvtIter)->GetChanID().GetGroup().c_str()) - 1;

		bool isFront = (*chanEvtIter)->GetChanID().HasTag("front");
		bool isBack = (*chanEvtIter)->GetChanID().HasTag("back");
		int chanOffset;
		if (isFront) {
			chanOffset = 0;
		} else if (isBack) {
			chanOffset = 1;
		} else if ((isFront && isBack) || (!isBack && !isFront) ) {
			cout<<"ERROR::BSMProcessor:PreProcess BOTHESSSSS ("<<(*chanEvtIter)->GetModuleNumber()<<" , " << (*chanEvtIter)->GetChannelNumber() << ") !"<<endl;
			return false;
		} else {
			chanOffset = -9999;
		}
		if (chanOffset == -9999 ){
			cout<<"ERROR::BSMProcessor:PreProcess Channel ("<<(*chanEvtIter)->GetModuleNumber()<<" , " << (*chanEvtIter)->GetChannelNumber() << ") found which doesnt have a front or back tag, or you didnt set the Ring right (xml subtype). This means the XML is not right so you need to fix it!!"<<endl;
			return false;
		}
		
		int GlobalChanID = 2*(segmentNum) + chanOffset;

		//THE SATURATE AND PILEUP CHECK SHOULD NOT BE PERFORMED ON LOGIC SIGNALS IN HISTORICAL DATA
		//THIS IS BECAUSE ALL BUT THE MTC SIGNAL TRIP ONE OF THESE FLAGS SO BE CAREFUL
		//YOU'VE BEEN WARNED.
		//-THOMAS RULAND 04/25/2022
		if( (*chanEvtIter)->IsSaturated() || (*chanEvtIter)->IsPileup() or (*chanEvtIter)->GetEnergy() > 30000 ){
			continue;
		} else {
			BSMSegMulti.at(GlobalChanID)++; // increment the multipliciy "map" based on GlobalMtasSegID

			BSMSegVec.at(segmentNum).gBSMSegID_ = segmentNum;
			if(isFront && BSMSegVec.at(segmentNum).segFront_ == nullptr){  
				if( (*chanEvtIter)->GetTimeSansCfd() < EarliestTime )
					EarliestTime = (*chanEvtIter)->GetTimeSansCfd(); 
				BSMSegVec.at(segmentNum).segFront_ = (*chanEvtIter);
				BSMSegVec.at(segmentNum).PixieRev = PixieRev;
			}
			//! Thomas Ruland Gets a gold star 
			else if (isBack && BSMSegVec.at(segmentNum).segBack_ == nullptr) { 
				if( (*chanEvtIter)->GetTimeSansCfd() < EarliestTime )
					EarliestTime = (*chanEvtIter)->GetTimeSansCfd(); 
				BSMSegVec.at(segmentNum).segBack_ = (*chanEvtIter);
				BSMSegVec.at(segmentNum).PixieRev = PixieRev;
			}
		}
	}

	if( not FoundFirst ){
		FoundFirst = true;
		CurrTime = EarliestTime;
	}else{
		PreviousTime = CurrTime;
		CurrTime = EarliestTime;
		plot(D_TDIFF_EVENTS,(CurrTime - PreviousTime));
	}

	//reset this during pre-process
	BSMTotal = make_pair(0,not HasZeroSuppression);
	int NumFire = 0;
	for( auto& segIter : BSMSegVec ){
		if( segIter.IsValidSegment() )
			++NumFire;
	}
	
	for( auto& segIter : BSMSegVec ){
		auto result = segIter.GetSegmentAverageEnergy();
		if( not result.second ){
			continue;
		}else{
			BSMTotal.first += result.first/static_cast<double>(NumFire);
			BSMTotal.second = true;	
		}
	}

	if( BSMTotal.second )
		plot(D_BSM_TOTAL,BSMTotal.first);

	EventData TotalData(EarliestTime,BSMTotal.first);
	TreeCorrelator::get()->place("BSM_Total")->activate(TotalData);

	return true;
}

bool BSMProcessor::Process(RawEvent &event) {
	if (!EventProcessor::Process(event))
		return false;

	if(TreeCorrelator::get()->checkPlace("MTAS_Total")){
		if(TreeCorrelator::get()->place("MTAS_Total")->status()){
			PlaceDetector* mtas_total = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_Total")); 
			//for( auto& e : mtas_total->info_ )
			//	cout << e.energy << " " << e.time << '\t';
			//cout << endl;
			if( BSMTotal.second and mtas_total->info_.size() > 0 ){
				plot(DD_BSM_MTAS_TOTAL,mtas_total->last().energy,BSMTotal.first);
				plot(D_BSM_MTAS_SUM,mtas_total->last().energy + BSMTotal.first);
			}
		}else{
			if( BSMTotal.second ){
					plot(D_BSM_ZERO_MTAS,BSMTotal.first);
					plot(DD_BSM_MTAS_TOTAL,0.0,BSMTotal.first);
					plot(D_BSM_MTAS_SUM,BSMTotal.first);
			}
		}
	}else{
		if( not StandAlone )
			throw "ERROR BSMProcessor::Process BSM isn't in standalone mode and the treecorrelator for MTAS_Total isn't setup";
	}

	EndProcess();
	return true;
}
