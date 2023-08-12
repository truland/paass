/** @file MtasProcessor.cpp
 * @brief  Basic MtasProcessor for MTAS at FRIB
 * @authors T.T. King, T. Ruland, B.C. Rasco
 * @date 03/25/2022
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
#include "MtasProcessor.hpp"
#include "DammPlotIds.hpp"

namespace dammIds {
	namespace mtas {
		//7200 is beginning
		const unsigned TOTAL_OFFSET = 0;
		const unsigned CENTER_OFFSET = 10;
		const unsigned INNER_OFFSET = 20;
		const unsigned MIDDLE_OFFSET = 30;
		const unsigned OUTER_OFFSET = 40;
		const unsigned DD_OFFSET = 50;
		const unsigned DD_LIGHT_SHARE_OFFSET = 60;
		const unsigned D_SINGLE_CENTER_OFFSET = 70;

		const unsigned INDIVIDUALS_OFFSET = 5; 
		const unsigned BETA_GATED_OFFSET = 100; 
		const unsigned D_MTAS_SUM_FB = 200;
		const unsigned D_MTAS_TDIFF_OFFSET = 300;
		const unsigned D_ONE_OFFS = 400; 
		const unsigned D_CONICIDENCE = 500;

		const unsigned D_MTAS_TOTAL = TOTAL_OFFSET;
		const unsigned D_MTAS_CENTER = CENTER_OFFSET;
		const unsigned D_MTAS_INNER = INNER_OFFSET;
		const unsigned D_MTAS_MIDDLE = MIDDLE_OFFSET;
		const unsigned D_MTAS_OUTER = OUTER_OFFSET;

		const unsigned D_MTAS_BETA_TOTAL = TOTAL_OFFSET + BETA_GATED_OFFSET;
		const unsigned D_MTAS_BETA_CENTER = CENTER_OFFSET + BETA_GATED_OFFSET;
		const unsigned D_MTAS_BETA_INNER = INNER_OFFSET + BETA_GATED_OFFSET;
		const unsigned D_MTAS_BETA_MIDDLE = MIDDLE_OFFSET + BETA_GATED_OFFSET;
		const unsigned D_MTAS_BETA_OUTER = OUTER_OFFSET + BETA_GATED_OFFSET;

		const unsigned D_MTAS_TOTAL_INDI = TOTAL_OFFSET + INDIVIDUALS_OFFSET;
		const unsigned D_MTAS_CENTER_INDI = CENTER_OFFSET + INDIVIDUALS_OFFSET;
		const unsigned D_MTAS_INNER_INDI = INNER_OFFSET + INDIVIDUALS_OFFSET;
		const unsigned D_MTAS_MIDDLE_INDI = MIDDLE_OFFSET + INDIVIDUALS_OFFSET;
		const unsigned D_MTAS_OUTER_INDI = OUTER_OFFSET + INDIVIDUALS_OFFSET;

		const unsigned DD_MTAS_CS_T = DD_OFFSET; //! MTAS  Center Segments vs Total Sum
		const unsigned DD_MTAS_CR_T = DD_OFFSET + 1; //! MTAS Center Ring Sum vs Total Sum
		const unsigned DD_MTAS_CS_CR = DD_OFFSET + 2; //! MTAS Center Segments vs Center Ring Sum
		const unsigned DD_MTAS_IMO_T = DD_OFFSET + 3;  //! MTAS I,M,O Segments vs Total
		const unsigned DD_MTAS_C1_POSITION = DD_OFFSET + 4; //!MTAS C1 energy vs position
		const unsigned DD_MTAS_C1_C2 = DD_LIGHT_SHARE_OFFSET + 0; //! MTAS C1 VS C2
		const unsigned DD_MTAS_C1_C6 = DD_LIGHT_SHARE_OFFSET + 1; //! MTAS C1 vs C2
		const unsigned DD_MTAS_C3_C2 = DD_LIGHT_SHARE_OFFSET + 2; //! MTAS C3 vs C2
		const unsigned DD_MTAS_C3_C4 = DD_LIGHT_SHARE_OFFSET + 3; //! MTAS C3 vs C4
		const unsigned DD_MTAS_C5_C6 = DD_LIGHT_SHARE_OFFSET + 4; //! MTAS C5 vs C6
		const unsigned DD_MTAS_C5_C4 = DD_LIGHT_SHARE_OFFSET + 5; //! MTAS C5 vs C4
		const unsigned DD_MTAS_C1_C3 = DD_LIGHT_SHARE_OFFSET + 6; //! MTAS C1 vs C3
		const unsigned DD_MTAS_C1_C4 = DD_LIGHT_SHARE_OFFSET + 7; //! MTAS C1 vs C4
		const unsigned DD_MTAS_C1_C5 = DD_LIGHT_SHARE_OFFSET + 8; //! MTAS C1 vs C5
	}
}
using namespace std;
using namespace dammIds::mtas;

void MtasProcessor::DeclarePlots(void){
	DeclareHistogram1D(D_MTAS_TOTAL,SE, "Mtas TOTAL");
	DeclareHistogram1D(D_MTAS_CENTER,SE, "Mtas CENTER");
	DeclareHistogram1D(D_MTAS_INNER,SE, "Mtas INNER");
	DeclareHistogram1D(D_MTAS_MIDDLE,SE, "Mtas MIDDLE");
	DeclareHistogram1D(D_MTAS_OUTER,SE, "Mtas OUTER");

	DeclareHistogram1D(D_MTAS_BETA_TOTAL,SE, "Mtas Beta TOTAL");
	DeclareHistogram1D(D_MTAS_BETA_CENTER,SE, "Mtas Beta CENTER");
	DeclareHistogram1D(D_MTAS_BETA_INNER,SE, "Mtas Beta INNER");
	DeclareHistogram1D(D_MTAS_BETA_MIDDLE,SE, "Mtas Beta MIDDLE");
	DeclareHistogram1D(D_MTAS_BETA_OUTER,SE, "Mtas Beta OUTER");

	DeclareHistogram1D(D_MTAS_TOTAL_INDI,SE, "Mtas TOTAL Individuals");
	DeclareHistogram1D(D_MTAS_CENTER_INDI,SE, "Mtas CENTER Individuals");
	DeclareHistogram1D(D_MTAS_INNER_INDI,SE, "Mtas INNER Individuals");
	DeclareHistogram1D(D_MTAS_MIDDLE_INDI,SE, "Mtas MIDDLE Individuals");
	DeclareHistogram1D(D_MTAS_OUTER_INDI,SE, "Mtas OUTER Individuals");

	for (unsigned i= 0; i < 6; ++i){ //! to_string(i+1) to account for 0 counting
		std::string title = "Mtas Center " + to_string(i+1) + " Sum F+B";
		DeclareHistogram1D(D_MTAS_SUM_FB + CENTER_OFFSET + i , SE, title.c_str());
		std::string title2 = "Mtas Center " + to_string(i+1) + " FB-TDIFF";
		DeclareHistogram1D(D_MTAS_TDIFF_OFFSET + CENTER_OFFSET + i , SE, title2.c_str());
		std::string title3 = "Mtas Center " + to_string(i+1) + " front when both";
		std::string title4 = "Mtas Center " + to_string(i+1) + " back when both";
		DeclareHistogram1D(D_CONICIDENCE + (2*i), SE, title3.c_str());
		DeclareHistogram1D(D_CONICIDENCE + (2*i + 1), SE, title4.c_str());
		std::string title5 = "Mtas Center " + to_string(i+1) + " single segment fire";
		DeclareHistogram1D(D_SINGLE_CENTER_OFFSET + i,SE,title5.c_str());
		std::string title6 = "Mtas Center " + to_string(i+1) + " no neighbors fire";
		DeclareHistogram1D(D_SINGLE_CENTER_OFFSET + i + 6,SE,title6.c_str());
	}
	for (unsigned i= 0; i < 6; ++i){
		std::string title = "Mtas Inner " + to_string(i+1) + " Sum F+B";
		DeclareHistogram1D(D_MTAS_SUM_FB + INNER_OFFSET + i , SE, title.c_str());
		std::string title2 = "Mtas Inner " + to_string(i+1) + " FB-TDIFF";
		DeclareHistogram1D(D_MTAS_TDIFF_OFFSET + INNER_OFFSET + i , SE, title2.c_str());
		std::string title3 = "Mtas Inner " + to_string(i+1) + " front when both";
		std::string title4 = "Mtas Inner " + to_string(i+1) + " back when both";
		DeclareHistogram1D(D_CONICIDENCE + 12 + (2*i), SE, title3.c_str());
		DeclareHistogram1D(D_CONICIDENCE + 12 + (2*i + 1) , SE, title4.c_str());
	}
	for (unsigned i= 0; i < 6; ++i){
		std::string title = "Mtas Middle " + to_string(i+1) + " Sum F+B";
		DeclareHistogram1D(D_MTAS_SUM_FB + MIDDLE_OFFSET + i , SE, title.c_str());
		std::string title2 = "Mtas Middle " + to_string(i+1) + " FB-TDIFF";
		DeclareHistogram1D(D_MTAS_TDIFF_OFFSET + MIDDLE_OFFSET + i , SE, title2.c_str());
		std::string title3 = "Mtas Middle " + to_string(i+1) + " front when both";
		std::string title4 = "Mtas Middle " + to_string(i+1) + " back when both";
		DeclareHistogram1D(D_CONICIDENCE + 24 + (2*i), SE, title3.c_str());
		DeclareHistogram1D(D_CONICIDENCE + 24 + (2*i + 1), SE, title4.c_str());
	}
	for (unsigned i= 0; i < 6; ++i){
		std::string title = "Mtas Outer " + to_string(i+1) + " Sum F+B";
		DeclareHistogram1D(D_MTAS_SUM_FB + OUTER_OFFSET + i , SE, title.c_str());
		std::string title2 = "Mtas Outer " + to_string(i+1) + " FB-TDIFF";
		DeclareHistogram1D(D_MTAS_TDIFF_OFFSET + OUTER_OFFSET + i , SE, title2.c_str());
		std::string title3 = "Mtas Outer " + to_string(i+1) + " front when both";
		std::string title4 = "Mtas Outer " + to_string(i+1) + " back when both";
		DeclareHistogram1D(D_CONICIDENCE + 36+ (2*i), SE, title3.c_str());
		DeclareHistogram1D(D_CONICIDENCE + 36+ (2*i + 1), SE, title4.c_str());
	}

	DeclareHistogram2D(DD_MTAS_CR_T,SD,SD,"Mtas Center Ring vs Total");
	DeclareHistogram2D(DD_MTAS_CS_CR,SD,SD,"Mtas Center Segment vs Center Ring");
	DeclareHistogram2D(DD_MTAS_CS_T,SD,SD,"Mtas Center Segment vs Total");
	DeclareHistogram2D(DD_MTAS_IMO_T,SD,SD,"Mtas I,M,O Segment vs Total");
	DeclareHistogram2D(DD_MTAS_C1_POSITION,SD,SD,"Mtas C1 vs Position");
	//DeclareHistogram2D(DD_MTAS_C1_C2,SD,SD,"Mtas C1 vs C2");
	//DeclareHistogram2D(DD_MTAS_C1_C6,SD,SD,"Mtas C1 vs C6");
	//DeclareHistogram2D(DD_MTAS_C3_C2,SD,SD,"Mtas C3 vs C2");
	//DeclareHistogram2D(DD_MTAS_C3_C4,SD,SD,"Mtas C3 vs C4");
	//DeclareHistogram2D(DD_MTAS_C5_C4,SD,SD,"Mtas C5 vs C4");
	//DeclareHistogram2D(DD_MTAS_C5_C6,SD,SD,"Mtas C5 vs C6");
	//DeclareHistogram2D(DD_MTAS_C1_C3,SD,SD,"Mtas C1 vs C3");
	//DeclareHistogram2D(DD_MTAS_C1_C4,SD,SD,"Mtas C1 vs C4");
	//DeclareHistogram2D(DD_MTAS_C1_C5,SD,SD,"Mtas C1 vs C5");

	DeclareHistogram2D(D_ONE_OFFS + 0, S6, S6, "Mtas Multiplicity x=ID y=multi");
	DeclareHistogram1D(D_ONE_OFFS + 1,SE,"Isomer Energy");
	DeclareHistogram2D(D_ONE_OFFS + 2,SD,SD,"Time since Last Beta Event vs Total");
	DeclareHistogram2D(D_ONE_OFFS + 3,SD,SD,"Isomeric BSM vs MTAS Total");
}


MtasProcessor::MtasProcessor(bool newcenter,bool zerosuppress,double thresh,int pmtthresh) : EventProcessor(OFFSET, RANGE, "MtasProcessor") {
	associatedTypes.insert("mtas");
	PixieRev = Globals::get()->GetPixieRevision();
	IsNewCenter = newcenter;
	HasZeroSuppression = zerosuppress;
	BetaThreshold = thresh;
	IsPrevBetaTriggered = false;
	NumPMTFire = pmtthresh; 
}

bool MtasProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::PreProcess(event))
		return false;

	static const auto &chanEvents = event.GetSummary("mtas", true)->GetList();

	MtasSegVec = vector<MtasSegment>(24, MtasSegment(HasZeroSuppression));
	vector<short> MtasSegMulti(48,0); // MTAS segment multiplicity "map"

	double EarliestTime = 1.0e99;
	for (auto chanEvtIter = chanEvents.begin(); chanEvtIter != chanEvents.end(); ++chanEvtIter){

		if (PixieRev == "F"){
			clockInSeconds = Globals::get()->GetClockInSeconds((*chanEvtIter)->GetChanID().GetModFreq());
		} else {
			clockInSeconds = Globals::get()->GetClockInSeconds();
		}
		//! needs try/catch for non numeric string in group
		int segmentNum = stoi((*chanEvtIter)->GetChanID().GetGroup().c_str());
		string Ring = StringManipulation::StringLower((*chanEvtIter)->GetChanID().GetSubtype());
		int RingOffset = -9999;
		if (strcmp(Ring.c_str(), "center") == 0) {
			RingOffset = -1;
		}
		if (strcmp(Ring.c_str(), "inner") == 0) {
			RingOffset = 5;
		}
		if (strcmp(Ring.c_str(), "middle") == 0) {
			RingOffset = 11;
		}
		if (strcmp(Ring.c_str(), "outer") == 0) {
			RingOffset = 17;
		}

		bool isFront = (*chanEvtIter)->GetChanID().HasTag("front");
		bool isBack = (*chanEvtIter)->GetChanID().HasTag("back");
		int chanOffset;
		if (isFront) {
			chanOffset = 0;
		} else if (isBack) {
			chanOffset = 1;
		} else if ((isFront && isBack) || (!isBack && !isFront) ) {
			cout<<"ERROR::MtasProcessor:PreProcess BOTHESSSSS ("<<(*chanEvtIter)->GetModuleNumber()<<" , " << (*chanEvtIter)->GetChannelNumber() << ") !"<<endl;
			return false;
		} else {
			chanOffset = -9999;
		}
		if (chanOffset == -9999 || RingOffset == -9999){
			cout<<"ERROR::MtasProcessor:PreProcess Channel ("<<(*chanEvtIter)->GetModuleNumber()<<" , " << (*chanEvtIter)->GetChannelNumber() << ") found which doesnt have a front or back tag, or you didnt set the Ring right (xml subtype). This means the XML is not right so you need to fix it!!"<<endl;
			return false;
		}

		int GlobalMtasSegID = RingOffset + segmentNum;
		int GlobalMtasChanID = (segmentNum + RingOffset) * 2 + chanOffset;

		//THE SATURATE AND PILEUP CHECK SHOULD NOT BE PERFORMED ON LOGIC SIGNALS IN HISTORICAL DATA
		//THIS IS BECAUSE ALL BUT THE MTC SIGNAL TRIP ONE OF THESE FLAGS SO BE CAREFUL
		//YOU'VE BEEN WARNED.
		//-THOMAS RULAND 04/25/2022
		if( (*chanEvtIter)->IsSaturated() || (*chanEvtIter)->IsPileup() or (*chanEvtIter)->GetEnergy() > 30000 ){
			continue;
		} else {
			MtasSegMulti.at(GlobalMtasChanID)++; // increment the multipliciy "map" based on GlobalMtasSegID

			MtasSegVec.at(GlobalMtasSegID).gMtasSegID_ = GlobalMtasSegID;
			MtasSegVec.at(GlobalMtasSegID).segRing_ = Ring;
			if(isFront && MtasSegVec.at(GlobalMtasSegID).segFront_ == nullptr){  
				if( (*chanEvtIter)->GetTimeSansCfd() < EarliestTime )
					EarliestTime = (*chanEvtIter)->GetTimeSansCfd(); 
				MtasSegVec.at(GlobalMtasSegID).segFront_ = (*chanEvtIter);
				MtasSegVec.at(GlobalMtasSegID).PixieRev = PixieRev;
			}
			//! Thomas Ruland Gets a gold star 
			else if (isBack && MtasSegVec.at(GlobalMtasSegID).segBack_ == nullptr) { 
				if( (*chanEvtIter)->GetTimeSansCfd() < EarliestTime )
					EarliestTime = (*chanEvtIter)->GetTimeSansCfd(); 
				MtasSegVec.at(GlobalMtasSegID).segBack_ = (*chanEvtIter);
				MtasSegVec.at(GlobalMtasSegID).PixieRev = PixieRev;
			}
		}
	}  //! end loop over chanEvents.

	//! begin loop over segments for sums
	centerSum = make_pair(0,not HasZeroSuppression);
	innerSum = make_pair(0,not HasZeroSuppression);
	middleSum = make_pair(0,not HasZeroSuppression);
	outerSum = make_pair(0,not HasZeroSuppression);
	totalSum = make_pair(0,not HasZeroSuppression);
	int NumCenter = 0;
	int CurrNumCenter = 0;
	for( auto ii = 0; ii < 6; ++ii ){
		if( MtasSegVec.at(ii).IsValidSegment() ){
			NumCenter++;
			CurrNumCenter += 2;
		}else if( MtasSegVec.at(ii).IsValidFront() || MtasSegVec.at(ii).IsValidBack() ){
			CurrNumCenter++;
		}
	}
	double segmentAvg = 0.0;
	double segTdiff = 0.0;
	double position = 0.0;
	double frontenergy = 0.0;
	double backenergy = 0.0;
	for (auto segIter = MtasSegVec.begin(); segIter != MtasSegVec.end(); ++segIter) {
		int segmentID = segIter->gMtasSegID_;
		auto result = segIter->GetSegmentAverageEnergy();
		if( not result.second ){
			continue;
		}else{
			totalSum.second = true;
			segmentAvg = result.first;
		}
		if (segmentID >= 0 && segmentID <= 5 ){
			centerSum.second = true;
			if( IsNewCenter ){
				totalSum.first += segmentAvg;
				centerSum.first += segmentAvg;
			}else{
				if( CurrNumCenter >= NumPMTFire ){
					totalSum.first += segmentAvg/static_cast<double>(NumCenter);
					centerSum.first += segmentAvg/static_cast<double>(NumCenter);
				}
			}
		} else if (segmentID >= 6 && segmentID <= 11 ){
			innerSum.second = true;
			totalSum.first += segmentAvg;
			innerSum.first += segmentAvg;
		} else if (segmentID >= 12 && segmentID <= 17 ){
			middleSum.second = true;
			totalSum.first += segmentAvg;
			middleSum.first += segmentAvg;
		} else if (segmentID >= 18 && segmentID <= 23 ){
			outerSum.second = true;
			totalSum.first += segmentAvg;
			outerSum.first += segmentAvg;
		}
	}
	//! loop over segments for ploting
	for (auto segIter = MtasSegVec.begin(); segIter != MtasSegVec.end(); ++segIter) {
		int segmentID = segIter->gMtasSegID_;
		auto segmentAvgresult = segIter->GetSegmentAverageEnergy();
		auto segTdiffresult = segIter->GetSegmentTdiffInNS();
		auto positionresult = segIter->GetSegmentPosition();
		auto FrontEnergyResult = segIter->GetFrontEnergy();
		auto BackEnergyResult = segIter->GetBackEnergy();
		if( (not segmentAvgresult.second) or (not segTdiffresult.second) or (not positionresult.second) or (not FrontEnergyResult.second) or (not BackEnergyResult.second)){
			continue;
		}else{
			segmentAvg = segmentAvgresult.first;
			segTdiff = segTdiffresult.first;
			position = ((SD/2)*(1.0 + positionresult.first));
			frontenergy = FrontEnergyResult.first;
			backenergy = BackEnergyResult.first;
		}
		//! Begin Root Output stuff.
		if (DetectorDriver::get()->GetSysRootOutput() && segIter->IsValidSegment()) {
			Mtasstruct.energy = segIter->GetSegmentAverageEnergy().first;
			Mtasstruct.fEnergy = segIter->GetFrontEnergy().first;
			Mtasstruct.bEnergy = segIter->GetBackEnergy().first;
			Mtasstruct.time = (segIter->GetFrontTimeInNS().first + segIter->GetBackTimeInNS().first)/2.0;
			Mtasstruct.tdiff = (segIter->GetFrontTimeInNS().first - segIter->GetBackTimeInNS().first);
			Mtasstruct.gSegmentID = segIter->gMtasSegID_;
			Mtasstruct.segmentNum = (segIter->gMtasSegID_+1)/2;
			Mtasstruct.Ring = segIter->segRing_;
			if (strcmp(segIter->segRing_.c_str(), "center") == 0 ) {
				Mtasstruct.RingNum = 1;
			} else if (strcmp(segIter->segRing_.c_str(), "inner") == 0 ) {
				Mtasstruct.RingNum = 2;
			} else if (strcmp(segIter->segRing_.c_str(), "middle") == 0 ) {
				Mtasstruct.RingNum = 3;
			} else if (strcmp(segIter->segRing_.c_str(), "outer") == 0 ) {
				Mtasstruct.RingNum = 4;
			}
			
			pixie_tree_event_->mtas_vec_.emplace_back(Mtasstruct);
			Mtasstruct = processor_struct::MTAS_DEFAULT_STRUCT;
		}
	
		if (segmentID >= 0 && segmentID <= 5 ){
			//D_MTAS_SUM_FB + 50 + OUTER_OFFSET + (2* i ),

			plot(D_MTAS_TDIFF_OFFSET + CENTER_OFFSET + segmentID, segTdiff + SE/2 );

			plot(D_CONICIDENCE+ (2 * segmentID), frontenergy);
			plot(D_CONICIDENCE+ (2 * segmentID + 1), backenergy);

			//! per segment plots
			plot(D_MTAS_SUM_FB + CENTER_OFFSET + segmentID, segmentAvg);
			plot(D_MTAS_CENTER_INDI, segmentAvg);
			plot(DD_MTAS_CS_T, totalSum.first, segmentAvg);
			plot(DD_MTAS_CS_CR,centerSum.first,segmentAvg);
			
			if( segmentID == 0 )
				plot(DD_MTAS_C1_POSITION,position,segmentAvg);
		} else if (segmentID >= 6 && segmentID <= 11 ){
			plot(D_MTAS_TDIFF_OFFSET + INNER_OFFSET + (segmentID - 6), segTdiff + SE / 2);

			plot(D_CONICIDENCE+ 12 + (2 * (segmentID-6)), frontenergy);
			plot(D_CONICIDENCE+ 12 + (2 * (segmentID-6) + 1), backenergy);

			//! per segment plots
			plot(D_MTAS_SUM_FB + INNER_OFFSET + (segmentID - 6), segmentAvg);  //! subtract 12 to put it back in segments 1 - 6
			plot(DD_MTAS_IMO_T, totalSum.first, segmentAvg);
		} else if (segmentID >= 12 && segmentID <= 17 ){
			plot(D_MTAS_TDIFF_OFFSET + MIDDLE_OFFSET +  (segmentID-12), segTdiff + SE/2 );

			plot(D_CONICIDENCE+ 24 + (2 * (segmentID - 12)), frontenergy);
			plot(D_CONICIDENCE+ 24 + (2 * (segmentID - 12) + 1), backenergy);

			//! per segment plots
			plot(D_MTAS_SUM_FB + MIDDLE_OFFSET + (segmentID-12), segmentAvg); //! subtract 12 to put it back in segments 1 - 6 
			plot(DD_MTAS_IMO_T,totalSum.first,segmentAvg);
		}else if (segmentID >= 18 && segmentID <= 23 ){
			plot(D_MTAS_TDIFF_OFFSET + OUTER_OFFSET +  (segmentID-18), segTdiff + SE/2 );

			plot(D_CONICIDENCE+ 36 + (2 * (segmentID-18)), frontenergy);
			plot(D_CONICIDENCE+ 36 + (2 * (segmentID-18) + 1), backenergy);

			//! per segment plots
			plot(D_MTAS_SUM_FB + OUTER_OFFSET + (segmentID-18), segmentAvg); //! subtract 12 to put it back in segments 1 - 6 
			plot(DD_MTAS_IMO_T,totalSum.first,segmentAvg);
		}
	}

	//for(unsigned int ii = 0; ii < 6; ++ii){
	//	bool isSingleSegmentFire = true;
	//	for(unsigned int jj = 0; jj < 6; ++jj){
	//		if( ii == jj )
	//			continue;
	//		if( MtasSegVec.at(jj).IsValidSegment() ){
	//			isSingleSegmentFire = false;
	//		}
	//	}
	//	if( isSingleSegmentFire ){
	//		plot(D_SINGLE_CENTER_OFFSET + ii,MtasSegVec.at(ii).GetSegmentAverageEnergy());
	//	}
	//	if( ii == 0 and (not MtasSegVec.at(1).IsValidSegment()) and (not MtasSegVec.at(5).IsValidSegment()) )
	//		plot(D_SINGLE_CENTER_OFFSET + ii + 6,MtasSegVec.at(0).GetSegmentAverageEnergy());
	//	if( ii == 1 and (not MtasSegVec.at(2).IsValidSegment()) and (not MtasSegVec.at(0).IsValidSegment()) )
	//		plot(D_SINGLE_CENTER_OFFSET + ii + 6,MtasSegVec.at(1).GetSegmentAverageEnergy());
	//	if( ii == 2 and (not MtasSegVec.at(3).IsValidSegment()) and (not MtasSegVec.at(1).IsValidSegment()) )
	//		plot(D_SINGLE_CENTER_OFFSET + ii + 6,MtasSegVec.at(2).GetSegmentAverageEnergy());
	//	if( ii == 3 and (not MtasSegVec.at(4).IsValidSegment()) and (not MtasSegVec.at(2).IsValidSegment()) )
	//		plot(D_SINGLE_CENTER_OFFSET + ii + 6,MtasSegVec.at(3).GetSegmentAverageEnergy());
	//	if( ii == 4 and (not MtasSegVec.at(5).IsValidSegment()) and (not MtasSegVec.at(3).IsValidSegment()) )
	//		plot(D_SINGLE_CENTER_OFFSET + ii + 6,MtasSegVec.at(4).GetSegmentAverageEnergy());
	//	if( ii == 5 and (not MtasSegVec.at(0).IsValidSegment()) and (not MtasSegVec.at(4).IsValidSegment()) )
	//		plot(D_SINGLE_CENTER_OFFSET + ii + 6,MtasSegVec.at(5).GetSegmentAverageEnergy());
	//}

	//! sum spectra
	if( totalSum.second )
		plot(D_MTAS_TOTAL, totalSum.first);
	if( centerSum.second )
		plot(D_MTAS_CENTER, centerSum.first);
	if( innerSum.second )
		plot(D_MTAS_INNER, innerSum.first);
	if( middleSum.second )
		plot(D_MTAS_MIDDLE, middleSum.first);
	if( outerSum.second )
		plot(D_MTAS_OUTER, outerSum.first);

	if( totalSum.second and centerSum.second )
		plot(DD_MTAS_CR_T, totalSum.first, centerSum.first);

	//plot(DD_MTAS_C1_C2,MtasSegVec.at(1).GetSegmentAverageEnergy(),MtasSegVec.at(0).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C1_C6,MtasSegVec.at(5).GetSegmentAverageEnergy(),MtasSegVec.at(0).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C3_C2,MtasSegVec.at(1).GetSegmentAverageEnergy(),MtasSegVec.at(2).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C3_C4,MtasSegVec.at(3).GetSegmentAverageEnergy(),MtasSegVec.at(2).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C5_C4,MtasSegVec.at(3).GetSegmentAverageEnergy(),MtasSegVec.at(4).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C5_C6,MtasSegVec.at(5).GetSegmentAverageEnergy(),MtasSegVec.at(4).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C1_C3,MtasSegVec.at(2).GetSegmentAverageEnergy(),MtasSegVec.at(0).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C1_C4,MtasSegVec.at(3).GetSegmentAverageEnergy(),MtasSegVec.at(0).GetSegmentAverageEnergy());
	//plot(DD_MTAS_C1_C5,MtasSegVec.at(4).GetSegmentAverageEnergy(),MtasSegVec.at(0).GetSegmentAverageEnergy());

	//! Plot Multis 
	for (unsigned it = 0; it < MtasSegMulti.size(); ++it){
		plot(D_ONE_OFFS + 0, it, MtasSegMulti.at(it));
	}

	//place MTAS_Total so other classes have access to it
	EventData TotalData(EarliestTime,totalSum.first);
	TreeCorrelator::get()->place("MTAS_Total")->activate(TotalData);
	EventData CenterSum(EarliestTime,centerSum.first);
	TreeCorrelator::get()->place("MTAS_CenterSum")->activate(CenterSum);
	EventData InnerSum(EarliestTime,innerSum.first);
	TreeCorrelator::get()->place("MTAS_InnerSum")->activate(InnerSum);
	EventData MiddleSum(EarliestTime,middleSum.first);
	TreeCorrelator::get()->place("MTAS_MiddleSum")->activate(MiddleSum);
	EventData OutterSum(EarliestTime,outerSum.first);
	TreeCorrelator::get()->place("MTAS_OutterSum")->activate(OutterSum);
	for( int ii = 0; ii < 6; ++ii ){
		EventData CurrCenter(EarliestTime,MtasSegVec.at(ii).GetSegmentAverageEnergy().first);
		EventData CurrInner(EarliestTime,MtasSegVec.at(ii+6).GetSegmentAverageEnergy().first);
		EventData CurrMiddle(EarliestTime,MtasSegVec.at(ii+12).GetSegmentAverageEnergy().first);
		EventData CurrOutter(EarliestTime,MtasSegVec.at(ii+18).GetSegmentAverageEnergy().first);
	
		TreeCorrelator::get()->place("MTAS_C"+to_string(ii))->activate(CurrCenter);
		TreeCorrelator::get()->place("MTAS_I"+to_string(ii))->activate(CurrInner);
		TreeCorrelator::get()->place("MTAS_M"+to_string(ii))->activate(CurrMiddle);
		TreeCorrelator::get()->place("MTAS_O"+to_string(ii))->activate(CurrOutter);
	}

	SignalTime = EarliestTime;

	return true;
}

bool MtasProcessor::Process(RawEvent &event) {
	if (!EventProcessor::Process(event))
		return false;

	bool isBeta = false;
	if(TreeCorrelator::get()->checkPlace("BSM_Total")){
		//This is the BSM beta gated stuff
		if(TreeCorrelator::get()->place("BSM_Total")->status()){
			PlaceDetector* bsm_total = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("BSM_Total")); 
		
			if( bsm_total->info_.size() > 0 and bsm_total->last().energy >= BetaThreshold ){
				isBeta = true;
			}
			
			if( isBeta ){
				IsPrevBetaTriggered = true;
				PrevBetaTime = SignalTime;
				PrevBetaEnergy = bsm_total->last().energy;
				if( totalSum.second )
					plot(D_MTAS_BETA_TOTAL,totalSum.first);
				if( centerSum.second )
					plot(D_MTAS_BETA_CENTER, centerSum.first);
				if( innerSum.second )
					plot(D_MTAS_BETA_INNER, innerSum.first);
				if( middleSum.second )
					plot(D_MTAS_BETA_MIDDLE, middleSum.first);
				if( outerSum.second )
					plot(D_MTAS_BETA_OUTER, outerSum.first);
			}
		}else{
			if( IsPrevBetaTriggered ){
				IsPrevBetaTriggered = false;
				plot(D_ONE_OFFS + 1,totalSum.first);
				plot(D_ONE_OFFS + 2,totalSum.first,(SignalTime - PrevBetaTime)*clockInSeconds*1.0e6);
				plot(D_ONE_OFFS + 3,totalSum.first,PrevBetaEnergy);
			}
		}
	}

	EndProcess();
	return true;
}
