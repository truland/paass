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
#include <vector>

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
		const unsigned POS_OFFSET = 100;
		
		const unsigned DD_OFFSET = 250;

		const unsigned D_BSM_TOTAL = TOTAL_OFFSET; 
		const unsigned D_BSM_MTAS_SUM = TOTAL_OFFSET + 1;
		const unsigned D_BSM_ZERO_MTAS = TOTAL_OFFSET + 2;
		const unsigned D_BSM_MTAS_GATES = TOTAL_OFFSET + 10;

		const unsigned D_BSM_POSITION = POS_OFFSET;

		const unsigned DD_BSM_MTAS_TOTAL = DD_OFFSET; 
		const unsigned DD_BSM_TOTAL_POS = DD_OFFSET + 1;
		const unsigned DD_BSM_TOTAL_POS_ZERO_MTAS = DD_OFFSET + 2;
		const unsigned DD_BSM_TOTAL_POS_MTAS_GATES = DD_OFFSET + 10;
		const unsigned DD_BSM_F_B = DD_OFFSET + 20;
		const unsigned DD_BSM_F_B_MTAS_GATES = DD_OFFSET + 30;
	}
}
using namespace std;
using namespace dammIds::bsm;

void BSMProcessor::DeclarePlots(void){
	DeclareHistogram1D(D_BSM_TOTAL,SE,"BSM Total");
	DeclareHistogram1D(D_BSM_MTAS_SUM,SE,"BSM Total + MTAS Total");
	DeclareHistogram1D(D_BSM_ZERO_MTAS,SE,"BSM Total No MTAS");
	for( unsigned int ii = 0; ii < NumGates; ++ii ){
		string hisname = "BSM MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
		DeclareHistogram1D(D_BSM_MTAS_GATES+ii,SE,hisname.c_str());
		DeclareHistogram2D(D_BSM_MTAS_GATES+(NumGates + ii),SD,S4,("Cal "+hisname).c_str());
	}

	DeclareHistogram1D(D_BSM_POSITION,SD,"BSM Position");

	DeclareHistogram2D(DD_BSM_MTAS_TOTAL,SC,SC,"BSM Total vs MTAS Total");
	DeclareHistogram2D(DD_BSM_F_B,SC,SC,"BSM Front Avg vs BSM Back Avg");
	DeclareHistogram2D(DD_BSM_F_B + 1,SD,SD,"BSM Avg vs BSM Sqrt"); 
	DeclareHistogram2D(DD_BSM_TOTAL_POS,SD,SC,"BSM Energy vs Position");
	DeclareHistogram2D(DD_BSM_TOTAL_POS_ZERO_MTAS,SD,SC,"BSM Energy vs Position No MTAS");
	for( unsigned int ii = 0; ii < NumGates; ++ii ){
		string hisname = "BSM Energy vs Position MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
		DeclareHistogram2D(DD_BSM_TOTAL_POS_MTAS_GATES+ii,SD,SC,hisname.c_str());
		hisname = "BSM Front vs Back MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
		DeclareHistogram2D(DD_BSM_F_B_MTAS_GATES+ii,SC,SC,hisname.c_str());
	}

}


BSMProcessor::BSMProcessor(int numsegments,bool zerosuppress,bool alone,vector<pair<double,double>> mtasgates,double thresh) : EventProcessor(OFFSET, RANGE, "BSMProcessor") {
	associatedTypes.insert("bsm");
	PixieRev = Globals::get()->GetPixieRevision();
	NumSegments = numsegments;
	HasZeroSuppression = zerosuppress;
	StandAlone = alone;
	MTASGates = mtasgates;
	if( MTASGates.size() > MaxGates )
		cout << "Cannot have more than " << MaxGates << " MTAS gates for the BSM. Only the first " << MaxGates << " will be used" << endl;
	NumGates = (MTASGates.size() <= MaxGates) ? MTASGates.size() : MaxGates;
	cout << "Using " << NumGates << " MTAS Gates for the BSM" << endl;
	for( unsigned int ii = 0; ii < NumGates; ++ii )
		cout << " Gate " << ii << " : [" << MTASGates.at(ii).first << "," << MTASGates.at(ii).second << "]" << endl;
	Threshold = thresh;
}

bool BSMProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::PreProcess(event))
		return false;

	static const auto &chanEvents = event.GetSummary("bsm", true)->GetList();
	//vector<BSMSegment> BSMSegVec(NumSegments,BSMSegment(HasZeroSuppression));
	BSMSegVec = vector<BSMSegment>(NumSegments,BSMSegment(HasZeroSuppression));
	vector<short> BSMSegMulti(2*NumSegments,0); // MTAS segment multiplicity "map"

	double EarliestTime = 1.0e99;
	double clockInSeconds;
	for (auto chanEvtIter = chanEvents.begin(); chanEvtIter != chanEvents.end(); ++chanEvtIter){
		int segmentNum = stoi((*chanEvtIter)->GetChanID().GetGroup().c_str()) - 1;

		if (PixieRev == "F"){
			clockInSeconds = Globals::get()->GetClockInSeconds((*chanEvtIter)->GetChanID().GetModFreq());
		} else {
			clockInSeconds = Globals::get()->GetClockInSeconds();
		}

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
			if( (*chanEvtIter)->GetCalibratedEnergy() > Threshold ){
				BSMSegMulti.at(GlobalChanID)++; 

				BSMSegVec.at(segmentNum).gBSMSegID_ = segmentNum;
				if(isFront && BSMSegVec.at(segmentNum).segFront_ == nullptr){  
					if( (*chanEvtIter)->GetTimeSansCfd() < EarliestTime )
						EarliestTime = (*chanEvtIter)->GetTimeSansCfd(); 
					BSMSegVec.at(segmentNum).segFront_ = (*chanEvtIter);
					BSMSegVec.at(segmentNum).PixieRev = PixieRev;
				}else if (isBack && BSMSegVec.at(segmentNum).segBack_ == nullptr) { 
					if( (*chanEvtIter)->GetTimeSansCfd() < EarliestTime )
						EarliestTime = (*chanEvtIter)->GetTimeSansCfd(); 
					BSMSegVec.at(segmentNum).segBack_ = (*chanEvtIter);
					BSMSegVec.at(segmentNum).PixieRev = PixieRev;
				}
			}
		}
	}

	//reset this during pre-process
	BSMTotal = make_pair(0,not HasZeroSuppression);
	FrontAvg = make_pair(0,not HasZeroSuppression);
	BackAvg = make_pair(0,not HasZeroSuppression);
	int NumFire = 0;
	for( auto& segIter : BSMSegVec ){
		if( segIter.IsValidSegment() )
			++NumFire;
	}
	
	for( auto& segIter : BSMSegVec ){
		auto result = segIter.GetSegmentAverageEnergy();
		auto frontresult = segIter.GetFrontEnergy();
		auto backresult = segIter.GetBackEnergy();
		if( not result.second ){
			continue;
		}else{
			BSMTotal.first += result.first/static_cast<double>(NumFire);
			BSMTotal.second = true;	
			FrontAvg.first += frontresult.first/static_cast<double>(NumFire);
			FrontAvg.second = true;
			BackAvg.first += backresult.first/static_cast<double>(NumFire);
			BackAvg.second = true;
		}
	}

	if( BSMTotal.second )
		plot(D_BSM_TOTAL,BSMTotal.first);
	if( FrontAvg.second or BackAvg.second ){
		plot(DD_BSM_F_B,FrontAvg.first,BackAvg.first);
		plot(DD_BSM_F_B + 1,BSMTotal.first,sqrt(FrontAvg.first*BackAvg.first));
		BSMPosition = (SD/2)*(1.0+((FrontAvg.first-BackAvg.first)/(FrontAvg.first+BackAvg.first)));
		plot(D_BSM_POSITION,(SD/2)*(1.0+((FrontAvg.first-BackAvg.first)/(FrontAvg.first+BackAvg.first))));
		plot(DD_BSM_TOTAL_POS,(SD/2)*(1.0+((FrontAvg.first-BackAvg.first)/(FrontAvg.first+BackAvg.first))),BSMTotal.first);
	}

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
				double MTASTotal = mtas_total->last().energy;
				plot(DD_BSM_MTAS_TOTAL,MTASTotal,BSMTotal.first);
				plot(D_BSM_MTAS_SUM,MTASTotal + BSMTotal.first);
				for( unsigned int ii = 0; ii < NumGates; ++ii ){
					if( MTASTotal >= MTASGates.at(ii).first and MTASTotal <= MTASGates.at(ii).second ){
						plot(D_BSM_MTAS_GATES+ii,BSMTotal.first);
						plot(DD_BSM_TOTAL_POS_MTAS_GATES+ii,BSMPosition,BSMTotal.first);
						for( size_t jj = 0; jj < BSMSegVec.size(); ++jj ){
							plot(D_BSM_MTAS_GATES+(NumGates + ii ),BSMSegVec.at(jj).GetFrontEnergy().first,2*jj);
							plot(D_BSM_MTAS_GATES+(NumGates + ii ),BSMSegVec.at(jj).GetBackEnergy().first,2*jj + 1);
						}
						plot(DD_BSM_F_B_MTAS_GATES+ii,FrontAvg.first,BackAvg.first);
					}	
				}
			}
		}else{
			if( BSMTotal.second ){
				plot(D_BSM_ZERO_MTAS,BSMTotal.first);
				plot(DD_BSM_TOTAL_POS_ZERO_MTAS,BSMPosition,BSMTotal.first);
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
