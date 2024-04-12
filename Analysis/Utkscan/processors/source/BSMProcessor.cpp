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
		const unsigned DD_COMPRESSED_OFFSET = 400;

		const unsigned D_BSM_TOTAL = TOTAL_OFFSET; 
		const unsigned D_BSM_MTAS_SUM = TOTAL_OFFSET + 1;
		const unsigned D_BSM_ZERO_MTAS = TOTAL_OFFSET + 2;
		const unsigned D_BSM_SATPIL_MTAS = TOTAL_OFFSET + 3;
        //sum for yap
		const unsigned D_BSM_SUM = TOTAL_OFFSET + 4;
		const unsigned D_BSM_SUM_MTAS_SUM = TOTAL_OFFSET + 5;
		const unsigned D_BSM_SUM_ZERO_MTAS = TOTAL_OFFSET + 6;
		const unsigned D_BSM_SUM_SATPIL_MTAS = TOTAL_OFFSET + 7;
		
        const unsigned D_BSM_FRONT = TOTAL_OFFSET + 8;
		const unsigned D_BSM_FRONT_MTAS_SUM = TOTAL_OFFSET + 9;
		const unsigned D_BSM_FRONT_ZERO_MTAS = TOTAL_OFFSET + 10;
		const unsigned D_BSM_FRONT_SATPIL_MTAS = TOTAL_OFFSET + 11;

		const unsigned D_BSM_BACK = TOTAL_OFFSET + 12;
		const unsigned D_BSM_BACK_MTAS_SUM = TOTAL_OFFSET + 13;
		const unsigned D_BSM_BACK_ZERO_MTAS = TOTAL_OFFSET + 14;
		const unsigned D_BSM_BACK_SATPIL_MTAS = TOTAL_OFFSET + 15;

		const unsigned D_BSM_MTAS_GATES = TOTAL_OFFSET + 20;

		const unsigned D_BSM_POSITION = POS_OFFSET;

		const unsigned DD_BSM_MTAS_TOTAL = DD_OFFSET; 
		const unsigned DD_BSM_TOTAL_POS = DD_OFFSET + 1;
		const unsigned DD_BSM_TOTAL_POS_ZERO_MTAS = DD_OFFSET + 2;
		const unsigned DD_BSM_MTAS_GATES = DD_OFFSET +3;
		const unsigned DD_BSM_TOTAL_POS_MTAS_GATES = DD_OFFSET + 10;
		const unsigned DD_BSM_F_B = DD_OFFSET + 20;
		const unsigned DD_BSM_TOTAL_AVG = DD_OFFSET + 21;
		const unsigned DD_BSM_F_B_MTAS_GATES = DD_OFFSET + 30;
		const unsigned DD_BSM_F_POS_MTAS_GATES = DD_OFFSET + 40;
		const unsigned DD_BSM_B_POS_MTAS_GATES = DD_OFFSET + 50;
		const unsigned DD_BSM_MTAS_CENTER_SUM = DD_OFFSET + 60;
		const unsigned DD_BSM_MTAS_CENTER_INDIVIDUAL = DD_OFFSET + 61;
		const unsigned DD_BSM_MTAS_CENTER_SUM_VETO_MIDDLE_OUTTER = DD_OFFSET + 62;
		const unsigned DD_BSM_MTAS_CENTER_INDIVIDUAL_VETO_MIDDLE_OUTTER = DD_OFFSET + 63;
		const unsigned DD_BSM_MTAS_CENTER_SUM_INNER_SUM_VETO_MIDDLE_OUTTER = DD_OFFSET + 64;
		const unsigned DD_BSM_MTAS_CENTER_SUM_VETO_CENTER_NEIGHBORS = DD_OFFSET + 65;
		const unsigned DD_BSM_MTAS_CENTER_INDIVIDUAL_VETO_CENTER_NEIGHBORS = DD_OFFSET + 66;
		const unsigned DD_BSM_MTAS_INNER_SUM = DD_OFFSET + 70;
		const unsigned DD_BSM_MTAS_INNER_INDIVIDUAL = DD_OFFSET + 71;
		const unsigned DD_BSM_MTAS_INNER_SUM_VETO_MIDDLE_OUTTER = DD_OFFSET + 72;
		const unsigned DD_BSM_MTAS_INNER_INDIVIDUAL_VETO_MIDDLE_OUTTER = DD_OFFSET + 73;
		const unsigned DD_BSM_MTAS_MIDDLE_SUM = DD_OFFSET + 80;
		const unsigned DD_BSM_MTAS_MIDDLE_INDIVIDUAL = DD_OFFSET + 81;
		const unsigned DD_BSM_MTAS_OUTTER_SUM = DD_OFFSET + 90;
		const unsigned DD_BSM_MTAS_OUTTER_INDIVIDUAL= DD_OFFSET + 91;

		const unsigned DD_BSM_MTAS_TOTAL_COMPRESSED = DD_COMPRESSED_OFFSET;
	}
}
using namespace std;
using namespace dammIds::bsm;

void BSMProcessor::DeclarePlots(void){
	DeclareHistogram1D(D_BSM_TOTAL,SE,"BSM Total");
	DeclareHistogram1D(D_BSM_MTAS_SUM,SE,"BSM Total + MTAS Total");
	DeclareHistogram1D(D_BSM_ZERO_MTAS,SE,"BSM Total No MTAS");
	DeclareHistogram1D(D_BSM_SATPIL_MTAS,SE,"BSM Total MTAS Saturate");
	
    DeclareHistogram1D(D_BSM_SUM,SE,"BSM Sum");
	DeclareHistogram1D(D_BSM_SUM_MTAS_SUM,SE,"BSM Sum + MTAS Total");
	DeclareHistogram1D(D_BSM_SUM_ZERO_MTAS,SE,"BSM Sum No MTAS");
	DeclareHistogram1D(D_BSM_SUM_SATPIL_MTAS,SE,"BSM Sum MTAS Saturate");
    
    DeclareHistogram1D(D_BSM_FRONT,SE,"BSM Front Only");
	DeclareHistogram1D(D_BSM_FRONT_MTAS_SUM,SE,"BSM Front Only + MTAS Total");
	DeclareHistogram1D(D_BSM_FRONT_ZERO_MTAS,SE,"BSM Front Only No MTAS");
	DeclareHistogram1D(D_BSM_FRONT_SATPIL_MTAS,SE,"BSM Front Only MTAS Saturate");

	DeclareHistogram1D(D_BSM_BACK,SE,"BSM Back Only");
	DeclareHistogram1D(D_BSM_BACK_MTAS_SUM,SE,"BSM Back Only + MTAS Total");
	DeclareHistogram1D(D_BSM_BACK_ZERO_MTAS,SE,"BSM Back Only No MTAS");
	DeclareHistogram1D(D_BSM_BACK_SATPIL_MTAS,SE,"BSM Back Only MTAS Saturate");


	for( unsigned int ii = 0; ii < NumGates; ++ii ){
		string hisname = "BSM MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
		DeclareHistogram1D(D_BSM_MTAS_GATES+ii,SE,hisname.c_str());
		//DeclareHistogram2D(D_BSM_MTAS_GATES+(NumGates + ii),PositionSize,S4,("Cal "+hisname).c_str());
	}

	DeclareHistogram1D(D_BSM_POSITION,PositionSize,"BSM Position");

	DeclareHistogram2D(DD_BSM_MTAS_GATES,S4,SC,"MTAS Gates");
	DeclareHistogram2D(DD_BSM_MTAS_TOTAL,SC,SC,"BSM Total vs MTAS Total");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_SUM,SC,SC,"BSM Total vs MTAS Center Sum");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_INDIVIDUAL,SC,SC,"BSM Total vs MTAS Center Individual");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_SUM_VETO_MIDDLE_OUTTER,SC,SC,"BSM Total vs MTAS Center Sum veto M+O");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_INDIVIDUAL_VETO_MIDDLE_OUTTER,SC,SC,"BSM Total vs MTAS Center Individual veto M+O");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_SUM_INNER_SUM_VETO_MIDDLE_OUTTER,SC,SC,"BSM Total vs MTAS Center+Inner Sum veto M+O");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_SUM_VETO_CENTER_NEIGHBORS,SC,SC,"BSM Total vs MTAS Center Sum veto C Neighbors");
	DeclareHistogram2D(DD_BSM_MTAS_CENTER_INDIVIDUAL_VETO_CENTER_NEIGHBORS,SC,SC,"BSM Total vs MTAS Center Individual veto C Neighbors");
	DeclareHistogram2D(DD_BSM_MTAS_INNER_SUM_VETO_MIDDLE_OUTTER,SC,SC,"BSM Total vs MTAS Inner Sum veto M+O");
	DeclareHistogram2D(DD_BSM_MTAS_INNER_INDIVIDUAL_VETO_MIDDLE_OUTTER,SC,SC,"BSM Total vs MTAS Inner Individual veto M+O");
	DeclareHistogram2D(DD_BSM_MTAS_INNER_SUM,SC,SC,"BSM Total vs MTAS Inner Sum");
	DeclareHistogram2D(DD_BSM_MTAS_INNER_INDIVIDUAL,SC,SC,"BSM Total vs MTAS Inner Individual");
	DeclareHistogram2D(DD_BSM_MTAS_MIDDLE_SUM,SC,SC,"BSM Total vs MTAS Middle Sum");
	DeclareHistogram2D(DD_BSM_MTAS_MIDDLE_INDIVIDUAL,SC,SC,"BSM Total vs MTAS Middle Individual");
	DeclareHistogram2D(DD_BSM_MTAS_OUTTER_SUM,SC,SC,"BSM Total vs MTAS Outter Sum");
	DeclareHistogram2D(DD_BSM_MTAS_OUTTER_INDIVIDUAL,SC,SC,"BSM Total vs MTAS Outter Individual");
	DeclareHistogram2D(DD_BSM_F_B,SC,SC,"BSM Front Avg vs BSM Back Avg");
	//DeclareHistogram2D(DD_BSM_TOTAL_AVG,SD,SD,"BSM Avg vs BSM Sqrt"); 
	//DeclareHistogram2D(DD_BSM_TOTAL_POS,PositionSize,SC,"BSM Energy vs Position");
	//DeclareHistogram2D(DD_BSM_TOTAL_POS_ZERO_MTAS,PositionSize,SC,"BSM Energy vs Position No MTAS");
	//for( unsigned int ii = 0; ii < NumGates; ++ii ){
	//	string hisname = "BSM Energy vs Position MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
	//	DeclareHistogram2D(DD_BSM_TOTAL_POS_MTAS_GATES+ii,PositionSize,SC,hisname.c_str());
	//	hisname = "BSM Front vs Back MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
	//	DeclareHistogram2D(DD_BSM_F_B_MTAS_GATES+ii,SC,SC,hisname.c_str());
	//	hisname = "BSM Front Energy vs Position MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
	//	DeclareHistogram2D(DD_BSM_F_POS_MTAS_GATES+ii,PositionSize,SC,hisname.c_str());
	//	hisname = "BSM Back Energy vs Position MTAS T["+to_string((int)MTASGates.at(ii).first)+","+to_string((int)MTASGates.at(ii).second)+"]";
	//	DeclareHistogram2D(DD_BSM_B_POS_MTAS_GATES+ii,PositionSize,SC,hisname.c_str());
	//}

	DeclareHistogram2D(DD_BSM_MTAS_TOTAL_COMPRESSED,SB,SB,"BSM Total/10 vs MTAS Total/10");
}


BSMProcessor::BSMProcessor(int numsegments,bool alone,vector<pair<double,double>> mtasgates,double thresh,double fc,double fs,double fm,double bc,double bs,double bm) : EventProcessor(OFFSET, RANGE, "BSMProcessor") {
	associatedTypes.insert("bsm");
	PixieRev = Globals::get()->GetPixieRevision();
	NumSegments = numsegments;
	StandAlone = alone;
	MTASGates = mtasgates;
	if( MTASGates.size() > MaxGates )
		cout << "Cannot have more than " << MaxGates << " MTAS gates for the BSM. Only the first " << MaxGates << " will be used" << endl;
	NumGates = (MTASGates.size() <= MaxGates) ? MTASGates.size() : MaxGates;
	cout << "       * Using " << NumGates << " MTAS Gates for the BSM" << endl;
	for( unsigned int ii = 0; ii < NumGates; ++ii )
		cout << "       * Gate " << ii << " : [" << MTASGates.at(ii).first << "," << MTASGates.at(ii).second << "]" << endl;
	Threshold = thresh;
	FrontCorrection.constant = fc;
	FrontCorrection.slope = fs;
	FrontCorrection.mean = fm;
	BackCorrection.constant = bc;
	BackCorrection.slope = bs;
	BackCorrection.mean = bm;
	PositionSize = SA;
}
		
bool BSMProcessor::PreProcess(RawEvent &event) {
	//cout << "INSIDE BSMProcessor::PreProcess" << endl;
	if (!EventProcessor::PreProcess(event))
		return false;

	const auto &chanEvents = event.GetSummary("bsm", true)->GetList();
	BSMSegVec = vector<BSMSegment>(NumSegments,BSMSegment());
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
	BSMTotal = 0;
	FrontAvg = 0;
	BackAvg = 0;
    BSMSum = 0;
    FrontSolo = 0;
    BackSolo = 0;
	int NumFire = 0;
	for( auto& segIter : BSMSegVec ){
        BSMSum = (segIter.GetFrontEnergy() + segIter.GetBackEnergy());
		if( segIter.IsValidSegment() ){
			++NumFire;
        }else{
            FrontSolo = segIter.GetFrontEnergy();
            BackSolo = segIter.GetBackEnergy();
            if( FrontSolo <= 0.0 ){
                plot(D_BSM_BACK,BackSolo);
            }else{
                plot(D_BSM_FRONT,FrontSolo);
            }       
        }
	}
    plot(D_BSM_SUM,BSMSum);
	if( NumFire > 0 ){
		for( auto& segIter : BSMSegVec ){
			auto result = segIter.GetSegmentAverageEnergy();
			auto frontresult = segIter.GetFrontEnergy();
			auto backresult = segIter.GetBackEnergy();
			//BSMTotal.first += result.first/static_cast<double>(NumFire);
			//BSMTotal.second = true;	
			FrontAvg += frontresult/static_cast<double>(NumFire);
			BackAvg += backresult/static_cast<double>(NumFire);

			//! Begin Root Output stuff.
			if (DetectorDriver::get()->GetSysRootOutput() && segIter.IsValidSegment()) {
				Bsmstruct.energy = segIter.GetSegmentAverageEnergy();
				Bsmstruct.fEnergy = segIter.GetFrontEnergy();
				Bsmstruct.bEnergy = segIter.GetBackEnergy();
				Bsmstruct.time = (segIter.GetFrontTimeInNS() + segIter.GetBackTimeInNS())/2.0;
				Bsmstruct.tdiff = (segIter.GetFrontTimeInNS() - segIter.GetBackTimeInNS());
				Bsmstruct.gSegmentID = segIter.gBSMSegID_;
				Bsmstruct.segmentNum = (segIter.gBSMSegID_+1)/2;
				Bsmstruct.ftrace = segIter.segFront_->GetTrace();
				Bsmstruct.btrace = segIter.segBack_->GetTrace();
				
				pixie_tree_event_->bsm_vec_.emplace_back(Bsmstruct);
				Bsmstruct = processor_struct::BSM_DEFAULT_STRUCT;
			}


		}
		BSMPosition = (PositionSize/2)*(1.0+((FrontAvg-BackAvg)/(FrontAvg+BackAvg)));
		FrontAvg = FrontCorrection.Correct(FrontAvg,BSMPosition);
		BackAvg = BackCorrection.Correct(BackAvg,BSMPosition);
		BSMTotal = (FrontAvg + BackAvg)/static_cast<double>(2*NumFire);

		plot(D_BSM_TOTAL,BSMTotal);
		plot(DD_BSM_F_B,FrontAvg,BackAvg);
		//plot(DD_BSM_TOTAL_AVG,BSMTotal,sqrt(FrontAvg*BackAvg));
		//BSMPosition = (PositionSize/2)*(1.0+((FrontAvg-BackAvg)/(FrontAvg+BackAvg)));
		plot(D_BSM_POSITION,(PositionSize/2)*(1.0+((FrontAvg-BackAvg)/(FrontAvg+BackAvg))));
		//plot(DD_BSM_TOTAL_POS,(PositionSize/2)*(1.0+((FrontAvg-BackAvg)/(FrontAvg+BackAvg))),BSMTotal);
	}

	EventData TotalData(EarliestTime,BSMTotal);
	TreeCorrelator::get()->place("BSM_Total")->activate(TotalData);

	//cout << "LEAVE BSMProcessor::PreProcess" << endl;
        EndProcess();
	return true;
}

bool BSMProcessor::Process(RawEvent &event) {
	//cout << "INSIDE BSMProcessor::Process" << endl;
	if (!EventProcessor::Process(event))
		return false;

	bool IsValid = true;
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_Total");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_CenterSum");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_InnerSum");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_MiddleSum");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_OutterSum");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_C0");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_C1");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_C2");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_C3");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_C4");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_C5");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_I0");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_I1");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_I2");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_I3");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_I4");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_I5");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_M0");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_M1");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_M2");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_M3");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_M4");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_M5");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_O0");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_O1");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_O2");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_O3");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_O4");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_O5");
	IsValid &= TreeCorrelator::get()->checkPlace("MTAS_SaturatePileup");

	if(IsValid){
		double MTASTotal = 0.0;
		double MTASCenter = 0.0;
		double MTASInner = 0.0;
		double MTASMiddle = 0.0;
		double MTASOutter = 0.0;
		vector<double> CenterSegment(6,0.0);
		vector<double> InnerSegment(6,0.0);
		vector<double> MiddleSegment(6,0.0);
		vector<double> OutterSegment(6,0.0);
		double MTASSatPile = -1.0;

		PlaceDetector* mtas_saturatepileup = nullptr;
		if(TreeCorrelator::get()->place("MTAS_SaturatePileup")->status()){
			mtas_saturatepileup = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_SaturatePileup")); 
			if( mtas_saturatepileup->info_.size() > 0 )
				MTASSatPile = mtas_saturatepileup->last().energy;
		}

		if( MTASSatPile < 0 ){
			PlaceDetector* mtas_total = nullptr;
			if(TreeCorrelator::get()->place("MTAS_Total")->status()){
				mtas_total = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_Total")); 
				if( mtas_total->info_.size() > 0 )
					MTASTotal = mtas_total->last().energy;
			}

			PlaceDetector* mtas_center = nullptr;
			if(TreeCorrelator::get()->place("MTAS_CenterSum")->status()){
				mtas_center = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_CenterSum")); 
				if( mtas_center->info_.size() > 0 )
					MTASCenter = mtas_center->last().energy;
			}

			PlaceDetector* mtas_inner = nullptr;
			if(TreeCorrelator::get()->place("MTAS_InnerSum")->status()){
				mtas_inner = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_InnerSum")); 
				if( mtas_inner->info_.size() > 0 )
					MTASInner = mtas_inner->last().energy;
			}

			PlaceDetector* mtas_middle = nullptr;
			if(TreeCorrelator::get()->place("MTAS_MiddleSum")->status()){
				mtas_middle = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_MiddleSum")); 
				if( mtas_middle->info_.size() > 0 )
					MTASMiddle = mtas_middle->last().energy;
			}

			PlaceDetector* mtas_outter = nullptr;
			if(TreeCorrelator::get()->place("MTAS_OutterSum")->status()){
				mtas_outter = dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_OutterSum")); 
				if( mtas_outter->info_.size() > 0 )
					MTASOutter = mtas_outter->last().energy;
			}

			vector<PlaceDetector*> mtas_center_individual(6,nullptr);
			vector<PlaceDetector*> mtas_inner_individual(6,nullptr);
			vector<PlaceDetector*> mtas_middle_individual(6,nullptr);
			vector<PlaceDetector*> mtas_outter_individual(6,nullptr);
			for( int ii = 0; ii < 6; ++ii ){
				if(TreeCorrelator::get()->place("MTAS_C"+to_string(ii))->status()){
					mtas_center_individual.at(ii) =  dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_C"+to_string(ii)));
					if( mtas_center_individual.at(ii)->info_.size() > 0 )
						CenterSegment.at(ii) = mtas_center_individual.at(ii)->last().energy;
				}

				if(TreeCorrelator::get()->place("MTAS_I"+to_string(ii))->status()){
					mtas_inner_individual.at(ii) =  dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_I"+to_string(ii)));
					if( mtas_inner_individual.at(ii)->info_.size() > 0 )
						InnerSegment.at(ii) = mtas_inner_individual.at(ii)->last().energy;
				}

				if(TreeCorrelator::get()->place("MTAS_M"+to_string(ii))->status()){
					mtas_middle_individual.at(ii) =  dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_M"+to_string(ii)));
					if( mtas_middle_individual.at(ii)->info_.size() > 0 )
						MiddleSegment.at(ii) = mtas_middle_individual.at(ii)->last().energy;
				}

				if(TreeCorrelator::get()->place("MTAS_O"+to_string(ii))->status()){
					mtas_outter_individual.at(ii) =  dynamic_cast<PlaceDetector*>(TreeCorrelator::get()->place("MTAS_O"+to_string(ii)));
					if( mtas_outter_individual.at(ii)->info_.size() > 0 )
						OutterSegment.at(ii) = mtas_outter_individual.at(ii)->last().energy;
				}

			}
		}else{
			plot(D_BSM_SATPIL_MTAS,BSMTotal);
            plot(D_BSM_SUM_SATPIL_MTAS,BSMSum);
            plot(D_BSM_FRONT_SATPIL_MTAS,FrontSolo);
            plot(D_BSM_BACK_SATPIL_MTAS,BackSolo);
		}


		bool NeighborsFire = false;
		for( int ii = 0; ii < 6; ++ii ){
			NeighborsFire |= ((CenterSegment.at(ii) > 0.0) & ((CenterSegment.at((ii+1)%6) > 0.0) | (CenterSegment.at((5+ii)%6) > 0.0 )) );
		}
		//this catches everything that is multiplicity 3 and less in center
		if( !NeighborsFire ){
			plot(DD_BSM_MTAS_CENTER_SUM_VETO_CENTER_NEIGHBORS,MTASCenter,BSMTotal);
			for( int ii = 0; ii < 6; ++ii ){
				plot(DD_BSM_MTAS_CENTER_INDIVIDUAL_VETO_CENTER_NEIGHBORS,CenterSegment.at(ii),BSMTotal);
			}
		}

		//these plots can fill with real zeros and fake zeros(events that triggered pixie but were below threshold/only one pmt in a segment)
		plot(DD_BSM_MTAS_CENTER_SUM,MTASCenter,BSMTotal);
		plot(DD_BSM_MTAS_INNER_SUM,MTASInner,BSMTotal);
		plot(DD_BSM_MTAS_MIDDLE_SUM,MTASMiddle,BSMTotal);
		plot(DD_BSM_MTAS_OUTTER_SUM,MTASOutter,BSMTotal);
		//these will fill with many zeros since many of the pmts wont actually fire
		//will need to talk with rasshole if this is actually what we want to do
		for( int ii = 0; ii < 6; ++ii ){
			plot( DD_BSM_MTAS_CENTER_INDIVIDUAL, CenterSegment.at(ii),BSMTotal);
			plot( DD_BSM_MTAS_INNER_INDIVIDUAL, InnerSegment.at(ii),BSMTotal);
			plot( DD_BSM_MTAS_MIDDLE_INDIVIDUAL, MiddleSegment.at(ii),BSMTotal);
			plot( DD_BSM_MTAS_OUTTER_INDIVIDUAL, OutterSegment.at(ii),BSMTotal);
		}

		plot(DD_BSM_MTAS_TOTAL,MTASTotal,BSMTotal);
		plot(DD_BSM_MTAS_TOTAL_COMPRESSED,MTASTotal/10.0,BSMTotal/10.0);
		plot(D_BSM_MTAS_SUM,MTASTotal + BSMTotal);
        plot(D_BSM_SUM_MTAS_SUM,MTASTotal + BSMSum);
        plot(D_BSM_FRONT_MTAS_SUM,MTASTotal + FrontSolo);
        plot(D_BSM_BACK_MTAS_SUM,MTASTotal + BackSolo);
		for( unsigned int ii = 0; ii < NumGates; ++ii ){
			//plot(DD_BSM_MTAS_GATES,ii,MTASGates.at(ii).first);
			//plot(DD_BSM_MTAS_GATES,ii,MTASGates.at(ii).second);
			//if( MTASTotal >= MTASGates.at(ii).first and MTASTotal <= MTASGates.at(ii).second ){
			//	plot(D_BSM_MTAS_GATES+ii,BSMTotal);
			//	plot(DD_BSM_TOTAL_POS_MTAS_GATES+ii,BSMPosition,BSMTotal);
			//	plot(DD_BSM_F_POS_MTAS_GATES+ii,BSMPosition,FrontAvg);
			//	plot(DD_BSM_B_POS_MTAS_GATES+ii,BSMPosition,BackAvg);
			//	for( size_t jj = 0; jj < BSMSegVec.size(); ++jj ){
			//		plot(D_BSM_MTAS_GATES+(NumGates + ii ),BSMSegVec.at(jj).GetFrontEnergy(),2*jj);
			//		plot(D_BSM_MTAS_GATES+(NumGates + ii ),BSMSegVec.at(jj).GetBackEnergy(),2*jj + 1);
			//	}
			//	plot(DD_BSM_F_B_MTAS_GATES+ii,FrontAvg,BackAvg);
			//}	
		}

		if( (MTASMiddle+MTASOutter) <= 0.0 ){
			plot(DD_BSM_MTAS_CENTER_SUM_VETO_MIDDLE_OUTTER,MTASCenter,BSMTotal);
			plot(DD_BSM_MTAS_CENTER_SUM_INNER_SUM_VETO_MIDDLE_OUTTER,MTASCenter+MTASInner,BSMTotal);
			plot(DD_BSM_MTAS_INNER_SUM_VETO_MIDDLE_OUTTER,MTASInner,BSMTotal);
			for( int ii = 0; ii < 6; ++ii ){
				plot(DD_BSM_MTAS_CENTER_INDIVIDUAL_VETO_MIDDLE_OUTTER,CenterSegment.at(ii),BSMTotal);
				plot(DD_BSM_MTAS_INNER_INDIVIDUAL_VETO_MIDDLE_OUTTER,InnerSegment.at(ii),BSMTotal);
			}
		}

		//nothing fired inside MTAS
		//these fill only when none of the MTAS pmts fire in the event otherwise mtas_total will exist as an object
		if( MTASTotal <= 0.0 and MTASSatPile < 0 ){
			plot(D_BSM_ZERO_MTAS,BSMTotal);
            plot(D_BSM_SUM_ZERO_MTAS,BSMSum);
            plot(D_BSM_FRONT_ZERO_MTAS,FrontSolo);
            plot(D_BSM_BACK_ZERO_MTAS,BackSolo);
			//plot(DD_BSM_TOTAL_POS_ZERO_MTAS,BSMPosition,BSMTotal);
		}
	}else{
		if( not StandAlone )
			throw "ERROR BSMProcessor::Process BSM isn't in standalone mode and the treecorrelator for MTAS_Total and others isn't setup";
	}

	//cout << "LEAVING BSMProcessor::Process" << endl;
	EndProcess();
	//cout << "LEAVING BSMProcessor::Process" << endl;
	return true;
}
