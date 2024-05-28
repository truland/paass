/** @file PMTASProcessor.cpp
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
#include <stdexcept>
#include <string>

#include "Globals.hpp"
#include "DetectorDriver.hpp"
#include "DetectorLibrary.hpp"
#include "HelperFunctions.hpp"
#include "RawEvent.hpp"
#include "PMTASProcessor.hpp"

using namespace std;

PMTASProcessor::PMTASProcessor() : EventProcessor() {
	associatedTypes.insert("PMTAS");
	
	Rev = Globals::get()->GetPixieRevision();
	std::string name = Globals::get()->GetOutputPath() + Globals::get()->GetOutputFileName() + "_PMTAS.root";

	Zero = std::vector<double>(6,0.0);
	IndZero = std::vector<double>(12,0.0);
	ZeroHits = std::vector<int>(12,0);
	ZeroBetaHits = std::vector<int>(2,0);
	ZeroTrace = {};
	ZeroBeta = std::vector<double>(2,0.0);

	Reset();

	outputfile = new TFile(name.c_str(),"RECREATE");
	outputtree = new TTree("Event","Output from PMTASProcess in utkscan. For ORNL Pixie setup");
	
	outputtree->Branch("BetaSumFrontBack",&BetaSumFrontBack);
	outputtree->Branch("BetaAverage",&BetaAverage);

	outputtree->Branch("MTASTotalSum",&MTASTotalSum);
	outputtree->Branch("MTASCenterSum",&MTASCenterSum);
	outputtree->Branch("MTASInnerSum",&MTASInnerSum);
	outputtree->Branch("MTASMiddleSum",&MTASMiddleSum);
	outputtree->Branch("MTASOutterSum",&MTASOutterSum);

	outputtree->Branch("MTASCenterSumFrontBack",&MTASCenterSumFrontBack);
	outputtree->Branch("MTASInnerSumFrontBack",&MTASInnerSumFrontBack);
	outputtree->Branch("MTASMiddleSumFrontBack",&MTASMiddleSumFrontBack);
	outputtree->Branch("MTASOutterSumFrontBack",&MTASOutterSumFrontBack);

	outputtree->Branch("MTASCenter",&MTASCenter);
	outputtree->Branch("MTASInner",&MTASInner);
	outputtree->Branch("MTASMiddle",&MTASMiddle);
	outputtree->Branch("MTASOutter",&MTASOutter);

	outputtree->Branch("MTASCenterTS",&MTASCenterTS);
	outputtree->Branch("MTASInnerTS",&MTASInnerTS);
	outputtree->Branch("MTASMiddleTS",&MTASMiddleTS);
	outputtree->Branch("MTASOutterTS",&MTASOutterTS);

	outputtree->Branch("Beta",&Beta);
	outputtree->Branch("BetaTS",&BetaTS);
	outputtree->Branch("BetaFrontTrace",&BetaFrontTrace);
	outputtree->Branch("BetaBackTrace",&BetaBackTrace);
}

bool PMTASProcessor::PreProcess(RawEvent &event) {
	if (!EventProcessor::PreProcess(event))
		return false;

	return true;
}

bool PMTASProcessor::Process(RawEvent &event) {
	if (!EventProcessor::Process(event))
		return false;
	
	Reset();

	const auto &Events = event.GetSummary("PMTAS", true)->GetList();
	for (auto it = Events.begin(); it != Events.end(); it++) {
		ChanNum = (*it)->GetChannelNumber();
		ModNum = (*it)->GetModuleNumber();
		subtype = (*it)->GetChanID().GetSubtype();
		if( subtype.compare("center") == 0 ){
			currtype = SUBTYPE::CENTER;
		}else if(subtype.compare("inner") == 0 ){
			currtype = SUBTYPE::INNER;
		}else if(subtype.compare("middle") == 0 ){
			currtype = SUBTYPE::MIDDLE;
		}else if(subtype.compare("outter") == 0 ){
			currtype = SUBTYPE::OUTTER;
		}else if(subtype.compare("yap") == 0 or subtype.compare("bsm") == 0 ){
			currtype = SUBTYPE::BETA;
		}else{
			currtype = SUBTYPE::UNKNOWN;
		}
		group = (*it)->GetChanID().GetGroup();
		position = std::stoi(group) - 1;

		DidSaturate = (*it)->IsSaturated();
		DidPileup = (*it)->IsPileup();

		IsFront = (*it)->GetChanID().HasTag("front");
		IsBack = (*it)->GetChanID().HasTag("back");

		detectorposition = 2*position + IsBack;
		if( IsFront and IsBack ){
			throw std::runtime_error("Incorrect config file as Module "+std::to_string(ModNum)+" channel "+std::to_string(ChanNum)+" is both a front and back");
		}

		if( DidSaturate or DidPileup ){
			if( currtype == CENTER or currtype == INNER or currtype == MIDDLE or currtype == OUTTER ){
				MTASSaturatePileup = true;
			}
			if( currtype == BETA ){
				BetaSaturatePileup = true;
			}
			continue;
		}

		if( currtype == SUBTYPE::CENTER ){
			if( !CenterHits[detectorposition] ){
				MTASCenter[detectorposition] += (*it)->GetCalibratedEnergy();
				MTASCenterTS[detectorposition] = (*it)->GetTimeSansCfd();
				CenterHits[detectorposition]++;
			}
		}
		
		if( currtype == SUBTYPE::INNER ){
			if( !InnerHits[detectorposition] ){
				MTASInner[detectorposition] += (*it)->GetCalibratedEnergy();
				MTASInnerTS[detectorposition] = (*it)->GetTimeSansCfd();
				InnerHits[detectorposition]++;
			}
		}

		if( currtype == SUBTYPE::MIDDLE ){
			if( !MiddleHits[detectorposition] ){
				MTASMiddle[detectorposition] += (*it)->GetCalibratedEnergy();
				MTASMiddleTS[detectorposition] = (*it)->GetTimeSansCfd();
				MiddleHits[detectorposition]++;
			}
		}

		if( currtype == SUBTYPE::OUTTER ){
			if( !OutterHits[detectorposition] ){
				MTASOutter[detectorposition] += (*it)->GetCalibratedEnergy();
				MTASOutterTS[detectorposition] = (*it)->GetTimeSansCfd();
				OutterHits[detectorposition]++;
			}
		}

		if( currtype == SUBTYPE::BETA ){
			if( !BetaHits[detectorposition] ){
				Beta[detectorposition] += (*it)->GetCalibratedEnergy();
				BetaTS[detectorposition] = (*it)->GetTimeSansCfd();
				BetaHits[detectorposition]++;
				if( IsFront ){
					BetaFrontTrace = (*it)->GetTrace();
				}else{
					BetaBackTrace = (*it)->GetTrace();
				}
			}
		}
	}

	if( not MTASSaturatePileup and not BetaSaturatePileup ){
		for( int ii = 0; ii < 6; ++ii ){
			if( CenterHits[2*ii] and CenterHits[2*ii + 1] ){
				MTASCenterSumFrontBack[ii] += (MTASCenter[2*ii] + MTASCenter[2*ii + 1])/2.0;
			}
			if( InnerHits[2*ii] and InnerHits[2*ii + 1] ){
				MTASInnerSumFrontBack[ii] += (MTASInner[2*ii] + MTASInner[2*ii + 1])/2.0;
			}
			if( MiddleHits[2*ii] and MiddleHits[2*ii + 1] ){
				MTASMiddleSumFrontBack[ii] += (MTASMiddle[2*ii] + MTASMiddle[2*ii + 1])/2.0;
			}
			if( OutterHits[2*ii] and OutterHits[2*ii + 1] ){
				MTASOutterSumFrontBack[ii] += (MTASOutter[2*ii] + MTASOutter[2*ii + 1])/2.0;
			}
		}
		if( BetaHits[0]  and BetaHits[1] ){
			BetaSumFrontBack += (Beta[0] + Beta[1]);
			BetaAverage += (Beta[0] + Beta[1])/2.0;
		}

		for( int ii = 0; ii < 6; ++ii ){
			MTASCenterSum += MTASCenterSumFrontBack[ii];
			MTASInnerSum += MTASInnerSumFrontBack[ii];
			MTASMiddleSum += MTASMiddleSumFrontBack[ii];
			MTASOutterSum += MTASOutterSumFrontBack[ii];
			
			MTASTotalSum += MTASCenterSumFrontBack[ii];
			MTASTotalSum += MTASInnerSumFrontBack[ii];
			MTASTotalSum += MTASMiddleSumFrontBack[ii];
			MTASTotalSum += MTASOutterSumFrontBack[ii];
		}

		outputtree->Fill();
	}
	EndProcess();
	return true;
}

void PMTASProcessor::Reset(){
	BetaHits = ZeroBetaHits;
	BetaFrontTrace = ZeroTrace;
	BetaBackTrace = ZeroTrace;
	BetaSumFrontBack = 0.0;
	BetaAverage = 0.0;
	Beta = ZeroBeta;
	BetaTS = ZeroBeta;

	BetaSaturatePileup = false;

	CenterHits = ZeroHits;
	InnerHits = ZeroHits;
	MiddleHits = ZeroHits;
	OutterHits = ZeroHits;

	MTASTotalSum = 0.0;
	MTASCenterSum = 0.0;
	MTASInnerSum = 0.0;
	MTASMiddleSum = 0.0;
	MTASOutterSum = 0.0;

	MTASSaturatePileup = false;

	MTASCenterSumFrontBack = Zero;
	MTASInnerSumFrontBack = Zero;
	MTASMiddleSumFrontBack = Zero;
	MTASOutterSumFrontBack = Zero;

	MTASCenter = IndZero;
	MTASInner = IndZero;
	MTASMiddle = IndZero;
	MTASOutter = IndZero;

	MTASCenterTS = IndZero;
	MTASInnerTS = IndZero;
	MTASMiddleTS = IndZero;
	MTASOutterTS = IndZero;

}
