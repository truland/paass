/**@file MMTASProcessor.hpp
 *@brief  Basic ROOT output. Fills root tree to mimic CAEN CoMPASS output. It has NO damm output
 *@authors T. Ruland 
 *@date 07/27/2021
 */
#ifndef PAASS_PMTASProcessor_H
#define PAASS_PMTASProcessor_H

#include "EventProcessor.hpp"
#include "RawEvent.hpp"

#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"

class PMTASProcessor : public EventProcessor {
	public:
		/**Constructor */
		PMTASProcessor();

		/** Deconstructor */
		~PMTASProcessor() {
			std::cout << "---Deleted PMTASProcessor---" << std::endl;
			outputfile = outputtree->GetCurrentFile();
			outputfile->Write(0,2,0);
			outputfile->Close();
		}

		/** Preprocess the event
		 * \param [in] event : the event to preprocess
		 * \return true if successful
		 */
		virtual bool PreProcess(RawEvent &event);

		/** Process the event
		 * \param [in] event : the event to process
		 * \return true if successful
		 */
		virtual bool Process(RawEvent &event);

		void Reset();

	private:
		enum SUBTYPE{
			CENTER,
			INNER,
			MIDDLE,
			OUTTER,
			BETA,
			UNKNOWN
		};

		std::string Rev;
		TFile* outputfile;
		TTree* outputtree;

		std::vector<double> Zero;
		std::vector<double> IndZero;
		std::vector<int> ZeroHits;
		std::vector<unsigned int> ZeroTrace;
		std::vector<int> ZeroBetaHits;
		std::vector<double> ZeroBeta;

		bool IsFront;
		bool IsBack;
		bool DidSaturate;
		bool DidPileup;
		int ChanNum;
		int ModNum;
		int position;
		int detectorposition;
		SUBTYPE currtype;
		std::string subtype;
		std::string group;

		std::vector<int> CenterHits;
		std::vector<int> InnerHits;
		std::vector<int> MiddleHits;
		std::vector<int> OutterHits;
		
		//YAP, etc
		std::vector<int> BetaHits;
		std::vector<double> Beta;
		double BetaSumFrontBack;
		double BetaAverage;
		std::vector<double> BetaTS;
		std::vector<unsigned int> BetaFrontTrace;
		std::vector<unsigned int> BetaBackTrace;
		bool BetaSaturatePileup;

		double MTASTotalSum;
		double MTASCenterSum;
		double MTASInnerSum;
		double MTASMiddleSum;
		double MTASOutterSum;

		bool MTASSaturatePileup;

		std::vector<double> MTASCenterSumFrontBack;
		std::vector<double> MTASInnerSumFrontBack;
		std::vector<double> MTASMiddleSumFrontBack;
		std::vector<double> MTASOutterSumFrontBack;

		std::vector<double> MTASCenter;
		std::vector<double> MTASInner;
		std::vector<double> MTASMiddle;
		std::vector<double> MTASOutter;


		std::vector<double> MTASCenterTS;
		std::vector<double> MTASInnerTS;
		std::vector<double> MTASMiddleTS;
		std::vector<double> MTASOutterTS;

};

#endif  //PAASS_PMTASProcessor_H
