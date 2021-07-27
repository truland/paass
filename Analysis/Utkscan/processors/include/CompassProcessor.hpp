/**@file CompassProcessor.hpp
 *@brief  Basic ROOT output. Fills root tree to mimic CAEN CoMPASS output. It has NO damm output
 *@authors T. Ruland 
 *@date 07/27/2021
 */
#ifndef PAASS_CompassProcessor_H
#define PAASS_CompassProcessor_H

#include "EventProcessor.hpp"
#include "PaassRootStruct.hpp"
#include "RawEvent.hpp"

#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TArrayS.h"

class CompassProcessor : public EventProcessor {
	public:
		/**Constructor */
		CompassProcessor();

		/** Deconstructor */
		~CompassProcessor() {
			//PixieFile = PTree->GetCurrentFile();
			//PixieFile->Write();
			//PixieFile->Close();
			outputfile = outputtree->GetCurrentFile();
			outputfile->Write();
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

	private:
		std::string Rev;
		TFile* outputfile;
		TTree* outputtree;
		unsigned long long timestamp;
		unsigned short board;
		unsigned short channel;
		unsigned short energy;
		unsigned short energyshort;
		unsigned int flags;
		TArrayS* samples;
};

#endif  //PAASS_CompassProcessor_H
