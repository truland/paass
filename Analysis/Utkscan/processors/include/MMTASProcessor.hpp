/**@file MMTASProcessor.hpp
 *@brief  Basic ROOT output. Fills root tree to mimic CAEN CoMPASS output. It has NO damm output
 *@authors T. Ruland 
 *@date 07/27/2021
 */
#ifndef PAASS_MMTASProcessor_H
#define PAASS_MMTASProcessor_H

#include "EventProcessor.hpp"
#include "MMTASStruct.hpp"
#include "RawEvent.hpp"

#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TArrayS.h"

class MMTASProcessor : public EventProcessor {
	public:
		/**Constructor */
		MMTASProcessor();

		/** Deconstructor */
		~MMTASProcessor() {
			std::cout << "---Deleted MMTASProcessor---" << std::endl;
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

	private:
		std::string Rev;
		TFile* outputfile;
		TTree* outputtree;

		MMTASStruct MMstruct;
		std::vector<MMTASStruct> mmtas_vec_;
};

#endif  //PAASS_MMTASProcessor_H
