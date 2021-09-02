/**@file MMTASProcessor.hpp
 *@brief  Basic ROOT output. Fills root tree to mimic CAEN CoMPASS output. It has NO damm output
 *@authors T. Ruland 
 *@date 07/27/2021
 */
#ifndef PAASS_MMTASProcessor_H
#define PAASS_MMTASProcessor_H

#include "EventProcessor.hpp"
#include "PaassRootStruct.hpp"
#include "RawEvent.hpp"

#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TArrayS.h"

struct MMTASSingleDetector{
	double finetimestamp;
	double timestamp;
	double rawenergy;
	double energy;
	double psd;
	std::vector<short> trace;

	MMTASSingleDetector(){
		finetimestamp = -1.0;
		timestamp = -1.0;
		rawenergy = 0.0;
		energy = 0.0;
		psd = 1.0;
		trace.push_back(-1)
	}

	MMTASSingleDetector(const MMTASSingleDetector& other){
		this->finetimestamp = other.finetimestamp;
		this->timestamp = other.timestamp;
		this->rawenergy = other.rawenergy;
		this->energy = other.energy;
		this->psd = other.psd;
		this->trace = other.trace;
	}

	const MMTASSingleDetector& operator=(const MMTASSingleDetector& other){
		this->finetimestamp = other.finetimestamp;
		this->timestamp = other.timestamp;
		this->rawenergy = other.rawenergy;
		this->energy = other.energy;
		this->psd = other.psd;
		this->trace = other.trace;

		return *this;
	}
};

class MMTASProcessor : public EventProcessor {
	public:
		/**Constructor */
		MMTASProcessor();

		/** Deconstructor */
		~MMTASProcessor() {
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
		std::map<std::string,MMTASSingleDetector> HitsMap;
		std::map<std::string,short> Count;
		std::map<short,std::string> DetNameMap;

		bool RecordTraces;
				
		const short NumLeftNaI = 12;
		const short NumRightNaI = 12;
		const short NumLeftBees = 1;
		const short NumRightBees = 1;

		void ClearHitsMap() noexcept;
};

#endif  //PAASS_MMTASProcessor_H
