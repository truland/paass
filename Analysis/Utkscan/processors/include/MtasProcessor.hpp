/** @file MtasProcessor.hpp
 * @brief  Basic MtasProcessor for MTAS at FRIB
 * @authors T.T. King, T. Ruland, B.C. Rasco
 * @date 03/25/2022
 */
#ifndef PAASS_MtasProcessor_H
#define PAASS_MtasProcessor_H

#include <utility>

#include "EventProcessor.hpp"
#include "PaassRootStruct.hpp"
#include "RawEvent.hpp"
#include "Globals.hpp"

class MtasSegment {
	public: 
		/** Constructor */
		MtasSegment(){
			segFront_ = nullptr;
			segBack_ = nullptr;
			gMtasSegID_ = -1;

		};
		/** Destructor */
		~MtasSegment() = default;
		bool IsValidSegment() const { return segBack_ != nullptr and segFront_ != nullptr; }
		std::pair<double,bool> GetSegmentPosition() const{
			if( IsValidSegment() ){
				return std::make_pair((segFront_->GetCalibratedEnergy() - segBack_->GetCalibratedEnergy())/(segFront_->GetCalibratedEnergy() + segBack_->GetCalibratedEnergy()),true);
			}else{
				return std::make_pair(0.0,false);
			}
		}
				
		std::pair<double,bool> GetSegmentAverageEnergy() const { 
			if( IsValidSegment() ){
				return std::make_pair((segFront_->GetCalibratedEnergy() + segBack_->GetCalibratedEnergy())/2.0,true); 
			}else{
				return std::make_pair(0.0,false);
			}
		}
		std::pair<double,bool> GetSegmentTdiffInNS() const { 
			double clockInSeconds;
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( IsValidSegment() ){
				return std::make_pair((segFront_->GetTimeSansCfd() - segBack_->GetTimeSansCfd()) * clockInSeconds * 1.0e9,true);
			}else{
				return std::make_pair(0.0,false);
			}
		}
		std::pair<double,bool> GetFrontEnergy() const{
			if( segFront_ == nullptr ){
				return std::make_pair(0.0,false);
			}else{
				return std::make_pair(segFront_->GetCalibratedEnergy(),true);
			}
		}
		std::pair<double,bool> GetBackEnergy() const{
			if( segBack_ == nullptr ){
				return std::make_pair(0.0,false);
			}else{
				return std::make_pair(segBack_->GetCalibratedEnergy(),true);
			}
		}

	public:
		int gMtasSegID_;
		ChanEvent* segFront_;
		ChanEvent* segBack_ ;
		std::string PixieRev;
};

class MtasProcessor : public EventProcessor {
	public:
		/**Constructor */
		MtasProcessor(bool);

		/** Deconstructor */
		~MtasProcessor() = default;

		/** Preprocess the event
		 * \param [in] event : the event to preprocess
		 * \return true if successful
		 */
		bool PreProcess(RawEvent &event);

		/** Process the event
		 * \param [in] event : the event to process
		 * \return true if successful
		 */
		bool Process(RawEvent &event);


		/** Declares the plots for the class */
		void DeclarePlots(void);

	private:
		//processor_struct::MTAS Mtasstruct;  //!<Root Struct
		std::string PixieRev; //! pixie revision
		bool IsNewCenter;

};

#endif  //PAASS_MtasProcessor_H
