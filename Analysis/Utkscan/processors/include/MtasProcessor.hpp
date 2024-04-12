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

#include <utility>

class MtasSegment {
	public: 
		/** Constructor */
		MtasSegment(bool zerosuppress){
			segFront_ = nullptr;
			segBack_ = nullptr;
			gMtasSegID_ = -1;
			HasZeroSuppression = zerosuppress;
		};
		/** Destructor */
		~MtasSegment() = default;
		bool IsValidSegment() const { return segBack_ != nullptr and segFront_ != nullptr; }
		bool IsValidFront() const { return segFront_ != nullptr; }
		bool IsValidBack() const { return segBack_ != nullptr; }
		std::pair<double,bool> GetSegmentPosition() const{
			if( IsValidSegment() ){
				return std::make_pair((segFront_->GetCalibratedEnergy() - segBack_->GetCalibratedEnergy())/(segFront_->GetCalibratedEnergy() + segBack_->GetCalibratedEnergy()),true);
			}else{
				return std::make_pair(0.0,not HasZeroSuppression);
			}
		}
				
		std::pair<double,bool> GetSegmentAverageEnergy() const { 
			if( IsValidSegment() ){
				return std::make_pair((segFront_->GetCalibratedEnergy() + segBack_->GetCalibratedEnergy())/2.0,true); 
			}else{
				return std::make_pair(0.0,not HasZeroSuppression);
			}
		}
		std::pair<double,bool> GetSegmentTdiffInNS() const { 
			double clockInSeconds;
			//std::cout << "Get Clock" << std::endl;
			if (PixieRev == "F"){
				//std::cout << "Rev F" << std::endl;
				//std::cout << segFront_->GetChanID().GetModFreq() << std::endl;
				//std::cout << "Rev F" << std::endl;
				if( segFront_ != nullptr ){
					clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
				}else{
					clockInSeconds = Globals::get()->GetClockInSeconds(segBack_->GetChanID().GetModFreq());
				}
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			//std::cout << "Finished Get Clock" << std::endl;
			if( IsValidSegment() ){
				//std::cout << "Valid Segment" << std::endl;
				return std::make_pair((segFront_->GetTimeSansCfd() - segBack_->GetTimeSansCfd()) * clockInSeconds * 1.0e9,true);
			}else{
				//std::cout << "InValid Segment" << std::endl;
				return std::make_pair(0.0,not HasZeroSuppression);
			}
		}
		std::pair<double,bool> GetFrontEnergy() const{
			if( segFront_ == nullptr ){
				return std::make_pair(0.0,not HasZeroSuppression);
			}else{
				return std::make_pair(segFront_->GetCalibratedEnergy(),true);
			}
		}
		std::pair<double,bool> GetBackEnergy() const{
			if( segBack_ == nullptr ){
				return std::make_pair(0.0,not HasZeroSuppression);
			}else{
				return std::make_pair(segBack_->GetCalibratedEnergy(),true);
			}
		}
	
		std::pair<double,bool> GetFrontTimeInNS() const{
			double clockInSeconds;
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( segFront_ == nullptr ){
				return std::make_pair(0.0,not HasZeroSuppression);
			}else{
				return std::make_pair(segFront_->GetTimeSansCfd() * clockInSeconds * 1.0e9,true);
			}
		}
		std::pair<double,bool> GetBackTimeInNS() const{
			double clockInSeconds;
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( segBack_ == nullptr ){
				return std::make_pair(0.0,not HasZeroSuppression);
			}else{
				return std::make_pair(segFront_->GetTimeSansCfd() * clockInSeconds * 1.0e9,true);
			}
		}


	public:
		int gMtasSegID_;
		std::string segRing_;
		ChanEvent* segFront_;
		ChanEvent* segBack_ ;
		std::string PixieRev;
		bool HasZeroSuppression;
};

class MtasProcessor : public EventProcessor {
	public:
		/**Constructor */
		MtasProcessor(bool,bool,double,int);

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
		processor_struct::MTAS Mtasstruct;  //!<Root Struct
		std::string PixieRev; //! pixie revision
		bool IsNewCenter;
		bool HasZeroSuppression;
		bool IsPrevBetaTriggered;
		double BetaThreshold;
		int NumPMTFire;

		double SignalTime;
		double PrevBetaEnergy;
		double PrevBetaTime;
		double clockInSeconds;

		bool MTASSaturatePileup;
	
		std::vector<MtasSegment> MtasSegVec;
		std::pair<double,bool> centerSum;
		std::pair<double,bool> innerSum;
		std::pair<double,bool> middleSum;
		std::pair<double,bool> outerSum;
		std::pair<double,bool> totalSum;
};

#endif  //PAASS_MtasProcessor_H
