/** @file BSMProcessor.hpp
 * @brief  Basic BSMProcessor
 * @author T. Ruland
 * @date 04/25/2022
 */
#ifndef PAASS_BSMProcessor_H
#define PAASS_BSMProcessor_H

#include <utility>

#include "EventProcessor.hpp"
#include "PaassRootStruct.hpp"
#include "RawEvent.hpp"
#include "Globals.hpp"

#include <vector>
#include <utility>
#include <functional>

class BSMSegment {
	public: 
		/** Constructor */
		BSMSegment(bool zerosuppress){
			segFront_ = nullptr;
			segBack_ = nullptr;
			gBSMSegID_ = -1;
			HasZeroSuppression = zerosuppress;

		};
		/** Destructor */
		~BSMSegment() = default;
		bool IsValidSegment() const { return segBack_ != nullptr and segFront_ != nullptr; }
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
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( IsValidSegment() ){
				return std::make_pair((segFront_->GetTimeSansCfd() - segBack_->GetTimeSansCfd()) * clockInSeconds * 1.0e9,true);
			}else{
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
				return std::make_pair(segBack_->GetTimeSansCfd() * clockInSeconds * 1.0e9,true);
			}
		}

		std::pair<double,bool> GetBackEnergy() const{
			if( segBack_ == nullptr ){
				return std::make_pair(0.0,not HasZeroSuppression);
			}else{
				return std::make_pair(segBack_->GetCalibratedEnergy(),true);
			}
		}

	public:
		int gBSMSegID_;
		ChanEvent* segFront_;
		ChanEvent* segBack_ ;
		std::string PixieRev;
		bool HasZeroSuppression;
};


class BSMProcessor : public EventProcessor {
	public:
		/**Constructor */
		BSMProcessor(int,bool,bool,std::vector<std::pair<double,double>>,double,double,double,double,double);

		/** Deconstructor */
		~BSMProcessor() = default;

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
		std::string PixieRev; //! pixie revision
		int NumSegments;
		bool HasZeroSuppression;
		bool StandAlone;
		double Threshold;
		double MeanEnergy;
		double a_0;
		double a_1;
		double a_2;

		std::pair<double,bool> BSMTotal;
		std::pair<double,bool> FrontAvg;
		std::pair<double,bool> BackAvg;
		std::vector<BSMSegment> BSMSegVec;
		double BSMPosition;

		std::vector<std::pair<double,double>> MTASGates;
		const unsigned int MaxGates = 10;
		unsigned int NumGates;
};

#endif  //PAASS_BSMProcessor_H
