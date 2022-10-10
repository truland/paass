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
#include <cmath>

class BSMSegment {
	public: 
		/** Constructor */
		BSMSegment(){
			segFront_ = nullptr;
			segBack_ = nullptr;
			gBSMSegID_ = -1;
		};
		/** Destructor */
		~BSMSegment() = default;
		bool IsValidSegment() const { return segBack_ != nullptr and segFront_ != nullptr; }
		double GetSegmentPosition() const{
			if( IsValidSegment() ){
				return (segFront_->GetCalibratedEnergy() - segBack_->GetCalibratedEnergy())/(segFront_->GetCalibratedEnergy() + segBack_->GetCalibratedEnergy());
			}else{
				return 0.0;
			}
		}
				
		double GetSegmentAverageEnergy() const { 
			if( IsValidSegment() ){
				return (segFront_->GetCalibratedEnergy() + segBack_->GetCalibratedEnergy())/2.0; 
			}else{
				return 0.0;
			}
		}
		double GetSegmentTdiffInNS() const { 
			double clockInSeconds;
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( IsValidSegment() ){
				return (segFront_->GetTimeSansCfd() - segBack_->GetTimeSansCfd()) * clockInSeconds * 1.0e9;
			}else{
				return 0.0;
			}
		}
		double GetFrontEnergy() const{
			if( segFront_ == nullptr ){
				return 0.0;
			}else{
				return segFront_->GetCalibratedEnergy();
			}
		}
		double GetFrontTimeInNS() const{
			double clockInSeconds;
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( segFront_ == nullptr ){
				return 0.0;
			}else{
				return segFront_->GetTimeSansCfd() * clockInSeconds * 1.0e9;
			}
		}
		double GetBackTimeInNS() const{
			double clockInSeconds;
			if (PixieRev == "F"){
				clockInSeconds = Globals::get()->GetClockInSeconds(segFront_->GetChanID().GetModFreq());
			} else {
				clockInSeconds = Globals::get()->GetClockInSeconds();
			}
			if( segBack_ == nullptr ){
				return 0.0;
			}else{
				return segBack_->GetTimeSansCfd() * clockInSeconds * 1.0e9;
			}
		}

		double GetBackEnergy() const{
			if( segBack_ == nullptr ){
				return 0.0;
			}else{
				return segBack_->GetCalibratedEnergy();
			}
		}

	public:
		int gBSMSegID_;
		ChanEvent* segFront_;
		ChanEvent* segBack_ ;
		std::string PixieRev;
		bool HasZeroSuppression;
};

struct BSMPositionCorrection{
	double constant = 0.0;
	double slope = 0.0;
	double mean = 1.0;
	
	double Correct(double erg,double pos){
		double val = mean/std::exp(constant + slope*pos);
		//std::cout << erg << '\t' << pos << '\t' << erg/val << std::endl;
		return erg*val;
	}
};

class BSMProcessor : public EventProcessor {
	public:
		/**Constructor */
		BSMProcessor(int,bool,std::vector<std::pair<double,double>>,double,double,double,double,double,double,double);

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
		bool StandAlone;
		double Threshold;
		double MeanEnergy;
		double a_0;
		double a_1;
		double a_2;

		double BSMTotal;
		double FrontAvg;
		double BackAvg;
		std::vector<BSMSegment> BSMSegVec;
		BSMPositionCorrection FrontCorrection;
		BSMPositionCorrection BackCorrection;
		double BSMPosition;

		std::vector<std::pair<double,double>> MTASGates;
		const unsigned int MaxGates = 10;
		unsigned int NumGates;
};

#endif  //PAASS_BSMProcessor_H
