/** \file GeCalibProcessor.cpp
 *
 * implementation for germanium processor for calibration and diagnostic
 */

//? Make a clover specific processor

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

#include <limits>
#include <cmath>
#include <cstdlib>

#include "Plots.hpp"
#include "PlotsRegister.hpp"
#include "DammPlotIds.hpp"

#include "Correlator.hpp"
#include "DetectorLibrary.hpp"
#include "GeProcessor.hpp"
#include "GeCalibProcessor.hpp"
#include "RawEvent.hpp"
#include "Messenger.hpp"
#include "Exceptions.hpp"

using namespace std;

namespace dammIds {
    namespace ge {
        namespace calib {
            const int D_E_SUM = 400;
            const int D_ENERGY_HIGHGAIN = 401;
            const int DD_CLOVER_ENERGY_RATIO = 402;
            const int D_E_CRYSTALX = 403;

            const int DD_E_DETX = 419;
            const int DD_E_DETX_BETA_GATED = 499;

            const int DD_EGAMMA__EBETA_ALL = 420;
            const int DD_EGAMMA__EBETA_B0 = 421;
            const int DD_EGAMMA__EBETA_B1 = 422;
            const int DD_EGAMMA__EBETA = 423;

            const int DD_TDIFF__EGAMMA_ALL = 460;
            const int DD_TDIFF__EGAMMA_B0 = 461;
            const int DD_TDIFF__EGAMMA_B1 = 462;
            const int DD_TDIFF__EGAMMA = 463;
        }
    }
}

GeCalibProcessor::GeCalibProcessor() : GeProcessor() {
}

/** Declare plots including many for decay/implant/neutron gated analysis  */
void GeCalibProcessor::DeclarePlots(void) 
{
    const int energyBins = SE;
    const int energyBins2 = SC;
    const int energyBins3 = S9;
    const int timeBins = S7;

    using namespace dammIds::ge;

    DetectorLibrary* modChan = DetectorLibrary::get();

    const set<int> &cloverLocations = modChan->GetLocations("ge", "clover_high");
    // could set it now but we'll iterate through the locations to set this
    unsigned int cloverChans = 0;

    for ( set<int>::const_iterator it = cloverLocations.begin();
	  it != cloverLocations.end(); it++) {
        leafToClover[*it] = int(cloverChans / 4); 
        cloverChans++;
    }

    if (cloverChans % chansPerClover != 0) {
        cout << " There does not appear to be the proper number of"
            << " channels per clover.\n Program terminating." << endl;
        exit(EXIT_FAILURE);	
    }

    if (cloverChans != 0) {
        numClovers = cloverChans / chansPerClover;
        Messenger m;
        m.start("Building clovers");

        stringstream ss;
        ss << "A total of " << cloverChans 
           << " clover channels were detected: ";
        int lastClover = numeric_limits<int>::min();
        for ( map<int, int>::const_iterator it = leafToClover.begin();
            it != leafToClover.end(); it++ ) {
            if (it->second != lastClover) {
                m.detail(ss.str());
                ss.str("");
                lastClover = it->second;
                ss << "Clover " << lastClover << " : ";
            } else {
                ss << ", ";
            }
            ss << setw(2) << it->first;
        }
        m.detail(ss.str());

        if (numClovers > dammIds::ge::MAX_CLOVERS) {
            m.fail();
            stringstream ss;
            ss << "Number of detected clovers is greater than defined" 
               << " MAX_CLOVERS = " << dammIds::ge::MAX_CLOVERS << "."
               << " See GeCalibProcessor.hpp for details.";
            throw GeneralException(ss.str());
        }   
        m.done();
    }


    DeclareHistogram1D(calib::D_E_SUM, energyBins, "Gamma energy");
    for (unsigned i = 0; i < cloverChans; ++i) {
        stringstream ss;
        ss << "Gamma energy crystal " << i << " Clover " << leafToClover[i];
        DeclareHistogram1D(calib::D_E_CRYSTALX + i,
                           energyBins, ss.str().c_str());
    }
    DeclareHistogram2D(calib::DD_E_DETX, energyBins, S5,
                       "Gamma E vs. crystal number ");
    DeclareHistogram2D(calib::DD_E_DETX_BETA_GATED, energyBins, S5,
                       "Gamma E vs. crystal number beta gated");

    DeclareHistogram2D(calib::DD_EGAMMA__EBETA_ALL,
                       energyBins2, energyBins2, 
                       "Gamma E, Beta E (all crystals)");
    DeclareHistogram2D(calib::DD_TDIFF__EGAMMA_ALL,
                       timeBins, energyBins, 
                       "dt gamma-beta, gamma E (all)");

    DeclareHistogram2D(calib::DD_EGAMMA__EBETA_B0,
                    energyBins2, energyBins2, "Gamma E, beta E; beta 0");
    DeclareHistogram2D(calib::DD_EGAMMA__EBETA_B1,
                    energyBins2, energyBins2, "Gamma E, beta E; beta 0");
    DeclareHistogram2D(calib::DD_TDIFF__EGAMMA_B0,
                       timeBins, energyBins, "dt gamma-beta, gamma E; beta 0");
    DeclareHistogram2D(calib::DD_TDIFF__EGAMMA_B1,
                       timeBins, energyBins, "dt gamma-beta, gamma E; beta 0");

    for (unsigned b = 0; b < 2; ++b) {
        stringstream ss;
        for (unsigned i = 0; i < cloverChans; ++i) {
            ss.str("");
            ss << "Gamma E, beta E; crystal " << i 
               << " beta " << b << " (energy/2)";
            DeclareHistogram2D(calib::DD_EGAMMA__EBETA + i * 2 + b,
                               energyBins2, energyBins3, ss.str().c_str());

            ss.str("");
            ss << "dt gamma-beta, gamma E; crystal " 
               << i << " beta " << b;
            DeclareHistogram2D(calib::DD_TDIFF__EGAMMA + i * 2 + b,
                               timeBins, energyBins, ss.str().c_str());
        }
    }

    DeclareHistogram1D(calib::D_ENERGY_HIGHGAIN, energyBins,
                       "Gamma singles, high gain");
    DeclareHistogram2D(calib::DD_CLOVER_ENERGY_RATIO, S4, S6,
                       "high/low energy ratio (x10)");
    DeclareHistogram1D(dammIds::ge::D_MULT, S3, 
                       "Gamma multiplicity");                  
}

bool GeCalibProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return false;

    using namespace dammIds::ge;

    // Clear all events stored in vectors from previous event
    geEvents_.clear();

    static const vector<ChanEvent*> &highEvents = event.GetSummary("ge:clover_high", true)->GetList();
    static const vector<ChanEvent*> &lowEvents  = event.GetSummary("ge:clover_low", true)->GetList();

    /** Only the high gain events are going to be used. The events where
     * low/high gain mismatches, saturation or pileup is marked are rejected
     */
    for (vector<ChanEvent*>::const_iterator itHigh = highEvents.begin();
	 itHigh != highEvents.end(); itHigh++) {
        int location = (*itHigh)->GetChanID().GetLocation();
        plot(calib::D_ENERGY_HIGHGAIN, (*itHigh)->GetCalEnergy());

        if ( (*itHigh)->IsSaturated() || (*itHigh)->IsPileup() )
            continue;


        // find the matching low gain event
        vector <ChanEvent*>::const_iterator itLow = lowEvents.begin();
        for (; itLow != lowEvents.end(); itLow++) {
            if ( (*itLow)->GetChanID().GetLocation() == location ) {
                break;
            }
        }
        if ( itLow != lowEvents.end() ) {
            double ratio = (*itHigh)->GetEnergy() / (*itLow)->GetEnergy();
            plot(calib::DD_CLOVER_ENERGY_RATIO, location, ratio * 10.);
            if ( (ratio < detectors::geLowRatio ||
                  ratio > detectors::geHighRatio) )
                continue;
        }
        geEvents_.push_back(*itHigh);
    }


    /** NOTE we do permanents changes to events here
     *  Necessary in order to set corrected time for use in correlator
     */

    /** 
     * Walk Correction turned off to investigate walk correction parameters
     *
    for (vector<ChanEvent*>::iterator it = geEvents_.begin(); 
	 it != geEvents_.end(); it++) {
        double energy = (*it)->GetCalEnergy();
        double time   = (*it)->GetTime() - WalkCorrection(energy);	
        (*it)->SetCorrectedTime(time);
    }
    */

    // now we sort the germanium events according to their corrected time
    sort(geEvents_.begin(), geEvents_.end(), CompareCorrectedTime);

    return true;
} 


bool GeCalibProcessor::Process(RawEvent &event) {
    using namespace dammIds::ge;

    if (!EventProcessor::Process(event))
        return false;

    bool hasBeta = TreeCorrelator::get()->place("Beta")->status();

    plot(D_MULT, geEvents_.size());

    for (vector<ChanEvent*>::iterator it = geEvents_.begin(); 
         it != geEvents_.end(); ++it) {
        ChanEvent* chan = *it;
        double gEnergy = chan->GetCalEnergy();	

        if (gEnergy < detectors::gammaThreshold)
            continue;

        double gTime = chan->GetCorrectedTime();
        int det = chan->GetChanID().GetLocation();

        plot(calib::D_E_CRYSTALX + det, gEnergy);
        plot(calib::D_E_SUM, gEnergy);
        plot(calib::DD_E_DETX, gEnergy, det);

        if (hasBeta) {
            plot(calib::DD_E_DETX_BETA_GATED, gEnergy, det);

            EventData beta = TreeCorrelator::get()->place("Beta")->last();
            double gb_dtime = (gTime - beta.time);
            double betaEnergy = beta.energy;
            int betaLocation = beta.location;
            int timeShift = 20;

            plot(calib::DD_TDIFF__EGAMMA_ALL,
                    (gb_dtime + timeShift), gEnergy);
            plot(calib::DD_EGAMMA__EBETA_ALL, gEnergy, betaEnergy);

            if (betaLocation == 0) {
                plot(calib::DD_EGAMMA__EBETA_B0,
                    gEnergy, betaEnergy);
                plot(calib::DD_TDIFF__EGAMMA_B0,
                    (gb_dtime + timeShift), gEnergy);
            } else if (betaLocation == 1) {
                plot(calib::DD_EGAMMA__EBETA_B1,
                    gEnergy, betaEnergy);
                plot(calib::DD_TDIFF__EGAMMA_B1,
                    (gb_dtime + timeShift), gEnergy);
            }

            plot(calib::DD_TDIFF__EGAMMA + 2 * det + betaLocation,
                 (gb_dtime + timeShift), gEnergy);

            plot(calib::DD_EGAMMA__EBETA + 2 * det + betaLocation,
                 gEnergy, betaEnergy / 2.0);
        }

    }

    EndProcess(); 
    return true;
}
