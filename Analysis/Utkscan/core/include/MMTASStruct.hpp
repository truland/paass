/**Created by T.T. King 06/22/2018
 * */
/** PixTreeEvent class added by R. Yokoyama 08/31/2018 **/
#ifndef MMTAS_MMTASSTRUC_HPP
#define MMTAS_MMTASSTRUC_HPP

#include <TObject.h>
#include <TString.h>

struct  MMTASStruct{
    double energy = -999;
    double rawEnergy = -999;
    double timeSansCfd = -999;
    double time = -999;
    int detNum = -999;   //the instance number of RD in the xml Map
    int modNum = -999;   // the physical module number
    int chanNum = -999;  // the physical channel number
    TString subtype = "";
    TString group = "";
    bool pileup = false;                   //Did pixie detect pileup in the event
    bool saturation = false;               //Did the trace go out of the ADC range
    std::vector<unsigned int> trace = {};  //The trace if present
    double baseline = -999;
    double stdBaseline = -999;
    double phase = -999;
    double tqdc = -999;                      //QDC from trace (requires the Waveform Analyzer)
    int maxPos = -999;                       //Max location in the trace  (requires the Waveform Analyzer)
    double maxVal = -999;                    //Max value in the trace (requires the Waveform Analyzer)
    double extMaxVal = -999;                 // Extrapolated Max value in the trace (requires the Waveform Analyzer)
    double highResTime = -999;               //High Resolution Time derived from the trace fitting (requires the Waveform and Fitting Analyzer)
    std::vector<unsigned int> qdcSums = {};  //output the onboard qdc sums if present
    bool hasValidTimingAnalysis = false;
    bool hasValidWaveformAnalysis = false;
};
static const MMTASStruct MMTAS_DEFAULT_STRUCT;

#endif  //MMTAS_PROCESSORSTRUC_HPP
