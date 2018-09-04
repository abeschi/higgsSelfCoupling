#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "interface/CMS_lumi.h"
//#include "interface/TreeUtils.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <fstream>

#include "TLorentzVector.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

#define MZ 91.187
#define PI 3.14159265359

#define bDiscriminantThresholdLooseCSV 0.5426
#define bDiscriminantThresholdMediumCSV 0.8484
#define bDiscriminantThresholdTightCSV 0.9535

#define bDiscriminantThresholdLooseDeep  0.1522
#define bDiscriminantThresholdMediumDeep 0.4941
#define bDiscriminantThresholdTightDeep  0.8001


float DeltaEta(const float& eta1, const float& eta2);
float DeltaPhi(const float& phi1, const float& phi2);
float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);

void  writeMap (std::vector<std::string> names, std::map<std::string,std::vector<float>> C1, std::string outputFile);

void MakePlot(TH1F**, TString title);
void MakePlot2(std::map<std::string,TH1F*>& histos, TString title);

#endif
