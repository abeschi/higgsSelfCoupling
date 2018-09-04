/*
g++ -Wall -o CutBasedDiphotons `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit -lRooStats  -lMathMore CutBasedDiphotons.cpp
*/

#ifndef CMS_LUMI_H
#include "interface/CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "interface/tdrstyle.h"
#endif

#ifndef ANALYSIS_UTILS_H
#include "interface/AnalysisUtils.h"
#endif

#ifndef TREE_UTILS_H
#include "interface/TreeUtils.h"
#endif

#include "CfgManager/interface/CfgManagerT.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLeaf.h"

#include <iostream>
#include <iomanip>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#else
// Refer to a class implemented in libRooFit to force its loading                                                                                                                                                   
// via the autoloader.                                                                                                                                                                                               
class Roo2DKeysPdf;
#endif

#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include <RooMCStudy.h>
#include "RooMsgService.h"
#include "RooCustomizer.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooClassFactory.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>

#include "Math/QuantFuncMathCore.h"
#include "Math/DistFunc.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h" 
  

using namespace std;
using namespace RooFit;




int main(int argc, char *argv[])
{
	if( argc < 2 )
	{
		std::cerr << ">>>>> MakeC1ForMassirozzo.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
		return -1;
	}

	writeExtraText = true;       // if extra text
	extraText  = "Preliminary Simulation";  // default extra text is "Preliminary"
	lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	// Parsing of the config file

	CfgManager opts;
	opts.ParseConfigFile(argv[1]);
  
	std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  	std::vector<std::string> C1 = opts.GetOpt<std::vector<std::string> >("C1Histo.C1Histo");

	std::string outputFile = opts.GetOpt<std::string>("Output.outputFile");
	std::string outputPlots = opts.GetOpt<std::string>("Output.outputPlots");
  
	std::vector<float>  nBins = opts.GetOpt<std::vector<float>>("Settings.nBins");
	std::vector<float>  ptBins_input = opts.GetOpt<std::vector<float>>("Settings.ptBins");
	int Entries = opts.GetOpt<int>("Settings.Entries");
  	bool doPlots = opts.GetOpt<int>("Settings.doPlots");
	std::vector<std::string> nameOfTags = opts.GetOpt< std::vector<std::string> >("Settings.nameOfTags");

	std::map<std::string,TChain*> trees;
	std::map<std::string,std::string> processes;
	std::map<std::string,TH1F*> pT_histo;

	std::map<std::string,TFile*> C1_file;
	std::map<std::string,TH1F*> C1_h;

	std::map<std::string,std::vector<float>> ptBins;

	int numberOfExpectedTagNames = 0;
	bool singleBinning = 0;

	if(nBins.size()!=input.size()/4 && nBins.size()!=1)
	{
		cerr << "Binning specified is not coherent with the number of input. Please specify the number of bin boundaries for each input file." << endl;
		return 0;
	}
	else if(nBins.size()==1 && input.size()/4!=1)
	{
		cout << "Only one bin size specified, assuming all the input files has the same bins." << endl;
		for(unsigned int i=1; i<input.size()/4; ++i)
			nBins.push_back(nBins.at(0));
		singleBinning = 1;
	}


	std::map<std::string,std::vector<float>> C1_Binned;

  	for(unsigned int n = 0; n < C1.size()/3; ++n)
  	{
		std::string FileName = C1.at(0+n*3);
    		std::string hName   = C1.at(1+n*3);
    		std::string label      = C1.at(2+n*3);

		if( C1_file[label] == 0 )
		{	C1_file[label] = new TFile((TString)FileName);
 			C1_h[label] = (TH1F*) C1_file[label] -> Get((TString)hName);
		}
  		std::cout << ">>> Getting C1 histo " << FileName << std::endl;
  	}

	int totalSize = 0;

  	for(unsigned int n = 0; n < input.size()/4; ++n)
  	{
		std::string inFileName = input.at(0+n*4);
    		std::string treeName   = input.at(1+n*4);
    		std::string label      = input.at(2+n*4);
    		std::string process    = input.at(3+n*4);

		if( trees[label] == 0 )
			trees[label] = new TChain((TString)treeName);
  		std::cout << ">>> Adding trees " << inFileName << " to chain " << label << std::endl;
    		trees[label] -> Add((TString)inFileName);
  		std::cout << ">>> Added trees " << inFileName << " to chain " << label << ", total entries; " << trees[label] -> GetEntries() << endl;
		if( pT_histo[label] == 0 )
			pT_histo[label] = new TH1F((TString)label, "; p_{T}^{H}; Entires", 10000, 0, 10000);
		if( processes[label] == "")
			processes[label] = process;
		auto it = ptBins.find(label);
		if( it==ptBins.end() )
		{
			for(int i=0; i<nBins[n]; ++i)
				ptBins[label].push_back(ptBins_input[totalSize+i]);
			if(!singleBinning)
				totalSize += nBins[n];
		}
  	}

	//Check binning in pT

	for(auto& it : ptBins)
	{
		std::string label = it.first;
		std::vector<float> bins = it.second;

		if(bins.size()==0)
		{	cerr << "No binning specified for one of the input in the config file, please check nBins and ptBins in the Setting section of the config file" << endl;
			return 0;
		}
		std::vector<int> indexesToRemove;
		//Pruning pTboundaries: remove pT bins higher than 10000 GeV or negative, remove doubled boundaries, inserting 0 and 10000 and lowest and higest boundary, necessary for the following code
		for(unsigned int i=0; i<bins.size(); ++i)
		{	
			if(bins.at(i)>=10000. || bins.at(i)<0.)
				indexesToRemove.push_back(i);
		}

		for(unsigned int i=0; i<indexesToRemove.size(); ++i)
		{	bins.at(indexesToRemove.at(i)) = 0;
			cout << "Find bin boundaries higher than 10000 GeV or negative, which is pretty unlikely. Removing them from the boundaries." << endl;
		}

		bins.push_back(10000);
		std::sort (bins.begin(), bins.end());
		bins.erase( unique( bins.begin(), bins.end() ), bins.end() ); //Remove doubled boundaries

		if(bins[0]!=0.)
			bins.insert(bins.begin(), 0.);

		it.second = bins;


		numberOfExpectedTagNames += bins.size()-1;
	}

  	// Starting computation for C1

	for(std::map<std::string,TChain*>::const_iterator treeIt = trees.begin(); treeIt != trees.end(); ++treeIt)
	{
		std::string label = treeIt -> first;
		TChain* tree = treeIt -> second;

		int nentries = 0;
		if(Entries<0)
			nentries = tree -> GetEntries();
		else
			nentries = min(Entries, (int)tree -> GetEntries());

		//Fill an histogram with the higgs pT
		for(int i=0; i<nentries; i++)
		{	
			if(i%10000==0) cout << "Processing tag " << label << ", event " << i << " out of " << nentries << "\r" << flush;
			tree -> GetEntry(i);
			std::map<std::string,TLeaf*> treeVars = initTreeVars(tree);

			TLorentzVector ph1, ph2;
			ph1.SetPtEtaPhiE(treeVars["dipho_leadPt"]->GetValue(0), treeVars["dipho_leadEta"]->GetValue(0), treeVars["dipho_leadPhi"]->GetValue(0), treeVars["dipho_leadEnergy"]->GetValue(0));
			ph2.SetPtEtaPhiE(treeVars["dipho_subleadPt"]->GetValue(0), treeVars["dipho_subleadEta"]->GetValue(0), treeVars["dipho_subleadPhi"]->GetValue(0), treeVars["dipho_subleadEnergy"]->GetValue(0));
			float pT = (ph1+ph2).Pt();
			pT_histo[label] -> Fill(pT, treeVars["weight"]->GetValue(0));
		}


		//Choose C1 to use
		TH1F* c;
		std::string proc = processes[label];
		cout << endl;

		if(C1_h.find(proc)!=C1_h.end())
			c = (TH1F*)C1_h[proc] -> Clone();

		else if(proc=="wh" && (C1_h.find("wmh")!=C1_h.end()) && (C1_h.find("wph")!=C1_h.end()) )
		{	c = (TH1F*)C1_h["wmh"]->Clone();
			c -> Add(C1_h["wph"], C1_h["wmh"], 0.333, 0.667);
		}
		else
		{
			c = new TH1F("c", "", 1, 0, 10000);
			c -> SetBinContent(1, 0);
			cout << "No C1 file found for the specified process, C1 will be assumed to be zero." << endl;
		}


		for(unsigned int bin=1; bin<ptBins[label].size(); ++bin)
		{	
			float tmp = 0.;
			float tmpNorm = 0.;
			for(int i=ptBins[label].at(bin-1)+1; i<ptBins[label].at(bin); ++i)
			{	
				tmp += pT_histo[label]->GetBinContent(i)*c->GetBinContent(c->FindBin(pT_histo[label]-> GetBinCenter(i)));

				tmpNorm += pT_histo[label]-> GetBinContent(i);
			}
			
			C1_Binned[label].push_back(tmp/tmpNorm);
		}

		cout << "Processed tag " << label << ", " << nentries << " events out of " << tree -> GetEntries() << " in the tree" << endl;
	}

	if(doPlots)
	{	TCanvas* c1 = new TCanvas();
		c1 -> cd();
		TLegend* leg = new TLegend(0.65, 0.70, 0.9, 0.85);

		for(auto it = C1_h.begin(); it != C1_h.end(); ++it)
		{
			std::string label = it -> first;
			TH1F* histo = it -> second;

			histo -> SetLineWidth(3);
			histo -> SetLineColor(std::distance(C1_h.begin(), it)+1);

			if(it==C1_h.begin())
				histo -> Draw("histo");
			else
				histo -> Draw("histo SAME");

			leg -> AddEntry(histo, (TString)label,  "l");
		}
		
		leg -> Draw("SAME");
		CMS_lumi(c1, 0, 0);

		c1 -> SaveAs("C1Tags.png");
		c1 -> SaveAs("C1Tags.pdf");

		TCanvas* c2 = new TCanvas();
		c2 -> cd();
		TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
		std::vector<float> maxima;

		for(auto it = pT_histo.begin(); it != pT_histo.end(); ++it)
		{	TH1F* histo = it -> second;
			maxima.push_back(histo->GetMaximum()/histo->Integral());
		}

		double maximum = *max_element(maxima.begin(), maxima.end());
		for(auto it = pT_histo.begin(); it != pT_histo.end(); ++it)
		{
			TH1F* histo = it -> second;
			std::string label = it -> first;
			histo -> SetLineWidth(3);
			histo -> SetLineColor(std::distance(pT_histo.begin(), it)+38);
			histo -> Rebin(10);

			if(it==pT_histo.begin())
			{	
				histo -> GetYaxis() -> SetRangeUser(0, 1.2*maximum);
				histo -> GetXaxis() -> SetRangeUser(0, 500);
				histo -> DrawNormalized("histo");
			}
			else
				histo -> DrawNormalized("histo SAME");
			leg -> AddEntry(histo, (TString)label,  "l");
		}

		CMS_lumi(c2, 0, 0);

		c2 -> SaveAs("PtTags.png");
		c2 -> SaveAs("PtTags.pdf");
		leg2 -> Draw("SAME");

		system(("mv *.png " + outputPlots).c_str());
		system(("mv *.pdf " + outputPlots).c_str());
	}

	if(numberOfExpectedTagNames==(int)nameOfTags.size() && outputFile!="")
	{
		writeMap(nameOfTags, C1_Binned, outputFile);
		cout << "C1 constant written in " << outputFile << endl;
	}

	else if (outputFile=="")
	{
		cerr << "No outputfile specified. No code snippet will be procuced, anyway you can read the value of C1 in the following printout" << endl;
	}
	else 
	{
		cerr << "Number of tags names does not correspond to the expected one. No code snippet will be procuced, anyway you can read the value of C1 in the following printout" << endl;
	}

	cout << endl << "|-------------------------------------------------------------------------------------------------|" << endl << endl << endl;


	for ( auto C1_it = C1_Binned.begin(), pt_it = ptBins.begin(); C1_it != C1_Binned.end(); ++C1_it , ++pt_it )
	{
		string label = C1_it -> first;
		std::vector<float> C1 = C1_it -> second;
		std::vector<float> pt = pt_it -> second;

		cout << "Tag: " << label << ";" << endl << "pT binning: ";
		for(unsigned int i=0; i<pt.size(); ++i)
		{	cout << pt.at(i);
			if(i!=pt.size()-1)
				cout  << ", ";
			else
				cout << ";";
		}

		cout << endl << "C1: ";
		for(unsigned int i=0; i<C1.size(); ++i)
		{	cout << C1.at(i);
			if(i!=C1.size()-1)
				cout  << ", ";
			else
				cout << ";";
		}

		cout << endl << endl;
	}

	cout << endl << "|-------------------------------------------------------------------------------------------------|" << endl << endl;

	cout << endl << endl;





	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosing our company. We hope to see you again on board of our aircrafts" << endl << endl;
}

























