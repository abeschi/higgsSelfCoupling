#include "interface/AnalysisUtils.h"
#include <assert.h> 


float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) + 
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}




void  writeMap (std::vector<std::string> names, std::map<std::string,std::vector<float>> C1, std::string outputFile)
{
	std::vector<float> C1_;
 
	for(auto it = C1.rbegin(); it != C1.rend(); ++it)
	{
		std::vector<float> tmp = it -> second;
		for(unsigned int i=0; i<tmp.size(); ++i)
		{	
			if(tmp.at(i)<0)
				C1_.push_back(0.);
			else
				C1_.push_back(tmp.at(i));
		}
	}

	assert(C1_.size()==names.size());

	std::ofstream myfile;
	myfile.open (outputFile);

	myfile << "C1Map = { \n";
	for(unsigned int i=0; i<names.size(); ++i)
		myfile << "        \"" << names.at(i) << "\":" << C1_.at(i) << ",\n";

	myfile << "}\n\n\n";

	myfile <<  "Processes = ["; 
	for(unsigned int i=0; i<names.size(); ++i)
	{	myfile << names.at(i);
		if(i!=names.size()-1)
			myfile << ", ";
		else
			myfile << "]\n";
	}

	myfile.close();

	return ;
}



void MakePlot(TH1F** histos, TString title)
{
  histos[0] -> SetLineWidth(3);			//ttH
  histos[0] -> SetLineColor(kRed + 1);
  histos[0] -> SetFillStyle(0);
  
  histos[1] -> SetLineWidth(3);			//ggH
  histos[1] -> SetLineColor(kGreen + 2);
  histos[1] -> SetFillStyle(0);
  
  histos[2] -> SetLineWidth(3);			//VBF
  histos[2] -> SetLineColor(kAzure);
  histos[2] -> SetFillStyle(0);
  
  histos[3] -> SetLineWidth(3);			//VH
  histos[3] -> SetLineColor(kViolet - 2);
  histos[3] -> SetFillStyle(0);
  
  histos[4] -> SetLineWidth(3);			//bbH
  histos[4] -> SetLineColor(kOrange);
  histos[4] -> SetFillStyle(0);
  
  histos[5] -> SetLineWidth(3);			//tHq
  histos[5] -> SetLineColor(kAzure + 8);
  histos[5] -> SetFillStyle(0);
  
  histos[6] -> SetLineWidth(3);			//tHW
  histos[6] -> SetLineColor(kViolet + 2);
  histos[6] -> SetFillStyle(0);
  
  
  
  histos[7] -> SetMarkerStyle(20);		//Data
  histos[7] -> SetMarkerSize(1);
  histos[7] -> SetMarkerColor(kBlack);
  histos[7] -> SetFillStyle(0);
  
  
  
  histos[8] -> SetLineWidth(1);			//Diphoton
  histos[8] -> SetFillColor(kAzure + 1);
  histos[8] -> SetFillStyle(1001);
  
  histos[9] -> SetLineWidth(1);			//Gamma + jets
  histos[9] -> SetFillStyle(1001);
  histos[9] -> SetFillColor(kYellow - 4);
  
  histos[10] -> SetLineWidth(1);			//QCD
  histos[10] -> SetFillColor(kTeal + 9);
  histos[10] -> SetFillStyle(1001);
  
  histos[11] -> SetLineWidth(1);			//ttGG
  histos[11] -> SetFillColor(kMagenta + 1);
  histos[11] -> SetFillStyle(1001);
  
  histos[12] -> SetLineWidth(1);			//ttGJets
  histos[12] -> SetFillColor(kMagenta + 1);
  histos[12] -> SetFillStyle(1001);
  
  histos[13] -> SetLineWidth(1);			//ttJets
  histos[13] -> SetFillColor(kMagenta + 1);
  histos[13] -> SetFillStyle(1001);
  
  
  TLegend* leg = new TLegend(0.65, 0.70, 0.9, 0.85);
  leg -> AddEntry(histos[0], "ttH", "l");
  leg -> AddEntry(histos[1], "ggH", "l");
  leg -> AddEntry(histos[2], "VBF", "l");
  leg -> AddEntry(histos[3], "VH", "l");
  leg -> AddEntry(histos[4], "bbH", "l");
  leg -> AddEntry(histos[5], "tHq", "l");
  leg -> AddEntry(histos[6], "tHW", "l");
  leg -> AddEntry(histos[7], "Data sidebands", "p");
  
  TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
  leg2 -> AddEntry(histos[7], "Data sidebands", "p");
  leg2 -> AddEntry(histos[8], "Diphotons", "f");
  leg2 -> AddEntry(histos[9], "Gamma + jets", "f");
  leg2 -> AddEntry(histos[10], "QCD", "f");
  leg2 -> AddEntry(histos[11], "ttGJets", "f");
  
  
  histos[11] -> Add(histos[12]);
  histos[11] -> Add(histos[13]);
  
  
  TCanvas* c = new TCanvas();
  c -> cd();
  
  float m = std::max(std::max(histos[0]->GetMaximum()/histos[0]->Integral(), histos[1]->GetMaximum()/histos[1]->Integral()), std::max(histos[2]->GetMaximum()/histos[2]->Integral(), histos[3]->GetMaximum()/histos[3]->Integral()));
  float m2 = std::max(std::max(histos[4]->GetMaximum()/histos[4]->Integral(), histos[5]->GetMaximum()/histos[5]->Integral()), std::max(histos[6]->GetMaximum()/histos[6]->Integral(), histos[7]->GetMaximum()/histos[7]->Integral()));
  m = std::max((double)m, (double)m2);
  
  TH1F* axis = new TH1F(*histos[0]);
  axis -> SetMarkerSize(0);
  axis -> SetLineWidth(0);
  axis -> GetYaxis() -> SetTitleOffset(1.5);
  axis -> GetYaxis() -> SetRangeUser(0, 1.1*m);
  
  axis -> Draw("histo");
  histos[1] -> DrawNormalized("histo SAME");
  histos[2] -> DrawNormalized("histo SAME");
  histos[3] -> DrawNormalized("histo SAME");
  histos[4] -> DrawNormalized("histo SAME");
  histos[5] -> DrawNormalized("histo SAME");
  histos[6] -> DrawNormalized("histo SAME");
  histos[0] -> DrawNormalized("histo SAME");
  histos[7] -> DrawNormalized("SAME E1");
  
  leg -> Draw("SAME");
  CMS_lumi(c, 0, 0);
  
  c -> SaveAs("c_" + title + "Signal.png");
  c -> SaveAs("c_" +title + "Signal.pdf");
  
  
  TCanvas* c2 = new TCanvas();
  c2 -> cd();
  
  histos[9] -> Add(histos[8]);
  histos[10] -> Add(histos[9]);
  histos[11] -> Add(histos[10]);
  
  
  double max2 = histos[7]->GetBinContent(histos[7]->GetMaximumBin());
  if(histos[11]->GetBinContent(histos[11]->GetMaximumBin()) > max2 )
    max2 = histos[11]->GetBinContent(histos[11]->GetMaximumBin());
  histos[11]-> GetYaxis() -> SetRangeUser(0, 1.1*max2);
  histos[11]-> GetYaxis() -> SetTitleSize(0.05);
  histos[11]-> GetYaxis() -> SetTitleFont(42);
  histos[11]-> GetYaxis() -> SetLabelSize(0.045);
  histos[11]-> GetYaxis() -> SetLabelFont(42);
  histos[11]-> GetXaxis() -> SetLabelSize(0);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.97);
  pad1->SetBottomMargin(0.035);
  pad1->Draw();
  pad1->cd();
  
  histos[11] -> Draw("histo");
  histos[10] -> Draw("histo SAME");
  histos[9] -> Draw("histo SAME");
  histos[8] -> Draw("histo SAME");
  histos[7] -> Draw("SAME E1");
  leg2 -> Draw("SAME");
  
  c2->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.10, 1, 0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->Draw();
  pad2->cd();
  
  TH1F *h = (TH1F*)histos[7]->Clone("h");
  h->SetLineColor(kBlack);
  h->SetMinimum(0.5);  // Define Y ..
  h->SetMaximum(1.5); // .. range
  h->Sumw2();
  h->SetStats(0);      // No statistics on lower plot
  h->Divide(histos[11]);
  h->SetMarkerStyle(21);
  h -> SetTitle("");
  h-> GetYaxis() -> SetTitle("Data/MC");
  
  
  // Y axis ratio plot settings
  h->GetYaxis()->SetNdivisions(-10);
  h->GetYaxis()->SetTitleSize(0.13);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.5);
  h->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
  h->GetYaxis()->SetLabelSize(0.12);
  
  // X axis ratio plot settings
  h->GetXaxis()->SetTitleSize(0.15);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
  h->GetXaxis()->SetLabelSize(0.12);
  
  h->Draw("EP");
  
  
  CMS_lumi(c2, 0, 10);
  
  c2 -> SaveAs("c_" +title + "Bkg.png");
  c2 -> SaveAs("c_" +title + "Bkg.pdf");
  
  return;
}



void MakePlot2(std::map<std::string,TH1F*>& histos, TString title)
{
  histos["ttH"] -> SetLineWidth(3);			//ttH
  histos["ttH"] -> SetLineColor(kRed + 1);
  histos["ttH"] -> SetFillStyle(0);
  
  histos["ggH"] -> SetLineWidth(3);			//ggH
  histos["ggH"] -> SetLineColor(kGreen + 2);
  histos["ggH"] -> SetFillStyle(0);
  
  histos["VBF"] -> SetLineWidth(3);			//VBF
  histos["VBF"] -> SetLineColor(kAzure);
  histos["VBF"] -> SetFillStyle(0);
  
  histos["VH"] -> SetLineWidth(3);			//VH
  histos["VH"] -> SetLineColor(kViolet - 2);
  histos["VH"] -> SetFillStyle(0);
  
  // histos["bbH"] -> SetLineWidth(3);			//bbH
  // histos["bbH"] -> SetLineColor(kOrange);
  // histos["bbH"] -> SetFillStyle(0);
  
  // histos["tHq"] -> SetLineWidth(3);			//tHq
  // histos["tHq"] -> SetLineColor(kAzure + 8);
  // histos["tHq"] -> SetFillStyle(0);
  
  // histos["tHW"] -> SetLineWidth(3);			//tHW
  // histos["tHW"] -> SetLineColor(kViolet + 2);
  // histos["tHW"] -> SetFillStyle(0);
  
  histos["data"] -> SetMarkerStyle(20);		//Data
  histos["data"] -> SetMarkerSize(1);
  histos["data"] -> SetMarkerColor(kBlack);
  histos["data"] -> SetFillStyle(0);
  
  /*
  histos[8] -> SetLineWidth(1);			//Diphoton
  histos[8] -> SetFillColor(kAzure + 1);
  histos[8] -> SetFillStyle(1001);
  
  histos[9] -> SetLineWidth(1);			//Gamma + jets
  histos[9] -> SetFillStyle(1001);
  histos[9] -> SetFillColor(kYellow - 4);
  
  histos[10] -> SetLineWidth(1);			//QCD
  histos[10] -> SetFillColor(kTeal + 9);
  histos[10] -> SetFillStyle(1001);
  
  histos[11] -> SetLineWidth(1);			//ttGG
  histos[11] -> SetFillColor(kMagenta + 1);
  histos[11] -> SetFillStyle(1001);
  
  histos[12] -> SetLineWidth(1);			//ttGJets
  histos[12] -> SetFillColor(kMagenta + 1);
  histos[12] -> SetFillStyle(1001);
  
  histos[13] -> SetLineWidth(1);			//ttJets
  histos[13] -> SetFillColor(kMagenta + 1);
  histos[13] -> SetFillStyle(1001);
  */
  
  histos[Form("CS_%s",title.Data())] -> SetLineWidth(3);			//CS
  histos[Form("CS_%s",title.Data())] -> SetLineColor(kAzure);
  histos[Form("CS_%s",title.Data())] -> SetFillStyle(1);
  histos[Form("CS_%s",title.Data())] -> SetFillColor(kAzure-9);
  
  
  TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
  leg2 -> AddEntry(histos["ttH"], "Signal", "l");
  leg2 -> AddEntry(histos["data"], "Data sidebands", "p");
  leg2 -> AddEntry(histos[Form("CS_%s",title.Data())], "Control Sample", "f");
  //	leg2 -> AddEntry(histos[9], "Gamma + jets", "f");
  //	leg2 -> AddEntry(histos[10], "QCD", "f");
  //	leg2 -> AddEntry(histos[11], "ttGJets", "f");
  
  
  // histos[11] -> Add(histos[12]);
  // histos[11] -> Add(histos[13]);
  
  TH1F* signal = new TH1F(*histos["ttH"]);
  signal -> Add(histos["ggH"]);
  signal -> Add(histos["VBF"]);
  signal -> Add(histos["VH"]);
  // signal -> Add(histos["bbH"]);
  // signal -> Add(histos["tHq"]);
  // signal -> Add(histos["tHW"]);
  
  std::cout << "Tag Purity: ttH " << std::setprecision(2) << histos["ttH"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(2) << histos["ttH"]->Integral() << ") events" << std::endl;
  std::cout << "Tag Purity: ggH " << std::setprecision(2) << histos["ggH"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["ggH"]->Integral() << ") events" << std::endl;
  std::cout << "Tag Purity: VBF " << std::setprecision(2) << histos["VBF"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["VBF"]->Integral() << ") events" << std::endl;
  std::cout << "Tag Purity: VH  " << std::setprecision(2) << histos["VH"] ->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["VH"] ->Integral() << ") events" << std::endl;
  // std::cout << "Tag Purity: bbH " << std::setprecision(2) << histos["bbH"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["bbH"]->Integral() << ") events" << std::endl;
  // std::cout << "Tag Purity: tHq " << std::setprecision(2) << histos["tHq"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["tHq"]->Integral() << ") events" << std::endl;
  // std::cout << "Tag Purity: tHW " << std::setprecision(2) << histos["tHW"]->Integral()/signal->Integral()*100 << "%, (" << std::setprecision(4) << histos["tHW"]->Integral() << ") events" << std::endl;
  
  // histos[9] -> Add(histos[8]);
  // histos[10] -> Add(histos[9]);
  // histos[11] -> Add(histos[10]);
  
  
  TCanvas* c = new TCanvas();
  c -> cd();
  
  TH1F* axis = new TH1F(*histos["ttH"]);
  axis -> SetMarkerSize(0);
  axis -> SetLineWidth(0);
  axis -> GetYaxis() -> SetTitleOffset(1.5);
  axis -> GetYaxis() -> SetRangeUser(0, 1.1*histos["data"]->GetMaximum());
  
  axis -> Draw("histo");
  
  //	histos[11] -> Draw("histo SAME");
  //	histos[10] -> Draw("histo SAME");
  //	histos[9] -> Draw("histo SAME");
  //	histos[8] -> Draw("histo SAME");
  histos[Form("CS_%s",title.Data())] -> Draw("histo,same");
  histos["data"] -> Draw("SAME E1");
  signal -> Draw("histo SAME");
  leg2 -> Draw("SAME");
  
  CMS_lumi(c, 0, 0);
  
  c -> SaveAs("c_" +title + "Signal.png");
  c -> SaveAs("c_" +title + "Signal.pdf");
  
  return;
}


