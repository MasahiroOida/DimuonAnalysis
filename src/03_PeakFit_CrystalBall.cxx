// include
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TObject.h>
#include <TStyle.h>
#include <TLegend.h>
#include <cmath>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <TKey.h>
#include <TParameter.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

Double_t myFreeFunction(Double_t *x, Double_t *par)
{
    if (0.5 <= x[0] && x[0] < 0.70) return par[0] * exp(-par[1] * x[0]);
    else if (0.70 <= x[0] && x[0] <= 0.86) return 0;
    else if (0.86 < x[0] && x[0] < 0.92) return par[0] * exp(-par[1] * x[0]);
    else if (0.92 <= x[0] && x[0] <= 1.15) return 0;
    else if (1.15 < x[0] && x[0] <= 1.3) return par[0] * exp(-par[1] * x[0]);
    return 0;
};

void Peakfit(TH1 *Sig, double nEvents = -1)
{
    if (Sig == nullptr) return;
    TH1::SetDefaultSumw2();
    
    TH1F *BGHistogram = (TH1F *)Sig->Clone("BGHistogram");
    TH1F *SigHistogram = (TH1F *)Sig->Clone("SigHistogram");
    
    double splitLow1 = 0.70; double splitHigh1 = 0.86;
    double splitLow2 = 0.92; double splitHigh2 = 1.15;

    for (int bin = 1; bin <= BGHistogram->GetNbinsX(); ++bin) {
        double binCenter = BGHistogram->GetBinCenter(bin);
        if ((binCenter >= splitLow1 && binCenter <= splitHigh1) || (binCenter >= splitLow2 && binCenter <= splitHigh2))
            BGHistogram->SetBinContent(bin, 0);
    }

    TF1 *fitFunction = new TF1("fitFunction", myFreeFunction, 0.6, 1.4, 2);
    fitFunction->SetParameters(500, 0.7);
    BGHistogram->Fit("fitFunction", "RQ");
    
    double exact5 = fitFunction->GetParameter(0);
    double exact6 = fitFunction->GetParameter(1);
    TF1 *BGFunction = new TF1("BGFunction", "[0]*exp(-[1]*x)", 0.6, 1.4);
    BGFunction->SetParameters(exact5, exact6);

    for (int bin = 1; bin <= Sig->GetNbinsX(); ++bin) {
        SigHistogram->SetBinContent(bin, Sig->GetBinContent(bin) - BGFunction->Eval(Sig->GetBinCenter(bin)));
    }

    // Individual Peak Fit Constraints: Alpha [3], n [4]
   // gStyle->SetOptStat(0);
   // gStyle->SetOptFit(0);

    // ... (rest of the initial code)
    TF1 *omegafitfunc = new TF1("omegafitfunc", "crystalball", splitLow1, splitHigh1);
    omegafitfunc->SetParameters(exact5 / 10, 0.78, 0.02, 1.5, 2.0);
    omegafitfunc->SetParLimits(1, 0.75, 0.85);
    omegafitfunc->SetParLimits(2, 0.01, 1.0); // Sigma >= 0.01
    omegafitfunc->SetParLimits(3, 0.5, 3.0);  // alpha > 1
    omegafitfunc->SetParLimits(4, 1, 10.0); // n > 0
    SigHistogram->Fit("omegafitfunc", "RQ");

    TF1 *phifitfunc = new TF1("phifitfunc", "crystalball", splitLow2, splitHigh2);
    phifitfunc->SetParameters(exact5 / 10, 1.02, 0.02, 1.5, 2.0);
    phifitfunc->SetParLimits(1, 0.9, 1.1);
    phifitfunc->SetParLimits(2, 0.01, 1.0); // Sigma >= 0.01
    phifitfunc->SetParLimits(3, 0.5, 3.0);  // alpha > 1
    phifitfunc->SetParLimits(4, 1, 10.0); // n > 0
    SigHistogram->Fit("phifitfunc", "RQ");

    // Total Fit Mapping:
    // BG: [0, 1]
    // Omega: [2, 3, 4, 5, 6] (Sigma is [4])
    // Phi: [7, 8, 9, 10, 11] (Sigma is [9])
    TF1 *totalfit = new TF1("totalfit", "[0]*exp(-[1]*x) + crystalball(2) + crystalball(7)", 0.6, 1.2);
    
    // Set initial values from sub-fits first
    totalfit->SetParameter(0, exact5);
    totalfit->SetParameter(1, exact6);
    for(int i=0; i<5; ++i) {
        totalfit->SetParameter(2+i, omegafitfunc->GetParameter(i));
        totalfit->SetParameter(7+i, phifitfunc->GetParameter(i));
    }

    // Now apply limits. If the current value is outside, ROOT will still warn, 
    // so we clamp the values to be safe.
    auto clamp = [](double val, double low, double high) {
        if (val <= low) return low + 0.0001;
        if (val >= high) return high - 0.0001;
        return val;
    };

    totalfit->SetParameter(3, clamp(totalfit->GetParameter(3), 0.75, 0.85));
    totalfit->SetParLimits(3, 0.75, 0.85);
    totalfit->SetParameter(8, clamp(totalfit->GetParameter(8), 0.9, 1.1));
    totalfit->SetParLimits(8, 0.9, 1.1);

    totalfit->SetParameter(4, clamp(totalfit->GetParameter(4), 0.01, 1.0));
    totalfit->SetParLimits(4, 0.01, 1.0);
    totalfit->SetParameter(9, clamp(totalfit->GetParameter(9), 0.01, 1.0));
    totalfit->SetParLimits(9, 0.01, 1.0);

    totalfit->SetParameter(5, clamp(totalfit->GetParameter(5), 0.5, 3.0));
    totalfit->SetParLimits(5, 0.5, 3.0);
    totalfit->SetParameter(6, clamp(totalfit->GetParameter(6), 1.0, 10.0));
    totalfit->SetParLimits(6, 1, 10.0);
    totalfit->SetParameter(10, clamp(totalfit->GetParameter(10), 0.5, 3.0));
    totalfit->SetParLimits(10, 0.5, 3.0);
    totalfit->SetParameter(11, clamp(totalfit->GetParameter(11), 1.0, 10.0));
    totalfit->SetParLimits(11, 1, 10.0);


    TFitResultPtr fitRes = Sig->Fit("totalfit", "RQS");
    if (fitRes.Get()) fitRes->Write("fitResult");
    
    // Save Peak-only histogram by subtracting final BG (raw counts)
    TH1F *SigOnly = (TH1F *)Sig->Clone("LikeSignSig_PeakOnly");
    SigOnly->SetTitle("Signal after subtracting final BG");
    SigOnly->GetListOfFunctions()->Clear(); // Remove any associated fit functions
    TF1 *fBG_final = new TF1("fBG_final", "[0]*exp(-[1]*x)", 0.6, 1.4);
    fBG_final->SetParameters(totalfit->GetParameter(0), totalfit->GetParameter(1));
    for (int bin = 1; bin <= SigOnly->GetNbinsX(); ++bin) {
        double binCenter = SigOnly->GetBinCenter(bin);
        SigOnly->SetBinContent(bin, SigOnly->GetBinContent(bin) - fBG_final->Eval(binCenter));
    }
    
    // Keep raw counts (do not normalize by event count and bin width)
    Sig->SetYTitle("Counts / (Mass bin)");
    SigOnly->SetYTitle("Counts / (Mass bin)");

    Sig->Write("LikeSignSig_fit");    totalfit->Write("totalfit");
    SigOnly->Write("LikeSignSig_PeakOnly");

    // Create a canvas to overlay Sig and PeakOnly (both already normalized)
    TCanvas *cOverlay = new TCanvas("Signal_Overlay_Canvas", "Overlay of LS Signal and PeakOnly", 800, 600);
    cOverlay->cd();
    
    Sig->SetLineColor(kBlack);
    Sig->SetMarkerColor(kBlack);
    Sig->SetMarkerStyle(20);
    Sig->SetMarkerSize(0.7);
    Sig->SetTitle("Signal Comparison;Mass (GeV/c^{2});dN/dm");
    Sig->GetXaxis()->SetRangeUser(0.6, 1.4);
    Sig->Draw("E");

    SigOnly->SetLineColor(kRed);
    SigOnly->SetMarkerColor(kRed);
    SigOnly->SetMarkerStyle(20);
    SigOnly->SetMarkerSize(0.7);
    SigOnly->Draw("E SAME");

    // --- Draw Fit Components ---
    TF1 *fBG_draw = new TF1("fBG_draw", "[0]*exp(-[1]*x)", 0.6, 1.2);
    fBG_draw->SetParameters(totalfit->GetParameter(0), totalfit->GetParameter(1));
    fBG_draw->SetLineColor(kGreen+2); fBG_draw->SetLineStyle(2); fBG_draw->SetLineWidth(2);

    TF1 *fOmega_draw = new TF1("fOmega_draw", "crystalball", 0.6, 1.2);
    for(int i=0; i<5; ++i) fOmega_draw->SetParameter(i, totalfit->GetParameter(2+i));
    fOmega_draw->SetLineColor(kBlue); fOmega_draw->SetLineStyle(2); fOmega_draw->SetLineWidth(2);

    TF1 *fPhi_draw = new TF1("fPhi_draw", "crystalball", 0.6, 1.2);
    for(int i=0; i<5; ++i) fPhi_draw->SetParameter(i, totalfit->GetParameter(7+i));
    fPhi_draw->SetLineColor(kRed); fPhi_draw->SetLineStyle(2); fPhi_draw->SetLineWidth(2);

    fBG_draw->Draw("SAME"); fOmega_draw->Draw("SAME"); fPhi_draw->Draw("SAME");
    totalfit->SetLineColor(kMagenta); totalfit->SetLineWidth(3);
    totalfit->Draw("SAME");

    TLegend *legOverlay = new TLegend(0.6, 0.6, 0.88, 0.88);
    legOverlay->AddEntry(Sig, "Data (incl. BG)", "lep");
    legOverlay->AddEntry(SigOnly, "Data (BG sub)", "lep");
    legOverlay->AddEntry(totalfit, "Total Fit", "l");
    legOverlay->AddEntry(fOmega_draw, "#omega Signal (CB)", "l");
    legOverlay->AddEntry(fPhi_draw, "#phi Signal (CB)", "l");
    legOverlay->AddEntry(fBG_draw, "Background (Exp)", "l");
    legOverlay->Draw();

    cOverlay->Write();
    delete cOverlay; delete legOverlay; delete fBG_draw; delete fOmega_draw; delete fPhi_draw;

    delete BGHistogram; delete SigHistogram; delete fitFunction; delete BGFunction;
    delete omegafitfunc; delete phifitfunc; delete totalfit; delete SigOnly; delete fBG_final;
}

void PeakFit_CrystalBall(TFile *input_Sig){
    TFile *outfile = new TFile("PeakFit_CrystalBall_Results.root", "RECREATE");
    TIter next(input_Sig->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (TString(key->GetClassName()) == "TDirectoryFile") {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TDirectory *outTopDir = outfile->mkdir(topDir->GetName());
            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt())) {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile") {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    TDirectory *outPtDir = outTopDir->mkdir(ptDir->GetName());
                    outPtDir->cd();
                    TH1F *LikeSignSig = (TH1F *)ptDir->Get("LikeSignSig");
                    TParameter<double> *pN = (TParameter<double> *)ptDir->Get("nEvents");
                    double nEvents = (pN) ? pN->GetVal() : -1;
                    if(LikeSignSig) Peakfit(LikeSignSig, nEvents);
                    TObject *p;
                    if((p = ptDir->Get("ptrange"))) p->Write("ptrange");
                    if((p = ptDir->Get("ptmin"))) p->Write("ptmin");
                    if((p = ptDir->Get("ptmax"))) p->Write("ptmax");
                    if((p = ptDir->Get("nEvents"))) p->Write("nEvents");
                }
            }
        }
    }
    outfile->Close();
    delete outfile;
    std::cout << "Finish Crystal Ball fitting with alpha > 1 and n > 0 constraints" << std::endl;
}
