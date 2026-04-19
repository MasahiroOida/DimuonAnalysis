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

void Peakfit(TH1 *Sig)
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
    totalfit->SetParameter(0, exact5);
    totalfit->SetParameter(1, exact6);
    totalfit->SetParLimits(1, 0.001, 10.0); // Force alpha > 0 in exp(-alpha*x) which means slope is negative

    for(int i=0; i<5; ++i) {
        totalfit->SetParameter(2+i, omegafitfunc->GetParameter(i));
        totalfit->SetParameter(7+i, phifitfunc->GetParameter(i));
    };

    // Apply limits for totalfit
    totalfit->SetParLimits(3, 0.75,0.85);
    totalfit->SetParLimits(8, 0.9,1.1);

    totalfit->SetParLimits(4, 0.01, 1.0);  // Omega sigma >= 0.01
    totalfit->SetParLimits(9, 0.01, 1.0);  // Phi sigma >= 0.01

    totalfit->SetParLimits(5, 0.5 , 3.0);  // Omega CB alpha > 1
    totalfit->SetParLimits(6, 1, 10.0);  // Omega CB n > 0
    totalfit->SetParLimits(10, 0.5, 3.0); // Phi CB alpha > 1
    totalfit->SetParLimits(11, 1, 10.0); // Phi CB n > 0


    TFitResultPtr fitRes = Sig->Fit("totalfit", "RQS");
    if (fitRes.Get()) fitRes->Write("fitResult");
    Sig->Write("LikeSignSig_fit");
    totalfit->Write("totalfit");

    delete BGHistogram; delete SigHistogram; delete fitFunction; delete BGFunction;
    delete omegafitfunc; delete phifitfunc; delete totalfit;
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
                    if(LikeSignSig) Peakfit(LikeSignSig);
                    TObject *p;
                    if((p = ptDir->Get("ptrange"))) p->Write("ptrange");
                    if((p = ptDir->Get("ptmin"))) p->Write("ptmin");
                    if((p = ptDir->Get("ptmax"))) p->Write("ptmax");
                }
            }
        }
    }
    outfile->Close();
    delete outfile;
    std::cout << "Finish Crystal Ball fitting with alpha > 1 and n > 0 constraints" << std::endl;
}
