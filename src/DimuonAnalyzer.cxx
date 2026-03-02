#include "DimuonAnalyzer.h"
#include <TMath.h>
#include <iostream>

// 2D to 1D projection
TH1D* DimuonAnalyzer::ProjectTo1D(TH2* hist2D, double ptMin, double ptMax, const std::string& name) {
    TAxis *yAxis = hist2D->GetYaxis();
    int binMin = yAxis->FindBin(ptMin);
    int binMax = yAxis->FindBin(ptMax) - 1; // Be careful with the upper bin limit
    TH1D* hist1D = hist2D->ProjectionX(name.c_str(), binMin, binMax);
    hist1D->Sumw2();
    return hist1D;
}

// Calculate Like-Sign background
std::tuple<TH1D*, TH1D*> DimuonAnalyzer::CalculateLikeSign(
    TH1D* SEPM, TH1D* SEPP, TH1D* SEMM, TH1D* MEPM, TH1D* MEPP, TH1D* MEMM, const std::string& suffix) {
    
    TH1D* histBG = (TH1D*)SEPM->Clone(("LikeSignBG_" + suffix).c_str());
    TH1D* histSig = (TH1D*)SEPM->Clone(("LikeSignSig_" + suffix).c_str());
    histBG->Reset();
    histSig->Reset();

    bool useMixed = (MEPM != nullptr && MEPP != nullptr && MEMM != nullptr);

    for (int bin = 1; bin <= SEPM->GetNbinsX(); ++bin) {
        double N_se_pp = SEPP->GetBinContent(bin);
        double N_se_mm = SEMM->GetBinContent(bin);
        double Rfactor = 1.0;

        if (useMixed) {
            double N_me_pm = MEPM->GetBinContent(bin);
            double N_me_pp = MEPP->GetBinContent(bin);
            double N_me_mm = MEMM->GetBinContent(bin);
            if (N_me_pp > 0 && N_me_mm > 0) {
                Rfactor = 0.5 * N_me_pm / TMath::Sqrt(N_me_pp * N_me_mm);
            }
        }
        
        // N_BG = 2 * R * sqrt(N_same_pp * N_same_mm)
        double bgVal = 2.0 * Rfactor * TMath::Sqrt(N_se_pp * N_se_mm);
        double sigVal = SEPM->GetBinContent(bin) - bgVal;

        // Calculate errors here (same as original code)
        histBG->SetBinContent(bin, bgVal);
        histSig->SetBinContent(bin, sigVal);
    }
    return {histBG, histSig};
}

// Peak fitting
std::tuple<TF1*, TF1*, TF1*, TF1*> DimuonAnalyzer::PerformPeakFit(TH1D* sigHist, const std::string& suffix) {
    // 1. Define and fit background function
    TF1* bgFunc = new TF1(("BGFit_" + suffix).c_str(), "[0]*exp(-[1]*x)", 0.5, 1.3);
    // Implement original discontinuous fit logic here
    
    // 2. Fit Omega and Phi
    TF1* totalFit = new TF1(("TotalFit_" + suffix).c_str(), 
        "[0]*exp(-[1]*x)+[2]*exp(-(x-[3])*(x-[3])/[4]/[4]/2)+[5]*exp(-(x-[6])*(x-[6])/[7]/[7]/2)", 0.5, 1.3);
    
    // Initial parameter settings
    // totalFit->SetParameters(...);
    sigHist->Fit(totalFit, "RQS");

    // 3. Create and return individual functions
    TF1* omegaFunc = new TF1(("OmegaFit_" + suffix).c_str(), "gaus", 0, 2);
    // omegaFunc->SetParameters(...);
    
    TF1* phiFunc = new TF1(("PhiFit_" + suffix).c_str(), "gaus", 0, 2);
    // phiFunc->SetParameters(...);

    return {totalFit, bgFunc, omegaFunc, phiFunc};
}
