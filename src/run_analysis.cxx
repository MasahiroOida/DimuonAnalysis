#include <TFile.h>
#include <TH2D.h>
#include "DimuonAnalyzer.h"

int main() {
    // Assuming 2D histograms are already created
    TFile* inFile = new TFile("Made2Dhist.root", "READ"); 
    TFile* outFile = new TFile("Processed_Histograms.root", "RECREATE");

    TH2D* h2_uls = (TH2D*)inFile->Get("Mass_pt_same_uls");
    TH2D* h2_lspp = (TH2D*)inFile->Get("Mass_pt_same_lspp");
    TH2D* h2_lsmm = (TH2D*)inFile->Get("Mass_pt_same_lsmm");
    // Get mixed events here as well...

    DimuonAnalyzer analyzer;

    // Example: analyze 0.0 < pT < 1.0 GeV/c
    double ptMin = 0.0, ptMax = 1.0;
    std::string ptName = "Pt_0.0_1.0";
    
    outFile->mkdir(ptName.c_str())->cd();
    
    // 1. Project to 1D
    TH1D* h1_uls = analyzer.ProjectTo1D(h2_uls, ptMin, ptMax, "SEPM");
    TH1D* h1_lspp = analyzer.ProjectTo1D(h2_lspp, ptMin, ptMax, "SEPP");
    TH1D* h1_lsmm = analyzer.ProjectTo1D(h2_lsmm, ptMin, ptMax, "SEMM");

    // 2. Apply Like-Sign method
    auto [histBG, histSig] = analyzer.CalculateLikeSign(h1_uls, h1_lspp, h1_lsmm, nullptr, nullptr, nullptr, ptName);

    // 3. Apply peak fitting
    auto [totalFit, bgFit, omegaFit, phiFit] = analyzer.PerformPeakFit(histSig, ptName);

    // 4. Save calculation results to file
    h1_uls->Write();
    histBG->Write();
    histSig->Write();
    totalFit->Write();
    bgFit->Write();
    omegaFit->Write();
    phiFit->Write();

    outFile->Close();
    inFile->Close();
    return 0;
}
