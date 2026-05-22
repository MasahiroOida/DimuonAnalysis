#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TObject.h>
#include <TStyle.h>
#include <TLegend.h>
#include <cmath>
#include "TMath.h"
#include <iostream>
#include <sstream>
#include <TStopwatch.h>
#include "myAnalysis.h"

void dimuon_macro_forEM_minpt_pol4()
{
    TH1::SetDefaultSumw2();
    
    TFile *inputResults_SE = new TFile("AnalysisResults.root");
    if (!inputResults_SE || inputResults_SE->IsZombie()) {
        std::cerr << "Error: Could not open AnalysisResults.root" << std::endl;
        return;
    }

    // Step 1: Make 1D histograms from HnSparse
    Making1DfromHnSparse(inputResults_SE);
    
    // Step 2: Like-Sign Method
    TFile *outfile_Made1D = new TFile("Made1Dhist.root");
    if (outfile_Made1D && !outfile_Made1D->IsZombie()) {
        LikeSignMethod(outfile_Made1D);
        outfile_Made1D->Close();
    }
    
    // Step 3: Peak Fitting with Crystal Ball (pol4 background)
    TFile *outfile_LikeSignMethod = new TFile("LikeSignMethod.root");
    if (outfile_LikeSignMethod && !outfile_LikeSignMethod->IsZombie()) {
        PeakFit_CrystalBall_pol4(outfile_LikeSignMethod);
        outfile_LikeSignMethod->Close();
    }
    
    // Step 4: Yield Calculation and Significance (minpt specific)
    TFile *outfile_PeakFit = new TFile("PeakFit_CrystalBall_pol4_Results.root");
    if (outfile_PeakFit && !outfile_PeakFit->IsZombie()) {
        YieldCalcuration_CrystalBall_pol4(outfile_PeakFit);
        Significance_CrystalBall_pol4_minpt(outfile_PeakFit);
        outfile_PeakFit->Close();
    }

    // Step 5: Scan (Optional, can be run separately)
    // _06_dR_minpt_Scan_2D_pol4();
}

int main()
{
    dimuon_macro_forEM_minpt_pol4();
    return 0;
}
