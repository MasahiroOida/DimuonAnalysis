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

void dimuon_macro_forEM_MFTnCluster()
{
    TH1::SetDefaultSumw2();
    
    TFile *inputResults_SE = new TFile("AnalysisResults.root");
    if (!inputResults_SE || inputResults_SE->IsZombie()) {
        std::cerr << "Error: Could not open AnalysisResults.root" << std::endl;
        return;
    }

    // Step 1: Make 1D histograms from HnSparse
    // (This processes all directories in AnalysisResults.root)
    Making1DfromHnSparse(inputResults_SE);
    
    // Step 2: Like-Sign Method
    TFile *outfile_Made1D = new TFile("Made1Dhist.root");
    if (outfile_Made1D && !outfile_Made1D->IsZombie()) {
        LikeSignMethod(outfile_Made1D);
        outfile_Made1D->Close();
    }
    
    // Step 3: Peak Fitting with Crystal Ball
    TFile *outfile_LikeSignMethod = new TFile("LikeSignMethod.root");
    if (outfile_LikeSignMethod && !outfile_LikeSignMethod->IsZombie()) {
        PeakFit_CrystalBall(outfile_LikeSignMethod);
        outfile_LikeSignMethod->Close();
    }
    
    // Step 4: Yield Calculation and Significance (MFTnCluster specific)
    TFile *outfile_PeakFit = new TFile("PeakFit_CrystalBall_Results.root");
    if (outfile_PeakFit && !outfile_PeakFit->IsZombie()) {
        YieldCalcuration_CrystalBall(outfile_PeakFit);
        Significance_MFTnCluster(outfile_PeakFit); // Use the new MFTnCluster significance script
        outfile_PeakFit->Close();
    }
}

int main()
{
    dimuon_macro_forEM_MFTnCluster();
    return 0;
}
