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
#include <cmath> // TMathとは別
#include "TMath.h"
#include <iostream>
#include <sstream>
#include <TStopwatch.h>
#include "myAnalysis.h"

void dimuon_macro_forEM_production()
{
    TH1::SetDefaultSumw2();
    // input Hyperloop result
    TFile *inputResults_SE = new TFile("AnalysisResults.root");
    if (inputResults_SE->IsZombie()) {
        std::cerr << "Error: AnalysisResults.root not found. Exiting." << std::endl;
        return;
    }

    Making1DfromHnSparse(inputResults_SE);
    TFile *outfile_Made1D = new TFile("Made1Dhist.root");

    LikeSignMethod(outfile_Made1D);
    TFile *outfile_LikeSignMethod = new TFile("LikeSignMethod.root");

    PeakFit_CrystalBall_pol4(outfile_LikeSignMethod);
    TFile *outfile_PeakFit = new TFile("PeakFit_CrystalBall_pol4_Results.root");

    YieldCalcuration_CrystalBall_pol4(outfile_PeakFit);
    TFile *outfile_Yield = new TFile("Yield_Results_CrystalBall_pol4.root");

    TFile *inputEfficiency = new TFile("AcceptanceEfficiency_Both.root");
    if (inputEfficiency->IsZombie()) {
        std::cout << "Warning: AcceptanceEfficiency_Both.root not found. Stopping analysis before correction." << std::endl;
        return;
    }

    correction(inputResults_SE, outfile_Yield, inputEfficiency, "CB");
}

int main()
{
    dimuon_macro_forEM_production();
    return 0;
}
