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

void dimuon_macro_forEM_production_CB_exp()
{
    TH1::SetDefaultSumw2();
    // input Hyperloop result
    TFile *inputResults_SE = new TFile("AnalysisResults.root");
    if (inputResults_SE->IsZombie()) {
        std::cerr << "Error: AnalysisResults.root not found. Exiting." << std::endl;
        return;
    }

    // 1Dヒストグラムの作成
    Making1DfromHnSparse(inputResults_SE);
    TFile *outfile_Made1D = new TFile("Made1Dhist.root");

    // Like-sign法による背景推定
    LikeSignMethod(outfile_Made1D);
    TFile *outfile_LikeSignMethod = new TFile("LikeSignMethod.root");

    // Peak Fit (CrystalBall + Exponential)
    PeakFit_CrystalBall(outfile_LikeSignMethod);
    TFile *outfile_PeakFit = new TFile("PeakFit_CrystalBall_Results.root");

    // Yield計算
    YieldCalcuration_CrystalBall(outfile_PeakFit);
    TFile *outfile_Yield = new TFile("Yield_Results_CrystalBall.root");

    // 効率補正
    TFile *inputEfficiency = new TFile("AcceptanceEfficiency_Both.root");
    if (inputEfficiency->IsZombie()) {
        std::cout << "Warning: AcceptanceEfficiency_Both.root not found. Stopping analysis before correction." << std::endl;
        return;
    }

    correction(inputResults_SE, outfile_Yield, inputEfficiency, "CB");
}

int main()
{
    dimuon_macro_forEM_production_CB_exp();
    return 0;
}
