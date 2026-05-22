#include <TFile.h>
#include <TH1.h>
#include <iostream>
#include "myAnalysis.h"

int main() {
    TH1::SetDefaultSumw2();
    
    // input Hyperloop result
    TFile *inputResults_mc = new TFile("AnalysisResults.root", "READ");
    if (!inputResults_mc || inputResults_mc->IsZombie()) {
        std::cerr << "Error: AnalysisResults.root not found. Exiting." << std::endl;
        return 1;
    }

    std::cout << "Starting MC analysis for Mass Check..." << std::endl;
    
    // 1Dヒストグラムの作成 (Made1Dhist_mc.root)
    Making1DfromHnSparse_mc(inputResults_mc);
    
    // 質量分布のチェックとCanvasの作成 (check_mc_mass.root)
    Check_mc_mass();
    
    inputResults_mc->Close();
    delete inputResults_mc;
    
    std::cout << "Analysis finished successfully." << std::endl;
    return 0;
}
