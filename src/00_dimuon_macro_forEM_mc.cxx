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

void dimuon_macro_forEM_mc()
{
    TH1::SetDefaultSumw2();
    // input Hyperloop result
    TFile *inputResults_mc = new TFile("AnalysisResults.root");
    if (inputResults_mc->IsZombie()) {
        std::cerr << "Error: AnalysisResults.root not found. Exiting." << std::endl;
        return;
    }
    Making1DfromHnSparse_mc(inputResults_mc);
    Check_mc_mass();
    CalAcceptance_Efficiency();
}

int main()
{
    dimuon_macro_forEM_mc();
    return 0;
}
