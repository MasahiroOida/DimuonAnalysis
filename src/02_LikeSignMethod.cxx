#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TObject.h>
#include <TStyle.h>
#include <TLegend.h>
#include <cmath>
#include "TMath.h"
#include <iostream>
#include "myAnalysis.h"
#include <TKey.h>
#include <TParameter.h>

const double EPSILON = 1e-7;

void LikeSignMethod_cal(TH1 *SEPM, TH1 *SEPP, TH1 *SEMM, TH1 *MEPM, TH1 *MEPP, TH1 *MEMM)
{
    if (!SEPM || !SEPP || !SEMM) return;
    TH1::SetDefaultSumw2();
    bool noRfactor = (MEPM == nullptr);
    
    TAxis *xAxis = SEPM->GetXaxis();
    int numBins = xAxis->GetNbins();
    TH1D* BG = (TH1D*)SEPM->Clone("BG");
    TH1D* Sig = (TH1D*)SEPM->Clone("Sig");
    TH1D* SEBG = (TH1D*)SEPM->Clone("SEBG");
    TH1D* Rfactor = (TH1D*)SEPM->Clone("Rfactor");
    BG->Reset(); Sig->Reset(); SEBG->Reset(); Rfactor->Reset();

    for (int bin = 1; bin <= numBins; ++bin)
    {
        double SEPMValue = SEPM->GetBinContent(bin);
        double SEPPValue = SEPP->GetBinContent(bin);
        double SEMMValue = SEMM->GetBinContent(bin);
        double SEPMValue_err = SEPM->GetBinError(bin);
        double SEPPValue_err = SEPP->GetBinError(bin);
        double SEMMValue_err = SEMM->GetBinError(bin);

        double MEPMValue = 0, MEPPValue = 0, MEMMValue = 0;
        double MEPMValue_err = 0, MEPPValue_err = 0, MEMMValue_err = 0;

        if (!noRfactor)
        {
            MEPMValue = MEPM->GetBinContent(bin);
            MEPPValue = MEPP->GetBinContent(bin);
            MEMMValue = MEMM->GetBinContent(bin);
            MEPMValue_err = MEPM->GetBinError(bin);
            MEPPValue_err = MEPP->GetBinError(bin);
            MEMMValue_err = MEMM->GetBinError(bin);
        }

        double SEsqrtValue = TMath::Sqrt(SEPPValue * SEMMValue);
        double RfactorValue = 1.0;
        
        if (!noRfactor)
        {
            double MEsqrtValue = TMath::Sqrt(MEPPValue * MEMMValue);
            RfactorValue = 0.5 * MEPMValue / std::max(MEsqrtValue, EPSILON);
        }
        
        double BGValue = 2 * RfactorValue * SEsqrtValue;
        double SigValue = SEPMValue - BGValue;

        // SEsqrt Error
        double SEsqrtValue_err = 0.5 * TMath::Sqrt( 
            (SEMMValue / std::max(SEPPValue, EPSILON)) * pow(SEPPValue_err, 2) + 
            (SEPPValue / std::max(SEMMValue, EPSILON)) * pow(SEMMValue_err, 2) 
        );

        // Rfactor Error
        double RfactorValue_err = 0;
        if (!noRfactor) {
            double term1 = pow(MEPMValue_err / std::max(MEPMValue, EPSILON), 2);
            double term2 = 0.25 * pow(MEPPValue_err / std::max(MEPPValue, EPSILON), 2);
            double term3 = 0.25 * pow(MEMMValue_err / std::max(MEMMValue, EPSILON), 2);
            RfactorValue_err = RfactorValue * TMath::Sqrt(term1 + term2 + term3);
        }

        double BGValue_err = BGValue * TMath::Sqrt(
            pow(RfactorValue_err / std::max(RfactorValue, EPSILON), 2) + 
            pow(SEsqrtValue_err / std::max(SEsqrtValue, EPSILON), 2)
        );

        double SigValue_err = TMath::Sqrt(pow(SEPMValue_err, 2) + pow(BGValue_err, 2));

        Sig->SetBinContent(bin, SigValue); Sig->SetBinError(bin, SigValue_err);
        BG->SetBinContent(bin, BGValue); BG->SetBinError(bin, BGValue_err);
        SEBG->SetBinContent(bin, SEsqrtValue); SEBG->SetBinError(bin, SEsqrtValue_err);
        Rfactor->SetBinContent(bin, RfactorValue); Rfactor->SetBinError(bin, RfactorValue_err);
    }

    SEPM->Write("SEPM"); SEPP->Write("SEPP"); SEMM->Write("SEMM"); SEBG->Write("SEBG");
    if (!noRfactor) { MEPM->Write("MEPM"); MEPP->Write("MEPP"); MEMM->Write("MEMM"); }
    Rfactor->Write("Rfactor"); BG->Write("LikeSignBG"); Sig->Write("LikeSignSig");

    // Create a canvas to overlay SEPM, BG, and Sig
    TCanvas *c1 = new TCanvas("LikeSign_Canvas", "LikeSign Method Results", 800, 600);
    c1->cd();
    
    SEPM->SetLineColor(kBlack);
    SEPM->SetMarkerColor(kBlack);
    SEPM->SetMarkerStyle(20);
    SEPM->SetMarkerSize(0.7);
    SEPM->GetXaxis()->SetRangeUser(0, 1.5);
    SEPM->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    SEPM->GetYaxis()->SetTitle("Counts");
    SEPM->Draw("E");

    BG->SetLineColor(kBlue);
    BG->SetMarkerColor(kBlue);
    BG->SetMarkerStyle(20);
    BG->SetMarkerSize(0.7);
    BG->Draw("E SAME");

    Sig->SetLineColor(kRed);
    Sig->SetMarkerColor(kRed);
    Sig->SetMarkerStyle(20);
    Sig->SetMarkerSize(0.7);
    Sig->Draw("E SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(SEPM, "SEPM", "lep");
    leg->AddEntry(BG, "Background (LS)", "lep");
    leg->AddEntry(Sig, "Signal (SEPM - BG)", "lep");
    leg->Draw();

    c1->Write();
    delete c1;

    // Create a canvas for Rfactor
    TCanvas *c2 = new TCanvas("Rfactor_Canvas", "R-factor", 800, 600);
    c2->cd();
    c2->SetGrid();
    
    Rfactor->SetLineColor(kBlue);
    Rfactor->SetMarkerColor(kBlue);
    Rfactor->SetMarkerStyle(20);
    Rfactor->SetMarkerSize(0.7);
    Rfactor->GetXaxis()->SetRangeUser(0, 1.5);
    Rfactor->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    Rfactor->GetYaxis()->SetTitle("Rfactor");
    Rfactor->Draw("E");
    
    c2->Write();
    delete c2;
    
    delete SEBG; delete Rfactor; delete BG; delete Sig;
}

void ProcessPtDirectory(TDirectory *inputDir, TDirectory *outputDir) {
    outputDir->cd();
    TH1F *SEPM = (TH1F *)inputDir->Get("SEPM");
    TH1F *SEPP = (TH1F *)inputDir->Get("SEPP");
    TH1F *SEMM = (TH1F *)inputDir->Get("SEMM");
    TObject *p;
    if((p = inputDir->Get("ptrange"))) p->Write("ptrange");
    if((p = inputDir->Get("ptmin"))) p->Write("ptmin");
    if((p = inputDir->Get("ptmax"))) p->Write("ptmax");
    if((p = inputDir->Get("nEvents"))) p->Write("nEvents");
    TH1F *MEPM = (TH1F *)inputDir->Get("MEPM");
    TH1F *MEPP = (TH1F *)inputDir->Get("MEPP");
    TH1F *MEMM = (TH1F *)inputDir->Get("MEMM");
    LikeSignMethod_cal(SEPM, SEPP, SEMM, MEPM, MEPP, MEMM);
}

void LikeSignMethod(TFile *input_SeparetePt) {
    TFile *outfile = new TFile("LikeSignMethod.root", "RECREATE");
    TIter next(input_SeparetePt->GetListOfKeys());
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
                    ProcessPtDirectory(ptDir, outPtDir);
                }
            }
        }
    }
    outfile->Close();
    delete outfile;
}
