#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <TString.h>

void CalculateRatio1D(TH1D* hNum, TH1D* hDen, const char* name, const char* title) {
    if (!hNum || !hDen) {
        std::cout << "Error: Histogram not found for " << name << std::endl;
        return;
    }

    TH1D* hRatio = (TH1D*)hNum->Clone(name);
    hRatio->SetTitle(title);
    hRatio->Reset();

    int nBins = hNum->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double num = hNum->GetBinContent(i);
        double den = hDen->GetBinContent(i);
        
        if (den > 0) {
            double ratio = num / den;
            if (ratio > 1.0) ratio = 1.0;
            double err = TMath::Sqrt(ratio * (1.0 - ratio) / den);
            
            hRatio->SetBinContent(i, ratio);
            hRatio->SetBinError(i, err);
        } else {
            hRatio->SetBinContent(i, 0);
            hRatio->SetBinError(i, 0);
        }
    }
    hRatio->Write();
}

void CalculateRatio2D(TH2D* hNum, TH2D* hDen, const char* name, const char* title) {
    if (!hNum || !hDen) {
        std::cout << "Error: 2D Histogram not found for " << name << std::endl;
        return;
    }

    TH2D* hRatio = (TH2D*)hNum->Clone(name);
    hRatio->SetTitle(title);
    hRatio->Reset();

    for (int ix = 1; ix <= hNum->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hNum->GetNbinsY(); ++iy) {
            double num = hNum->GetBinContent(ix, iy);
            double den = hDen->GetBinContent(ix, iy);
            
            if (den > 0) {
                double ratio = num / den;
                if (ratio > 1.0) ratio = 1.0;
                double err = TMath::Sqrt(ratio * (1.0 - ratio) / den);
                
                hRatio->SetBinContent(ix, iy, ratio);
                hRatio->SetBinError(ix, iy, err);
            } else {
                hRatio->SetBinContent(ix, iy, 0);
                hRatio->SetBinError(ix, iy, 0);
            }
        }
    }
    hRatio->Write();
}

void CalAcceptance_Efficiency() {
    TFile *fileIn = new TFile("Made1Dhist_mc.root", "READ");
    if (!fileIn || fileIn->IsZombie()) {
        std::cout << "Error: Could not open Made1Dhist_mc.root" << std::endl;
        return;
    }

    TFile *fileOut = new TFile("AcceptanceEfficiency.root", "RECREATE");

    std::vector<TString> particles = {"Omega", "Phi"};
    std::vector<TString> vars = {"Pt", "Y", "PtY"};

    for (const auto& particle : particles) {
        TDirectory *partDir = fileOut->mkdir(particle);
        partDir->cd();

        for (const auto& var : vars) {
            TString pathGen = Form("Generated_All/%s/%s_uls", particle.Data(), var.Data());
            TString pathAcc = Form("Generated_Acc/%s/%s_uls", particle.Data(), var.Data());
            TString pathRec = Form("Reconstructed/%s/%s_uls", particle.Data(), var.Data());

            if (var == "PtY") {
                TH2D* hGen = (TH2D*)fileIn->Get(pathGen);
                TH2D* hAcc = (TH2D*)fileIn->Get(pathAcc);
                TH2D* hRec = (TH2D*)fileIn->Get(pathRec);

                // Acceptance (Acc / Gen)
                if (hGen && hAcc) {
                    TString name = Form("Acceptance_%s_%s", particle.Data(), var.Data());
                    TString title = Form("Acceptance (Acc/Gen) - %s %s", particle.Data(), var.Data());
                    CalculateRatio2D(hAcc, hGen, name.Data(), title.Data());
                }
                // Efficiency (Reco / Acc)
                if (hAcc && hRec) {
                    TString name = Form("Efficiency_%s_%s", particle.Data(), var.Data());
                    TString title = Form("Efficiency (Reco/Acc) - %s %s", particle.Data(), var.Data());
                    CalculateRatio2D(hRec, hAcc, name.Data(), title.Data());
                }
                // Acceptance * Efficiency (Reco / Gen)
                if (hGen && hRec) {
                    TString name = Form("AccEff_%s_%s", particle.Data(), var.Data());
                    TString title = Form("Acceptance * Efficiency (Reco/Gen) - %s %s", particle.Data(), var.Data());
                    CalculateRatio2D(hRec, hGen, name.Data(), title.Data());
                }
            } else {
                TH1D* hGen = (TH1D*)fileIn->Get(pathGen);
                TH1D* hAcc = (TH1D*)fileIn->Get(pathAcc);
                TH1D* hRec = (TH1D*)fileIn->Get(pathRec);

                // Acceptance (Acc / Gen)
                if (hGen && hAcc) {
                    TString name = Form("Acceptance_%s_%s", particle.Data(), var.Data());
                    TString title = Form("Acceptance (Acc/Gen) - %s %s", particle.Data(), var.Data());
                    CalculateRatio1D(hAcc, hGen, name.Data(), title.Data());
                }
                // Efficiency (Reco / Acc)
                if (hAcc && hRec) {
                    TString name = Form("Efficiency_%s_%s", particle.Data(), var.Data());
                    TString title = Form("Efficiency (Reco/Acc) - %s %s", particle.Data(), var.Data());
                    CalculateRatio1D(hRec, hAcc, name.Data(), title.Data());
                }
                // Acceptance * Efficiency (Reco / Gen)
                if (hGen && hRec) {
                    TString name = Form("AccEff_%s_%s", particle.Data(), var.Data());
                    TString title = Form("Acceptance * Efficiency (Reco/Gen) - %s %s", particle.Data(), var.Data());
                    CalculateRatio1D(hRec, hGen, name.Data(), title.Data());
                }
            }
        }
    }

    fileOut->Close();
    fileIn->Close();
    std::cout << "Acceptance and Efficiency calculation finished. Output: AcceptanceEfficiency.root" << std::endl;
}
