#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <iostream>
#include <vector>
#include <TString.h>
#include <TDirectory.h>

void SetPlotStyle(TH1D* h, int color, int marker, double size, TString xTitle) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(size);
    h->SetStats(0);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
}

double GetYieldFromFit(TH1D* h, TString particle, double &yieldErr, TFile* fitFile, TString nameTag) {
    if (!h || h->GetEntries() < 10) {
        yieldErr = 0;
        return 0;
    }

    TF1 *fitFunc = nullptr;
    if (particle.Contains("Omega", TString::kIgnoreCase)) {
        fitFunc = new TF1(Form("fCB_%s", nameTag.Data()), "crystalball", 0.70, 0.86);
        fitFunc->SetParameters(h->GetMaximum(), 0.78, 0.02, 1.5, 2.0);
        fitFunc->SetParLimits(1, 0.75, 0.85);
    } else if (particle.Contains("Phi", TString::kIgnoreCase)) {
        fitFunc = new TF1(Form("fCB_%s", nameTag.Data()), "crystalball", 0.92, 1.15);
        fitFunc->SetParameters(h->GetMaximum(), 1.02, 0.02, 1.5, 2.0);
        fitFunc->SetParLimits(1, 0.9, 1.1);
    } else {
        yieldErr = 0;
        return 0;
    }

    fitFunc->SetParLimits(2, 0.005, 0.05); // sigma
    fitFunc->SetParLimits(3, 0.5, 5.0);    // alpha
    fitFunc->SetParLimits(4, 1.0, 10.0);   // n

    TFitResultPtr r = h->Fit(fitFunc, "SRQ");
    
    double yield = 0;
    if (int(r) == 0) {
        double mean = fitFunc->GetParameter(1);
        double sigma = fitFunc->GetParameter(2);
        double binWidth = h->GetBinWidth(1);
        
        yield = fitFunc->Integral(mean - 3 * sigma, mean + 3 * sigma) / binWidth;
        yieldErr = fitFunc->IntegralError(mean - 3 * sigma, mean + 3 * sigma) / binWidth;
        
        if (fitFile) {
            fitFile->cd();
            TCanvas *c = new TCanvas(Form("c_%s", nameTag.Data()), nameTag, 800, 600);
            h->Draw();
            fitFunc->Draw("same");
            c->Write();
            delete c;
        }
    } else {
        std::cout << "Warning: Fit failed for " << nameTag << std::endl;
        yieldErr = 0;
        yield = 0;
    }

    delete fitFunc;
    return yield;
}

void CalculateRatio1D(TH1D* hNum, TH1D* hDen, const char* name, const char* title) {
    if (!hNum || !hDen) return;

    TH1D* hRatio = (TH1D*)hNum->Clone(name);
    hRatio->SetTitle(title);
    hRatio->Reset();

    for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
        double num = hNum->GetBinContent(i);
        double den = hDen->GetBinContent(i);
        double numErr = hNum->GetBinError(i);
        double denErr = hDen->GetBinError(i);
        
        if (den > 0) {
            double ratio = num / den;
            if (ratio > 1.2) ratio = 1.0; 
            
            double err = 0;
            if (num <= den) {
                err = TMath::Sqrt(ratio * (1.0 - ratio) / den);
            } else {
                err = ratio * TMath::Sqrt(TMath::Power(numErr/num, 2) + TMath::Power(denErr/den, 2));
            }
            
            hRatio->SetBinContent(i, ratio);
            hRatio->SetBinError(i, err);
        }
    }
    hRatio->Write();
}

void CalAcceptance_Efficiency() {
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);

    TFile *fileIn = new TFile("Made1Dhist_mc.root", "READ");
    if (!fileIn || fileIn->IsZombie()) {
        std::cout << "Error: Could not open Made1Dhist_mc.root" << std::endl;
        return;
    }

    TFile *fileOut = new TFile("AcceptanceEfficiency.root", "RECREATE");
    TFile *fileFits = new TFile("AcceptanceEfficiency_Fits.root", "RECREATE");

    std::vector<TString> particles = {"Omega", "Phi"};
    std::vector<TString> categories = {"Generated_All", "Generated_Acc", "Reconstructed"};
    
    std::vector<double> pt_bin_edges;
    for (double pt = 0.0; pt <= 10.0; pt += 0.5) pt_bin_edges.push_back(pt);
    std::vector<double> y_bin_edges;
    for (double y = -4.0; y <= -2.5; y += 0.25) y_bin_edges.push_back(y);

    for (const auto& particle : particles) {
        TDirectory *partDir = fileOut->mkdir(particle);
        partDir->cd();

        TH1D* hGenPt = new TH1D(Form("hGenPt_%s", particle.Data()), "Gen Pt", pt_bin_edges.size()-1, &pt_bin_edges[0]);
        TH1D* hAccPt = new TH1D(Form("hAccPt_%s", particle.Data()), "Acc Pt", pt_bin_edges.size()-1, &pt_bin_edges[0]);
        TH1D* hRecPt = new TH1D(Form("hRecPt_%s", particle.Data()), "Rec Pt", pt_bin_edges.size()-1, &pt_bin_edges[0]);

        TH1D* hGenY = new TH1D(Form("hGenY_%s", particle.Data()), "Gen Y", y_bin_edges.size()-1, &y_bin_edges[0]);
        TH1D* hAccY = new TH1D(Form("hAccY_%s", particle.Data()), "Acc Y", y_bin_edges.size()-1, &y_bin_edges[0]);
        TH1D* hRecY = new TH1D(Form("hRecY_%s", particle.Data()), "Rec Y", y_bin_edges.size()-1, &y_bin_edges[0]);

        for (size_t i = 0; i < pt_bin_edges.size() - 1; ++i) {
            TString ptRange = Form("Pt_%.1fto%.1f", pt_bin_edges[i], pt_bin_edges[i+1]);
            for (const auto& cat : categories) {
                TString path = Form("%s/%s/%s/SEPM", cat.Data(), particle.Data(), ptRange.Data());
                TH1D* hMass = (TH1D*)fileIn->Get(path);
                if (!hMass) continue;

                double yield = 0, yieldErr = 0;
                if (cat.Contains("Generated")) {
                    yield = hMass->IntegralAndError(1, hMass->GetNbinsX(), yieldErr);
                } else {
                    yield = GetYieldFromFit(hMass, particle, yieldErr, fileFits, Form("%s_%s_%s", particle.Data(), cat.Data(), ptRange.Data()));
                }
                
                if (cat == "Generated_All") { hGenPt->SetBinContent(i+1, yield); hGenPt->SetBinError(i+1, yieldErr); }
                else if (cat == "Generated_Acc") { hAccPt->SetBinContent(i+1, yield); hAccPt->SetBinError(i+1, yieldErr); }
                else if (cat == "Reconstructed") { hRecPt->SetBinContent(i+1, yield); hRecPt->SetBinError(i+1, yieldErr); }
            }
        }

        for (size_t i = 0; i < y_bin_edges.size() - 1; ++i) {
            TString yRange = Form("Y_%.2fto%.2f", y_bin_edges[i], y_bin_edges[i+1]);
            for (const auto& cat : categories) {
                TString path = Form("%s/%s/%s/SEPM", cat.Data(), particle.Data(), yRange.Data());
                TH1D* hMass = (TH1D*)fileIn->Get(path);
                if (!hMass) continue;

                double yield = 0, yieldErr = 0;
                if (cat.Contains("Generated")) {
                    yield = hMass->IntegralAndError(1, hMass->GetNbinsX(), yieldErr);
                } else {
                    yield = GetYieldFromFit(hMass, particle, yieldErr, fileFits, Form("%s_%s_%s", particle.Data(), cat.Data(), yRange.Data()));
                }
                
                if (cat == "Generated_All") { hGenY->SetBinContent(i+1, yield); hGenY->SetBinError(i+1, yieldErr); }
                else if (cat == "Generated_Acc") { hAccY->SetBinContent(i+1, yield); hAccY->SetBinError(i+1, yieldErr); }
                else if (cat == "Reconstructed") { hRecY->SetBinContent(i+1, yield); hRecY->SetBinError(i+1, yieldErr); }
            }
        }

        partDir->cd();
        hGenPt->Write(); hAccPt->Write(); hRecPt->Write();
        hGenY->Write(); hAccY->Write(); hRecY->Write();

        CalculateRatio1D(hAccPt, hGenPt, Form("Acceptance_%s_Pt", particle.Data()), "Acceptance vs Pt");
        CalculateRatio1D(hRecPt, hAccPt, Form("Efficiency_%s_Pt", particle.Data()), "Efficiency vs Pt");
        CalculateRatio1D(hRecPt, hGenPt, Form("AccEff_%s_Pt", particle.Data()), "Acc*Eff vs Pt");

        CalculateRatio1D(hAccY, hGenY, Form("Acceptance_%s_Y", particle.Data()), "Acceptance vs Y");
        CalculateRatio1D(hRecY, hAccY, Form("Efficiency_%s_Y", particle.Data()), "Efficiency vs Y");
        CalculateRatio1D(hRecY, hGenY, Form("AccEff_%s_Y", particle.Data()), "Acc*Eff vs Y");
    }

    // --- Plotting Section ---
    TString variables[] = {"Pt", "Y"};
    TString xTitles[] = {"p_{T} (GeV/c)", "y"};

    for (int i = 0; i < 2; ++i) {
        TString var = variables[i];
        TString xTitle = xTitles[i];

        TCanvas *cCounts = new TCanvas(Form("cCounts_%s", var.Data()), Form("Counts vs %s", var.Data()), 1200, 600);
        cCounts->Divide(2, 1);
        for (int p = 0; p < 2; ++p) {
            cCounts->cd(p + 1);
            gPad->SetLogy();
            gPad->SetTicks(1, 1);
            TH1D *hGen = (TH1D*)fileOut->Get(Form("%s/hGen%s_%s", particles[p].Data(), var.Data(), particles[p].Data()));
            TH1D *hAcc = (TH1D*)fileOut->Get(Form("%s/hAcc%s_%s", particles[p].Data(), var.Data(), particles[p].Data()));
            TH1D *hRec = (TH1D*)fileOut->Get(Form("%s/hRec%s_%s", particles[p].Data(), var.Data(), particles[p].Data()));
            if (hGen && hAcc && hRec) {
                SetPlotStyle(hGen, kBlack, 20, 0.6, xTitle);
                SetPlotStyle(hAcc, kBlue, 20, 0.6, xTitle);
                SetPlotStyle(hRec, kGreen+2, 20, 0.6, xTitle);
                hGen->GetYaxis()->SetTitle("Yield");
                hGen->Draw("P"); hAcc->Draw("P same"); hRec->Draw("P same");
                TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                TString displayName = (particles[p] == "Omega") ? "#omega #rightarrow #mu#mu" : "#phi #rightarrow #mu#mu";
                leg->SetHeader(Form("#bf{%s}", displayName.Data()));
                leg->AddEntry(hGen, "Generated", "p"); leg->AddEntry(hAcc, "Accepted", "p"); leg->AddEntry(hRec, "Reconstructed", "p");
                leg->Draw();
            }
        }
        fileOut->cd();
        cCounts->Write();

        TCanvas *cRatios = new TCanvas(Form("cRatios_%s", var.Data()), Form("Ratios vs %s", var.Data()), 1200, 600);
        cRatios->Divide(2, 1);
        for (int p = 0; p < 2; ++p) {
            cRatios->cd(p + 1);
            gPad->SetLogy();
            gPad->SetTicks(1, 1);
            TH1D *hA = (TH1D*)fileOut->Get(Form("%s/Acceptance_%s_%s", particles[p].Data(), particles[p].Data(), var.Data()));
            TH1D *hE = (TH1D*)fileOut->Get(Form("%s/Efficiency_%s_%s", particles[p].Data(), particles[p].Data(), var.Data()));
            TH1D *hAE = (TH1D*)fileOut->Get(Form("%s/AccEff_%s_%s", particles[p].Data(), particles[p].Data(), var.Data()));
            if (hA && hE && hAE) {
                SetPlotStyle(hA, kBlack, 20, 0.6, xTitle);
                SetPlotStyle(hE, kBlue, 20, 0.6, xTitle);
                SetPlotStyle(hAE, kGreen+2, 20, 0.6, xTitle);
                hA->GetYaxis()->SetRangeUser(0.001, 1.0);
                hA->GetYaxis()->SetTitle("Acceptance / Efficiency");
                hA->Draw("P"); hE->Draw("P same"); hAE->Draw("P same");
                TLegend *leg = new TLegend(0.55, 0.15, 0.88, 0.35);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                TString displayName = (particles[p] == "Omega") ? "#omega #rightarrow #mu#mu" : "#phi #rightarrow #mu#mu";
                leg->SetHeader(Form("#bf{%s}", displayName.Data()));
                leg->AddEntry(hA, "Acceptance", "p"); leg->AddEntry(hE, "Efficiency", "p"); leg->AddEntry(hAE, "Acc * Eff", "p");
                leg->Draw();
            }
        }
        fileOut->cd();
        cRatios->Write();
    }

    fileOut->Close();
    fileFits->Close();
    fileIn->Close();
    std::cout << "Acceptance and Efficiency calculation finished. Canvases saved in AcceptanceEfficiency.root." << std::endl;
}
