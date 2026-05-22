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
#include <TKey.h>

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

// --- Crystal Ball Fit ---
double GetYieldFromFitCB(TH1D* h, TString particle, double &yieldErr, TFile* fitFile, TString nameTag) {
    if (!h || h->GetEntries() < 10) {
        yieldErr = 0;
        return 0;
    }

    TF1 *fitFunc = nullptr;
    if (particle.Contains("Omega", TString::kIgnoreCase)) {
        fitFunc = new TF1(Form("fCB_%s", nameTag.Data()), "crystalball", 0.70, 0.86);
        fitFunc->SetParLimits(1, 0.75, 0.85);
        fitFunc->SetParLimits(2, 0.005, 0.05); // sigma
        fitFunc->SetParLimits(3, 0.5, 5.0);    // alpha
        fitFunc->SetParLimits(4, 1.0, 10.0);   // n
        fitFunc->SetParameters(h->GetMaximum(), 0.78, 0.02, 1.5, 2.0);
    } else if (particle.Contains("Phi", TString::kIgnoreCase)) {
        fitFunc = new TF1(Form("fCB_%s", nameTag.Data()), "crystalball", 0.92, 1.15);
        fitFunc->SetParLimits(1, 0.9, 1.1);
        fitFunc->SetParLimits(2, 0.005, 0.05); // sigma
        fitFunc->SetParLimits(3, 0.5, 5.0);    // alpha
        fitFunc->SetParLimits(4, 1.0, 10.0);   // n
        fitFunc->SetParameters(h->GetMaximum(), 1.02, 0.02, 1.5, 2.0);
    } else {
        yieldErr = 0; return 0;
    }

    TFitResultPtr r = h->Fit(fitFunc, "SRQ");
    double yield = 0; yieldErr = 0;
    if (int(r) == 0 || int(r) == 2) {
        double mean = fitFunc->GetParameter(1);
        double sigma = fitFunc->GetParameter(2);
        int binLow = h->FindBin(mean - 3 * sigma);
        int binHigh = h->FindBin(mean + 3 * sigma);
        double sum = 0, sumErr2 = 0;
        for (int i = binLow; i <= binHigh; ++i) {
            sum += h->GetBinContent(i);
            double err = h->GetBinError(i);
            sumErr2 += err * err;
        }
        yield = sum; yieldErr = TMath::Sqrt(sumErr2);
        if (fitFile) {
            fitFile->cd();
            TCanvas *c = new TCanvas(Form("cCB_%s", nameTag.Data()), nameTag + "_CB", 800, 600);
            h->Draw(); fitFunc->Draw("same"); c->Write(); delete c;
        }
    }
    delete fitFunc;
    return yield;
}

// --- Gaussian Fit ---
double GetYieldFromFitGauss(TH1D* h, TString particle, double &yieldErr, TFile* fitFile, TString nameTag) {
    if (!h || h->GetEntries() < 10) {
        yieldErr = 0; return 0;
    }

    TF1 *fitFunc = nullptr;
    if (particle.Contains("Omega", TString::kIgnoreCase)) {
        fitFunc = new TF1(Form("fGaus_%s", nameTag.Data()), "gaus", 0.70, 0.86);
        fitFunc->SetParLimits(1, 0.75, 0.85);
        fitFunc->SetParLimits(2, 0.005, 0.05); // sigma
        fitFunc->SetParameters(h->GetMaximum(), 0.78, 0.02);
    } else if (particle.Contains("Phi", TString::kIgnoreCase)) {
        fitFunc = new TF1(Form("fGaus_%s", nameTag.Data()), "gaus", 0.92, 1.15);
        fitFunc->SetParLimits(1, 0.9, 1.1);
        fitFunc->SetParLimits(2, 0.005, 0.05); // sigma
        fitFunc->SetParameters(h->GetMaximum(), 1.02, 0.02);
    } else {
        yieldErr = 0; return 0;
    }

    TFitResultPtr r = h->Fit(fitFunc, "SRQ");
    double yield = 0; yieldErr = 0;
    if (int(r) == 0 || int(r) == 2) {
        double mean = fitFunc->GetParameter(1);
        double sigma = fitFunc->GetParameter(2);
        int binLow = h->FindBin(mean - 3 * sigma);
        int binHigh = h->FindBin(mean + 3 * sigma);
        double sum = 0, sumErr2 = 0;
        for (int i = binLow; i <= binHigh; ++i) {
            sum += h->GetBinContent(i);
            double err = h->GetBinError(i);
            sumErr2 += err * err;
        }
        yield = sum; yieldErr = TMath::Sqrt(sumErr2);
        if (fitFile) {
            fitFile->cd();
            TCanvas *c = new TCanvas(Form("cGaus_%s", nameTag.Data()), nameTag + "_Gauss", 800, 600);
            h->Draw(); fitFunc->Draw("same"); c->Write(); delete c;
        }
    }
    delete fitFunc;
    return yield;
}

void CalculateRatio1D(TH1D* hNum, TH1D* hDen, const char* name, const char* title) {
    if (!hNum || !hDen) return;
    TH1D* hRatio = (TH1D*)hNum->Clone(name);
    hRatio->SetTitle(title); hRatio->Reset();
    for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
        double num = hNum->GetBinContent(i); double den = hDen->GetBinContent(i);
        double numErr = hNum->GetBinError(i); double denErr = hDen->GetBinError(i);
        if (den > 0) {
            double ratio = num / den; if (ratio > 1.2) ratio = 1.0; 
            double err = (num <= den) ? TMath::Sqrt(ratio * (1.0 - ratio) / den) : ratio * TMath::Sqrt(TMath::Power(numErr/num, 2) + TMath::Power(denErr/den, 2));
            hRatio->SetBinContent(i, ratio); hRatio->SetBinError(i, err);
        }
    }
    hRatio->Write();
}

void ProcessAnalysisForTopDir(TDirectory* topDirIn, TDirectory* topDirOut, TFile* fileFits, TString fitType) {
    std::vector<TString> particles = {"Omega", "Phi"};
    std::vector<TString> categories = {"Generated_All", "Generated_Acc", "Reconstructed"};
    std::vector<double> pt_bin_edges = {0, 0.5, 0.75, 1, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10};
    std::vector<double> y_bin_edges; for (double y = -4.0; y <= -2.5; y += 0.25) y_bin_edges.push_back(y);

    for (const auto& particle : particles) {
        TDirectory *partDir = topDirOut->mkdir(Form("%s/%s", fitType.Data(), particle.Data()), "", true);
        partDir->cd();

        TH1D* hGenPt = new TH1D(Form("hGenPt_%s", particle.Data()), "Gen Pt", pt_bin_edges.size()-1, &pt_bin_edges[0]);
        TH1D* hAccPt = new TH1D(Form("hAccPt_%s", particle.Data()), "Acc Pt", pt_bin_edges.size()-1, &pt_bin_edges[0]);
        TH1D* hRecPt = new TH1D(Form("hRecPt_%s", particle.Data()), "Rec Pt", pt_bin_edges.size()-1, &pt_bin_edges[0]);
        TH1D* hGenY = new TH1D(Form("hGenY_%s", particle.Data()), "Gen Y", y_bin_edges.size()-1, &y_bin_edges[0]);
        TH1D* hAccY = new TH1D(Form("hAccY_%s", particle.Data()), "Acc Y", y_bin_edges.size()-1, &y_bin_edges[0]);
        TH1D* hRecY = new TH1D(Form("hRecY_%s", particle.Data()), "Rec Y", y_bin_edges.size()-1, &y_bin_edges[0]);

        for (size_t i = 0; i < pt_bin_edges.size() - 1; ++i) {
            TString ptRange = Form("Pt_%.2fto%.2f", pt_bin_edges[i], pt_bin_edges[i+1]);
            for (const auto& cat : categories) {
                TString path = Form("%s/%s/%s/SEPM", cat.Data(), particle.Data(), ptRange.Data());
                TH1D* hMass = (TH1D*)topDirIn->Get(path); if (!hMass) continue;
                double yield = 0, yieldErr = 0;
                if (cat.Contains("Generated")) yield = hMass->IntegralAndError(1, hMass->GetNbinsX(), yieldErr);
                else {
                    if (fitType == "CB") yield = GetYieldFromFitCB(hMass, particle, yieldErr, fileFits, Form("%s_CB_%s_%s_%s", topDirIn->GetName(), particle.Data(), cat.Data(), ptRange.Data()));
                    else if (fitType == "Gauss") yield = GetYieldFromFitGauss(hMass, particle, yieldErr, fileFits, Form("%s_Gauss_%s_%s_%s", topDirIn->GetName(), particle.Data(), cat.Data(), ptRange.Data()));
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
                TH1D* hMass = (TH1D*)topDirIn->Get(path); if (!hMass) continue;
                double yield = 0, yieldErr = 0;
                if (cat.Contains("Generated")) yield = hMass->IntegralAndError(1, hMass->GetNbinsX(), yieldErr);
                else {
                    if (fitType == "CB") yield = GetYieldFromFitCB(hMass, particle, yieldErr, fileFits, Form("%s_CB_%s_%s_%s", topDirIn->GetName(), particle.Data(), cat.Data(), yRange.Data()));
                    else if (fitType == "Gauss") yield = GetYieldFromFitGauss(hMass, particle, yieldErr, fileFits, Form("%s_Gauss_%s_%s_%s", topDirIn->GetName(), particle.Data(), cat.Data(), yRange.Data()));
                }
                if (cat == "Generated_All") { hGenY->SetBinContent(i+1, yield); hGenY->SetBinError(i+1, yieldErr); }
                else if (cat == "Generated_Acc") { hAccY->SetBinContent(i+1, yield); hAccY->SetBinError(i+1, yieldErr); }
                else if (cat == "Reconstructed") { hRecY->SetBinContent(i+1, yield); hRecY->SetBinError(i+1, yieldErr); }
            }
        }
        partDir->cd();
        hGenPt->Write(); hAccPt->Write(); hRecPt->Write(); hGenY->Write(); hAccY->Write(); hRecY->Write();
        CalculateRatio1D(hAccPt, hGenPt, Form("Acceptance_%s_Pt", particle.Data()), "Acceptance vs Pt");
        CalculateRatio1D(hRecPt, hAccPt, Form("Efficiency_%s_Pt", particle.Data()), "Efficiency vs Pt");
        CalculateRatio1D(hRecPt, hGenPt, Form("AccEff_%s_Pt", particle.Data()), "Acc*Eff vs Pt");
        CalculateRatio1D(hAccY, hGenY, Form("Acceptance_%s_Y", particle.Data()), "Acceptance vs Y");
        CalculateRatio1D(hRecY, hAccY, Form("Efficiency_%s_Y", particle.Data()), "Efficiency vs Y");
        CalculateRatio1D(hRecY, hGenY, Form("AccEff_%s_Y", particle.Data()), "Acc*Eff vs Y");

        // --- Integrated (0-10 GeV) ---
        double genAll = hGenPt->Integral();
        double accAll = hAccPt->Integral();
        double recAll = hRecPt->Integral();
        if (genAll > 0 && accAll > 0) {
            double acc = accAll / genAll;
            double eff = recAll / accAll;
            double acceff = recAll / genAll;
            std::cout << "    [Acc/Eff 0-10 GeV] " << particle << " (" << fitType << "): "
                      << "Acc=" << acc << ", Eff=" << eff << ", Acc*Eff=" << acceff << std::endl;
        }
    }
}

void CalAcceptance_Efficiency() {
    gStyle->SetOptFit(1111); gStyle->SetOptTitle(0); gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.12);
    TFile *fileIn = new TFile("Made1Dhist_mc.root", "READ"); if (!fileIn || fileIn->IsZombie()) return;
    TFile *fileOut = new TFile("AcceptanceEfficiency_Both.root", "RECREATE");
    TFile *fileFits = new TFile("AcceptanceEfficiency_Fits_Both.root", "RECREATE");

    TIter nextKey(fileIn->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextKey())) {
        if (TString(key->GetClassName()) != "TDirectoryFile") continue;
        TString topDirName = key->GetName();
        if (!topDirName.BeginsWith("dimuon-mc")) continue;
        TDirectory *topDirIn = (TDirectory*)fileIn->Get(topDirName);
        TDirectory *topDirOut = fileOut->mkdir(topDirName);
        std::cout << "Calculating Efficiency for Top Directory: " << topDirName << std::endl;
        ProcessAnalysisForTopDir(topDirIn, topDirOut, fileFits, "CB");
        ProcessAnalysisForTopDir(topDirIn, topDirOut, fileFits, "Gauss");
    }
    fileOut->Close(); fileFits->Close(); fileIn->Close();
    std::cout << "Acceptance and Efficiency calculation finished." << std::endl;
}
