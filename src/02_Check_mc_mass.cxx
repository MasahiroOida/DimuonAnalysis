#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TKey.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TParameter.h>
#include <vector>
#include <iostream>
#include "myAnalysis.h"

void Check_mc_mass() {
    gStyle->SetOptFit(1111);
    TFile *file = TFile::Open("Made1Dhist_mc.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open Made1Dhist_mc.root" << std::endl;
        return;
    }

    // 出力用のROOTファイルを作成
    TFile *outFile = new TFile("check_mc_mass.root", "RECREATE");

    TIter nextCat(file->GetListOfKeys());
    TKey *catKey;

    std::vector<TString> distNames = {"SEPM", "Pt", "Y"};
    std::vector<int> colors = {1, 632, 600, 418, 807, 616}; // 6つのビン用

    // Categories (Generated_All, Generated_Acc, Reconstructed) をループ
    while ((catKey = (TKey*)nextCat())) {
        if (TString(catKey->GetClassName()) != "TDirectoryFile") continue;
        
        TString catName = catKey->GetName();
        TDirectory *catDir = (TDirectory*)file->Get(catName);
        std::cout << "Processing category: " << catName << std::endl;

        TIter nextPart(catDir->GetListOfKeys());
        TKey *partKey;

        // Omega, Phi などの粒子ディレクトリをループ
        while ((partKey = (TKey*)nextPart())) {
            if (TString(partKey->GetClassName()) != "TDirectoryFile") continue;
            
            TString partName = partKey->GetName();
            TDirectory *partDir = (TDirectory*)catDir->Get(partName);
            std::cout << "  Processing particle: " << partName << std::endl;

            for (const auto& distName : distNames) {
                TString canvasName = Form("c_%s_%s_%s", catName.Data(), partName.Data(), distName.Data());
                TCanvas *c = new TCanvas(canvasName, canvasName, 1200, 800);
                c->SetGrid();
                
                TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
                leg->SetHeader(Form("%s %s %s distributions", catName.Data(), partName.Data(), distName.Data()));
                leg->SetTextSize(0.025);

                TIter nextPt(partDir->GetListOfKeys());
                TKey *ptKey;
                std::vector<TH1D*> histograms;
                double globalMax = 0;

                // フィット結果格納用
                std::vector<double> vec_pt, vec_pt_err;
                std::vector<double> vec_mean, vec_mean_err;
                std::vector<double> vec_sigma, vec_sigma_err;
                std::vector<double> vec_alpha, vec_alpha_err;
                std::vector<double> vec_n, vec_n_err;

                while ((ptKey = (TKey*)nextPt())) {
                    if (TString(ptKey->GetClassName()) != "TDirectoryFile") continue;
                    
                    TString ptDirName = ptKey->GetName();
                    TDirectory *ptDir = (TDirectory*)partDir->Get(ptDirName);
                    TH1D *h = (TH1D*)ptDir->Get(distName);

                    if (!h || h->GetEntries() == 0) continue;

                    h->SetDirectory(0);
                    h->SetName(ptDirName);
                    
                    double currentMax = h->GetMaximum();
                    if (currentMax > globalMax) globalMax = currentMax;
                    
                    // SEPM (Invariant Mass) の場合にフィッティングを行う
                    if (distName == "SEPM") {
                        TF1 *fitFunc = nullptr;
                        if (partName.Contains("Omega", TString::kIgnoreCase)) {
                            fitFunc = new TF1(Form("fit_%s_%s", h->GetName(), partName.Data()), "crystalball", 0.70, 0.86);
                            fitFunc->SetParameters(h->GetMaximum() / 2.0, 0.78, 0.02, 1.5, 2.0);
                            fitFunc->SetParLimits(1, 0.75, 0.85);
                            fitFunc->SetParLimits(2, 0.01, 1.0);
                            fitFunc->SetParLimits(3, 0.5, 3.0);
                            fitFunc->SetParLimits(4, 1, 10.0);
                        } else if (partName.Contains("Phi", TString::kIgnoreCase)) {
                            fitFunc = new TF1(Form("fit_%s_%s", h->GetName(), partName.Data()), "crystalball", 0.92, 1.15);
                            fitFunc->SetParameters(h->GetMaximum() / 2.0, 1.02, 0.02, 1.5, 2.0);
                            fitFunc->SetParLimits(1, 0.9, 1.1);
                            fitFunc->SetParLimits(2, 0.01, 1.0);
                            fitFunc->SetParLimits(3, 0.5, 3.0);
                            fitFunc->SetParLimits(4, 1, 10.0);
                        }

                        if (fitFunc) {
                            h->Fit(fitFunc, "RQ");

                            TParameter<double> *pMin = (TParameter<double>*)ptDir->Get("ptmin");
                            TParameter<double> *pMax = (TParameter<double>*)ptDir->Get("ptmax");
                            double val_ptmin = -1, val_ptmax = -1;

                            if (pMin && pMax) {
                                val_ptmin = pMin->GetVal();
                                val_ptmax = pMax->GetVal();
                            } else {
                                // Fallback: parse from ptDirName (e.g., "Pt_0.0to0.5")
                                sscanf(ptDirName.Data(), "Pt_%lfto%lf", &val_ptmin, &val_ptmax);
                            }

                            if (val_ptmin != -1 && val_ptmax != -1) {
                                if (val_ptmax - val_ptmin > 9.0) continue; // 0-10 GeV/c などの全領域ビンを除外
                                vec_pt.push_back((val_ptmin + val_ptmax) / 2.0);
                                vec_pt_err.push_back((val_ptmax - val_ptmin) / 2.0);
                                vec_mean.push_back(fitFunc->GetParameter(1));
                                vec_mean_err.push_back(fitFunc->GetParError(1));
                                vec_sigma.push_back(fitFunc->GetParameter(2));
                                vec_sigma_err.push_back(fitFunc->GetParError(2));
                                vec_alpha.push_back(fitFunc->GetParameter(3));
                                vec_alpha_err.push_back(fitFunc->GetParError(3));
                                vec_n.push_back(fitFunc->GetParameter(4));
                                vec_n_err.push_back(fitFunc->GetParError(4));
                            }

                            // 個別の結果を保存するためのキャンバス作成
                            TString safeCatName = catName; safeCatName.ReplaceAll("/", "_");
                            TString individualCanvasName = Form("fit_%s_%s_%s_%s", safeCatName.Data(), partName.Data(), distName.Data(), h->GetName());
                            TCanvas *cIndiv = new TCanvas(individualCanvasName, individualCanvasName, 800, 600);
                            h->SetStats(1);
                            h->Draw("PE");
                            cIndiv->Modified();
                            cIndiv->Update();
                            outFile->cd();
                            cIndiv->Write();
                            delete cIndiv;
                            h->SetStats(0);
                        }
                    }
                    histograms.push_back(h);
                }

                if (histograms.empty()) {
                    delete c;
                    continue;
                }

                c->cd();
                for (size_t i = 0; i < histograms.size(); ++i) {
                    TH1D *h = histograms[i];
                    int color = colors[i % colors.size()];

                    h->SetMarkerStyle(20);
                    h->SetMarkerSize(0.8);
                    h->SetMarkerColor(color);
                    h->SetLineColor(color);
                    h->SetLineWidth(2);
                    
                    if (i == 0) {
                        TString title = Form("%s %s %s;%s;Counts", catName.Data(), partName.Data(), distName.Data(), distName.Data());
                        h->SetTitle(title);
                        h->GetYaxis()->SetRangeUser(0, globalMax * 1.3);
                        h->Draw("PE");
                    } else {
                        h->Draw("PE SAME");
                    }
                    leg->AddEntry(h, h->GetName(), "pe");
                }

                leg->Draw();
                c->Modified();
                c->Update();
                
                outFile->cd();
                c->Write();
                delete c;

                // フィット結果のサマリーキャンバス作成（個別に保存）
                if (distName == "SEPM" && !vec_pt.empty()) {
                    TString safeCatName = catName; safeCatName.ReplaceAll("/", "_");
                    
                    auto saveSummary = [&](std::vector<double>& y, std::vector<double>& yerr, const char* type, const char* title, const char* ytitle) {
                        TString gName = Form("g_%s_%s_%s", type, safeCatName.Data(), partName.Data());
                        TString cName = Form("Summary_%s_%s_%s", type, safeCatName.Data(), partName.Data());
                        TCanvas *cTemp = new TCanvas(cName, cName, 800, 600);
                        TGraphErrors *gr = new TGraphErrors(vec_pt.size(), &vec_pt[0], &y[0], &vec_pt_err[0], &yerr[0]);
                        gr->SetName(gName);
                        gr->SetTitle(Form("%s;p_{T} (GeV/c);%s", title, ytitle));
                        gr->SetMarkerStyle(20);
                        gr->SetMarkerColor(kBlue);
                        gr->SetLineColor(kBlue);
                        gr->Draw("APE");
                        outFile->cd();
                        cTemp->Write();
                        gr->Write();
                        delete cTemp;
                        delete gr;
                    };

                    saveSummary(vec_mean, vec_mean_err, "Mean", "Mean Mass", "Mean (GeV/c^{2})");
                    saveSummary(vec_sigma, vec_sigma_err, "Sigma", "Mass Width (Sigma)", "#sigma (GeV/c^{2})");
                    saveSummary(vec_alpha, vec_alpha_err, "Alpha", "CB Alpha", "#alpha");
                    saveSummary(vec_n, vec_n_err, "n", "CB n", "n");
                }
            }
        }
    }

    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
    std::cout << "Successfully saved canvases to check_mc_mass.root" << std::endl;
}
