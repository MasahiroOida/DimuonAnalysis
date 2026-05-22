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
#include <map>
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

    TIter nextKey(file->GetListOfKeys());
    TKey *key;

    std::vector<TString> distNames = {"SEPM", "Pt", "Y"};
    std::vector<int> colors = {1, 632, 600, 418, 807, 616}; // 6つのビン用

    while ((key = (TKey*)nextKey())) {
        if (TString(key->GetClassName()) != "TDirectoryFile") continue;
        TString topDirName = key->GetName();
        if (!topDirName.BeginsWith("dimuon-mc")) continue;
        
        TDirectory *topDirInput = (TDirectory*)file->Get(topDirName);
        std::cout << "Checking MC Mass for Top Directory: " << topDirName << std::endl;
        
        TDirectory *outTopDir = outFile->mkdir(topDirName);

        std::vector<TString> categories = {"Generated_All", "Generated_Acc", "Reconstructed"};
        for (const auto& catName : categories) {
            TDirectory *catDir = (TDirectory*)topDirInput->Get(catName);
            if (!catDir) continue;
            std::cout << "  Processing category: " << catName << std::endl;

            std::map<TString, std::vector<TString>> binNamesPerDist;

            for (const auto& partName : {"Omega", "Phi"}) {
                TDirectory *partDir = (TDirectory*)catDir->Get(partName);
                if (!partDir) continue;

                for (const auto& distName : distNames) {
                    TIter nextPt(partDir->GetListOfKeys());
                    TKey *ptKey;

                    std::vector<double> vec_pt, vec_pt_err;
                    std::vector<double> vec_mean, vec_mean_err;
                    std::vector<double> vec_sigma, vec_sigma_err;
                    std::vector<double> vec_alpha, vec_alpha_err;
                    std::vector<double> vec_n, vec_n_err;

                    while ((ptKey = (TKey*)nextPt())) {
                        if (TString(ptKey->GetClassName()) != "TDirectoryFile") continue;
                        TString ptDirName = ptKey->GetName();

                        if (distName == "SEPM" && ptDirName.BeginsWith("Y_")) continue;
                        if (distName == "Pt" && ptDirName.BeginsWith("Y_")) continue;
                        if (distName == "Y" && ptDirName.BeginsWith("Pt_")) continue;

                        bool found = false;
                        for(const auto& b : binNamesPerDist[distName]) if(b == ptDirName) found = true;
                        if(!found) binNamesPerDist[distName].push_back(ptDirName);

                        TDirectory *ptDir = (TDirectory*)partDir->Get(ptDirName);
                        TH1D *h = (TH1D*)ptDir->Get(distName);
                        if (!h || h->GetEntries() == 0) continue;

                        TH1D *hClone = (TH1D*)h->Clone(Form("%s_%s_%s_%s_%s", topDirName.Data(), catName.Data(), partName, ptDirName.Data(), distName.Data()));
                        hClone->SetDirectory(0);

                        if (distName == "SEPM") {
                            TF1 *fitFunc = nullptr;
                            if (TString(partName).Contains("Omega", TString::kIgnoreCase)) {
                                fitFunc = new TF1(Form("fit_%s_%s", hClone->GetName(), partName), "crystalball", 0.70, 0.86);
                                fitFunc->SetParameters(hClone->GetMaximum() / 2.0, 0.78, 0.02, 1.5, 2.0);
                                fitFunc->SetParLimits(1, 0.75, 0.85);
                            } else if (TString(partName).Contains("Phi", TString::kIgnoreCase)) {
                                fitFunc = new TF1(Form("fit_%s_%s", hClone->GetName(), partName), "crystalball", 0.92, 1.15);
                                fitFunc->SetParameters(hClone->GetMaximum() / 2.0, 1.02, 0.02, 1.5, 2.0);
                                fitFunc->SetParLimits(1, 0.9, 1.1);
                            }

                            if (fitFunc) {
                                hClone->Fit(fitFunc, "RQ");
                                double val_ptmin = -1, val_ptmax = -1;
                                TParameter<double> *pMin = (TParameter<double>*)ptDir->Get("ptmin");
                                TParameter<double> *pMax = (TParameter<double>*)ptDir->Get("ptmax");
                                if (pMin && pMax) {
                                    val_ptmin = pMin->GetVal();
                                    val_ptmax = pMax->GetVal();
                                } else {
                                    sscanf(ptDirName.Data(), "Pt_%lfto%lf", &val_ptmin, &val_ptmax);
                                }
                                if (val_ptmin != -1 && val_ptmax != -1 && (val_ptmax - val_ptmin < 9.0)) {
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
                                delete fitFunc;
                            }
                        }

                        TDirectory *outCatDir = outTopDir->GetDirectory(catName);
                        if (!outCatDir) outCatDir = outTopDir->mkdir(catName);
                        TDirectory *outPartDir = outCatDir->GetDirectory(partName);
                        if (!outPartDir) outPartDir = outCatDir->mkdir(partName);
                        TDirectory *outPtDir = outPartDir->GetDirectory(ptDirName);
                        if (!outPtDir) outPtDir = outPartDir->mkdir(ptDirName);
                        
                        outPtDir->cd();
                        hClone->GetListOfFunctions()->Clear();
                        hClone->Write(distName + "_fit");

                        TParameter<double> *pN = (TParameter<double> *)ptDir->Get("nEvents");
                        if (pN && pN->GetVal() > 0) {
                            TH1D *hNorm = (TH1D*)hClone->Clone(TString(hClone->GetName()) + "_norm");
                            // Keep raw counts
                            hNorm->SetYTitle("Counts / (Mass bin)");
                            hNorm->Write(distName + "_norm");
                            delete hNorm;
                        }
                        delete hClone;
                    }

                    if (distName == "SEPM" && !vec_pt.empty()) {
                        TDirectory *outPartDir = outTopDir->GetDirectory(Form("%s/%s", catName.Data(), partName));
                        outPartDir->cd();
                        TGraphErrors *grMean = new TGraphErrors(vec_pt.size(), &vec_pt[0], &vec_mean[0], &vec_pt_err[0], &vec_mean_err[0]);
                        grMean->SetName("g_Mean"); grMean->Write(); delete grMean;
                        TGraphErrors *grSigma = new TGraphErrors(vec_pt.size(), &vec_pt[0], &vec_sigma[0], &vec_pt_err[0], &vec_sigma_err[0]);
                        grSigma->SetName("g_Sigma"); grSigma->Write(); delete grSigma;
                    }
                }
            }

            // --- Combined Omega + Phi Plots ---
            for (const auto& distName : distNames) {
                for (const auto& ptDirName : binNamesPerDist[distName]) {
                    TH1D *hOmega = (TH1D*)topDirInput->Get(Form("%s/Omega/%s/%s", catName.Data(), ptDirName.Data(), distName.Data()));
                    TH1D *hPhi = (TH1D*)topDirInput->Get(Form("%s/Phi/%s/%s", catName.Data(), ptDirName.Data(), distName.Data()));
                    if (!hOmega && !hPhi) continue;

                    TCanvas *cCombined = new TCanvas(Form("c_%s_%s_%s", catName.Data(), ptDirName.Data(), distName.Data()), "Omega and Phi", 800, 600);
                    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
                    
                    double maxVal = 0;
                    if (hOmega) {
                        hOmega->SetLineColor(kBlue); hOmega->SetMarkerColor(kBlue); hOmega->SetMarkerStyle(20);
                        maxVal = std::max(maxVal, hOmega->GetMaximum());
                    }
                    if (hPhi) {
                        hPhi->SetLineColor(kRed); hPhi->SetMarkerColor(kRed); hPhi->SetMarkerStyle(24);
                        maxVal = std::max(maxVal, hPhi->GetMaximum());
                    }

                    if (hOmega) {
                        hOmega->GetYaxis()->SetRangeUser(0, maxVal * 1.2);
                        hOmega->Draw("PE");
                        leg->AddEntry(hOmega, "Omega", "pe");
                        if (hPhi) hPhi->Draw("PE SAME");
                    } else if (hPhi) {
                        hPhi->GetYaxis()->SetRangeUser(0, maxVal * 1.2);
                        hPhi->Draw("PE");
                    }
                    if (hPhi) leg->AddEntry(hPhi, "Phi", "pe");
                    
                    // Add Sum and Normalize by nEvents
                    TH1D *hSum = nullptr;
                    if (hOmega && hPhi) {
                        hSum = (TH1D*)hPhi->Clone(Form("Sum_%s_%s", ptDirName.Data(), distName.Data()));
                        hSum->SetDirectory(0);
                        hSum->Scale(5.72e-4); hSum->Add(hOmega, 14.8e-5);

                        // Normalize by nEvents
                        double nEvents = -1;
                        TDirectory *ptDirRef = (TDirectory*)topDirInput->Get(Form("%s/Omega/%s", catName.Data(), ptDirName.Data()));
                        if (!ptDirRef) ptDirRef = (TDirectory*)topDirInput->Get(Form("%s/Phi/%s", catName.Data(), ptDirName.Data()));
                        if (ptDirRef) {
                            TParameter<double> *pN = (TParameter<double> *)ptDirRef->Get("nEvents");
                            if (pN) nEvents = pN->GetVal();
                        }

                        if (nEvents > 0) {
                            // Keep raw counts
                            hSum->SetYTitle("Counts / (Mass bin)");
                        }

                        hSum->SetLineColor(kRed); hSum->SetMarkerColor(kRed); hSum->SetMarkerStyle(20); hSum->SetMarkerSize(0.6);
                        hSum->Draw("PE SAME");
                        leg->AddEntry(hSum, "Sum (Scaled/Norm)", "pe");
                    }

                    leg->Draw();
                    TDirectory *outCatDir = outTopDir->GetDirectory(catName);
                    TDirectory *outSumDir = outCatDir->GetDirectory("Sum_Omega_Phi");
                    if (!outSumDir) outSumDir = outCatDir->mkdir("Sum_Omega_Phi");
                    TDirectory *outPtDir = outSumDir->GetDirectory(ptDirName);
                    if (!outPtDir) outPtDir = outSumDir->mkdir(ptDirName);
                    outPtDir->cd();
                    cCombined->Write();
                    if (hSum) hSum->Write(distName + "_sum");
                    
                    delete cCombined;
                    if (hSum) delete hSum;
                }
            }
        }
    }

    outFile->Close();
    file->Close();
    std::cout << "Successfully saved results to check_mc_mass.root" << std::endl;
}
