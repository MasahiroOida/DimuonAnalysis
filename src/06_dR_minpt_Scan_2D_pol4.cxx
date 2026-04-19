#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

void _06_dR_minpt_Scan_2D_pol4() {
    // 1. Settings
    std::vector<TString> dR_dirs = {"dR_005","dR_007","dR_009","dR_011","dR_013","dR_015", "dR_017", "dR_019","dR_021"};
    std::vector<double> minpt_cuts = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<int> colors = {kRed, kOrange+1, kYellow+1, kSpring+2, kGreen+2, kCyan+1, kAzure-2, kBlue, kBlack};

    // Scan for available pT bins
    std::vector<TString> pt_bins;
    TString firstFname = dR_dirs[0] + "/Significance_Results_pol4_minpt.root";
    TFile *fFirst = TFile::Open(firstFname);
    if (!fFirst || fFirst->IsZombie()) fFirst = TFile::Open("../" + firstFname);
    if (fFirst && !fFirst->IsZombie()) {
        TString firstMinptDir = Form("dimuon_minpt_%.1f", minpt_cuts[0]);
        TDirectory *d = (TDirectory*)fFirst->Get(firstMinptDir);
        if (d) {
            TIter next(d->GetListOfKeys());
            TKey *key;
            while ((key = (TKey*)next())) {
                if (TString(key->GetClassName()) == "TDirectoryFile") {
                    TString name = key->GetName();
                    if (name.BeginsWith("Pt_")) pt_bins.push_back(name);
                }
            }
        }
        fFirst->Close();
    }

    if (pt_bins.empty()) {
        std::cout << "Error: No Pt bins found in " << firstFname << std::endl;
        // Fallback or exit
        pt_bins = {"Pt_0to30"}; 
    }

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".2f");

    TFile *outFile = new TFile("dR_minpt_Scan_2D_pol4.root", "RECREATE");

    struct PtData {
        double pt, pt_err;
        double val, err;
    };
    std::map<double, std::map<TString, std::map<TString, std::vector<PtData>>>> globalPtData;

    for (const auto& pt_bin : pt_bins) {
        TString pt_title = pt_bin;
        pt_title.ReplaceAll("Pt_", "p_{T}: "); pt_title.ReplaceAll("to", "-");
        
        TDirectory *ptDir = outFile->mkdir(pt_bin);
        ptDir->cd();

        TH2D *h2_sig_omega = new TH2D("h2_sig_omega", Form("Significance Map (#omega, %s, pol4);p_{T,#mu} (GeV/c);dR cut", pt_title.Data()), minpt_cuts.size(), 0.5, minpt_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_sig_phi = new TH2D("h2_sig_phi", Form("Significance Map (#phi, %s, pol4);p_{T,#mu} (GeV/c);dR cut", pt_title.Data()), minpt_cuts.size(), 0.5, minpt_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        
        for (int i = 0; i < (int)minpt_cuts.size(); ++i) {
            TString label = Form("p_{T,#mu} < %.1f", minpt_cuts[i]);
            h2_sig_omega->GetXaxis()->SetBinLabel(i+1, label); h2_sig_phi->GetXaxis()->SetBinLabel(i+1, label);
        }
        for (int j = 0; j < (int)dR_dirs.size(); ++j) {
            TString label = dR_dirs[j]; label.ReplaceAll("dR_0", "dR < 0.");
            h2_sig_omega->GetYaxis()->SetBinLabel(j+1, label); h2_sig_phi->GetYaxis()->SetBinLabel(j+1, label);
        }

        std::vector<std::vector<double>> vSigO(dR_dirs.size(), std::vector<double>(minpt_cuts.size(), 0));
        std::vector<std::vector<double>> vSigOErr(dR_dirs.size(), std::vector<double>(minpt_cuts.size(), 0));
        std::vector<std::vector<double>> vSigP(dR_dirs.size(), std::vector<double>(minpt_cuts.size(), 0));
        std::vector<std::vector<double>> vSigPErr(dR_dirs.size(), std::vector<double>(minpt_cuts.size(), 0));

        for (int j = 0; j < (int)dR_dirs.size(); ++j) {
            TString dR_dir = dR_dirs[j];
            TString fname = dR_dir + "/Significance_Results_pol4_minpt.root";
            TFile *fIn = TFile::Open(fname);
            if (!fIn || fIn->IsZombie()) fIn = TFile::Open("../" + fname);
            if (!fIn || fIn->IsZombie()) continue;

            for (int i = 0; i < (int)minpt_cuts.size(); ++i) {
                TString minptDirName = Form("dimuon_minpt_%.1f", minpt_cuts[i]);
                TString basePath = minptDirName + "/" + pt_bin + "/Significance_";

                TParameter<double> *pSigO = (TParameter<double>*)fIn->Get(basePath + "omega");
                TParameter<double> *pSigO_err = (TParameter<double>*)fIn->Get(basePath + "omega_err");
                TParameter<double> *pSigP = (TParameter<double>*)fIn->Get(basePath + "phi");
                TParameter<double> *pSigP_err = (TParameter<double>*)fIn->Get(basePath + "phi_err");

                if (pSigO) {
                    h2_sig_omega->SetBinContent(i+1, j+1, pSigO->GetVal());
                    if (pSigO_err) h2_sig_omega->SetBinError(i+1, j+1, pSigO_err->GetVal());
                    vSigO[j][i] = pSigO->GetVal(); vSigOErr[j][i] = pSigO_err ? pSigO_err->GetVal() : 0;
                }
                if (pSigP) {
                    h2_sig_phi->SetBinContent(i+1, j+1, pSigP->GetVal());
                    if (pSigP_err) h2_sig_phi->SetBinError(i+1, j+1, pSigP_err->GetVal());
                    vSigP[j][i] = pSigP->GetVal(); vSigPErr[j][i] = pSigP_err ? pSigP_err->GetVal() : 0;
                }
            }
            fIn->Close();
        }

        ptDir->cd();
        auto DrawWithValues = [&](TH2D* h, TString name) {
            TCanvas *c = new TCanvas(Form("c_%s_%s_pol4", name.Data(), pt_bin.Data()), h->GetTitle(), 1200, 900);
            c->SetRightMargin(0.15); c->SetLeftMargin(0.15); c->SetBottomMargin(0.15);
            h->SetMarkerSize(1.5); h->Draw("COLZ TEXT");
            c->Write(); delete c;
        };
        DrawWithValues(h2_sig_omega, "h2_sig_omega"); DrawWithValues(h2_sig_phi, "h2_sig_phi");
    }

    outFile->Close();
    std::cout << "Finish! Result saved in dR_minpt_Scan_2D_pol4.root" << std::endl;
}
