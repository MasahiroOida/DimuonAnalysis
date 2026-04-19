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

void _06_dR_MCHnCluster_Scan_2D() {
    std::vector<TString> dR_dirs = {"dR_005","dR_007","dR_009","dR_011","dR_013","dR_015", "dR_017", "dR_019","dR_021"};
    std::vector<double> cluster_cuts = {2, 3, 4, 5, 6}; // Typical MCH cluster cuts
    std::vector<int> colors = {kRed, kOrange+1, kYellow+1, kSpring+2, kGreen+2, kCyan+1, kAzure-2, kBlue, kBlack};

    std::vector<TString> pt_bins;
    TString firstFname = dR_dirs[0] + "/Significance_Results_MCHnCluster.root";
    TFile *fFirst = TFile::Open(firstFname);
    if (!fFirst || fFirst->IsZombie()) fFirst = TFile::Open("../" + firstFname);
    if (fFirst && !fFirst->IsZombie()) {
        TString firstClusterDir = "";
        TIter nextKey(fFirst->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)nextKey())) {
            if (TString(key->GetClassName()) == "TDirectoryFile") {
                TString name = key->GetName();
                if (name.Contains("cluster", TString::kIgnoreCase)) {
                    firstClusterDir = name;
                    break;
                }
            }
        }
        
        if (firstClusterDir != "") {
            TDirectory *d = (TDirectory*)fFirst->Get(firstClusterDir);
            TIter next(d->GetListOfKeys());
            TKey *k;
            while ((k = (TKey*)next())) {
                if (TString(k->GetClassName()) == "TDirectoryFile") {
                    TString name = k->GetName();
                    if (name.BeginsWith("Pt_")) pt_bins.push_back(name);
                }
            }
        }
        fFirst->Close();
    }

    if (pt_bins.empty()) {
        std::cout << "Error: No Pt bins found in " << firstFname << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    TFile *outFile = new TFile("dR_MCHnCluster_Scan_2D.root", "RECREATE");

    for (const auto& pt_bin : pt_bins) {
        TString pt_title = pt_bin;
        pt_title.ReplaceAll("Pt_", "p_{T}: "); pt_title.ReplaceAll("to", "-");
        
        TDirectory *ptDir = outFile->mkdir(pt_bin);
        ptDir->cd();

        TH2D *h2_sig_omega = new TH2D("h2_sig_omega", Form("Significance (#omega, %s);MCHnCluster cut;dR cut", pt_title.Data()), cluster_cuts.size(), 0.5, cluster_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_sig_phi = new TH2D("h2_sig_phi", Form("Significance (#phi, %s);MCHnCluster cut;dR cut", pt_title.Data()), cluster_cuts.size(), 0.5, cluster_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);

        for (int i = 0; i < (int)cluster_cuts.size(); ++i) {
            h2_sig_omega->GetXaxis()->SetBinLabel(i+1, Form(">= %.0f", cluster_cuts[i]));
            h2_sig_phi->GetXaxis()->SetBinLabel(i+1, Form(">= %.0f", cluster_cuts[i]));
        }
        for (int j = 0; j < (int)dR_dirs.size(); ++j) {
            TString label = dR_dirs[j]; label.ReplaceAll("dR_0", "dR < 0.");
            h2_sig_omega->GetYaxis()->SetBinLabel(j+1, label);
            h2_sig_phi->GetYaxis()->SetBinLabel(j+1, label);
        }

        std::vector<std::vector<double>> vSigO(dR_dirs.size(), std::vector<double>(cluster_cuts.size(), 0));
        std::vector<std::vector<double>> vSigOErr(dR_dirs.size(), std::vector<double>(cluster_cuts.size(), 0));
        std::vector<std::vector<double>> vSigP(dR_dirs.size(), std::vector<double>(cluster_cuts.size(), 0));
        std::vector<std::vector<double>> vSigPErr(dR_dirs.size(), std::vector<double>(cluster_cuts.size(), 0));

        for (int j = 0; j < (int)dR_dirs.size(); ++j) {
            TString fname = dR_dirs[j] + "/Significance_Results_MCHnCluster.root";
            TFile *fIn = TFile::Open(fname);
            if (!fIn || fIn->IsZombie()) fIn = TFile::Open("../" + fname);
            if (!fIn || fIn->IsZombie()) continue;

            for (int i = 0; i < (int)cluster_cuts.size(); ++i) {
                // Find directory for this cluster cut
                TString targetDir = "";
                TIter nextK(fIn->GetListOfKeys());
                TKey *tk;
                while((tk = (TKey*)nextK())){
                    TString kname = tk->GetName();
                    // Matching cluster_N format
                    if((kname.Contains("cluster", TString::kIgnoreCase) || kname.Contains("mch", TString::kIgnoreCase)) && 
                       (kname.EndsWith(Form("_%.0f", cluster_cuts[i])) || kname.Contains(Form("_%.0f_", cluster_cuts[i])))){
                        targetDir = kname;
                        break;
                    }
                }
                if(targetDir == "") continue;

                TString basePath = targetDir + "/" + pt_bin + "/";
                TParameter<double> *pSigO = (TParameter<double>*)fIn->Get(basePath + "Significance_omega");
                TParameter<double> *pSigO_err = (TParameter<double>*)fIn->Get(basePath + "Significance_omega_err");
                TParameter<double> *pSigP = (TParameter<double>*)fIn->Get(basePath + "Significance_phi");
                TParameter<double> *pSigP_err = (TParameter<double>*)fIn->Get(basePath + "Significance_phi_err");

                if (pSigO) {
                    h2_sig_omega->SetBinContent(i+1, j+1, pSigO->GetVal());
                    if (pSigO_err) h2_sig_omega->SetBinError(i+1, j+1, pSigO_err->GetVal());
                    vSigO[j][i] = pSigO->GetVal();
                    vSigOErr[j][i] = pSigO_err ? pSigO_err->GetVal() : 0;
                }
                if (pSigP) {
                    h2_sig_phi->SetBinContent(i+1, j+1, pSigP->GetVal());
                    if (pSigP_err) h2_sig_phi->SetBinError(i+1, j+1, pSigP_err->GetVal());
                    vSigP[j][i] = pSigP->GetVal();
                    vSigPErr[j][i] = pSigP_err ? pSigP_err->GetVal() : 0;
                }
            }
            fIn->Close();
        }

        ptDir->cd();
        h2_sig_omega->Write();
        h2_sig_phi->Write();
    }

    outFile->Close();
    std::cout << "Finish! Result saved in dR_MCHnCluster_Scan_2D.root" << std::endl;
}
