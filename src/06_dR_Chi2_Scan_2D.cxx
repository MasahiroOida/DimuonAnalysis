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

void _06_dR_Chi2_Scan_2D() {
    // 1. 設定
    std::vector<TString> dR_dirs = {"dR_005","dR_007","dR_009","dR_011","dR_013","dR_015", "dR_017", "dR_019","dR_021"};
    std::vector<double> chi2_cuts = {10, 20, 30, 40, 50, 60, 70, 80};
    // Spectral color palette (color-pencil order)
    std::vector<int> colors = {kRed, kOrange+1, kYellow+1, kSpring+2, kGreen+2, kCyan+1, kAzure-2, kBlue, kBlack};

    // Scan for available pT bins from the first available file
    std::vector<TString> pt_bins;
    TString firstFname = dR_dirs[0] + "/Significance_Results.root";
    TFile *fFirst = TFile::Open(firstFname);
    if (!fFirst || fFirst->IsZombie()) fFirst = TFile::Open("../" + firstFname);
    if (fFirst && !fFirst->IsZombie()) {
        TString firstChi2Dir = Form("dimuon_Chi2_%.0f", chi2_cuts[0]);
        TDirectory *d = (TDirectory*)fFirst->Get(firstChi2Dir);
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
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".2f");

    TFile *outFile = new TFile("dR_Chi2_Scan_2D.root", "RECREATE");

    // Structure to collect data for pT-dependence plots
    struct PtData {
        double pt, pt_err;
        double val, err;
    };
    // map[chi2_cut][dR_dir][variable_name] -> vector<PtData>
    std::map<double, std::map<TString, std::map<TString, std::vector<PtData>>>> globalPtData;

    for (const auto& pt_bin : pt_bins) {
        TString pt_title = pt_bin;
        pt_title.ReplaceAll("Pt_", "p_{T}: "); pt_title.ReplaceAll("to", "-");
        
        TDirectory *ptDir = outFile->mkdir(pt_bin);
        ptDir->cd();

        // 2D Maps (Keep existing ones)
        TH2D *h2_sig_omega = new TH2D("h2_sig_omega", Form("Significance Map (#omega, %s);#chi^{2} cut;dR cut", pt_title.Data()), chi2_cuts.size(), 0.5, chi2_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_sig_phi = new TH2D("h2_sig_phi", Form("Significance Map (#phi, %s);#chi^{2} cut;dR cut", pt_title.Data()), chi2_cuts.size(), 0.5, chi2_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_mean_omega = new TH2D("h2_mean_omega", Form("Mean Map (#omega, %s);#chi^{2} cut;dR cut", pt_title.Data()), chi2_cuts.size(), 0.5, chi2_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_mean_phi = new TH2D("h2_mean_phi", Form("Mean Map (#phi, %s);#chi^{2} cut;dR cut", pt_title.Data()), chi2_cuts.size(), 0.5, chi2_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_sigma_omega = new TH2D("h2_sigma_omega", Form("Width Map (#omega, %s);#chi^{2} cut;dR cut", pt_title.Data()), chi2_cuts.size(), 0.5, chi2_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);
        TH2D *h2_sigma_phi = new TH2D("h2_sigma_phi", Form("Width Map (#phi, %s);#chi^{2} cut;dR cut", pt_title.Data()), chi2_cuts.size(), 0.5, chi2_cuts.size() + 0.5, dR_dirs.size(), 0.5, dR_dirs.size() + 0.5);

        for (int i = 0; i < (int)chi2_cuts.size(); ++i) {
            TString label = Form("#chi^{2} < %.0f", chi2_cuts[i]);
            h2_sig_omega->GetXaxis()->SetBinLabel(i+1, label); h2_sig_phi->GetXaxis()->SetBinLabel(i+1, label);
            h2_mean_omega->GetXaxis()->SetBinLabel(i+1, label); h2_mean_phi->GetXaxis()->SetBinLabel(i+1, label);
            h2_sigma_omega->GetXaxis()->SetBinLabel(i+1, label); h2_sigma_phi->GetXaxis()->SetBinLabel(i+1, label);
        }
        for (int j = 0; j < (int)dR_dirs.size(); ++j) {
            TString label = dR_dirs[j]; label.ReplaceAll("dR_0", "dR < 0.");
            h2_sig_omega->GetYaxis()->SetBinLabel(j+1, label); h2_sig_phi->GetYaxis()->SetBinLabel(j+1, label);
            h2_mean_omega->GetYaxis()->SetBinLabel(j+1, label); h2_mean_phi->GetYaxis()->SetBinLabel(j+1, label);
            h2_sigma_omega->GetYaxis()->SetBinLabel(j+1, label); h2_sigma_phi->GetYaxis()->SetBinLabel(j+1, label);
        }

        // Data containers for 1D plots: [dR_idx][chi2_idx]
        std::vector<std::vector<double>> vSigO(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vSigOErr(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vSigP(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vSigPErr(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        
        std::vector<std::vector<double>> vMeanO(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vMeanOErr(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vMeanP(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vMeanPErr(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));

        std::vector<std::vector<double>> vSigmaO(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vSigmaOErr(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vSigmaP(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));
        std::vector<std::vector<double>> vSigmaPErr(dR_dirs.size(), std::vector<double>(chi2_cuts.size(), 0));

        // 3. データの読み込み
        for (int j = 0; j < (int)dR_dirs.size(); ++j) {
            TString dR_dir = dR_dirs[j];
            TString fname = dR_dir + "/Significance_Results.root";
            TFile *fIn = TFile::Open(fname);
            if (!fIn || fIn->IsZombie()) fIn = TFile::Open("../" + fname);
            if (!fIn || fIn->IsZombie()) continue;

            for (int i = 0; i < (int)chi2_cuts.size(); ++i) {
                TString chi2DirName = Form("dimuon_Chi2_%.0f", chi2_cuts[i]);
                TString basePath = chi2DirName + "/" + pt_bin + "/Significance_";
                TString meanPath = chi2DirName + "/" + pt_bin + "/Mean_";
                TString sigmaPath = chi2DirName + "/" + pt_bin + "/Sigma_";

                TParameter<double> *pSigO = (TParameter<double>*)fIn->Get(basePath + "omega");
                TParameter<double> *pSigO_err = (TParameter<double>*)fIn->Get(basePath + "omega_err");
                TParameter<double> *pSigP = (TParameter<double>*)fIn->Get(basePath + "phi");
                TParameter<double> *pSigP_err = (TParameter<double>*)fIn->Get(basePath + "phi_err");
                TParameter<double> *pMeanO = (TParameter<double>*)fIn->Get(meanPath + "omega");
                TParameter<double> *pMeanO_err = (TParameter<double>*)fIn->Get(meanPath + "omega_err");
                TParameter<double> *pMeanP = (TParameter<double>*)fIn->Get(meanPath + "phi");
                TParameter<double> *pMeanP_err = (TParameter<double>*)fIn->Get(meanPath + "phi_err");
                TParameter<double> *pSigmaO = (TParameter<double>*)fIn->Get(sigmaPath + "omega");
                TParameter<double> *pSigmaO_err = (TParameter<double>*)fIn->Get(sigmaPath + "omega_err");
                TParameter<double> *pSigmaP = (TParameter<double>*)fIn->Get(sigmaPath + "phi");
                TParameter<double> *pSigmaP_err = (TParameter<double>*)fIn->Get(sigmaPath + "phi_err");

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
                if (pMeanO) {
                    h2_mean_omega->SetBinContent(i+1, j+1, pMeanO->GetVal());
                    if (pMeanO_err) h2_mean_omega->SetBinError(i+1, j+1, pMeanO_err->GetVal());
                    vMeanO[j][i] = pMeanO->GetVal(); vMeanOErr[j][i] = pMeanO_err ? pMeanO_err->GetVal() : 0;
                }
                if (pMeanP) {
                    h2_mean_phi->SetBinContent(i+1, j+1, pMeanP->GetVal());
                    if (pMeanP_err) h2_mean_phi->SetBinError(i+1, j+1, pMeanP_err->GetVal());
                    vMeanP[j][i] = pMeanP->GetVal(); vMeanPErr[j][i] = pMeanP_err ? pMeanP_err->GetVal() : 0;
                }
                if (pSigmaO) {
                    h2_sigma_omega->SetBinContent(i+1, j+1, pSigmaO->GetVal());
                    if (pSigmaO_err) h2_sigma_omega->SetBinError(i+1, j+1, pSigmaO_err->GetVal());
                    vSigmaO[j][i] = pSigmaO->GetVal(); vSigmaOErr[j][i] = pSigmaO_err ? pSigmaO_err->GetVal() : 0;
                }
                if (pSigmaP) {
                    h2_sigma_phi->SetBinContent(i+1, j+1, pSigmaP->GetVal());
                    if (pSigmaP_err) h2_sigma_phi->SetBinError(i+1, j+1, pSigmaP_err->GetVal());
                    vSigmaP[j][i] = pSigmaP->GetVal(); vSigmaPErr[j][i] = pSigmaP_err ? pSigmaP_err->GetVal() : 0;
                }

                // Collect data for global Pt dependence plots
                TParameter<double> *pPtMin = (TParameter<double>*)fIn->Get(chi2DirName + "/" + pt_bin + "/ptmin");
                TParameter<double> *pPtMax = (TParameter<double>*)fIn->Get(chi2DirName + "/" + pt_bin + "/ptmax");
                if (pPtMin && pPtMax) {
                    double pt_val = (pPtMin->GetVal() + pPtMax->GetVal()) / 2.0;
                    double pt_err = (pPtMax->GetVal() - pPtMin->GetVal()) / 2.0;
                    double chi2 = chi2_cuts[i];
                    TString dR = dR_dirs[j];

                    if (pMeanO)  globalPtData[chi2][dR]["mean_omega"].push_back({pt_val, pt_err, pMeanO->GetVal(), pMeanO_err ? pMeanO_err->GetVal() : 0});
                    if (pMeanP)  globalPtData[chi2][dR]["mean_phi"].push_back({pt_val, pt_err, pMeanP->GetVal(), pMeanP_err ? pMeanP_err->GetVal() : 0});
                    if (pSigmaO) globalPtData[chi2][dR]["sigma_omega"].push_back({pt_val, pt_err, pSigmaO->GetVal(), pSigmaO_err ? pSigmaO_err->GetVal() : 0});
                    if (pSigmaP) globalPtData[chi2][dR]["sigma_phi"].push_back({pt_val, pt_err, pSigmaP->GetVal(), pSigmaP_err ? pSigmaP_err->GetVal() : 0});
                    if (pSigO)   globalPtData[chi2][dR]["sig_omega"].push_back({pt_val, pt_err, pSigO->GetVal(), pSigO_err ? pSigO_err->GetVal() : 0});
                    if (pSigP)   globalPtData[chi2][dR]["sig_phi"].push_back({pt_val, pt_err, pSigP->GetVal(), pSigP_err ? pSigP_err->GetVal() : 0});
                }
            }
            fIn->Close();
        }

        // 4. TH2D Draw
        auto DrawWithValues = [&](TH2D* h, TString name, TString format = ".2f") {
            TCanvas *c = new TCanvas(Form("c_%s_%s", name.Data(), pt_bin.Data()), h->GetTitle(), 1200, 900);
            c->SetRightMargin(0.15); c->SetLeftMargin(0.15); c->SetBottomMargin(0.15);
            gStyle->SetPaintTextFormat(format.Data());
            h->SetMarkerSize(1.5); h->Draw("COLZ TEXT");
            c->Write(); delete c;
        };
        ptDir->cd();
        DrawWithValues(h2_sig_omega, "h2_sig_omega", ".2f"); DrawWithValues(h2_sig_phi, "h2_sig_phi", ".2f");
        DrawWithValues(h2_mean_omega, "h2_mean_omega", ".4f"); DrawWithValues(h2_mean_phi, "h2_mean_phi", ".4f");
        DrawWithValues(h2_sigma_omega, "h2_sigma_omega", ".4f"); DrawWithValues(h2_sigma_phi, "h2_sigma_phi", ".4f");

        // 5. Comparison 1D plots (Value vs Chi2, colors for dR)
        auto DrawComparison = [&](TString name, TString title, std::vector<std::vector<double>>& vy, std::vector<std::vector<double>>& vyer) {
            TCanvas *cComp = new TCanvas("cComp_"+name+"_"+pt_bin, title, 1000, 700);
            cComp->SetGrid();
            TLegend *leg = new TLegend(0.7, 0.6, 0.88, 0.88);
            double maxY = -1e9, minY = 1e9;
            for(auto& row : vy) for(auto& val : row) {
                if(val > maxY) maxY = val;
                if(val < minY && val != 0) minY = val;
            }
            if (minY > maxY) { minY = 0; maxY = 1; }

            for (int j = 0; j < (int)dR_dirs.size(); ++j) {
                TGraphErrors *gr = new TGraphErrors(chi2_cuts.size(), &chi2_cuts[0], &vy[j][0], 0, &vyer[j][0]);
                int color = colors[j % colors.size()];
                gr->SetMarkerStyle(20); gr->SetMarkerColor(color); gr->SetLineColor(color);
                if (j == 0) {
                    gr->SetTitle(title + " vs #chi^{2} cut (" + pt_title + ");#chi^{2} cut;" + title);
                    
                    gr->Draw("APE");

                    // Set Y-axis range based on the plot name after Draw
                    if (name.Contains("sig")) gr->GetYaxis()->SetRangeUser(0, 200);
                    else if (name.Contains("mean_omega")) gr->GetYaxis()->SetRangeUser(0.76, 0.79);
                    else if (name.Contains("mean_phi")) gr->GetYaxis()->SetRangeUser(1.01, 1.025);
                    else if (name.Contains("sigma")) gr->GetYaxis()->SetRangeUser(0.015, 0.035);
                    else gr->GetYaxis()->SetRangeUser(minY * 0.9, maxY * 1.1);

                    cComp->Update();
                } else {
                    gr->Draw("PE SAME");
                }
                TString label = dR_dirs[j]; label.ReplaceAll("dR_0", "dR < 0.");
                leg->AddEntry(gr, label, "lp");
            }
            leg->Draw(); cComp->Write(); delete cComp;
        };
        DrawComparison("sig_omega", "Significance (#omega)", vSigO, vSigOErr);
        DrawComparison("sig_phi", "Significance (#phi)", vSigP, vSigPErr);
        DrawComparison("mean_omega", "Mean (#omega)", vMeanO, vMeanOErr);
        DrawComparison("mean_phi", "Mean (#phi)", vMeanP, vMeanPErr);
        DrawComparison("sigma_omega", "Sigma (#omega)", vSigmaO, vSigmaOErr);
        DrawComparison("sigma_phi", "Sigma (#phi)", vSigmaP, vSigmaPErr);

        h2_sig_omega->Write(); h2_sig_phi->Write();
        h2_mean_omega->Write(); h2_mean_phi->Write();
        h2_sigma_omega->Write(); h2_sigma_phi->Write();
    }

    // 6. Create Pt-dependence plots for each chi2 cut
    TDirectory *ptDepDir = outFile->mkdir("Pt_Dependence");
    for (double chi2 : chi2_cuts) {
        TDirectory *chi2Dir = ptDepDir->mkdir(Form("Chi2_%.0f", chi2));
        chi2Dir->cd();

        std::vector<TString> variables = {"sig_omega", "sig_phi", "mean_omega", "mean_phi", "sigma_omega", "sigma_phi"};
        std::map<TString, TString> titles = {
            {"sig_omega", "Significance (#omega)"}, {"sig_phi", "Significance (#phi)"},
            {"mean_omega", "Mean Mass (#omega)"}, {"mean_phi", "Mean Mass (#phi)"},
            {"sigma_omega", "Mass Width (#omega)"}, {"sigma_phi", "Mass Width (#phi)"}
        };

        for (const auto& var : variables) {
            TCanvas *cPt = new TCanvas(Form("c_pt_%s_chi2%.0f", var.Data(), chi2), titles[var], 1000, 700);
            cPt->SetGrid();
            TLegend *leg = new TLegend(0.7, 0.7, 0.88, 0.88);
            leg->SetHeader(Form("#chi^{2} < %.0f", chi2));

            bool first = true;
            for (int j = 0; j < (int)dR_dirs.size(); ++j) {
                TString dR = dR_dirs[j];
                auto& points = globalPtData[chi2][dR][var];
                if (points.empty()) continue;

                // Sort points by pt
                std::sort(points.begin(), points.end(), [](const PtData& a, const PtData& b) { return a.pt < b.pt; });

                int n = points.size();
                std::vector<double> x(n), xe(n), y(n), ye(n);
                for (int i = 0; i < n; ++i) {
                    x[i] = points[i].pt; xe[i] = points[i].pt_err;
                    y[i] = points[i].val; ye[i] = points[i].err;
                }

                TGraphErrors *gr = new TGraphErrors(n, &x[0], &y[0], &xe[0], &ye[0]);
                int color = colors[j % colors.size()];
                gr->SetMarkerStyle(20); gr->SetMarkerColor(color); gr->SetLineColor(color);
                
                if (first) {
                    gr->SetTitle(titles[var] + " vs p_{T} (#chi^{2} < " + Form("%.0f", chi2) + ");p_{T} (GeV/c);" + titles[var]);
                    gr->Draw("APE");
                    
                    // Apply requested Y-axis ranges after Draw
                    if (var.Contains("sig")) gr->GetYaxis()->SetRangeUser(0, 200);
                    else if (var.Contains("mean_omega")) gr->GetYaxis()->SetRangeUser(0.76, 0.79);
                    else if (var.Contains("mean_phi")) gr->GetYaxis()->SetRangeUser(1.01, 1.025);
                    else if (var.Contains("sigma")) gr->GetYaxis()->SetRangeUser(0.015, 0.035);
                    
                    cPt->Update();
                    first = false;
                } else {
                    gr->Draw("PE SAME");
                }
                TString label = dR; label.ReplaceAll("dR_0", "dR < 0.");
                leg->AddEntry(gr, label, "lp");
            }
            leg->Draw();
            cPt->Write();
            delete cPt;
        }
    }

    outFile->Close();
    std::cout << "\nFinish! Result saved in dR_Chi2_Scan_2D.root with comparison plots." << std::endl;
}
