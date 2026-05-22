#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TMath.h>
#include <TParameter.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TGraphErrors.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <map>
#include <algorithm>
#include <TObjString.h>

const double EPSILON = 1e-7;

double GetIntegralError(TF1 *f, double low, double high, TFitResult *res, const std::vector<int>& indices) {
    if (!res || !res->IsValid()) {
        return 0.0;
    }
    
    TMatrixDSym cov = res->GetCovarianceMatrix();
    int nRows = cov.GetNrows();
    if (nRows == 0) {
        return 0.0;
    }
    
    int n = indices.size();
    for(int i=0; i<n; ++i) {
        if (indices[i] < 0 || indices[i] >= nRows) {
            return 0.0;
        }
    }
    
    TMatrixDSym subCov(n);
    for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j) {
            subCov(i, j) = cov(indices[i], indices[j]);
        }
    }
    
    double err = f->IntegralError(low, high, f->GetParameters(), subCov.GetMatrixArray(), 1e-2);
    if (std::isnan(err) || std::isinf(err)) return 0.0;
    return err;
}

double CalcSignificanceError(double S, double S_err, double T, double T_err, double CovST) {
    double safeT = std::max(T, EPSILON);
    double termS = pow(S_err, 2) / safeT;
    double termT = pow(S, 2) * pow(T_err, 2) / (4.0 * pow(safeT, 3));
    double termCov = - S * CovST / pow(safeT, 2);
    return sqrt(std::max(0.0, termS + termT + termCov));
}

void Significance_CrystalBall(TFile *input_FitResult) {
    if (!input_FitResult || input_FitResult->IsZombie()) return;
    TFile *outfile = new TFile("Significance_Results.root", "RECREATE");
    
    std::set<std::string> validIDs;
    TIter scanIter(input_FitResult->GetListOfKeys());
    TKey *scanKey;
    while ((scanKey = (TKey *)scanIter())) {
        if (TString(scanKey->GetClassName()) == "TDirectoryFile") {
            TString dName = scanKey->GetName();
            if (dName.Contains("chi2", TString::kIgnoreCase)) {
                TObjArray *tokens = dName.Tokenize("_");
                for (int i = 0; i < tokens->GetEntries(); ++i) {
                    TString tok = tokens->At(i)->GetName();
                    TString tokLow = tok; tokLow.ToLower();
                    if (tokLow.BeginsWith("id")) {
                        validIDs.insert(tokLow.Data());
                        break;
                    }
                }
                delete tokens;
            }
        }
    }

    std::cout << "--- Valid IDs Found ---" << std::endl;
    for (const auto& vid : validIDs) std::cout << " " << vid;
    std::cout << "\n-----------------------" << std::endl;

    struct SummaryData {
        TString dirName;
        TString ptName;
        double chi2Val;
        double ptMin, ptMax;
        double o_sig, o_sig_err, p_sig, p_sig_err;
        double o_m, o_m_err, o_s, o_s_err;
        double p_m, p_m_err, p_s, p_s_err;
        TH1F* hSig;
    };
    std::vector<SummaryData> summaryDataList;

    TIter next(input_FitResult->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (TString(key->GetClassName()) == "TDirectoryFile") {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TString topName = topDir->GetName();

            bool process = false;
            if (topName.Contains("chi2", TString::kIgnoreCase)) {
                process = true; 
            } else if (topName.Contains("id", TString::kIgnoreCase)) {
                TObjArray *tokens = topName.Tokenize("_");
                for (int i = 0; i < tokens->GetEntries(); ++i) {
                    TString tok = tokens->At(i)->GetName();
                    TString tokLow = tok; tokLow.ToLower();
                    if (tokLow.BeginsWith("id")) {
                        if (validIDs.count(tokLow.Data())) process = true;
                        break;
                    }
                }
                delete tokens;
            }

            if (!process) {
                delete topDir;
                continue;
            }
            std::cout << "Processing directory: " << topName << std::endl;

            TDirectory *outTopDir = outfile->mkdir(topName);
            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt())) {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile") {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    TString ptName = ptDir->GetName();
                    TDirectory *outPtDir = outTopDir->mkdir(ptName);
                    outPtDir->cd();

                    TH1F *hSig = (TH1F *)ptDir->Get("LikeSignSig_fit");
                    TF1 *totalfit = (TF1 *)ptDir->Get("totalfit");
                    TFitResult *fitRes = (TFitResult*)ptDir->Get("fitResult");
                    if (!hSig || !totalfit || !fitRes) continue;

                    TParameter<double> *pPtMin = (TParameter<double>*)ptDir->Get("ptmin");
                    TParameter<double> *pPtMax = (TParameter<double>*)ptDir->Get("ptmax");
                    double ptMin = pPtMin ? pPtMin->GetVal() : -1;
                    double ptMax = pPtMax ? pPtMax->GetVal() : -1;

                    double binWidth = hSig->GetBinWidth(1);
                    double o_m = totalfit->GetParameter(3); double o_s = totalfit->GetParameter(4);
                    double p_m = totalfit->GetParameter(8); double p_s = totalfit->GetParameter(9);

                    bool physicallyReasonable = (o_m > 0.6 && o_m < 0.9 && p_m > 0.9 && p_m < 1.2 && o_s > 0 && p_s > 0);
                    if (!physicallyReasonable) {
                        continue;
                    }

                    double bg_p0 = totalfit->GetParameter(0);
                    double bg_p1 = totalfit->GetParameter(1);

                    TF1 *fBG = new TF1("fBG", "[0]*exp(-[1]*x)", 0.6, 1.3);
                    fBG->SetParameters(bg_p0, bg_p1);

                    TH1F *hPureSig = (TH1F*)hSig->Clone("hPureSig");
                    TH1F *hResBG = (TH1F*)hSig->Clone("hResBG");
                    hPureSig->Reset(); hResBG->Reset();

                    for (int bin=1; bin<=hSig->GetNbinsX(); ++bin) {
                        double x = hSig->GetBinCenter(bin);
                        double val_total = hSig->GetBinContent(bin);
                        double err_total = hSig->GetBinError(bin);
                        double val_bg = fBG->Eval(x);
                        hResBG->SetBinContent(bin, val_bg);
                        hResBG->SetBinError(bin, 0); 
                        hPureSig->SetBinContent(bin, val_total - val_bg);
                        hPureSig->SetBinError(bin, err_total); 
                    }

                    double o_low = o_m - 3.0*o_s; double o_high = o_m + 3.0*o_s;
                    double p_low = p_m - 3.0*p_s; double p_high = p_m + 3.0*p_s;

                    double S_o, S_o_err, B_o, B_o_err;
                    double S_p, S_p_err, B_p, B_p_err;

                    int o_bin_low = hPureSig->FindBin(o_low); int o_bin_high = hPureSig->FindBin(o_high);
                    int p_bin_low = hPureSig->FindBin(p_low); int p_bin_high = hPureSig->FindBin(p_high);

                    S_o = hPureSig->IntegralAndError(o_bin_low, o_bin_high, S_o_err);
                    B_o = hResBG->IntegralAndError(o_bin_low, o_bin_high, B_o_err);
                    S_p = hPureSig->IntegralAndError(p_bin_low, p_bin_high, S_p_err);
                    B_p = hResBG->IntegralAndError(p_bin_low, p_bin_high, B_p_err);

                    double T_o = S_o + B_o;
                    double T_p = S_p + B_p;
                    double T_o_err = sqrt(pow(S_o_err, 2) + pow(B_o_err, 2));
                    double T_p_err = sqrt(pow(S_p_err, 2) + pow(B_p_err, 2));

                    if (T_o <= 0 || T_p <= 0 || S_o < 0 || S_p < 0) {
                        delete fBG; delete hPureSig; delete hResBG;
                        continue;
                    }

                    double sig_o = S_o / sqrt(T_o);
                    double sig_p = S_p / sqrt(T_p);
                    double sig_o_err = CalcSignificanceError(S_o, S_o_err, T_o, T_o_err, 0.0);
                    double sig_p_err = CalcSignificanceError(S_p, S_p_err, T_p, T_p_err, 0.0);

                    delete fBG; delete hPureSig; delete hResBG;

                    (new TParameter<double>("Significance_omega", sig_o))->Write();
                    (new TParameter<double>("Significance_omega_err", sig_o_err))->Write();
                    (new TParameter<double>("Significance_phi", sig_p))->Write();
                    (new TParameter<double>("Significance_phi_err", sig_p_err))->Write();
                    (new TParameter<double>("Mean_omega", o_m))->Write();
                    (new TParameter<double>("Mean_omega_err", totalfit->GetParError(3)))->Write();
                    (new TParameter<double>("Mean_phi", p_m))->Write();
                    (new TParameter<double>("Mean_phi_err", totalfit->GetParError(8)))->Write();
                    (new TParameter<double>("Sigma_omega", o_s))->Write();
                    (new TParameter<double>("Sigma_omega_err", totalfit->GetParError(4)))->Write();
                    (new TParameter<double>("Sigma_phi", p_s))->Write();
                    (new TParameter<double>("Sigma_phi_err", totalfit->GetParError(9)))->Write();

                    SummaryData data;
                    data.dirName = topName;
                    data.ptName = ptName;
                    data.ptMin = ptMin;
                    data.ptMax = ptMax;
                    data.chi2Val = 999999.0;
                    TObjArray *tokens = topName.Tokenize("_");
                    for (int i = 0; i < tokens->GetEntries(); ++i) {
                        TString tok = tokens->At(i)->GetName();
                        TString tokLow = tok; tokLow.ToLower();
                        if (tokLow == "chi2" && i + 1 < tokens->GetEntries()) {
                            data.chi2Val = TString(tokens->At(i+1)->GetName()).Atof();
                            break;
                        }
                    }
                    delete tokens;

                    data.o_sig = sig_o; data.o_sig_err = sig_o_err;
                    data.p_sig = sig_p; data.p_sig_err = sig_p_err;
                    data.o_m = o_m; data.o_m_err = totalfit->GetParError(3);
                    data.o_s = o_s; data.o_s_err = totalfit->GetParError(4);
                    data.p_m = p_m; data.p_m_err = totalfit->GetParError(8);
                    data.p_s = p_s; data.p_s_err = totalfit->GetParError(9);

                    TH1F *h = (TH1F*)hSig->Clone(Form("hSig_%s_%s", topName.Data(), ptName.Data()));
                    h->SetTitle(Form("%s (%s)", topName.Data(), ptName.Data())); h->SetDirectory(0);
                    h->GetListOfFunctions()->Clear();
                    data.hSig = h;

                    summaryDataList.push_back(data);
                }
            }
        }
    }

    if (!summaryDataList.empty()) {
        outfile->cd();
        gStyle->SetOptStat(0);

        auto formatLabel = [](TString dirName) -> TString {
            TString dirNameLow = dirName; dirNameLow.ToLower();
            if (!dirNameLow.Contains("chi2")) return "nocut";
            TObjArray *tokens = dirName.Tokenize("_");
            TString label = "nocut";
            for (int i = 0; i < tokens->GetEntries(); ++i) {
                TString tok = tokens->At(i)->GetName();
                TString tokLow = tok; tokLow.ToLower();
                if (tokLow == "chi2" && i + 1 < tokens->GetEntries()) {
                    TString chi2Val = tokens->At(i+1)->GetName();
                    label = Form("#chi^{2} < %s", chi2Val.Data());
                    break;
                }
            }
            delete tokens;
            return label;
        };

        // Group by Pt
        std::map<TString, std::vector<SummaryData>> dataByPt;
        std::map<TString, std::vector<SummaryData>> dataByID;
        for (const auto& d : summaryDataList) {
            dataByPt[d.ptName].push_back(d);
            dataByID[d.dirName].push_back(d);
        }

        // 1. Summary plots for each Pt bin (Significance vs Chi2)
        for (auto& entry : dataByPt) {
            TString ptName = entry.first;
            auto& list = entry.second;
            std::sort(list.begin(), list.end(), [](const SummaryData& a, const SummaryData& b) {
                return a.chi2Val < b.chi2Val;
            });

            TDirectory *ptSummaryDir = outfile->mkdir(Form("Summary_%s", ptName.Data()));
            ptSummaryDir->cd();
            int nData = list.size();

            auto saveSummary = [&](TString name, TString title, auto getVal, auto getErr, int color) {
                TCanvas *c = new TCanvas("c_"+name, title, 800, 600);
                c->SetGrid();
                TH1F *h = new TH1F("h_"+name, title+";;Value", nData, 0, nData);
                for (int i=0; i<nData; ++i) {
                    h->SetBinContent(i+1, getVal(list[i]));
                    h->SetBinError(i+1, getErr(list[i]));
                    h->GetXaxis()->SetBinLabel(i+1, formatLabel(list[i].dirName));
                }
                h->SetMarkerStyle(20); h->SetMarkerSize(1.2); h->SetMarkerColor(color); h->SetLineColor(color);
                
                // Set Y-axis range based on the plot name
                if (name.Contains("Significance")) h->GetYaxis()->SetRangeUser(0, 200);
                else if (name.Contains("Mean_omega")) h->GetYaxis()->SetRangeUser(0.76, 0.79);
                else if (name.Contains("Mean_phi")) h->GetYaxis()->SetRangeUser(1.01, 1.025);
                else if (name.Contains("Sigma")) h->GetYaxis()->SetRangeUser(0.015, 0.035);

                h->Draw("E1P"); h->LabelsOption("v");
                c->Write(); 
                delete c; delete h;
            };

            saveSummary("Significance_omega", "Significance (#omega) "+ptName, [](const SummaryData& d){return d.o_sig;}, [](const SummaryData& d){return d.o_sig_err;}, kBlue);
            saveSummary("Significance_phi", "Significance (#phi) "+ptName, [](const SummaryData& d){return d.p_sig;}, [](const SummaryData& d){return d.p_sig_err;}, kRed);
            saveSummary("Mean_omega", "Mass Position (#omega) "+ptName, [](const SummaryData& d){return d.o_m;}, [](const SummaryData& d){return d.o_m_err;}, kBlue+2);
            saveSummary("Mean_phi", "Mass Position (#phi) "+ptName, [](const SummaryData& d){return d.p_m;}, [](const SummaryData& d){return d.p_m_err;}, kRed+2);
            saveSummary("Sigma_omega", "Mass Width (#omega) "+ptName, [](const SummaryData& d){return d.o_s;}, [](const SummaryData& d){return d.o_s_err;}, kCyan+2);
            saveSummary("Sigma_phi", "Mass Width (#phi) "+ptName, [](const SummaryData& d){return d.p_s;}, [](const SummaryData& d){return d.p_s_err;}, kOrange+2);

            TCanvas *cH = new TCanvas("cHistComp", "Hist Comparison " + ptName, 800, 600); cH->SetGrid();
            TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.88);
            std::vector<int> colors = {kRed, kOrange+1, kYellow+1, kSpring+2, kGreen+2, kCyan+1, kAzure-2, kBlue, kBlack};
            double maxVal = 0;
            for (size_t i=0; i<list.size(); ++i) {
                if (!list[i].hSig) continue;
                if (list[i].hSig->GetMaximum() > maxVal) maxVal = list[i].hSig->GetMaximum();
            }

            for (size_t i=0; i<list.size(); ++i) {
                if (!list[i].hSig) continue;
                int color = colors[i % colors.size()];
                list[i].hSig->SetLineColor(color); list[i].hSig->SetMarkerColor(color);
                list[i].hSig->SetMarkerStyle(20); list[i].hSig->SetMarkerSize(0.7);
                if (i == 0) {
                    list[i].hSig->GetXaxis()->SetRangeUser(0.0, 1.5);
                    list[i].hSig->GetYaxis()->SetRangeUser(0, maxVal*1.2);
                    list[i].hSig->Draw("E1P");
                } else {
                    list[i].hSig->Draw("E1P SAME");
                }
                leg->AddEntry(list[i].hSig, formatLabel(list[i].dirName), "lp");
            }
            leg->Draw(); cH->Write(); 
            delete cH; delete leg;
        }

        // 2. Pt-dependence plots for each ID (Significance vs Pt)
        for (auto& entry : dataByID) {
            TString idName = entry.first;
            auto& list = entry.second;
            // Filter out inclusive pt and sort by ptMin
            std::vector<SummaryData> ptList;
            for (const auto& d : list) {
                if (d.ptName != "Pt_0to30") ptList.push_back(d);
            }
            if (ptList.empty()) continue;

            std::sort(ptList.begin(), ptList.end(), [](const SummaryData& a, const SummaryData& b) {
                return a.ptMin < b.ptMin;
            });

            outfile->cd(idName);
            int nPt = ptList.size();
            std::vector<double> x(nPt), xerr(nPt), y_osig(nPt), yerr_osig(nPt), y_psig(nPt), yerr_psig(nPt);
            std::vector<double> y_om(nPt), yerr_om(nPt), y_pm(nPt), yerr_pm(nPt);
            std::vector<double> y_os(nPt), yerr_os(nPt), y_ps(nPt), yerr_ps(nPt);

            for (int i=0; i<nPt; ++i) {
                x[i] = (ptList[i].ptMin + ptList[i].ptMax) / 2.0;
                xerr[i] = (ptList[i].ptMax - ptList[i].ptMin) / 2.0;
                y_osig[i] = ptList[i].o_sig; yerr_osig[i] = ptList[i].o_sig_err;
                y_psig[i] = ptList[i].p_sig; yerr_psig[i] = ptList[i].p_sig_err;
                y_om[i] = ptList[i].o_m; yerr_om[i] = ptList[i].o_m_err;
                y_pm[i] = ptList[i].p_m; yerr_pm[i] = ptList[i].p_m_err;
                y_os[i] = ptList[i].o_s; yerr_os[i] = ptList[i].o_s_err;
                y_ps[i] = ptList[i].p_s; yerr_ps[i] = ptList[i].p_s_err;
            }

            auto saveGraph = [&](TString name, TString title, std::vector<double>& vy, std::vector<double>& vyer, int color) {
                TGraphErrors *gr = new TGraphErrors(nPt, &x[0], &vy[0], &xerr[0], &vyer[0]);
                gr->SetName("gr_"+name); gr->SetTitle(title + " vs p_{T};p_{T} (GeV/c);Value");
                gr->SetMarkerStyle(20); gr->SetMarkerSize(1.2); gr->SetMarkerColor(color); gr->SetLineColor(color);
                
                gr->Write();
                
                // Create a canvas to save the plot without connecting lines
                TCanvas *cGr = new TCanvas("c_gr_"+name, title, 800, 600);
                cGr->SetGrid();
                gr->Draw("APE"); // A: Axis, P: Points, E: Errors (No 'L' means no connecting lines)

                // Set Y-axis range based on the plot name after Draw
                if (name.Contains("Significance")) gr->GetYaxis()->SetRangeUser(0, 200);
                else if (name.Contains("Mean_omega")) gr->GetYaxis()->SetRangeUser(0.76, 0.79);
                else if (name.Contains("Mean_phi")) gr->GetYaxis()->SetRangeUser(1.01, 1.025);
                else if (name.Contains("Sigma")) gr->GetYaxis()->SetRangeUser(0.015, 0.035);

                cGr->Update();
                cGr->Write();
                
                delete cGr;
                delete gr;
            };

            saveGraph("Significance_omega_pt", "Significance (#omega)", y_osig, yerr_osig, kBlue);
            saveGraph("Significance_phi_pt", "Significance (#phi)", y_psig, yerr_psig, kRed);
            saveGraph("Mean_omega_pt", "Mass Position (#omega)", y_om, yerr_om, kBlue+2);
            saveGraph("Mean_phi_pt", "Mass Position (#phi)", y_pm, yerr_pm, kRed+2);
            saveGraph("Sigma_omega_pt", "Mass Width (#omega)", y_os, yerr_os, kCyan+2);
            saveGraph("Sigma_phi_pt", "Mass Width (#phi)", y_ps, yerr_ps, kOrange+2);
        }

        for (auto& d : summaryDataList) if (d.hSig) delete d.hSig;
    }
    outfile->Close();
    std::cout << "Finish summary plots generation with pT dependence" << std::endl;
}
