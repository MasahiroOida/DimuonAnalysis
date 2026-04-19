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

double CalcSignificanceErrorPol4_minpt(double S, double S_err, double T, double T_err, double CovST) {
    double safeT = std::max(T, EPSILON);
    double termS = pow(S_err, 2) / safeT;
    double termT = pow(S, 2) * pow(T_err, 2) / (4.0 * pow(safeT, 3));
    double termCov = - S * CovST / pow(safeT, 2);
    return sqrt(std::max(0.0, termS + termT + termCov));
}

void Significance_CrystalBall_pol4_minpt(TFile *input_FitResult) {
    if (!input_FitResult || input_FitResult->IsZombie()) return;
    TFile *outfile = new TFile("Significance_Results_pol4_minpt.root", "RECREATE");
    
    std::set<std::string> validIDs;
    TIter scanIter(input_FitResult->GetListOfKeys());
    TKey *scanKey;
    while ((scanKey = (TKey *)scanIter())) {
        if (TString(scanKey->GetClassName()) == "TDirectoryFile") {
            TString dName = scanKey->GetName();
            if (dName.Contains("minpt", TString::kIgnoreCase)) {
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

    struct SummaryData {
        TString dirName;
        TString ptName;
        double minptVal;
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
            if (topName.Contains("minpt", TString::kIgnoreCase)) {
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

                    // pol4 background case: Omega (5-9), Phi (10-14)
                    double o_m = totalfit->GetParameter(6); double o_s = totalfit->GetParameter(7);
                    double p_m = totalfit->GetParameter(11); double p_s = totalfit->GetParameter(12);

                    bool physicallyReasonable = (o_m > 0.6 && o_m < 0.9 && p_m > 0.9 && p_m < 1.2 && o_s > 0 && p_s > 0);
                    if (!physicallyReasonable) continue;

                    TF1 *fBG = new TF1("fBG", "pol4", 0.6, 1.3);
                    for(int i=0; i<5; ++i) fBG->SetParameter(i, totalfit->GetParameter(i));

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
                    double sig_o_err = CalcSignificanceErrorPol4_minpt(S_o, S_o_err, T_o, T_o_err, 0.0);
                    double sig_p_err = CalcSignificanceErrorPol4_minpt(S_p, S_p_err, T_p, T_p_err, 0.0);

                    delete fBG; delete hPureSig; delete hResBG;

                    (new TParameter<double>("Significance_omega", sig_o))->Write();
                    (new TParameter<double>("Significance_omega_err", sig_o_err))->Write();
                    (new TParameter<double>("Significance_phi", sig_p))->Write();
                    (new TParameter<double>("Significance_phi_err", sig_p_err))->Write();
                    (new TParameter<double>("Mean_omega", o_m))->Write();
                    (new TParameter<double>("Mean_omega_err", totalfit->GetParError(6)))->Write();
                    (new TParameter<double>("Mean_phi", p_m))->Write();
                    (new TParameter<double>("Mean_phi_err", totalfit->GetParError(11)))->Write();
                    (new TParameter<double>("Sigma_omega", o_s))->Write();
                    (new TParameter<double>("Sigma_omega_err", totalfit->GetParError(7)))->Write();
                    (new TParameter<double>("Sigma_phi", p_s))->Write();
                    (new TParameter<double>("Sigma_phi_err", totalfit->GetParError(12)))->Write();

                    SummaryData data;
                    data.dirName = topName;
                    data.ptName = ptName;
                    data.ptMin = ptMin;
                    data.ptMax = ptMax;
                    data.minptVal = 999999.0;
                    TObjArray *tokens = topName.Tokenize("_");
                    for (int i = 0; i < tokens->GetEntries(); ++i) {
                        TString tok = tokens->At(i)->GetName();
                        TString tokLow = tok; tokLow.ToLower();
                        if (tokLow == "minpt" && i + 1 < tokens->GetEntries()) {
                            TString valStr = tokens->At(i+1)->GetName();
                            double val = valStr.Atof();
                            if (valStr.Length() == 2 && valStr.IsDigit()) val /= 10.0;
                            data.minptVal = val;
                            break;
                        }
                    }
                    delete tokens;

                    data.o_sig = sig_o; data.o_sig_err = sig_o_err;
                    data.p_sig = sig_p; data.p_sig_err = sig_p_err;
                    data.o_m = o_m; data.o_m_err = totalfit->GetParError(6);
                    data.o_s = o_s; data.o_s_err = totalfit->GetParError(7);
                    data.p_m = p_m; data.p_m_err = totalfit->GetParError(11);
                    data.p_s = p_s; data.p_s_err = totalfit->GetParError(12);

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
            if (!dirNameLow.Contains("minpt")) return "nocut";
            TObjArray *tokens = dirName.Tokenize("_");
            TString label = "nocut";
            for (int i = 0; i < tokens->GetEntries(); ++i) {
                TString tok = tokens->At(i)->GetName();
                TString tokLow = tok; tokLow.ToLower();
                if (tokLow == "minpt" && i + 1 < tokens->GetEntries()) {
                    TString valStr = tokens->At(i+1)->GetName();
                    double val = valStr.Atof();
                    if (valStr.Length() == 2 && valStr.IsDigit()) val /= 10.0;
                    label = Form("p_{T,#mu} < %.1f", val);
                    break;
                }
            }
            delete tokens;
            return label;
        };
        std::map<TString, std::vector<SummaryData>> dataByPt;
        for (const auto& d : summaryDataList) dataByPt[d.ptName].push_back(d);

        for (auto& entry : dataByPt) {
            TString ptName = entry.first;
            auto& list = entry.second;
            std::sort(list.begin(), list.end(), [](const SummaryData& a, const SummaryData& b) {
                return a.minptVal < b.minptVal;
            });
            TDirectory *ptSummaryDir = outfile->mkdir(Form("Summary_%s", ptName.Data()));
            ptSummaryDir->cd();
            int nData = list.size();
            auto saveSummary = [&](TString name, TString title, auto getVal, auto getErr, int color) {
                TCanvas *c = new TCanvas("c_"+name, title, 800, 600);
                c->SetGrid();
                TH1F *h = new TH1F("h_"+name, title+";p_{T,#mu} (GeV/c);Value", nData, 0, nData);
                for (int i=0; i<nData; ++i) {
                    h->SetBinContent(i+1, getVal(list[i]));
                    h->SetBinError(i+1, getErr(list[i]));
                    h->GetXaxis()->SetBinLabel(i+1, formatLabel(list[i].dirName));
                }
                h->SetMarkerStyle(20); h->SetMarkerSize(1.2); h->SetMarkerColor(color); h->SetLineColor(color);
                if (name.Contains("Significance")) h->GetYaxis()->SetRangeUser(0, 200);
                h->Draw("E1P");
                c->Write(); delete c; delete h;
            };
            saveSummary("Significance_omega", "Significance (#omega) "+ptName, [](const SummaryData& d){return d.o_sig;}, [](const SummaryData& d){return d.o_sig_err;}, kBlue);
            saveSummary("Significance_phi", "Significance (#phi) "+ptName, [](const SummaryData& d){return d.p_sig;}, [](const SummaryData& d){return d.p_sig_err;}, kRed);
            saveSummary("Mean_omega", "Mass Position (#omega) "+ptName, [](const SummaryData& d){return d.o_m;}, [](const SummaryData& d){return d.o_m_err;}, kBlue+2);
            saveSummary("Mean_phi", "Mass Position (#phi) "+ptName, [](const SummaryData& d){return d.p_m;}, [](const SummaryData& d){return d.p_m_err;}, kRed+2);
            saveSummary("Sigma_omega", "Mass Width (#omega) "+ptName, [](const SummaryData& d){return d.o_s;}, [](const SummaryData& d){return d.o_s_err;}, kCyan+2);
            saveSummary("Sigma_phi", "Mass Width (#phi) "+ptName, [](const SummaryData& d){return d.p_s;}, [](const SummaryData& d){return d.p_s_err;}, kOrange+2);

            // Overlay Canvas for Mass Histograms
            TCanvas *cOverlay = new TCanvas("c_Overlay_Mass_"+ptName, "Overlay Mass Histograms "+ptName, 1000, 800);
            cOverlay->SetGrid();
            TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
            leg->SetBorderSize(0); leg->SetFillStyle(0);

            double maxVal = 0;
            for (int i=0; i<nData; ++i) {
                if (list[i].hSig->GetMaximum() > maxVal) maxVal = list[i].hSig->GetMaximum();
            }

            int colors[] = {kRed, kOrange+1, kYellow+1, kSpring+2, kGreen+2, kCyan+1, kAzure-2, kBlue, kBlack};
            int nColors = sizeof(colors)/sizeof(int);

            for (int i=0; i<nData; ++i) {
                TH1F *h = (TH1F*)list[i].hSig->Clone();
                h->GetListOfFunctions()->Clear();
                h->SetLineColor(colors[i % nColors]);
                h->SetMarkerColor(colors[i % nColors]);
                h->SetMarkerStyle(20 + (i / nColors));
                h->SetMarkerSize(0.8);
                h->SetLineWidth(2);
                h->GetYaxis()->SetRangeUser(0, maxVal * 1.3);
                h->GetXaxis()->SetRangeUser(0,1.5);
                h->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
                h->GetYaxis()->SetTitle("Counts");

                TString label = formatLabel(list[i].dirName);
                leg->AddEntry(h, label, "lp");

                if (i == 0) h->Draw("E1P");
                else h->Draw("E1P SAME");
            }
            leg->Draw();
            cOverlay->Write();
            delete cOverlay;
        }
    }
    outfile->Close();
    std::cout << "Finish Significance calculation for minpt (pol4) with overlay plots" << std::endl;
}
