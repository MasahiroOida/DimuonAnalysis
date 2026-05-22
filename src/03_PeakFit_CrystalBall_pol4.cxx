// include
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TObject.h>
#include <TStyle.h>
#include <TLegend.h>
#include <cmath>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <TKey.h>
#include <TParameter.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>

Double_t myFreeFunctionPol4(Double_t *x, Double_t *par)
{
    // pol4: par[0] + par[1]*x + par[2]*x^2 + par[3]*x^3 + par[4]*x^4
    double val = par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4);
    if (0.5 <= x[0] && x[0] < 0.70) return val;
    else if (0.70 <= x[0] && x[0] <= 0.86) return 0;
    else if (0.86 < x[0] && x[0] < 0.92) return val;
    else if (0.92 <= x[0] && x[0] <= 1.15) return 0;
    else if (1.15 < x[0] && x[0] <= 1.3) return val;
    return 0;
};

void Peakfit(TH1 *Sig, double nEvents = -1)
{
    if (Sig == nullptr) return;
    TH1::SetDefaultSumw2();
    
    TH1F *BGHistogram = (TH1F *)Sig->Clone("BGHistogram");
    TH1F *SigHistogram = (TH1F *)Sig->Clone("SigHistogram");
    
    double splitLow1 = 0.70; double splitHigh1 = 0.86;
    double splitLow2 = 0.92; double splitHigh2 = 1.15;

    for (int bin = 1; bin <= BGHistogram->GetNbinsX(); ++bin) {
        double binCenter = BGHistogram->GetBinCenter(bin);
        if ((binCenter >= splitLow1 && binCenter <= splitHigh1) || (binCenter >= splitLow2 && binCenter <= splitHigh2))
            BGHistogram->SetBinContent(bin, 0);
    }

    // pol4 has 5 parameters
    TF1 *fitFunction = new TF1("fitFunction", myFreeFunctionPol4, 0.6, 1.4, 5);
    fitFunction->SetParameters(100, -50, 20, -5, 1); // Rough initial values for pol4
    BGHistogram->Fit("fitFunction", "RQ");
    
    TF1 *BGFunction = new TF1("BGFunction", "pol4", 0.6, 1.4);
    for(int i=0; i<5; ++i) BGFunction->SetParameter(i, fitFunction->GetParameter(i));

    for (int bin = 1; bin <= Sig->GetNbinsX(); ++bin) {
        SigHistogram->SetBinContent(bin, Sig->GetBinContent(bin) - BGFunction->Eval(Sig->GetBinCenter(bin)));
    }

    TF1 *omegafitfunc = new TF1("omegafitfunc", "crystalball", splitLow1, splitHigh1);
    omegafitfunc->SetParameters(SigHistogram->GetMaximum() / 2, 0.78, 0.02, 1.5, 2.0);
    omegafitfunc->SetParLimits(1, 0.75, 0.85);
    omegafitfunc->SetParLimits(2, 0.01, 1.0);
    omegafitfunc->SetParLimits(3, 0.5, 3.0);
    omegafitfunc->SetParLimits(4, 1, 10.0);
    SigHistogram->Fit("omegafitfunc", "RQ");

    TF1 *phifitfunc = new TF1("phifitfunc", "crystalball", splitLow2, splitHigh2);
    phifitfunc->SetParameters(SigHistogram->GetMaximum() / 2, 1.02, 0.02, 1.5, 2.0);
    phifitfunc->SetParLimits(1, 0.9, 1.1);
    phifitfunc->SetParLimits(2, 0.01, 1.0);
    phifitfunc->SetParLimits(3, 0.5, 3.0);
    phifitfunc->SetParLimits(4, 1, 10.0);
    SigHistogram->Fit("phifitfunc", "RQ");

    // Total Fit Mapping:
    // BG (pol4): [0, 1, 2, 3, 4]
    // Omega (CB): [5, 6, 7, 8, 9]
    // Phi (CB): [10, 11, 12, 13, 14]
    TF1 *totalfit = new TF1("totalfit", "pol4(0) + crystalball(5) + crystalball(10)", 0.6, 1.2);

    for(int i=0; i<5; ++i) totalfit->SetParameter(i, fitFunction->GetParameter(i));
    for(int i=0; i<5; ++i) {
        totalfit->SetParameter(5+i, omegafitfunc->GetParameter(i));
        totalfit->SetParameter(10+i, phifitfunc->GetParameter(i));
    };

    auto clamp = [](double val, double low, double high) {
        if (val <= low) return low + 0.0001;
        if (val >= high) return high - 0.0001;
        return val;
    };

    // Apply limits for totalfit AFTER setting initial values and clamping
    totalfit->SetParameter(6, clamp(totalfit->GetParameter(6), 0.75, 0.85));
    totalfit->SetParLimits(6, 0.75, 0.85); // Omega Mean
    totalfit->SetParameter(11, clamp(totalfit->GetParameter(11), 0.9, 1.1));
    totalfit->SetParLimits(11, 0.9, 1.1);   // Phi Mean

    totalfit->SetParameter(7, clamp(totalfit->GetParameter(7), 0.01, 1.0));
    totalfit->SetParLimits(7, 0.01, 1.0);  // Omega Sigma
    totalfit->SetParameter(12, clamp(totalfit->GetParameter(12), 0.01, 1.0));
    totalfit->SetParLimits(12, 0.01, 1.0); // Phi Sigma

    totalfit->SetParameter(8, clamp(totalfit->GetParameter(8), 0.5, 3.0));
    totalfit->SetParLimits(8, 0.5, 3.0);   // Omega CB alpha
    totalfit->SetParameter(9, clamp(totalfit->GetParameter(9), 1.0, 10.0));
    totalfit->SetParLimits(9, 1, 10.0);    // Omega CB n
    totalfit->SetParameter(13, clamp(totalfit->GetParameter(13), 0.5, 3.0));
    totalfit->SetParLimits(13, 0.5, 3.0);  // Phi CB alpha
    totalfit->SetParameter(14, clamp(totalfit->GetParameter(14), 1.0, 10.0));
    totalfit->SetParLimits(14, 1, 10.0);   // Phi CB n

    TFitResultPtr fitRes = Sig->Fit("totalfit", "RQS");
    if (fitRes.Get()) fitRes->Write("fitResult");
    
    // Save Peak-only histogram by subtracting final BG (raw counts)
    TH1F *SigOnly = (TH1F *)Sig->Clone("LikeSignSig_PeakOnly");
    SigOnly->SetTitle("Signal after subtracting final BG");
    SigOnly->GetListOfFunctions()->Clear(); // Remove any associated fit functions
    TF1 *fBG_final = new TF1("fBG_final", "pol4", 0.6, 1.4);
    for(int i=0; i<5; ++i) fBG_final->SetParameter(i, totalfit->GetParameter(i));
    for (int bin = 1; bin <= SigOnly->GetNbinsX(); ++bin) {
        double binCenter = SigOnly->GetBinCenter(bin);
        SigOnly->SetBinContent(bin, SigOnly->GetBinContent(bin) - fBG_final->Eval(binCenter));
    }
    
    // Keep raw counts (do not normalize by event count and bin width)
    Sig->SetYTitle("Counts / (Mass bin)");
    SigOnly->SetYTitle("Counts / (Mass bin)");
    
    Sig->Write("LikeSignSig_fit");
    totalfit->Write("totalfit");
    SigOnly->Write("LikeSignSig_PeakOnly");

    // Create a canvas to overlay Sig and PeakOnly (both already normalized)
    TCanvas *cOverlay = new TCanvas("Signal_Overlay_Canvas", "Overlay of LS Signal and PeakOnly", 800, 600);
    cOverlay->cd();
    
    Sig->SetLineColor(kBlack);
    Sig->SetMarkerColor(kBlack);
    Sig->SetMarkerStyle(20);
    Sig->SetMarkerSize(0.7);
    Sig->SetTitle("Signal Comparison;Mass (GeV/c^{2});dN/dm");
    Sig->GetXaxis()->SetRangeUser(0.6, 1.4);
    Sig->Draw("E");

    SigOnly->SetLineColor(kRed);
    SigOnly->SetMarkerColor(kRed);
    SigOnly->SetMarkerStyle(20);
    SigOnly->SetMarkerSize(0.7);
    SigOnly->Draw("E SAME");

    // --- Draw Fit Components ---
    TF1 *fBG_draw = new TF1("fBG_draw", "pol4", 0.6, 1.2);
    for(int i=0; i<5; ++i) fBG_draw->SetParameter(i, totalfit->GetParameter(i));
    fBG_draw->SetLineColor(kGreen+2); fBG_draw->SetLineStyle(2); fBG_draw->SetLineWidth(2);

    TF1 *fOmega_draw = new TF1("fOmega_draw", "crystalball", 0.6, 1.2);
    for(int i=0; i<5; ++i) fOmega_draw->SetParameter(i, totalfit->GetParameter(5+i));
    fOmega_draw->SetLineColor(kBlue); fOmega_draw->SetLineStyle(2); fOmega_draw->SetLineWidth(2);

    TF1 *fPhi_draw = new TF1("fPhi_draw", "crystalball", 0.6, 1.2);
    for(int i=0; i<5; ++i) fPhi_draw->SetParameter(i, totalfit->GetParameter(10+i));
    fPhi_draw->SetLineColor(kRed); fPhi_draw->SetLineStyle(2); fPhi_draw->SetLineWidth(2);

    fBG_draw->Draw("SAME"); fOmega_draw->Draw("SAME"); fPhi_draw->Draw("SAME");
    totalfit->SetLineColor(kMagenta); totalfit->SetLineWidth(3);
    totalfit->Draw("SAME");

    TLegend *legOverlay = new TLegend(0.6, 0.6, 0.88, 0.88);
    legOverlay->AddEntry(Sig, "Data (incl. BG)", "lep");
    legOverlay->AddEntry(SigOnly, "Data (BG sub)", "lep");
    legOverlay->AddEntry(totalfit, "Total Fit", "l");
    legOverlay->AddEntry(fOmega_draw, "#omega Signal (CB)", "l");
    legOverlay->AddEntry(fPhi_draw, "#phi Signal (CB)", "l");
    legOverlay->AddEntry(fBG_draw, "Background (Pol4)", "l");
    legOverlay->Draw();

    cOverlay->Write();
    delete cOverlay; delete legOverlay; delete fBG_draw; delete fOmega_draw; delete fPhi_draw;

    delete BGHistogram; delete SigHistogram; delete fitFunction; delete BGFunction;
    delete omegafitfunc; delete phifitfunc; delete totalfit; delete SigOnly; delete fBG_final;
}

void PeakFit_CrystalBall_pol4(TFile *input_Sig){
    TFile *outfile = new TFile("PeakFit_CrystalBall_pol4_Results.root", "RECREATE");
    TIter next(input_Sig->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (TString(key->GetClassName()) == "TDirectoryFile") {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TString topDirName = topDir->GetName();
            TDirectory *outTopDir = outfile->mkdir(topDirName);

            // フィット結果格納用ベクトル
            std::vector<double> v_pt, v_pt_err;
            std::vector<double> v_om_mean, v_om_mean_err, v_om_sigma, v_om_sigma_err, v_om_alpha, v_om_alpha_err, v_om_n, v_om_n_err;
            std::vector<double> v_ph_mean, v_ph_mean_err, v_ph_sigma, v_ph_sigma_err, v_ph_alpha, v_ph_alpha_err, v_ph_n, v_ph_n_err;

            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt())) {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile") {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    TString ptDirName = ptDir->GetName();
                    TDirectory *outPtDir = outTopDir->mkdir(ptDirName);
                    outPtDir->cd();
                    TH1F *LikeSignSig = (TH1F *)ptDir->Get("LikeSignSig");
                    TParameter<double> *pN = (TParameter<double> *)ptDir->Get("nEvents");
                    double nEvents = (pN) ? pN->GetVal() : -1;
                    if(LikeSignSig) Peakfit(LikeSignSig, nEvents);
                    
                    // フィット結果の抽出 (Peakfit内でoutPtDirにtotalfitが保存されていることを前提)
                    TF1 *totalfit = (TF1*)outPtDir->Get("totalfit");
                    TParameter<double> *pMin = (TParameter<double>*)ptDir->Get("ptmin");
                    TParameter<double> *pMax = (TParameter<double>*)ptDir->Get("ptmax");
                    if (pN) pN->Write("nEvents");

                    if (totalfit && pMin && pMax) {
                        // Extract fit results for all bins
                        // Omega (Indices 5-9 in totalfit: pol4(0-4) + CB_om(5-9) + CB_ph(10-14))
                        
                        // Only add to summary vectors if it's not a total pt bin (e.g., range < 9.0)
                        if (pMax->GetVal() - pMin->GetVal() < 9.0) {
                            double pt_center = (pMin->GetVal() + pMax->GetVal()) / 2.0;
                            double pt_halfwidth = (pMax->GetVal() - pMin->GetVal()) / 2.0;
                            v_pt.push_back(pt_center);
                            v_pt_err.push_back(pt_halfwidth);

                            v_om_mean.push_back(totalfit->GetParameter(6));
                            v_om_mean_err.push_back(totalfit->GetParError(6));
                            v_om_sigma.push_back(totalfit->GetParameter(7));
                            v_om_sigma_err.push_back(totalfit->GetParError(7));
                            v_om_alpha.push_back(totalfit->GetParameter(8));
                            v_om_alpha_err.push_back(totalfit->GetParError(8));
                            v_om_n.push_back(totalfit->GetParameter(9));
                            v_om_n_err.push_back(totalfit->GetParError(9));

                            // Phi (Indices 11-15 in totalfit)
                            v_ph_mean.push_back(totalfit->GetParameter(11));
                            v_ph_mean_err.push_back(totalfit->GetParError(11));
                            v_ph_sigma.push_back(totalfit->GetParameter(12));
                            v_ph_sigma_err.push_back(totalfit->GetParError(12));
                            v_ph_alpha.push_back(totalfit->GetParameter(13));
                            v_ph_alpha_err.push_back(totalfit->GetParError(13));
                            v_ph_n.push_back(totalfit->GetParameter(14));
                            v_ph_n_err.push_back(totalfit->GetParError(14));
                        }
                    }

                    TObject *p;
                    if((p = ptDir->Get("ptrange"))) p->Write("ptrange");
                    if((p = ptDir->Get("ptmin"))) p->Write("ptmin");
                    if((p = ptDir->Get("ptmax"))) p->Write("ptmax");
                }
            }

            // サマリー図の作成
            if (!v_pt.empty()) {
                outTopDir->cd();
                auto saveSummary = [&](std::vector<double>& y, std::vector<double>& yerr, const char* part, const char* type, const char* title, const char* ytitle) {
                    TString baseName = Form("%s_%s_%s", part, type, topDirName.Data());
                    TString cName = "Summary_" + baseName;
                    TString gName = "Graph_" + baseName;

                    TCanvas *c = new TCanvas(cName, cName, 800, 600);
                    TGraphErrors *gr = new TGraphErrors(v_pt.size(), &v_pt[0], &y[0], &v_pt_err[0], &yerr[0]);
                    gr->SetName(gName);
                    gr->SetTitle(Form("%s (%s);p_{T} (GeV/c);%s", title, part, ytitle));
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerColor(kRed);
                    gr->SetLineColor(kRed);
                    gr->Draw("APE");
                    c->Write();
                    gr->Write();
                    delete c;
                    delete gr;
                };

                saveSummary(v_om_mean, v_om_mean_err, "Omega", "Mean", "Mean Mass", "Mean (GeV/c^{2})");
                saveSummary(v_om_sigma, v_om_sigma_err, "Omega", "Sigma", "Mass Width (Sigma)", "#sigma (GeV/c^{2})");
                saveSummary(v_om_alpha, v_om_alpha_err, "Omega", "Alpha", "CB Alpha", "#alpha");
                saveSummary(v_om_n, v_om_n_err, "Omega", "n", "CB n", "n");

                saveSummary(v_ph_mean, v_ph_mean_err, "Phi", "Mean", "Mean Mass", "Mean (GeV/c^{2})");
                saveSummary(v_ph_sigma, v_ph_sigma_err, "Phi", "Sigma", "Mass Width (Sigma)", "#sigma (GeV/c^{2})");
                saveSummary(v_ph_alpha, v_ph_alpha_err, "Phi", "Alpha", "CB Alpha", "#alpha");
                saveSummary(v_ph_n, v_ph_n_err, "Phi", "n", "CB n", "n");
            }
        }
    }
    outfile->Close();
    delete outfile;
    std::cout << "Finish Crystal Ball fitting with pol4 background and saved summaries" << std::endl;
}
