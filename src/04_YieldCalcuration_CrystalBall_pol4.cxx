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
#include <TVectorF.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>

struct Yield {
    double Omega_yield;
    double Omega_yield_err;
    double Phi_yield;
    double Phi_yield_err;
};

Yield yieldCal_CrystalBall_pol4(TH1 *Sig, TF1 *Totalfun, TFitResult *fitRes, double ptrange) {
    if (!Sig || !Totalfun || !fitRes) return {0,0,0,0};
    double binWidth = Sig->GetBinWidth(1);
    
    // Omega (Parameters 5-9 in Totalfun for pol4 background case)
    TF1 *omega_cb = new TF1("omega_cb", "crystalball", 0.6, 1.3);
    for(int i=0; i<5; ++i) omega_cb->SetParameter(i, Totalfun->GetParameter(5+i));
    
    // Phi (Parameters 10-14 in Totalfun for pol4 background case)
    TF1 *phi_cb = new TF1("phi_cb", "crystalball", 0.6, 1.3);
    for(int i=0; i<5; ++i) phi_cb->SetParameter(i, Totalfun->GetParameter(10+i));

    double Omega_mass = Totalfun->GetParameter(6);
    double Omega_width = Totalfun->GetParameter(7);
    double Phi_mass = Totalfun->GetParameter(11);
    double Phi_width = Totalfun->GetParameter(12);

    double Omega_Down = Omega_mass - 3 * Omega_width;
    double Omega_UP = Omega_mass + 3 * Omega_width;
    double Phi_Down = Phi_mass - 3 * Phi_width;
    double Phi_UP = Phi_mass + 3 * Phi_width;

    double Y_omega = omega_cb->Integral(Omega_Down, Omega_UP) / binWidth;
    double Y_phi = phi_cb->Integral(Phi_Down, Phi_UP) / binWidth;
    
    // Correct Error Propagation using Covariance Matrix
    TMatrixDSym cov = fitRes->GetCovarianceMatrix();
    TMatrixDSym subCovOmega(5);
    TMatrixDSym subCovPhi(5);
    for(int i=0; i<5; ++i) {
        for(int j=0; j<5; ++j) {
            subCovOmega(i, j) = cov(i+5, j+5);
            subCovPhi(i, j) = cov(i+10, j+10);
        }
    }
    
    double Y_omega_err = omega_cb->IntegralError(Omega_Down, Omega_UP, omega_cb->GetParameters(), subCovOmega.GetMatrixArray()) / binWidth;
    double Y_phi_err = phi_cb->IntegralError(Phi_Down, Phi_UP, phi_cb->GetParameters(), subCovPhi.GetMatrixArray()) / binWidth;

    delete omega_cb; delete phi_cb;
    return {Y_omega / ptrange, Y_omega_err / ptrange, Y_phi / ptrange, Y_phi_err / ptrange};
}

void YieldCalcuration_CrystalBall_pol4(TFile *input_FitResult) {
    TFile *outfile = new TFile("Yield_Results_CrystalBall_pol4.root", "RECREATE");
    TIter next(input_FitResult->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (TString(key->GetClassName()) == "TDirectoryFile") {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TDirectory *outTopDir = outfile->mkdir(topDir->GetName());
            std::vector<double> pt_vec, pt_err_vec, omega_yield_vec, omega_yield_err_vec, phi_yield_vec, phi_yield_err_vec;
            std::vector<double> omega_mean_vec, omega_mean_err_vec, omega_sigma_vec, omega_sigma_err_vec;
            std::vector<double> omega_alpha_vec, omega_alpha_err_vec, omega_n_vec, omega_n_err_vec;
            std::vector<double> phi_mean_vec, phi_mean_err_vec, phi_sigma_vec, phi_sigma_err_vec;
            std::vector<double> phi_alpha_vec, phi_alpha_err_vec, phi_n_vec, phi_n_err_vec;
            
            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt())) {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile") {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    if (TString(ptDir->GetName()) == "Pt_0to30" || TString(ptDir->GetName()) == "Pt_0to10") continue;
                    TDirectory *outPtDir = outTopDir->mkdir(ptDir->GetName());

                    TH1F *LikeSignSig = (TH1F *)ptDir->Get("LikeSignSig_fit");
                    TF1 *totalfit = (TF1 *)ptDir->Get("totalfit");
                    TFitResult *fitRes = (TFitResult*)ptDir->Get("fitResult");
                    TParameter<double>* p_range = (TParameter<double>*)ptDir->Get("ptrange");
                    TParameter<double>* p_min = (TParameter<double>*)ptDir->Get("ptmin");
                    TParameter<double>* p_max = (TParameter<double>*)ptDir->Get("ptmax");

                    if(LikeSignSig && totalfit && fitRes && p_range) {
                        Yield y = yieldCal_CrystalBall_pol4(LikeSignSig, totalfit, fitRes, p_range->GetVal());
                        pt_vec.push_back(0.5 * (p_max->GetVal() + p_min->GetVal()));
                        pt_err_vec.push_back(0.5 * p_range->GetVal());
                        omega_yield_vec.push_back(y.Omega_yield);
                        omega_yield_err_vec.push_back(y.Omega_yield_err);
                        phi_yield_vec.push_back(y.Phi_yield);
                        phi_yield_err_vec.push_back(y.Phi_yield_err);

                        // Omega parameters (Pars 5-9)
                        omega_mean_vec.push_back(totalfit->GetParameter(6));
                        omega_mean_err_vec.push_back(totalfit->GetParError(6));
                        omega_sigma_vec.push_back(totalfit->GetParameter(7));
                        omega_sigma_err_vec.push_back(totalfit->GetParError(7));
                        omega_alpha_vec.push_back(totalfit->GetParameter(8));
                        omega_alpha_err_vec.push_back(totalfit->GetParError(8));
                        omega_n_vec.push_back(totalfit->GetParameter(9));
                        omega_n_err_vec.push_back(totalfit->GetParError(9));

                        // Phi parameters (Pars 10-14)
                        phi_mean_vec.push_back(totalfit->GetParameter(11));
                        phi_mean_err_vec.push_back(totalfit->GetParError(11));
                        phi_sigma_vec.push_back(totalfit->GetParameter(12));
                        phi_sigma_err_vec.push_back(totalfit->GetParError(12));
                        phi_alpha_vec.push_back(totalfit->GetParameter(13));
                        phi_alpha_err_vec.push_back(totalfit->GetParError(13));
                        phi_n_vec.push_back(totalfit->GetParameter(14));
                        phi_n_err_vec.push_back(totalfit->GetParError(14));

                        // --- Add Canvas for Fit Visualization ---
                        outPtDir->cd();
                        TCanvas *cFit = new TCanvas("cFitResult", Form("Fit Result %s", ptDir->GetName()), 800, 600);
                        gStyle->SetOptFit(1111);
                        LikeSignSig->SetMarkerStyle(20); LikeSignSig->SetMarkerSize(0.8);
                        LikeSignSig->Draw("E1P");
                        LikeSignSig->GetXaxis()->SetRangeUser(0,1.5);

                        TF1 *fBG = new TF1("fBG_draw", "pol4", 0.6, 1.2);
                        for(int i=0; i<5; ++i) fBG->SetParameter(i, totalfit->GetParameter(i));
                        fBG->SetLineColor(kGreen+2); fBG->SetLineStyle(2);

                        TF1 *fOmega = new TF1("fOmega_draw", "crystalball", 0.6, 1.2);
                        for(int i=0; i<5; ++i) fOmega->SetParameter(i, totalfit->GetParameter(5+i));
                        fOmega->SetLineColor(kBlue); fOmega->SetLineStyle(2);

                        TF1 *fPhi = new TF1("fPhi_draw", "crystalball", 0.6, 1.2);
                        for(int i=0; i<5; ++i) fPhi->SetParameter(i, totalfit->GetParameter(10+i));
                        fPhi->SetLineColor(kRed); fPhi->SetLineStyle(2);

                        fBG->Draw("SAME"); fOmega->Draw("SAME"); fPhi->Draw("SAME");
                        totalfit->SetLineColor(kMagenta); totalfit->SetLineWidth(3);
                        totalfit->Draw("SAME");

                        TLegend *leg = new TLegend(0.6, 0.4, 0.88, 0.65);
                        leg->AddEntry(LikeSignSig, "Data (LS Subtracted)", "lp");
                        leg->AddEntry(totalfit, "Total Fit", "l");
                        leg->AddEntry(fOmega, "#omega Signal (CB)", "l");
                        leg->AddEntry(fPhi, "#phi Signal (CB)", "l");
                        leg->AddEntry(fBG, "Background (Pol4)", "l");
                        leg->Draw();

                        cFit->Write();
                        delete cFit; delete fBG; delete fOmega; delete fPhi;
                    }
                }
            }
            outTopDir->cd();
            if (!pt_vec.empty()) {
                TGraphErrors *g_omega = new TGraphErrors(pt_vec.size(), &pt_vec[0], &omega_yield_vec[0], &pt_err_vec[0], &omega_yield_err_vec[0]);
                g_omega->SetMarkerColor(kBlack);
                g_omega->SetLineColor(kBlack);
                g_omega->Write("omega_yield_ptspectrum");
                TGraphErrors *g_phi = new TGraphErrors(pt_vec.size(), &pt_vec[0], &phi_yield_vec[0], &pt_err_vec[0], &phi_yield_err_vec[0]);
                g_phi->SetMarkerColor(kBlack);
                g_phi->SetLineColor(kBlack);
                g_phi->Write("phi_yield_ptspectrum");

                auto saveParamGraph = [&](std::vector<double>& y, std::vector<double>& yerr, const char* name, const char* title, const char* ytitle, int color) {
                    TGraphErrors *gr = new TGraphErrors(pt_vec.size(), &pt_vec[0], &y[0], &pt_err_vec[0], &yerr[0]);
                    gr->SetName(name);
                    gr->SetTitle(Form("%s;p_{T} (GeV/c);%s", title, ytitle));
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerColor(color);
                    gr->SetLineColor(color);
                    gr->Write();
                    delete gr;
                };

                saveParamGraph(omega_mean_vec, omega_mean_err_vec, "omega_mean_ptspectrum", "Omega Mean Mass", "Mean (GeV/c^{2})", kBlue);
                saveParamGraph(omega_sigma_vec, omega_sigma_err_vec, "omega_sigma_ptspectrum", "Omega Sigma", "#sigma (GeV/c^{2})", kBlue);
                saveParamGraph(omega_alpha_vec, omega_alpha_err_vec, "omega_alpha_ptspectrum", "Omega Alpha", "#alpha", kBlue);
                saveParamGraph(omega_n_vec, omega_n_err_vec, "omega_n_ptspectrum", "Omega n", "n", kBlue);

                saveParamGraph(phi_mean_vec, phi_mean_err_vec, "phi_mean_ptspectrum", "Phi Mean Mass", "Mean (GeV/c^{2})", kRed);
                saveParamGraph(phi_sigma_vec, phi_sigma_err_vec, "phi_sigma_ptspectrum", "Phi Sigma", "#sigma (GeV/c^{2})", kRed);
                saveParamGraph(phi_alpha_vec, phi_alpha_err_vec, "phi_alpha_ptspectrum", "Phi Alpha", "#alpha", kRed);
                saveParamGraph(phi_n_vec, phi_n_err_vec, "phi_n_ptspectrum", "Phi n", "n", kRed);

                // --- Add Canvas for Omega Spectrum ---
                TCanvas *c_omega = new TCanvas("omega_yield_canvas", "Omega Yield Spectrum", 800, 600);
                c_omega->SetLogy();
                c_omega->SetGrid();
                g_omega->SetMarkerStyle(20);
                g_omega->SetMarkerSize(1.0);
                g_omega->SetMarkerColor(kBlack);
                g_omega->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                g_omega->GetYaxis()->SetTitle("raw yield");
                g_omega->Draw("AP"); // AP: Axis and Points (no lines)
                c_omega->Write();
                delete c_omega;

                // --- Add Canvas for Phi Spectrum ---
                TCanvas *c_phi = new TCanvas("phi_yield_canvas", "Phi Yield Spectrum", 800, 600);
                c_phi->SetLogy();
                c_phi->SetGrid();
                g_phi->SetMarkerStyle(20);
                g_phi->SetMarkerSize(1.0);
                g_phi->SetMarkerColor(kBlack);
                g_phi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                g_phi->GetYaxis()->SetTitle("raw yield");
                g_phi->Draw("AP"); // AP: Axis and Points (no lines)
                c_phi->Write();
                delete c_phi;
            }
        }
    }
    outfile->Close();
    delete outfile;
    std::cout << "Finish Yield calculation (pol4) with correct error propagation" << std::endl;
}
