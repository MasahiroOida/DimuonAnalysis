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

Yield yieldCal_CrystalBall(TH1 *Sig, TF1 *Totalfun, TFitResult *fitRes, double ptrange) {
    if (!Sig || !Totalfun || !fitRes) return {0,0,0,0};
    double binWidth = Sig->GetBinWidth(1);
    
    // Omega (Parameters 2-6 in Totalfun)
    TF1 *omega_cb = new TF1("omega_cb", "crystalball", 0.6, 1.3);
    for(int i=0; i<5; ++i) omega_cb->SetParameter(i, Totalfun->GetParameter(2+i));
    
    // Phi (Parameters 7-11 in Totalfun)
    TF1 *phi_cb = new TF1("phi_cb", "crystalball", 0.6, 1.3);
    for(int i=0; i<5; ++i) phi_cb->SetParameter(i, Totalfun->GetParameter(7+i));

    double Omega_mass = Totalfun->GetParameter(3);
    double Omega_width = Totalfun->GetParameter(4);
    double Phi_mass = Totalfun->GetParameter(8);
    double Phi_width = Totalfun->GetParameter(9);

    double Omega_Down = Omega_mass - 3 * Omega_width;
    double Omega_UP = Omega_mass + 3 * Omega_width;
    double Phi_Down = Phi_mass - 3 * Phi_width;
    double Phi_UP = Phi_mass + 3 * Phi_width;

    double Y_omega = omega_cb->Integral(Omega_Down, Omega_UP) / binWidth;
    double Y_phi = phi_cb->Integral(Phi_Down, Phi_UP) / binWidth;
    
    // Correct Error Propagation using Covariance Matrix
    // Extract sub-matrices for Omega (2-6) and Phi (7-11)
    TMatrixDSym cov = fitRes->GetCovarianceMatrix();
    TMatrixDSym subCovOmega(5);
    TMatrixDSym subCovPhi(5);
    for(int i=0; i<5; ++i) {
        for(int j=0; j<5; ++j) {
            subCovOmega(i, j) = cov(i+2, j+2);
            subCovPhi(i, j) = cov(i+7, j+7);
        }
    }
    
    double Y_omega_err = omega_cb->IntegralError(Omega_Down, Omega_UP, omega_cb->GetParameters(), subCovOmega.GetMatrixArray()) / binWidth;
    double Y_phi_err = phi_cb->IntegralError(Phi_Down, Phi_UP, phi_cb->GetParameters(), subCovPhi.GetMatrixArray()) / binWidth;

    delete omega_cb; delete phi_cb;
    return {Y_omega / ptrange, Y_omega_err / ptrange, Y_phi / ptrange, Y_phi_err / ptrange};
}

void YieldCalcuration_CrystalBall(TFile *input_FitResult) {
    TFile *outfile = new TFile("Yield_Results_CrystalBall.root", "RECREATE");
    TIter next(input_FitResult->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (TString(key->GetClassName()) == "TDirectoryFile") {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TDirectory *outTopDir = outfile->mkdir(topDir->GetName());
            std::vector<double> pt_vec, pt_err_vec, omega_yield_vec, omega_yield_err_vec, phi_yield_vec, phi_yield_err_vec;
            
            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt())) {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile") {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    TDirectory *outPtDir = outTopDir->mkdir(ptDir->GetName());

                    TH1F *LikeSignSig = (TH1F *)ptDir->Get("LikeSignSig_fit");
                    TF1 *totalfit = (TF1 *)ptDir->Get("totalfit");
                    TFitResult *fitRes = (TFitResult*)ptDir->Get("fitResult");
                    TParameter<double>* p_range = (TParameter<double>*)ptDir->Get("ptrange");
                    TParameter<double>* p_min = (TParameter<double>*)ptDir->Get("ptmin");
                    TParameter<double>* p_max = (TParameter<double>*)ptDir->Get("ptmax");

                    if(LikeSignSig && totalfit && fitRes && p_range && p_min && p_max) {
                        Yield y = yieldCal_CrystalBall(LikeSignSig, totalfit, fitRes, p_range->GetVal());
                        
                        // Only add to spectrum vectors if it's not a total pt bin (e.g., range < 9.0)
                        if (p_range->GetVal() < 9.0) {
                            pt_vec.push_back(0.5 * (p_max->GetVal() + p_min->GetVal()));
                            pt_err_vec.push_back(0.5 * p_range->GetVal());
                            omega_yield_vec.push_back(y.Omega_yield);
                            omega_yield_err_vec.push_back(y.Omega_yield_err);
                            phi_yield_vec.push_back(y.Phi_yield);
                            phi_yield_err_vec.push_back(y.Phi_yield_err);
                        }

                        // --- Add Canvas for Fit Visualization ---
                        outPtDir->cd();
                        // Save calculated Yields
                        (new TParameter<double>("omega_yield", y.Omega_yield))->Write();
                        (new TParameter<double>("omega_yield_err", y.Omega_yield_err))->Write();
                        (new TParameter<double>("phi_yield", y.Phi_yield))->Write();
                        (new TParameter<double>("phi_yield_err", y.Phi_yield_err))->Write();
                        p_range->Write("ptrange");
                        p_min->Write("ptmin");
                        p_max->Write("ptmax");

                        TCanvas *cFit = new TCanvas("cFitResult", Form("Fit Result %s", ptDir->GetName()), 800, 600);
                        gStyle->SetOptFit(1111);
                        LikeSignSig->SetMarkerStyle(20); LikeSignSig->SetMarkerSize(0.8);
                        LikeSignSig->Draw("E1P");

                        TF1 *fBG = new TF1("fBG_draw", "[0]*exp(-[1]*x)", 0.6, 1.2);
                        fBG->SetParameters(totalfit->GetParameter(0), totalfit->GetParameter(1));
                        fBG->SetLineColor(kGreen+2); fBG->SetLineStyle(2);

                        TF1 *fOmega = new TF1("fOmega_draw", "crystalball", 0.6, 1.2);
                        for(int i=0; i<5; ++i) fOmega->SetParameter(i, totalfit->GetParameter(2+i));
                        fOmega->SetLineColor(kBlue); fOmega->SetLineStyle(2);

                        TF1 *fPhi = new TF1("fPhi_draw", "crystalball", 0.6, 1.2);
                        for(int i=0; i<5; ++i) fPhi->SetParameter(i, totalfit->GetParameter(7+i));
                        fPhi->SetLineColor(kRed); fPhi->SetLineStyle(2);

                        fBG->Draw("SAME"); fOmega->Draw("SAME"); fPhi->Draw("SAME");
                        totalfit->SetLineColor(kMagenta); totalfit->SetLineWidth(3);
                        totalfit->Draw("SAME");

                        TLegend *leg = new TLegend(0.6, 0.4, 0.88, 0.65);
                        leg->AddEntry(LikeSignSig, "Data (LS Subtracted)", "lp");
                        leg->AddEntry(totalfit, "Total Fit", "l");
                        leg->AddEntry(fOmega, "#omega Signal (CB)", "l");
                        leg->AddEntry(fPhi, "#phi Signal (CB)", "l");
                        leg->AddEntry(fBG, "Background (Exp)", "l");
                        leg->Draw();

                        cFit->Write();
                        delete cFit; delete fBG; delete fOmega; delete fPhi;
                    }
                }
            }
            outTopDir->cd();
            if (!pt_vec.empty()) {
                (new TGraphErrors(pt_vec.size(), &pt_vec[0], &omega_yield_vec[0], &pt_err_vec[0], &omega_yield_err_vec[0]))->Write("omega_yield_ptspectrum");
                (new TGraphErrors(pt_vec.size(), &pt_vec[0], &phi_yield_vec[0], &pt_err_vec[0], &phi_yield_err_vec[0]))->Write("phi_yield_ptspectrum");
            }
        }
    }
    outfile->Close();
    delete outfile;
    std::cout << "Finish Yield calculation with correct error propagation" << std::endl;
}
