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

Yield yieldCal_CrystalBall(TH1 *Sig, TF1 *Totalfun) {
    if (!Sig || !Totalfun) return {0,0,0,0};
    
    double Omega_mass = Totalfun->GetParameter(3);
    double Omega_width = Totalfun->GetParameter(4);
    double Phi_mass = Totalfun->GetParameter(8);
    double Phi_width = Totalfun->GetParameter(9);

    double Omega_Down = Omega_mass - 3 * Omega_width;
    double Omega_UP = Omega_mass + 3 * Omega_width;
    double Phi_Down = Phi_mass - 3 * Phi_width;
    double Phi_UP = Phi_mass + 3 * Phi_width;

    TF1 *fBG = new TF1("fBG", "[0]*exp(-[1]*x)", 0.6, 1.4);
    fBG->SetParameters(Totalfun->GetParameter(0), Totalfun->GetParameter(1));

    double Omega_yield = 0;
    double Omega_err2 = 0;
    for (int i = Sig->FindBin(Omega_Down); i <= Sig->FindBin(Omega_UP); ++i) {
        double bin_center = Sig->GetBinCenter(i);
        double bin_content = Sig->GetBinContent(i);
        double bin_err = Sig->GetBinError(i);
        double bg_val = fBG->Eval(bin_center);
        
        Omega_yield += (bin_content - bg_val);
        Omega_err2 += bin_err * bin_err;
    }
    double Omega_yield_err = TMath::Sqrt(Omega_err2);

    double Phi_yield = 0;
    double Phi_err2 = 0;
    for (int i = Sig->FindBin(Phi_Down); i <= Sig->FindBin(Phi_UP); ++i) {
        double bin_center = Sig->GetBinCenter(i);
        double bin_content = Sig->GetBinContent(i);
        double bin_err = Sig->GetBinError(i);
        double bg_val = fBG->Eval(bin_center);
        
        Phi_yield += (bin_content - bg_val);
        Phi_err2 += bin_err * bin_err;
    }
    double Phi_yield_err = TMath::Sqrt(Phi_err2);

    delete fBG;

    return {Omega_yield, Omega_yield_err, Phi_yield, Phi_yield_err};
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

                    if(LikeSignSig && totalfit && p_range && p_min && p_max) {
                        Yield y = yieldCal_CrystalBall(LikeSignSig, totalfit);
                        
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
