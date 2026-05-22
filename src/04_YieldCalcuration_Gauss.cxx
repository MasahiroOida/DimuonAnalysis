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

struct Yield
{
    double Omega_yield;
    double Omega_yield_err;
    double Phi_yield;
    double Phi_yield_err;
};

// input (Sig, totalfit function) output (omega,phiの収量＋エラー)
Yield yieldCal_Gauss(TH1 *Sig, TF1 *Totalfun)
{
    if (!Sig || !Totalfun) return {0,0,0,0};

    double Omega_mass = Totalfun->GetParameter(3);
    double Omega_width = Totalfun->GetParameter(4);
    double Phi_mass = Totalfun->GetParameter(6);
    double Phi_width = Totalfun->GetParameter(7);

    double Omega_Down = Omega_mass - 3 * Omega_width;
    double Omega_UP = Omega_mass + 3 * Omega_width;
    double Phi_Down = Phi_mass - 3 * Phi_width;
    double Phi_UP = Phi_mass + 3 * Phi_width;

    TF1 *fBG = new TF1("fBG", "[0]*exp(-[1]*x)", 0.5, 1.3);
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

void YieldCalcuration_Gauss(TFile *input_FitResult){
    TFile *outfile_Yield = new TFile("Yield_Results.root", "RECREATE");
    TIter next(input_FitResult->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
        if (TString(key->GetClassName()) == "TDirectoryFile")
        {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TDirectory *outTopDir = outfile_Yield->mkdir(topDir->GetName());
            
            std::vector<double> v_pt, v_pt_err, v_omega_yield, v_omega_yield_err, v_phi_yield, v_phi_yield_err;

            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt()))
            {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile")
                {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    TDirectory *outPtDir = outTopDir->mkdir(ptDir->GetName());
                    outfile_Yield->cd();
                    outPtDir->cd();
                    
                    TH1F *LikeSignSig = (TH1F *)ptDir->Get("LikesignSig_fit");
                    TF1 *totalfit = (TF1 *)ptDir->Get("totalfit");
                    TParameter<double>* param_ptrange = (TParameter<double>*)ptDir->Get("ptrange");
                    TParameter<double>* param_ptmin = (TParameter<double>*)ptDir->Get("ptmin");
                    TParameter<double>* param_ptmax = (TParameter<double>*)ptDir->Get("ptmax");

                    if(LikeSignSig && totalfit && param_ptrange && param_ptmin && param_ptmax) {
                        Yield yield = yieldCal_Gauss(LikeSignSig, totalfit);
                        
                        // Save results to directory
                        (new TParameter<double>("omega_yield", yield.Omega_yield))->Write();
                        (new TParameter<double>("omega_yield_err", yield.Omega_yield_err))->Write();
                        (new TParameter<double>("phi_yield", yield.Phi_yield))->Write();
                        (new TParameter<double>("phi_yield_err", yield.Phi_yield_err))->Write();
                        param_ptrange->Write("ptrange");
                        param_ptmin -> Write("ptmin");
                        param_ptmax -> Write("ptmax");

                        // Only add to spectrum if not total pt bin
                        if (param_ptrange->GetVal() < 9.0) {
                            v_pt.push_back(0.5 * (param_ptmax->GetVal() + param_ptmin->GetVal()));
                            v_pt_err.push_back(0.5 * param_ptrange->GetVal());
                            v_omega_yield.push_back(yield.Omega_yield);
                            v_omega_yield_err.push_back(yield.Omega_yield_err);
                            v_phi_yield.push_back(yield.Phi_yield);
                            v_phi_yield_err.push_back(yield.Phi_yield_err);
                        }
                    }
                }
            }

            // サマリーの作成
            if (!v_pt.empty()) {
                outTopDir->cd();
                TGraphErrors *omega_yield_ptspectrum = new TGraphErrors(v_pt.size(), &v_pt[0], &v_omega_yield[0], &v_pt_err[0], &v_omega_yield_err[0]);
                TGraphErrors *phi_yield_ptspectrum = new TGraphErrors(v_pt.size(), &v_pt[0], &v_phi_yield[0], &v_pt_err[0], &v_phi_yield_err[0]);

                TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
                omega_yield_ptspectrum->Draw("P");
                omega_yield_ptspectrum->GetXaxis()->SetTitle("pt (GeV/#it{c})");
                omega_yield_ptspectrum->GetYaxis()->SetTitle("N");
                omega_yield_ptspectrum->SetMarkerSize(0.7);
                omega_yield_ptspectrum->SetMarkerStyle(20);
                omega_yield_ptspectrum->SetTitle(Form("#omega yield (%s)", topDir->GetName()));
                omega_yield_ptspectrum->GetXaxis()->SetRangeUser(0, 10);
                omega_yield_ptspectrum->GetYaxis()->SetRangeUser(1, 100000);
                gPad->Update();
                canvas -> Write("omega_yield_ptspectrum_canvas");
                omega_yield_ptspectrum->Write("omega_yield_ptspectrum");

                TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
                phi_yield_ptspectrum->Draw("P");
                phi_yield_ptspectrum->GetXaxis()->SetTitle("pt (GeV/#it{c})");
                phi_yield_ptspectrum->GetYaxis()->SetTitle("N");
                phi_yield_ptspectrum->SetTitle(Form("#phi yield (%s)", topDir->GetName()));
                phi_yield_ptspectrum->SetMarkerSize(0.7);
                phi_yield_ptspectrum->SetMarkerStyle(20);
                phi_yield_ptspectrum->GetXaxis()->SetRangeUser(0, 10);
                phi_yield_ptspectrum->GetYaxis()->SetRangeUser(1, 100000);
                gPad->Update();
                canvas2 -> Write("phi_yield_ptspectrum_canvas");
                phi_yield_ptspectrum->Write("phi_yield_ptspectrum");
                
                delete omega_yield_ptspectrum;
                delete phi_yield_ptspectrum;
                delete canvas;
                delete canvas2;
            }
        }
    }
    outfile_Yield->Close();
    std::cout << "finish calculating yield" << std::endl;
    delete outfile_Yield;
}
