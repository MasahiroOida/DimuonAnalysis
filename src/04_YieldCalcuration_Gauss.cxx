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

struct Yield
{
    double Omega_yield;
    double Omega_yield_err;
    double Phi_yield;
    double Phi_yield_err;
};

// input (Sig, totalfit function) output (omega,phiの収量＋エラー)
Yield yieldCal_Gauss(TH1 *Sig, TF1 *Totalfun, double ptrange)
{
    TH1::SetDefaultSumw2();
    // Resofun = [0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)+[3]*exp(-(x-[4])*(x-[4])/[5]/[5]/2)
    // Totalfun = [0]*exp(-[1]*x)+[2]*exp(-(x-[3])*(x-[3])/[4]/[4]/2)+[5]*exp(-(x-[6])*(x-[6])/[7]/[7]/2)
    // 同じ設定のhistを定義する
    TAxis *xAxis = Sig->GetXaxis();
    int numBins = xAxis->GetNbins();

    TH1D *Peakhist = (TH1D *)Sig->Clone("Peakhist");
    Peakhist->Reset();
    TH1D *BGhist = (TH1D *)Sig->Clone("BGhist");
    BGhist->Reset();
    //TH1F *Peakhist = new TH1F("Peakhist", "Peakhist", numBins, xMin, xMax);
    //TH1F *BGhist = new TH1F("BGhist", "BGhist", numBins, xMin, xMax);

    // Totalfit functionからBG部分のfunctionを抽出し、BGfuncからBGのヒストグラムを構成→SighistからBGhistを引くことでPeakのみのヒストグラムを作る
    TF1 *BGfunction = new TF1("BGfunction", "[0]*exp(-[1]*x)", 0.5, 1.3);
    TF1 *Peakfunction = new TF1("Peakfunction", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)+[3]*exp(-(x-[4])*(x-[4])/[5]/[5]/2)", 0.5, 1.3);
    double par1 = Totalfun->GetParameter(0);
    double par2 = Totalfun->GetParameter(1);
    BGfunction->SetParameters(par1, par2);

    for (int i = 1; i < numBins; ++i)
    {
        double Binrange = xAxis->GetBinWidth(i);
        double MassDown = i * Binrange;
        double MassUP = (i + 1) * Binrange;
        double MeanMass = (MassUP + MassDown) / 2;
        double BGfunValue = BGfunction->Eval(MeanMass);
        BGhist->SetBinContent(i + 1, BGfunValue);
    };
    Peakhist->Add(Sig, BGhist, 1, -1);

    // Resofunの3sigmaの範囲を決める(Omega,Phiの両方を)
    double Omega_Norm = Totalfun->GetParameter(2);
    double Omega_mass = abs(Totalfun->GetParameter(3));
    double Omega_mass_err = abs(Totalfun->GetParError(3));
    double Omega_width = abs(Totalfun->GetParameter(4));
    double Omega_width_err = abs(Totalfun->GetParError(4));
    double Phi_Norm = Totalfun->GetParameter(5);
    double Phi_mass = abs(Totalfun->GetParameter(6));
    double Phi_mass_err = abs(Totalfun->GetParError(6));
    double Phi_width = abs(Totalfun->GetParameter(7));
    double Phi_width_err = abs(Totalfun->GetParError(7));

    // Peak部分のフィット関数にパラメーターに入力
    Peakfunction->SetParameters(Omega_Norm, Omega_mass, Omega_width, Phi_Norm, Phi_mass, Phi_width);

    double Omega_Down = Omega_width * 3 * (-1) + Omega_mass;
    double Omega_UP = Omega_width * 3 + Omega_mass;
    double Phi_Down = Phi_width * 3 * (-1) + Phi_mass;
    double Phi_UP = Phi_width * 3 + Phi_mass;

    // 3sigamaの範囲に対応するbin数を探す
    int Omega_UPbin = Peakhist->FindBin(Omega_UP);
    int Omega_Downbin = Peakhist->FindBin(Omega_Down);
    int Phi_UPbin = Peakhist->FindBin(Phi_UP);
    int Phi_Downbin = Peakhist->FindBin(Phi_Down);

    // 指定した領域の各binのエントリーを積分する
    double Omega_Peakentry = 0;
    double Omega_error = 0;
    double bin_entry = 0;
    double bin_error = 0;

    for (int i = Omega_Downbin; i <= Omega_UPbin; ++i)
    {
        bin_entry = Peakhist->GetBinContent(i);
        Omega_Peakentry = bin_entry + Omega_Peakentry;
        bin_entry = 0;
        bin_error = Peakhist->GetBinError(i);
        Omega_error = bin_error * bin_error + Omega_error;
        bin_error = 0;
    }; // Omega_Peakentry,Omega_errorを計算
    Omega_error = TMath::Sqrt(Omega_error); // error^2をerror

    double Phi_Peakentry = 0;
    double Phi_error = 0;
    for (int i = Phi_Downbin; i <= Phi_UPbin; ++i)
    {
        bin_entry = Peakhist->GetBinContent(i);
        Phi_Peakentry = bin_entry + Phi_Peakentry;
        bin_entry = 0;
        bin_error = Peakhist->GetBinError(i);
        Phi_error = bin_error * bin_error + Phi_error;
        bin_error = 0;
    }; // Phi_Peakentry,Phi_errorを計算
    Phi_error = TMath::Sqrt(Phi_error); // error^2をerror

    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 800, 600);
    canvas1->SetFillColor(kWhite); // キャンバスの背景色を変更
    canvas1->SetGrid();            // グリッドを表示
    Peakhist->Draw();
    Peakfunction->Draw("same");
    Peakhist->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    Peakhist->GetYaxis()->SetTitle("Count");
    Peakhist->SetTitle("Signal Histogram");
    canvas1->Write("Signalhisto");
    delete Peakhist;
    delete BGhist;
    delete canvas1;
    delete Peakfunction;
    return {Omega_Peakentry / ptrange ,Omega_error / ptrange, Phi_Peakentry / ptrange, Phi_error / ptrange};
}

void YieldCalcuration_Gauss(TFile *input_FitResult){
    TFile *outfile_Yield = new TFile("Yield_Results.root", "RECREATE");
    std::vector<std::string> directoryNames;
    TIter next(input_FitResult->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
        // オブジェクトのクラス名を取得
        std::string className = key->GetClassName();

        // TDirectoryの場合、その名前をベクトルに追加
        if (className == "TDirectoryFile")
        {
            TDirectory *dir = (TDirectory *)key->ReadObj();
            directoryNames.push_back(dir->GetName());
        }
    }
    
    std::vector<double> v_pt, v_pt_err, v_omega_yield, v_omega_yield_err, v_phi_yield, v_phi_yield_err;

    for (const auto &name : directoryNames)
    { 
        TDirectoryFile *dir = (TDirectoryFile *)input_FitResult->Get(name.c_str());
        TDirectory *newDir = outfile_Yield->mkdir(name.c_str());
        outfile_Yield->cd();
        newDir->cd();
        TH1F *LikeSignSig = (TH1F *)dir->Get("LikesignSig_fit");

        TF1 *totalfit = (TF1 *)dir->Get("totalfit");
        TParameter<double>* param_ptrange = (TParameter<double>*)dir->Get("ptrange");
        TParameter<double>* param_ptmin = (TParameter<double>*)dir->Get("ptmin");
        TParameter<double>* param_ptmax = (TParameter<double>*)dir->Get("ptmax");

        if(LikeSignSig && totalfit && param_ptrange && param_ptmin && param_ptmax) {
            Yield yield = yieldCal_Gauss(LikeSignSig, totalfit, param_ptrange->GetVal());
            
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
    outfile_Yield->cd();
    if (!v_pt.empty()) {
        TGraphErrors *omega_yield_ptspectrum = new TGraphErrors(v_pt.size(), &v_pt[0], &v_omega_yield[0], &v_pt_err[0], &v_omega_yield_err[0]);
        TGraphErrors *phi_yield_ptspectrum = new TGraphErrors(v_pt.size(), &v_pt[0], &v_phi_yield[0], &v_pt_err[0], &v_phi_yield_err[0]);

        TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
        omega_yield_ptspectrum->Draw("P");
        omega_yield_ptspectrum->GetXaxis()->SetTitle("pt (GeV/#it{c})");
        omega_yield_ptspectrum->GetYaxis()->SetTitle("N");
        omega_yield_ptspectrum->SetMarkerSize(0.7);
        omega_yield_ptspectrum->SetMarkerStyle(20);
        omega_yield_ptspectrum->SetTitle("#omega yield");
        omega_yield_ptspectrum->GetXaxis()->SetRangeUser(0, 10);
        omega_yield_ptspectrum->GetYaxis()->SetRangeUser(1, 100000);
        gPad->Update();
        canvas -> Write("omega_yield_ptspectrum_canvas");
        omega_yield_ptspectrum->Write("omega_yield_ptspectrum");

        TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
        phi_yield_ptspectrum->Draw("P");
        phi_yield_ptspectrum->GetXaxis()->SetTitle("pt (GeV/#it{c})");
        phi_yield_ptspectrum->GetYaxis()->SetTitle("N");
        phi_yield_ptspectrum->SetTitle("#phi yield");
        phi_yield_ptspectrum->SetMarkerSize(0.7);
        phi_yield_ptspectrum->SetMarkerStyle(20);
        phi_yield_ptspectrum->GetXaxis()->SetRangeUser(0, 10);
        phi_yield_ptspectrum->GetYaxis()->SetRangeUser(1, 100000);
        gPad->Update();
        canvas2 -> Write("phi_yield_ptspectrum_canvas");
        phi_yield_ptspectrum->Write("phi_yield_ptspectrum");
    }
    outfile_Yield->Close();
    std::cout << "finish calculating yield" << std::endl;
    delete outfile_Yield;
    delete omega_yield_ptspectrum;
    delete phi_yield_ptspectrum;
}
