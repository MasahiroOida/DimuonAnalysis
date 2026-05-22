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

struct PeakfitResult
{
    TF1 *totalfitfunction;
};

static Double_t myFreeFunction(Double_t *x, Double_t *par)
{
    // 不連続な定義域を持つ例
    if (0.5 <= x[0] && x[0] < 0.70)
    {
        return par[0] * exp(-par[1] * x[0]);
    }
    else if (0.70 <= x[0] && x[0] <= 0.86)
    {
        return 0;
    }
    else if (0.86 < x[0] && x[0] < 0.92)
    {
        return par[0] * exp(-par[1] * x[0]);
    }
    else if (0.92 <= x[0] && x[0] <= 1.15)
    {
        return 0;
    }
    else if (1.15 < x[0] && x[0] <= 1.3)
    {
        return par[0] * exp(-par[1] * x[0]);
    }
    return 0;
};

static void Peakfit(TH1 *Sig, double nEvents = -1)
{
    TH1::SetDefaultSumw2();
    // Clone the Histogram
    if (Sig == nullptr)
    {
        std::cout << "Sig = nullptr. Skipping fit." << std::endl;
        return;
    }
    TH1F *BGHistogram = (TH1F *)Sig->Clone("BGHistogram");
    TH1F *SigHistogram = (TH1F *)Sig->Clone("SigHistogram");
    // Define split position
    double splitLow1 = 0.70;
    double splitHigh1 = 0.86;
    double splitLow2 = 0.92;
    double splitHigh2 = 1.15;

    // Sprit BG from LikesignSignal (Enter the BG to Clone histogram)
    for (int bin = 1; bin <= BGHistogram->GetNbinsX(); ++bin)
    {
        double binCenter = BGHistogram->GetBinCenter(bin);
        if ((binCenter >= splitLow1 && binCenter <= splitHigh1) || (binCenter >= splitLow2 && binCenter <= splitHigh2))
        {
            BGHistogram->SetBinContent(bin, 0);
        }
    }

    // Subtract expected BG
    int binCount = Sig->GetNbinsX();
    double binWidth = Sig->GetBinWidth(1);
    // Define Fit BG Function and Fit to BG(fit function -> [0]*exp(-[1]*x))
    // ===========
    TF1 *fitFunction = new TF1("fitFunction", myFreeFunction, 0.5, 1.3, 3);
    fitFunction->SetParameters(500, 0.7);
    BGHistogram->Fit("fitFunction", "RQ");
    double exact5 = fitFunction->GetParameter(0);
    double exact6 = fitFunction->GetParameter(1);
    TF1 *BGFunction = new TF1("BGFunction", "[0]*exp(-[1]*x)", 0.5, 1.3);
    BGFunction->SetParameters(exact5, exact6);
    // ============

    // Get peak histogram
    for (int bin = 1; bin <= binCount; ++bin)
    {
        double binContent = Sig->GetBinContent(bin);
        double binCenter = Sig->GetBinCenter(bin);
        double funcValue = BGFunction->Eval(binCenter);
        double NewSig = binContent - funcValue;
        SigHistogram->SetBinContent(bin, NewSig);
    }
    // Define Fit Peak Function and Fit to Peak (fit function -> [0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2))
    TF1 *omegafitfunc = new TF1("omegafitfunc", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)", splitLow1, splitHigh1);
    omegafitfunc->SetParameters(exact5 / 10, 0.78, 0.02);
    SigHistogram->Fit("omegafitfunc", "RQ");
    // Defin Fit Peak Function and Fit to Peak

    // Define Fit Peak Function and Fit to Peak (fit function -> [0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2))
    TF1 *phifitfunc = new TF1("phifitfunc", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)", splitLow2, splitHigh2);
    phifitfunc->SetParameters(exact5 / 10, 1.02, 0.02);
    SigHistogram->Fit("phifitfunc", "RQ");

    // 全体のシグナルに対して同時フィット
    double param0 = BGFunction->GetParameter(0);
    double param1 = BGFunction->GetParameter(1);
    double param2 = omegafitfunc->GetParameter(0);
    double param3 = omegafitfunc->GetParameter(1);
    double param4 = omegafitfunc->GetParameter(2);
    double param5 = phifitfunc->GetParameter(0);
    double param6 = phifitfunc->GetParameter(1);
    double param7 = phifitfunc->GetParameter(2);

    TF1 *totalfit = new TF1("totalfit", "[0]*exp(-[1]*x)+[2]*exp(-(x-[3])*(x-[3])/[4]/[4]/2)+[5]*exp(-(x-[6])*(x-[6])/[7]/[7]/2)", 0.5, 1.3);
    
    // Set parameter names first
    totalfit->SetParName(0, "N_{BG}");            // p0 → Constant
    totalfit->SetParName(1, "#alpha");            // p1 → Mean
    totalfit->SetParName(2, "N_{#omega}");        // p0 → Constant
    totalfit->SetParName(3, "#omega mean mass");  // p1 → Mean
    totalfit->SetParName(4, "#omega mass width"); // p2 → Sigma
    totalfit->SetParName(5, "N_{#phi}");          // p0 → Constant
    totalfit->SetParName(6, "#phi mean mass");    // p1 → Mean
    totalfit->SetParName(7, "#phi mass width");   // p2 → Sigma

    totalfit->SetParameters(param0, param1, param2, param3, param4, param5, param6, param7);

    auto clamp = [](double val, double low, double high) {
        if (val <= low) return low + 0.0001;
        if (val >= high) return high - 0.0001;
        return val;
    };

    // Apply limits AFTER setting parameters and clamping
    totalfit->SetParameter(3, clamp(totalfit->GetParameter(3), 0.75, 0.85));
    totalfit->SetParLimits(3, 0.75, 0.85);
    totalfit->SetParameter(4, clamp(totalfit->GetParameter(4), 0.01, 0.5));
    totalfit->SetParLimits(4, 0.01, 0.5);
    totalfit->SetParameter(6, clamp(totalfit->GetParameter(6), 0.9, 1.1));
    totalfit->SetParLimits(6, 0.9, 1.1);
    totalfit->SetParameter(7, clamp(totalfit->GetParameter(7), 0.01, 0.5));
    totalfit->SetParLimits(7, 0.01, 0.5);
    
    double exact1 = totalfit->GetParameter(3);
    double exact2 = totalfit->GetParameter(4);
    double exact3 = totalfit->GetParameter(6);
    double exact4 = totalfit->GetParameter(7);
    TFitResultPtr fitResult = Sig->Fit("totalfit", "RQS");
    
    // Save Peak-only histogram by subtracting final BG (raw counts)
    TH1F *SigOnly = (TH1F *)Sig->Clone("LikeSignSig_PeakOnly");
    SigOnly->SetTitle("Signal after subtracting final BG");
    SigOnly->GetListOfFunctions()->Clear(); // Remove any associated fit functions
    TF1 *fBG_final = new TF1("fBG_final", "[0]*exp(-[1]*x)", 0.5, 1.3);
    fBG_final->SetParameters(totalfit->GetParameter(0), totalfit->GetParameter(1));
    for (int bin = 1; bin <= SigOnly->GetNbinsX(); ++bin) {
        double binCenter = SigOnly->GetBinCenter(bin);
        SigOnly->SetBinContent(bin, SigOnly->GetBinContent(bin) - fBG_final->Eval(binCenter));
    }

    double Norm = totalfit->GetParameter(2);
    double Norm2 = totalfit->GetParameter(5);
    double Norm3 = totalfit->GetParameter(0);
    double Norm4 = totalfit->GetParameter(1);

    // Keep raw counts (do not normalize by event count and bin width)
    Sig->SetYTitle("Counts / (Mass bin)");
    SigOnly->SetYTitle("Counts / (Mass bin)");
    SigHistogram->SetYTitle("Counts / (Mass bin)");
    BGHistogram->SetYTitle("Counts / (Mass bin)");

    TF1 *func5 = new TF1("func5", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)+[3]*exp(-(x-[4])*(x-[4])/[5]/[5]/2)", 0.5, 1.3);
    func5->SetParameters(Norm, exact1, exact2, Norm2, exact3, exact4);

    TF1 *omegafunc = new TF1("omegafunc", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)", 0, 2);
    omegafunc->SetParameters(Norm, exact1, exact2);
    TF1 *phifunc = new TF1("phifitfunc", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)", 0, 2);
    phifunc->SetParameters(Norm2, exact3, exact4);
    TF1 *BGfunc = new TF1("BGFunction", "[0]*exp(-[1]*x)", 0.5, 1.3);
    BGfunc->SetParameters(Norm3, Norm4);

    // Define Color
    Sig->SetLineColor(kBlack);
    omegafitfunc->SetLineColor(kBlack);
    SigHistogram->SetLineColor(kBlue);
    BGHistogram->SetLineColor(kBlue);
    // Define TCanvas
    // gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1);
    gStyle->SetStatColor(0);      // 背景色を透明にする
    gStyle->SetStatBorderSize(0); // 枠線を消す
    gStyle->SetStatX(0.87);       // X座標（右端が1.0、左端が0.0）
    gStyle->SetStatY(0.87);       // Y座標（上端が1.0、下端が0.0）
    gStyle->SetStatW(0.2);        // 幅
    gStyle->SetStatH(0.15);       // 高さ

    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 800, 600);
    TLegend *leg1 = new TLegend(0.65, 0.25, 0.8, 0.45);
    canvas1->SetFillColor(kWhite); // キャンバスの背景色を変更
    canvas1->SetGrid();            // グリッドを表示
    canvas1->SetLeftMargin(0.15);
    leg1->SetBorderSize(0);  // 枠線を消す
    leg1->SetFillStyle(0);   // 背景透明
    leg1->SetTextSize(0.04); // 文字サイズ調整
    // canvas1 -> SetLogy();             // y軸を対数スケールに変更
    Sig->Draw();
    Sig->GetXaxis()->SetTitle("Dimuon invariant mass GeV/c^{2}");
    Sig->GetYaxis()->SetTitle("dN/dm");
    Sig->SetTitle("");
    Sig->GetXaxis()->SetRangeUser(0, 2);
    Sig->GetYaxis()->SetRangeUser(0, Norm3 * 2.5); 
    Sig->SetMarkerSize(0.7);
    Sig->SetMarkerStyle(8);
    totalfit->Draw("same");
    BGfunc->Draw("same");
    omegafunc->Draw("same");
    phifunc->Draw("same");
    leg1->AddEntry(Sig, "Dimuon Signal", "lep"); // "lep" は点と線を表示
    leg1->AddEntry(totalfit, "total fit", "l");
    leg1->AddEntry(BGfunc, "BG fit", "l");
    leg1->AddEntry(omegafunc, "#omega peak", "l"); // "l" は線を表示
    leg1->AddEntry(phifunc, "#phi peak", "l");
    leg1->Draw("same");
    totalfit->SetNpx(10000);
    totalfit->SetLineColor(kRed);
    BGfunc->SetNpx(10000);
    BGfunc->SetLineColor(kBlack);
    BGfunc->SetLineStyle(2);
    omegafunc->SetNpx(10000);
    omegafunc->SetLineColor(kGreen);
    phifunc->SetNpx(10000);
    phifunc->SetLineColor(kBlue);
    canvas1->Update();

    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
    canvas2->SetFillColor(kWhite); // キャンバスの背景色を変更
    canvas2->SetGrid();            // グリッドを表示
    SigHistogram->Draw();
    SigHistogram->GetXaxis()->SetTitle("Dimuon invariant mass GeV/c^{2}");
    SigHistogram->GetYaxis()->SetTitle("dN/dm");
    SigHistogram->SetTitle("Likesign method Signal");
    SigHistogram->GetXaxis()->SetRangeUser(0, 2);
    omegafitfunc->Draw("same");

    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 800, 600);
    canvas3->SetFillColor(kWhite); // キャンバスの背景色を変更
    canvas3->SetGrid();            // グリッドを表示
    BGHistogram->Draw();
    BGHistogram->GetXaxis()->SetTitle("Dimuon invariant mass GeV/c^{2}");
    BGHistogram->GetYaxis()->SetTitle("dN/dm");
    BGHistogram->SetTitle("Likesign method Signal");
    BGHistogram->GetXaxis()->SetRangeUser(0, 2);
    BGFunction->Draw("same");

    canvas1->Write("Total Fit");
    canvas2->Write("Peak Fit");
    canvas3->Write("BG Fit");
    Sig->Write("LikesignSig_fit");
    BGHistogram->Write("BGHistogram");
    SigHistogram->Write("SigHistogram");

    omegafitfunc->Write("omegafitfunc");
    phifitfunc->Write("phifitfunc");
    totalfit->Write("totalfit");
    fitFunction->Write();
    BGFunction->Write();

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
    TF1 *fBG_draw = new TF1("fBG_draw", "[0]*exp(-[1]*x)", 0.6, 1.2);
    fBG_draw->SetParameters(totalfit->GetParameter(0), totalfit->GetParameter(1));
    fBG_draw->SetLineColor(kGreen+2); fBG_draw->SetLineStyle(2); fBG_draw->SetLineWidth(2);

    TF1 *fOmega_draw = new TF1("fOmega_draw", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)", 0.6, 1.2);
    fOmega_draw->SetParameters(totalfit->GetParameter(2), totalfit->GetParameter(3), totalfit->GetParameter(4));
    fOmega_draw->SetLineColor(kBlue); fOmega_draw->SetLineStyle(2); fOmega_draw->SetLineWidth(2);

    TF1 *fPhi_draw = new TF1("fPhi_draw", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2]/2)", 0.6, 1.2);
    fPhi_draw->SetParameters(totalfit->GetParameter(5), totalfit->GetParameter(6), totalfit->GetParameter(7));
    fPhi_draw->SetLineColor(kRed); fPhi_draw->SetLineStyle(2); fPhi_draw->SetLineWidth(2);

    fBG_draw->Draw("SAME"); fOmega_draw->Draw("SAME"); fPhi_draw->Draw("SAME");
    totalfit->SetLineColor(kMagenta); totalfit->SetLineWidth(3);
    totalfit->Draw("SAME");

    TLegend *legOverlay = new TLegend(0.6, 0.6, 0.88, 0.88);
    legOverlay->AddEntry(Sig, "Data (incl. BG)", "lep");
    legOverlay->AddEntry(SigOnly, "Data (BG sub)", "lep");
    legOverlay->AddEntry(totalfit, "Total Fit", "l");
    legOverlay->AddEntry(fOmega_draw, "#omega Signal (Gauss)", "l");
    legOverlay->AddEntry(fPhi_draw, "#phi Signal (Gauss)", "l");
    legOverlay->AddEntry(fBG_draw, "Background (Exp)", "l");
    legOverlay->Draw();

    cOverlay->Write();
    delete cOverlay; delete legOverlay; delete fBG_draw; delete fOmega_draw; delete fPhi_draw;

    delete BGHistogram;
    delete SigHistogram;
    delete fitFunction;
    delete BGFunction;
    delete omegafitfunc;
    delete omegafunc;
    delete phifitfunc;
    delete phifunc;
    delete totalfit;
    delete canvas1;
    delete canvas2;
    delete canvas3;
    delete SigOnly;
    delete fBG_final;
}

void PeakFit_Gauss(TFile *input_Sig){
    TFile *outfile_PeakFit = new TFile("PeakFit_Gauss_Results.root", "RECREATE");
    TIter next(input_Sig->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
        if (TString(key->GetClassName()) == "TDirectoryFile")
        {
            TDirectory *topDir = (TDirectory *)key->ReadObj();
            TDirectory *outTopDir = outfile_PeakFit->mkdir(topDir->GetName());
            
            TIter nextPt(topDir->GetListOfKeys());
            TKey *keyPt;
            while ((keyPt = (TKey *)nextPt()))
            {
                if (TString(keyPt->GetClassName()) == "TDirectoryFile")
                {
                    TDirectory *ptDir = (TDirectory *)keyPt->ReadObj();
                    TDirectory *outPtDir = outTopDir->mkdir(ptDir->GetName());
                    outfile_PeakFit->cd();
                    outPtDir->cd();
                    
                    TH1F *LikeSignSig = (TH1F *)ptDir->Get("LikeSignSig");
                    TParameter<double> *pN = (TParameter<double> *)ptDir->Get("nEvents");
                    double nEvents = (pN) ? pN->GetVal() : -1;
                    Peakfit(LikeSignSig, nEvents);
                    
                    TParameter<double>* param_ptrange = (TParameter<double>*)ptDir->Get("ptrange");
                    if (param_ptrange) param_ptrange->Write("ptrange");
                    TParameter<double>* param_ptmin = (TParameter<double>*)ptDir->Get("ptmin");
                    if (param_ptmin) param_ptmin -> Write("ptmin");
                    TParameter<double>* param_ptmax = (TParameter<double>*)ptDir->Get("ptmax");
                    if (param_ptmax) param_ptmax -> Write("ptmax"); 
                    if (pN) pN->Write("nEvents");
                }
            }
        }
    }
    outfile_PeakFit -> Close();
    delete outfile_PeakFit;
    std::cout << "Finish to fit using Gauss function" << std::endl;
}
