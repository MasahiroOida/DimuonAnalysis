#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

void SetStyle(TH1D* h, int color, int marker, double size, TString xTitle) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(size);
    h->SetStats(0);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
}

void PlotAcceptanceEfficiency() {
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);
    
    TFile *fileIn = new TFile("AcceptanceEfficiency.root", "READ");
    if (!fileIn || fileIn->IsZombie()) {
        printf("Error: Could not open AcceptanceEfficiency.root\n");
        return;
    }

    TString variables[] = {"Pt", "Y"};
    TString xTitles[] = {"p_{T} (GeV/c)", "y"};

    for (int i = 0; i < 2; ++i) {
        TString var = variables[i];
        TString xTitle = xTitles[i];

        // --- 1. Counts Canvas ---
        TCanvas *cCounts = new TCanvas(Form("cCounts_%s", var.Data()), Form("Counts vs %s", var.Data()), 1200, 600);
        cCounts->Divide(2, 1);

        TString particles[] = {"Omega", "Phi"};
        for (int p = 0; p < 2; ++p) {
            cCounts->cd(p + 1);
            gPad->SetLogy();
            gPad->SetTicks(1, 1);
            
            TH1D *hGen = (TH1D*)fileIn->Get(Form("%s/hGen%s_%s", particles[p].Data(), var.Data(), particles[p].Data()));
            TH1D *hAcc = (TH1D*)fileIn->Get(Form("%s/hAcc%s_%s", particles[p].Data(), var.Data(), particles[p].Data()));
            TH1D *hRec = (TH1D*)fileIn->Get(Form("%s/hRec%s_%s", particles[p].Data(), var.Data(), particles[p].Data()));

            if (hGen && hAcc && hRec) {
                SetStyle(hGen, kBlack, 20, 0.6, xTitle);
                SetStyle(hAcc, kBlue, 20, 0.6, xTitle);
                SetStyle(hRec, kGreen+2, 20, 0.6, xTitle);
                hGen->GetYaxis()->SetTitle("Yield");
                
                hGen->Draw("P");
                hAcc->Draw("P same");
                hRec->Draw("P same");

                TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetHeader(Form("#bf{%s}", particles[p].Data()));
                leg->AddEntry(hGen, "Generated", "p");
                leg->AddEntry(hAcc, "Accepted", "p");
                leg->AddEntry(hRec, "Reconstructed", "p");
                leg->Draw();
            }
        }
        cCounts->SaveAs(Form("Combined_Counts_%s.png", var.Data()));

        // --- 2. Ratios Canvas ---
        TCanvas *cRatios = new TCanvas(Form("cRatios_%s", var.Data()), Form("Ratios vs %s", var.Data()), 1200, 600);
        cRatios->Divide(2, 1);

        for (int p = 0; p < 2; ++p) {
            cRatios->cd(p + 1);
            gPad->SetLogy();
            gPad->SetTicks(1, 1);

            TH1D *hAcc = (TH1D*)fileIn->Get(Form("%s/Acceptance_%s_%s", particles[p].Data(), particles[p].Data(), var.Data()));
            TH1D *hEff = (TH1D*)fileIn->Get(Form("%s/Efficiency_%s_%s", particles[p].Data(), particles[p].Data(), var.Data()));
            TH1D *hAE  = (TH1D*)fileIn->Get(Form("%s/AccEff_%s_%s", particles[p].Data(), particles[p].Data(), var.Data()));

            if (hAcc && hEff && hAE) {
                SetStyle(hAcc, kBlack, 20, 0.6, xTitle);
                SetStyle(hEff, kBlue, 20, 0.6, xTitle);
                SetStyle(hAE, kGreen+2, 20, 0.6, xTitle);
                hAcc->GetYaxis()->SetRangeUser(0.001, 1.0);
                hAcc->GetYaxis()->SetTitle("Acceptance / Efficiency");

                hAcc->Draw("P");
                hEff->Draw("P same");
                hAE->Draw("P same");

                TLegend *leg = new TLegend(0.55, 0.15, 0.88, 0.35);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetHeader(Form("#bf{%s}", particles[p].Data()));
                leg->AddEntry(hAcc, "Acceptance", "p");
                leg->AddEntry(hEff, "Efficiency", "p");
                leg->AddEntry(hAE, "Acc * Eff", "p");
                leg->Draw();
            }
        }
        cRatios->SaveAs(Form("Combined_Ratios_%s.png", var.Data()));
    }
}
