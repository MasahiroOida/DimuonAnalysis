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

void correction(TFile *input_Lumi, TFile *input_Yield, TFile *input_efficiency) {
    if (!input_Lumi || !input_Yield || !input_efficiency) {
        std::cerr << "Error: One of the input files is null." << std::endl;
        return;
    }

    TFile *output_file = new TFile("efficiency_corrected.root", "RECREATE");

    // --- 1. Efficiency Histograms (Acc x Eff) ---
    TH1D *Phi_eff_pt_orig = dynamic_cast<TH1D*>(input_efficiency->Get("Phi/AccEff_Phi_Pt"));
    TH1D *Omega_eff_pt_orig = dynamic_cast<TH1D*>(input_efficiency->Get("Omega/AccEff_Omega_Pt"));
    
    if (!Phi_eff_pt_orig || !Omega_eff_pt_orig) {
        std::cerr << "Error: Efficiency histograms not found in AcceptanceEfficiency.root" << std::endl;
        return;
    }

    output_file->cd();
    Phi_eff_pt_orig->Write("Phi_eff_pt_original");
    Omega_eff_pt_orig->Write("Omega_eff_pt_original");

    // --- 2. Find Top Directory in Yield File ---
    TIter nextTop(input_Yield->GetListOfKeys());
    TKey *keyTop;
    TString topDirName = "";
    while ((keyTop = (TKey *)nextTop())) {
        TString className = keyTop->GetClassName();
        if (className == "TDirectoryFile" || className == "TDirectory") {
            topDirName = keyTop->GetName();
            break;
        }
    }

    if (topDirName == "") {
        std::cerr << "Error: No TDirectory found in input_Yield file." << std::endl;
        return;
    }

    TDirectory *topDir = (TDirectory*)input_Yield->Get(topDirName);
    
    // --- 3. Raw Yields (Before Correction) ---
    TGraphErrors *g_omega_yield = (TGraphErrors*)topDir->Get("omega_yield_ptspectrum");
    TGraphErrors *g_phi_yield = (TGraphErrors*)topDir->Get("phi_yield_ptspectrum");
    
    if (!g_omega_yield || !g_phi_yield) {
        std::cerr << "Error: Yield graphs not found in " << topDirName << std::endl;
        return;
    }

    // --- 4. Cross Section Calculation (Luminosity Based) ---
    float BR_omega = 7.4e-5, BR_phi = 2.86e-4;
    // float BR_omega = 0.5, BR_phi = 0.5;
    TH1D *hColl = (TH1D *)input_Lumi->Get("dimuon/Event/norm/hCollisionCounter");
    if (!hColl) { 
        std::cerr << "Error: hCollisionCounter not found." << std::endl; 
        return; 
    }
    
    // Integrated Luminosity calculation
    // 0.0594e6 is assumed to be the pp inelastic cross section in microbarn (59.4 mb)
    double nEvents = hColl->GetBinContent(1);
    double sigma_inel_mb = 59.4; 
    double L_int_mb_inv = nEvents / sigma_inel_mb; 
    
    std::cout << "Luminosity normalization:" << std::endl;
    std::cout << " - Events: " << nEvents << std::endl;
    std::cout << " - Sigma_inel: " << sigma_inel_mb << " mb" << std::endl;
    std::cout << " - L_int: " << L_int_mb_inv << " mb^-1" << std::endl;

    std::vector<float> o_xs, o_xs_err, p_xs, p_xs_err;
    std::vector<float> final_pt_vec, final_pt_err_vec;
    std::vector<float> o_eff_vec, p_eff_vec, o_eff_err_vec, p_eff_err_vec;

    int nPoints = g_omega_yield->GetN();
    for (int i = 0; i < nPoints; i++) {
        double pt, y_o, y_p;
        g_omega_yield->GetPoint(i, pt, y_o);
        g_phi_yield->GetPoint(i, pt, y_p);

        double y_o_err = g_omega_yield->GetErrorY(i);
        double y_p_err = g_phi_yield->GetErrorY(i);
        double pt_err = g_omega_yield->GetErrorX(i);

        int bin_o = Omega_eff_pt_orig->FindBin(pt);
        int bin_p = Phi_eff_pt_orig->FindBin(pt);

        double eff_o = Omega_eff_pt_orig->GetBinContent(bin_o);
        double eff_o_err = Omega_eff_pt_orig->GetBinError(bin_o);
        double eff_p = Phi_eff_pt_orig->GetBinContent(bin_p);
        double eff_p_err = Phi_eff_pt_orig->GetBinError(bin_p);

        o_eff_vec.push_back(eff_o);
        o_eff_err_vec.push_back(eff_o_err);
        p_eff_vec.push_back(eff_p);
        p_eff_err_vec.push_back(eff_p_err);

        // Cross Section [mb/GeV/c] = Corrected Yield / L_int
        // Omega cross section calculation
        if (eff_o > 0) {
            o_xs.push_back(y_o / L_int_mb_inv / BR_omega / eff_o);
            o_xs_err.push_back(y_o_err / L_int_mb_inv / BR_omega / eff_o);
        } else {
            o_xs.push_back(0); o_xs_err.push_back(0);
        }

        if (eff_p > 0) {
            p_xs.push_back(y_p / L_int_mb_inv / BR_phi / eff_p);
            p_xs_err.push_back(y_p_err / L_int_mb_inv / BR_phi / eff_p);
        } else {
            p_xs.push_back(0); p_xs_err.push_back(0);
        }
        
        final_pt_vec.push_back(pt);
        final_pt_err_vec.push_back(pt_err);
    }

    TGraphErrors *g_omega_eff = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), o_eff_vec.data(), final_pt_err_vec.data(), o_eff_err_vec.data());
    TGraphErrors *g_phi_eff = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), p_eff_vec.data(), final_pt_err_vec.data(), p_eff_err_vec.data());
    
    TGraphErrors *g_omega_xs = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), o_xs.data(), final_pt_err_vec.data(), o_xs_err.data());
    TGraphErrors *g_phi_xs = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), p_xs.data(), final_pt_err_vec.data(), p_xs_err.data());
    
    output_file->cd();
    g_omega_eff->Write("Omega_efficiency_rebinned");
    g_phi_eff->Write("Phi_efficiency_rebinned");
    g_omega_xs->Write("Omega_cross_section");
    g_phi_xs->Write("Phi_cross_section");

    // --- 5. Visualization Canvases ---
    gStyle->SetOptStat(0);
    
    // --- Raw Yields ---
    TCanvas *c_Omega_Yield = new TCanvas("c_Omega_Yield", "Omega Raw Yield", 800, 600);
    gPad->SetLogy();
    g_omega_yield->SetTitle("Omega Raw Yield;p_{T} (GeV/c);dN/dp_{T} (raw)");
    g_omega_yield->Draw("APE");
    c_Omega_Yield->Write();

    TCanvas *c_Phi_Yield = new TCanvas("c_Phi_Yield", "Phi Raw Yield", 800, 600);
    gPad->SetLogy();
    g_phi_yield->SetTitle("Phi Raw Yield;p_{T} (GeV/c);dN/dp_{T} (raw)");
    g_phi_yield->Draw("APE");
    c_Phi_Yield->Write();

    // --- Efficiencies ---
    TCanvas *c_Omega_Eff = new TCanvas("c_Omega_Eff", "Omega Efficiency", 800, 600);
    g_omega_eff->SetTitle("Omega Efficiency (Rebinned);p_{T} (GeV/c);Acc #times Eff");
    g_omega_eff->Draw("APE");
    c_Omega_Eff->Write();

    TCanvas *c_Phi_Eff = new TCanvas("c_Phi_Eff", "Phi Efficiency", 800, 600);
    g_phi_eff->SetTitle("Phi Efficiency (Rebinned);p_{T} (GeV/c);Acc #times Eff");
    g_phi_eff->Draw("APE");
    c_Phi_Eff->Write();

    // --- Cross Sections ---
    TCanvas *c_Omega_XS = new TCanvas("c_Omega_XS", "Omega Cross Section", 800, 600);
    gPad->SetLogy();
    g_omega_xs->SetTitle("Omega Differential Cross Section;p_{T} (GeV/c);#frac{Yield}{L_{int} #cdot BR #cdot (Acc#timesEff)} (mb/GeV/c)");
    g_omega_xs->Draw("APE");
    c_Omega_XS->Write();

    TCanvas *c_Phi_XS = new TCanvas("c_Phi_XS", "Phi Cross Section", 800, 600);
    gPad->SetLogy();
    g_phi_xs->SetTitle("Phi Differential Cross Section;p_{T} (GeV/c);#frac{Yield}{L_{int} #cdot BR #cdot (Acc#timesEff)} (mb/GeV/c)");
    g_phi_xs->Draw("APE");
    c_Phi_XS->Write();

    output_file->Close();
    std::cout << "Correction and cross section calculation done. Results saved to 'efficiency_corrected.root'" << std::endl;
}
