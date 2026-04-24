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
    
    // Integrated Luminosity calculation using TVX counter
    TH1D *hTVX = (TH1D *)input_Lumi->Get("dimuon/BC/hTVXCounter");
    if (!hTVX) { 
        std::cerr << "Error: hTVXCounter not found in dimuon/BC/." << std::endl; 
        return; 
    }

    // --- TVX Efficiency Correction Calculation ---
    TFile *file_vtx_eff = TFile::Open("vertexZ_efficiency.root");
    if (!file_vtx_eff || file_vtx_eff->IsZombie()) {
        std::cerr << "Error: Could not open vertexZ_efficiency.root" << std::endl;
        return;
    }
    TH1D *hVtxZEff = (TH1D*)file_vtx_eff->Get("hVertexZEfficiency");
    if (hVtxZEff) hVtxZEff->SetDirectory(0);
    file_vtx_eff->Close();

    TH1D *hZvtx = (TH1D*)input_Lumi->Get("dimuon/Event/norm/hZvtx");

    if (!hVtxZEff || !hZvtx) {
        std::cerr << "Error: hVertexZEfficiency or hZvtx not found." << std::endl;
        return;
    }

    double sum_num = 0;
    double sum_den = 0;
    TH1D *hZvtx_corrected = (TH1D*)hZvtx->Clone("hZvtx_corrected");
    hZvtx_corrected->Reset();
    hZvtx_corrected->SetTitle("VtxZ distribution before and after efficiency correction;VtxZ (cm);Counts");

    for (int i = 1; i <= hZvtx->GetNbinsX(); ++i) {
        double z = hZvtx->GetXaxis()->GetBinCenter(i);
        double val_vtx = hZvtx->GetBinContent(i);
        double err_vtx = hZvtx->GetBinError(i);
        double eff_z = hVtxZEff->GetBinContent(hVtxZEff->FindBin(z));
        
        if (eff_z > 1e-6) { // Avoid division by zero
            double term = val_vtx / eff_z;
            double term_err = err_vtx / eff_z;
            sum_den += term;
            hZvtx_corrected->SetBinContent(i, term);
            hZvtx_corrected->SetBinError(i, term_err);
            if (z >= -10.0 && z <= 10.0) {
                sum_num += term;
            }
        }
    }
    double epsilon_tvx = (sum_den > 0) ? (sum_num / sum_den) : 0;
    std::cout << "TVX Efficiency (epsilon_tvx): " << epsilon_tvx << std::endl;

    // --- Plot VtxZ Comparison ---
    TCanvas *c_VtxZ_Correction = new TCanvas("c_VtxZ_Correction", "VtxZ Efficiency Correction Comparison", 800, 600);
    hZvtx_corrected->SetLineColor(kRed);
    hZvtx_corrected->SetMarkerColor(kRed);
    hZvtx->SetLineColor(kBlack);
    hZvtx->SetMarkerColor(kBlack);

    hZvtx_corrected->Draw("PE");
    hZvtx->Draw("PE SAME");

    TLegend *legVtx = new TLegend(0.6, 0.7, 0.85, 0.85);
    legVtx->AddEntry(hZvtx, "Original VtxZ", "lp");
    legVtx->AddEntry(hZvtx_corrected, "Corrected VtxZ (VtxZ/#epsilon(z))", "lp");
    legVtx->Draw();
    
    output_file->cd();
    c_VtxZ_Correction->Write();
    hZvtx_corrected->Write();

    TH1D *hColl = (TH1D *)input_Lumi->Get("dimuon/Event/norm/hCollisionCounter");
    
    // TVX events for normalization are typically in specific bins. 
    // Using bin at 5.0 to represent the 4.5 to 5.5 range.
    double nEvents = hTVX->GetBinContent(hTVX->FindBin(5.0));
    //double nEvents = hColl->GetBinContent(hColl->FindBin(21));
    double sigma_tvx_mb = 54.0; // TVX cross section for pp 13.6 TeV [mb]
    double L_int_mb_inv = (nEvents * epsilon_tvx) / sigma_tvx_mb; 
    
    std::cout << "Luminosity normalization (TVX based):" << std::endl;
    std::cout << " - TVX Events (total): " << nEvents << std::endl;
    std::cout << " - TVX Efficiency (epsilon_tvx): " << epsilon_tvx << std::endl;
    std::cout << " - TVX Events (corrected): " << nEvents * epsilon_tvx << std::endl;
    std::cout << " - Sigma_TVX: " << sigma_tvx_mb << " mb" << std::endl;
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
