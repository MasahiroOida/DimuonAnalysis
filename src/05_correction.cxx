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

void correction(TFile *input_Lumi, TFile *input_Yield, TFile *input_efficiency, TString fitType) {
    if (!input_Lumi || !input_Yield || !input_efficiency) {
        std::cerr << "Error: One of the input files is null." << std::endl;
        return;
    }

    TFile *output_file = new TFile(Form("efficiency_corrected_%s.root", fitType.Data()), "RECREATE");

    // --- 1. Compute Shared Corrections (Luminosity, VtxZ, etc.) ---
    float BR_omega = 7.4e-5, BR_phi = 2.86e-4;
    
    // Integrated Luminosity calculation using TVX counter
    TH1D *hTVX = (TH1D *)input_Lumi->Get("dimuon/BC/hTVXCounter");
    if (!hTVX) { 
        std::cerr << "Error: hTVXCounter not found in dimuon/BC/." << std::endl; 
        return; 
    }

    // --- TVX Efficiency Correction Calculation ---
    TFile *file_signal_loss = TFile::Open("Signal_loss.root");
    if (!file_signal_loss || file_signal_loss->IsZombie()) {
        std::cerr << "Error: Could not open Signal_loss.root" << std::endl;
        return;
    }
    TH1D *hEffPt_Omega_SL = (TH1D*)file_signal_loss->Get("Omega/hLossPt_Total_SignalLoss_Omega");
    TH1D *hEffPt_Phi_SL = (TH1D*)file_signal_loss->Get("Phi/hLossPt_Total_SignalLoss_Phi");
    if (!hEffPt_Omega_SL || !hEffPt_Phi_SL) {
        std::cerr << "Error: Omega/hLossPt_Total_SignalLoss_Omega or Phi/hLossPt_Total_SignalLoss_Phi not found in Signal_loss.root" << std::endl;
        return;
    }
    hEffPt_Omega_SL->SetDirectory(0);
    hEffPt_Phi_SL->SetDirectory(0);
    file_signal_loss->Close();

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
        
        if (eff_z > 1e-6) {
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

    double epsilon_evsel = 1.0;
    double nEvents = hTVX->GetBinContent(6);
    double sigma_tvx_mb = 54.0;
    double L_int_mb_inv = (nEvents / epsilon_evsel * epsilon_tvx) / sigma_tvx_mb; 
    
    std::cout << "L_int: " << L_int_mb_inv << " mb^-1" << std::endl;

    // --- 2. Loop over Top Directories in Yield File ---
    TIter nextTop(input_Yield->GetListOfKeys());
    TKey *keyTop;
    while ((keyTop = (TKey *)nextTop())) {
        if (TString(keyTop->GetClassName()) != "TDirectoryFile") continue;
        TString topDirName = keyTop->GetName();
        TDirectory *topDirYield = (TDirectory*)keyTop->ReadObj();
        
        std::cout << "Processing Yield Directory: " << topDirName << std::endl;

        // --- 3. Find Matching Efficiency Directory ---
        // Mapping: dimuon -> dimuon-mc, dimuon_xxx -> dimuon-mc_xxx
        TString topDirNameEff = topDirName;
        if (topDirName == "dimuon") {
            topDirNameEff = "dimuon-mc";
        } else if (topDirName.BeginsWith("dimuon_")) {
            topDirNameEff = topDirName;
            topDirNameEff.ReplaceAll("dimuon_", "dimuon-mc_");
        } else {
            // Default fallback if naming convention is different
            topDirNameEff = topDirName + "-mc";
        }

        TDirectory *topDirEff = (TDirectory*)input_efficiency->Get(topDirNameEff);
        if (!topDirEff) {
            // Fallback: try to find any directory that contains topDirName as a prefix
            TIter nextEff(input_efficiency->GetListOfKeys());
            TKey *keyEff;
            while ((keyEff = (TKey*)nextEff())) {
                if (TString(keyEff->GetClassName()) == "TDirectoryFile" && TString(keyEff->GetName()).Contains(topDirName)) {
                    topDirEff = (TDirectory*)keyEff->ReadObj();
                    topDirNameEff = keyEff->GetName();
                    break;
                }
            }
        }

        if (!topDirEff) {
            std::cerr << "Warning: Matching efficiency directory for " << topDirName << " not found. Skipping." << std::endl;
            continue;
        }
        std::cout << "  Using efficiency from: " << topDirNameEff << std::endl;

        TDirectory *outTopDir = output_file->mkdir(topDirName);
        outTopDir->cd();

        // --- 4. Get Yields and Efficiencies ---
        TGraphErrors *g_omega_yield = (TGraphErrors*)topDirYield->Get("omega_yield_ptspectrum");
        TGraphErrors *g_phi_yield = (TGraphErrors*)topDirYield->Get("phi_yield_ptspectrum");
        
        TString pathOmegaEff = Form("%s/Omega/AccEff_Omega_Pt", fitType.Data());
        TString pathPhiEff = Form("%s/Phi/AccEff_Phi_Pt", fitType.Data());
        TH1D *Omega_eff_pt_orig = (TH1D*)topDirEff->Get(pathOmegaEff);
        TH1D *Phi_eff_pt_orig = (TH1D*)topDirEff->Get(pathPhiEff);

        if (!g_omega_yield || !g_phi_yield || !Omega_eff_pt_orig || !Phi_eff_pt_orig) {
            std::cerr << "  Warning: Missing yields or efficiency histograms in " << topDirName << ". Skipping." << std::endl;
            continue;
        }

        Omega_eff_pt_orig->Write("Omega_eff_pt_original");
        Phi_eff_pt_orig->Write("Phi_eff_pt_original");

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

            double eff_sl_o = hEffPt_Omega_SL->GetBinContent(hEffPt_Omega_SL->FindBin(pt));
            double eff_sl_p = hEffPt_Phi_SL->GetBinContent(hEffPt_Phi_SL->FindBin(pt));

            o_eff_vec.push_back(eff_o * eff_sl_o);
            o_eff_err_vec.push_back(eff_o_err * eff_sl_o);
            p_eff_vec.push_back(eff_p * eff_sl_p);
            p_eff_err_vec.push_back(eff_p_err * eff_sl_p);

            double pt_width = pt_err * 2.0;

            if (eff_o > 0 && eff_sl_o > 0 && pt_width > 0) {
                o_xs.push_back(y_o / L_int_mb_inv / BR_omega / (eff_o * eff_sl_o) / pt_width);
                o_xs_err.push_back(y_o_err / L_int_mb_inv / BR_omega / (eff_o * eff_sl_o) / pt_width);
            } else { o_xs.push_back(0); o_xs_err.push_back(0); }

            if (eff_p > 0 && eff_sl_p > 0 && pt_width > 0) {
                p_xs.push_back(y_p / L_int_mb_inv / BR_phi / (eff_p * eff_sl_p) / pt_width);
                p_xs_err.push_back(y_p_err / L_int_mb_inv / BR_phi / (eff_p * eff_sl_p) / pt_width);
            } else { p_xs.push_back(0); p_xs_err.push_back(0); }
            
            final_pt_vec.push_back(pt);
            final_pt_err_vec.push_back(pt_err);
        }

        TGraphErrors *g_omega_eff = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), o_eff_vec.data(), final_pt_err_vec.data(), o_eff_err_vec.data());
        TGraphErrors *g_phi_eff = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), p_eff_vec.data(), final_pt_err_vec.data(), p_eff_err_vec.data());
        TGraphErrors *g_omega_xs = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), o_xs.data(), final_pt_err_vec.data(), o_xs_err.data());
        TGraphErrors *g_phi_xs = new TGraphErrors(final_pt_vec.size(), final_pt_vec.data(), p_xs.data(), final_pt_err_vec.data(), p_xs_err.data());
        
        g_omega_eff->Write("Omega_efficiency_rebinned");
        g_phi_eff->Write("Phi_efficiency_rebinned");
        g_omega_xs->Write("Omega_cross_section");
        g_phi_xs->Write("Phi_cross_section");

        TCanvas *c_Omega_XS = new TCanvas("c_Omega_XS", "Omega Cross Section", 800, 600);
        gPad->SetLogy(); g_omega_xs->SetTitle("Omega Differential Cross Section;p_{T} (GeV/c);#frac{Yield}{L_{int} #cdot BR #cdot (Acc#timesEff)} (mb/GeV/c)");
        g_omega_xs->Draw("APE"); c_Omega_XS->Write();

        TCanvas *c_Phi_XS = new TCanvas("c_Phi_XS", "Phi Cross Section", 800, 600);
        gPad->SetLogy(); g_phi_xs->SetTitle("Phi Differential Cross Section;p_{T} (GeV/c);#frac{Yield}{L_{int} #cdot BR #cdot (Acc#timesEff)} (mb/GeV/c)");
        g_phi_xs->Draw("APE"); c_Phi_XS->Write();
        
        delete g_omega_eff; delete g_phi_eff; delete g_omega_xs; delete g_phi_xs;
        delete c_Omega_XS; delete c_Phi_XS;
    }

    output_file->cd();
    hZvtx_corrected->Write();
    output_file->Close();
    std::cout << "Correction and cross section calculation done. Results saved to 'efficiency_corrected_" << fitType << ".root'" << std::endl;
}
