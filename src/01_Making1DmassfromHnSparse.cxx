#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TObject.h>
#include <TStyle.h>
#include <TLegend.h>
#include <cmath>
#include "TMath.h"
#include <iostream>
#include <TStopwatch.h>
#include <sstream>
#include "myAnalysis.h"
#include <THashList.h>
#include <TParameter.h>
#include <THnSparse.h>
#include <TObjString.h>

void Making2DmassfromHnSparse(TFile *inputROOT, TDirectory *outfileDir, TString dirName) {
    outfileDir->cd();
    std::vector<TString> signs = {"uls", "lspp", "lsmm"};

    for (int i = 0; i < signs.size(); ++i) {
        TString sign = signs[i];
        
        TString path_same = Form("%s/Pair/same/%s/hs", dirName.Data(), sign.Data());
        TString path_mix = Form("%s/Pair/mix/%s/hs", dirName.Data(), sign.Data());
        
        THnSparse* hn_same = (THnSparse*)inputROOT->Get(path_same);
        THnSparse* hn_mix = (THnSparse*)inputROOT->Get(path_mix);

        if (!hn_same || !hn_mix) {
            std::cout << "Warning: Could find HnSparse at " << dirName << " for sign " << sign << std::endl;
            continue;
        }
        
        // 誤差計算を明示的に有効化
        hn_same->CalculateErrors();
        hn_mix->CalculateErrors();

        int xAxisDim = 0; // mass
        int yAxisDim = 1; // pT
        
        TH2D* Mass_pt_same = (TH2D*)hn_same->Projection(yAxisDim, xAxisDim);
        TH2D* Mass_pt_mix = (TH2D*)hn_mix->Projection(yAxisDim, xAxisDim);
        
        Mass_pt_same->Sumw2();
        Mass_pt_mix->Sumw2();
        
        Mass_pt_same->SetName(Form("Mass_pt_same_%s", sign.Data()));
        Mass_pt_same->Write(Form("Mass_pt_same_%s",sign.Data()));
        Mass_pt_mix->SetName(Form("Mass_pt_mix_%s", sign.Data()));
        Mass_pt_mix->Write(Form("Mass_pt_mix_%s",sign.Data()));
        delete Mass_pt_same;
        delete Mass_pt_mix;
    }
}

void Projection(double Downpt, double Uppt, TH2 *SEPM2D, TH2 *SEPP2D, TH2 *SEMM2D, TH2 *MEPM2D, TH2 *MEPP2D, TH2 *MEMM2D)
{
    TAxis *yAxis = SEPM2D->GetYaxis();
    int Downbin = yAxis->FindBin(Downpt + 0.0001);
    int Upbin = yAxis->FindBin(Uppt - 0.0001);

    auto manualProjection = [&](TH2* h2, const char* name) -> TH1D* {
        if (!h2) return nullptr;
        // 構造をコピーするために一度ProjectionXを呼ぶが、中身は手動で詰め直す
        TH1D* h1 = h2->ProjectionX(Form("%s_temp", name), Downbin, Upbin);
        h1->SetName(name);
        h1->SetTitle(name);
        h1->Reset(); // 内容と誤差をリセット
        h1->Sumw2();

        for (int ix = 0; ix <= h2->GetNbinsX() + 1; ++ix) {
            double sumContent = 0;
            double sumError2 = 0;
            for (int iy = Downbin; iy <= Upbin; ++iy) {
                sumContent += h2->GetBinContent(ix, iy);
                double err = h2->GetBinError(ix, iy);
                sumError2 += err * err;
            }
            h1->SetBinContent(ix, sumContent);
            h1->SetBinError(ix, TMath::Sqrt(sumError2));
        }
        return h1;
    };

    TH1D *SEPM = manualProjection(SEPM2D, "SEPM");
    TH1D *SEPP = manualProjection(SEPP2D, "SEPP");
    TH1D *SEMM = manualProjection(SEMM2D, "SEMM");

    if (SEPM) SEPM->Write();
    if (SEPP) SEPP->Write();
    if (SEMM) SEMM->Write();

    if (MEPM2D != nullptr)
    {
        TH1D *MEPM = manualProjection(MEPM2D, "MEPM");
        TH1D *MEPP = manualProjection(MEPP2D, "MEPP");
        TH1D *MEMM = manualProjection(MEMM2D, "MEMM");
        if (MEPM) MEPM->Write();
        if (MEPP) MEPP->Write();
        if (MEMM) MEMM->Write();
        delete MEPM;
        delete MEPP;
        delete MEMM;
    }
    delete SEPM;
    delete SEPP;
    delete SEMM;
}

void make1Dfrom2Dhist(TDirectory *inputDir, TDirectory *outfile_projection, double nEvents = -1){
    TH2D* SEPM2D = (TH2D*)inputDir->Get("Mass_pt_same_uls");
    TH2D* SEPP2D = (TH2D*)inputDir->Get("Mass_pt_same_lspp");
    TH2D* SEMM2D = (TH2D*)inputDir->Get("Mass_pt_same_lsmm");
    TH2D* MEPM2D = (TH2D*)inputDir->Get("Mass_pt_mix_uls");
    TH2D* MEPP2D = (TH2D*)inputDir->Get("Mass_pt_mix_lspp");
    TH2D* MEMM2D = (TH2D*)inputDir->Get("Mass_pt_mix_lsmm");

    if (!SEPM2D) {
        std::cout << "Warning: SEPM2D (Mass_pt_same_uls) not found in directory " << inputDir->GetName() << std::endl;
        return;
    }

    TAxis *ptAxis = SEPM2D->GetYaxis();
    double ptMin_all = ptAxis->GetXmin();
    double ptMax_all = ptAxis->GetXmax();

    // All pt
    TString allPtDirName = Form("Pt_%.1fto%.1f", ptMin_all, ptMax_all);
    TDirectory *allpt_Dir = outfile_projection->mkdir(allPtDirName);
    allpt_Dir->cd();
    (new TParameter<double>("ptrange", ptMax_all - ptMin_all))->Write();
    (new TParameter<double>("ptmin", ptMin_all))->Write();
    (new TParameter<double>("ptmax", ptMax_all))->Write();
    if (nEvents > 0) (new TParameter<double>("nEvents", nEvents))->Write();
    Projection(ptMin_all, ptMax_all, SEPM2D, SEPP2D, SEMM2D, MEPM2D, MEPP2D, MEMM2D);    

    std::vector<double> pt_edges = {0, 0.5, 0.75, 1, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10};

    for (size_t i = 0; i < pt_edges.size() - 1; ++i)
    {
        double pt_min = pt_edges[i];
        double pt_max = pt_edges[i+1];
        double ptrange = pt_max - pt_min;
        TString dirName = Form("Pt_%.2fto%.2f", pt_min, pt_max);
        TDirectory *newDir = outfile_projection->mkdir(dirName);
        newDir->cd();
        Projection(pt_min, pt_max, SEPM2D, SEPP2D, SEMM2D, MEPM2D, MEPP2D, MEMM2D);
        (new TParameter<double>("ptrange", ptrange))->Write();
        (new TParameter<double>("ptmin", pt_min))->Write();
        (new TParameter<double>("ptmax", pt_max))->Write();
        if (nEvents > 0) (new TParameter<double>("nEvents", nEvents))->Write();
    }

    delete SEPM2D;
    delete SEPP2D;
    delete SEMM2D;
    delete MEPM2D;
    delete MEPP2D;
    delete MEMM2D;
}

void Making1DfromHnSparse(TFile *inputAnalysisResult_EM){
    TFile *outfile_2D = new TFile("Made2Dhist.root","RECREATE");
    TFile *outfile_1D = new TFile("Made1Dhist.root","RECREATE");
    
    // Extract event count (from dimuon/Event/after/hZVtx entries as requested)
    double nEvents = -1;
    TH1D *hZVtx = (TH1D *)inputAnalysisResult_EM->Get("dimuon/Event/after/hZvtx");
    if (hZVtx) {
        nEvents = hZVtx->GetEntries();
        std::cout << "Extracted nEvents (from hZVtx): " << nEvents << std::endl;
    } else {
        TH1D *hTVX = (TH1D *)inputAnalysisResult_EM->Get("dimuon/BC/hTVXCounter");
        if (hTVX) {
            nEvents = hTVX->GetBinContent(6);
            std::cout << "Extracted nEvents (from hTVX bin 6): " << nEvents << std::endl;
        } else {
            std::cerr << "Warning: Neither hZVtx nor hTVXCounter found. nEvents will not be saved." << std::endl;
        }
    }

    TIter next(inputAnalysisResult_EM->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next())) {
        if (TString(key->GetClassName()) == "TDirectoryFile") {
            TString dirName = key->GetName();
            TDirectory *dir = (TDirectory*)inputAnalysisResult_EM->Get(dirName);
            if (dir->Get("Pair")) {
                std::cout << "Processing directory: " << dirName << std::endl;
                TDirectory *dir2D = outfile_2D->mkdir(dirName);
                TDirectory *dir1D = outfile_1D->mkdir(dirName);
                
                Making2DmassfromHnSparse(inputAnalysisResult_EM, dir2D, dirName);
                make1Dfrom2Dhist(dir2D, dir1D, nEvents);
            }
        }
    }
    
    outfile_2D->Close();
    outfile_1D->Close();
    delete outfile_2D;
    delete outfile_1D;
}

void Making2DmassfromHnSparse_mc(TDirectory *inputROOT, TDirectory *outfileDir, TString dirName, TDirectory *outfile1DDir = nullptr) {
    outfileDir->cd();
    
    bool isGenerated = dirName.Contains("Generated");
    std::vector<TString> signs;
    if (isGenerated) {
        signs = {"uls"}; // Generatedの場合は1回だけ実行（名前は互換性のためにulsとする）
    } else {
        signs = {"uls", "lspp", "lsmm"};
    }

    for (int i = 0; i < signs.size(); ++i) {
        TString sign = signs[i];
        TString path;
        if (isGenerated) {
            path = dirName; // Generatedの場合は直下のディレクトリを見る
        } else {
            path = Form("%s/%s", dirName.Data(), sign.Data());
        }
        
        TObject* obj = inputROOT->Get(path);

        if (!obj) {
            std::cout << "Warning: Could not find object at " << path << std::endl;
            continue;
        }

        THnSparse* hn = nullptr;
        if (obj->InheritsFrom(THnSparse::Class())) {
            hn = (THnSparse*)obj;
        } else if (obj->InheritsFrom(TDirectory::Class())) {
            TDirectory* subdir = (TDirectory*)obj;
            // "hs" という名前のオブジェクトを直接探す
            TObject* hsObj = subdir->Get("hs");
            if (hsObj && hsObj->InheritsFrom(THnSparse::Class())) {
                hn = (THnSparse*)hsObj;
            } else {
                // 見つからない場合はリストから探す
                TIter next(subdir->GetListOfKeys());
                TKey* key;
                while ((key = (TKey*)next())) {
                    TObject* subobj = subdir->Get(key->GetName());
                    if (subobj && subobj->InheritsFrom(THnSparse::Class())) {
                        hn = (THnSparse*)subobj;
                        break;
                    }
                }
            }
        }

        if (!hn) {
            std::cout << "Warning: No THnSparse found at " << path << std::endl;
            continue;
        }

        hn->CalculateErrors(); // 誤差計算を有効化

        int xAxisDim = 0; // mass
        int yAxisDim = 1; // pT
        
        // 軸の数を確認してクラッシュを防止
        if (hn->GetNdimensions() < 2) {
            std::cout << "Error: THnSparse at " << path << " has less than 2 dimensions!" << std::endl;
            continue;
        }

        outfileDir->cd();
        TH2D* Mass_pt = (TH2D*)hn->Projection(yAxisDim, xAxisDim);
        if (!Mass_pt) {
            std::cout << "Error: Projection failed for " << path << std::endl;
            continue;
        }
        Mass_pt->Sumw2();
        
        Mass_pt->SetName(Form("Mass_pt_same_%s", sign.Data()));
        Mass_pt->Write();

        // SEPM (uls) の Pt と Y を追加で出力
        if (sign == "uls" && outfile1DDir) {
            outfile1DDir->cd();
            TH1D* hPt = (TH1D*)hn->Projection(1);
            hPt->SetName("Pt_uls");
            hPt->SetTitle("p_{T} distribution (uls)");
            hPt->Write();
            delete hPt;

            if (hn->GetNdimensions() > 2) {
                TH1D* hY = (TH1D*)hn->Projection(2);
                hY->SetName("Y_uls");
                hY->SetTitle("y distribution (uls)");
                hY->Write();
                delete hY;

                // 2D Pt-Y projection
                TH2D* hPtY = (TH2D*)hn->Projection(2, 1); // yDim=2 (Y), xDim=1 (Pt)
                hPtY->SetName("PtY_uls");
                hPtY->SetTitle("p_{T} vs y (uls)");
                hPtY->Write();
                delete hPtY;
            }
        }

        delete Mass_pt;
    }
}

void make1Dfrom2Dhist_mc(TDirectory *inputDir, TDirectory *outfile_projection, THnSparse *hn, double nEvents = -1) {
    if (!hn) return;

    int ptDim = 1;
    int massDim = 0;
    int yDim = 2;

    // Y bins
    std::vector<std::pair<double, double>> y_bins;
    for (double y = -4.0; y < -2.5; y += 0.25) {
        y_bins.push_back({y, y + 0.25});
    }

    TAxis *ptAxis = hn->GetAxis(ptDim);

    // All pt (0-10 GeV)
    {
        double pt_min = 0.0;
        double pt_max = 10.0;
        TString dirName = Form("Pt_%.2fto%.2f", pt_min, pt_max);
        TDirectory *newDir = outfile_projection->mkdir(dirName);
        newDir->cd();

        int binMin = ptAxis->FindBin(pt_min + 1e-6);
        int binMax = ptAxis->FindBin(pt_max - 1e-6);
        ptAxis->SetRange(binMin, binMax);

        TH1D *hMass = (TH1D*)hn->Projection(massDim);
        hMass->SetName("SEPM");
        hMass->SetTitle(Form("Mass distribution (%.2f < p_{T} < %.2f)", pt_min, pt_max));
        hMass->Write();
        delete hMass;

        TH1D *hPt = (TH1D*)hn->Projection(ptDim);
        hPt->SetName("Pt");
        hPt->Write();
        delete hPt;

        (new TParameter<double>("ptmin", pt_min))->Write();
        (new TParameter<double>("ptmax", pt_max))->Write();
        if (nEvents > 0) (new TParameter<double>("nEvents", nEvents))->Write();
        
        ptAxis->SetRange(0, 0);
    }

    std::vector<double> pt_edges = {0, 0.5, 0.75, 1, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10};

    // Process Pt bins
    for (size_t i = 0; i < pt_edges.size() - 1; ++i) {
        double pt_min = pt_edges[i];
        double pt_max = pt_edges[i+1];
        TString dirName = Form("Pt_%.2fto%.2f", pt_min, pt_max);
        TDirectory *newDir = outfile_projection->mkdir(dirName);
        newDir->cd();

        int binMin = ptAxis->FindBin(pt_min + 1e-6);
        int binMax = ptAxis->FindBin(pt_max - 1e-6);
        ptAxis->SetRange(binMin, binMax);

        TH1D *hMass = (TH1D*)hn->Projection(massDim);
        hMass->SetName("SEPM");
        hMass->SetTitle(Form("Mass distribution (%.2f < p_{T} < %.2f)", pt_min, pt_max));
        hMass->Write();
        delete hMass;

        TH1D *hPt = (TH1D*)hn->Projection(ptDim);
        hPt->SetName("Pt");
        hPt->Write();
        delete hPt;

        (new TParameter<double>("ptmin", pt_min))->Write();
        (new TParameter<double>("ptmax", pt_max))->Write();
        if (nEvents > 0) (new TParameter<double>("nEvents", nEvents))->Write();
        
        ptAxis->SetRange(0, 0);
    }

    // Process Y bins
    for (const auto& bin : y_bins) {
        double y_min = bin.first;
        double y_max = bin.second;
        TString dirName = Form("Y_%.2fto%.2f", y_min, y_max);
        TDirectory *newDir = outfile_projection->mkdir(dirName);
        newDir->cd();

        TAxis *yAxis = hn->GetAxis(yDim);
        int binMin = yAxis->FindBin(y_min + 1e-6);
        int binMax = yAxis->FindBin(y_max - 1e-6);
        yAxis->SetRange(binMin, binMax);

        TH1D *hMass = (TH1D*)hn->Projection(massDim);
        hMass->SetName("SEPM");
        hMass->SetTitle(Form("Mass distribution (%.2f < y < %.2f)", y_min, y_max));
        hMass->Write();
        delete hMass;

        TH1D *hY = (TH1D*)hn->Projection(yDim);
        hY->SetName("Y");
        hY->Write();
        delete hY;

        (new TParameter<double>("ymin", y_min))->Write();
        (new TParameter<double>("ymax", y_max))->Write();
        if (nEvents > 0) (new TParameter<double>("nEvents", nEvents))->Write();
        
        yAxis->SetRange(0, 0);
    }
}

void Making1DfromHnSparse_mc(TFile *inputAnalysisResult_EM_mc){
    TFile *outfile_2D = new TFile("Made2Dhist_mc.root","RECREATE");
    TFile *outfile_1D = new TFile("Made1Dhist_mc.root","RECREATE");
    
    TIter nextKey(inputAnalysisResult_EM_mc->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextKey())) {
        TString topDirName = key->GetName();
        if (!topDirName.BeginsWith("dimuon-mc")) continue;
        
        TDirectory *topDirInput = (TDirectory*)key->ReadObj();
        std::cout << "Processing MC Top Directory: " << topDirName << std::endl;

        // Extract event count from MC input (Event/after/hZVtx entries)
        double nEvents = -1;
        TH1D *hZVtx = (TH1D *)topDirInput->Get("Event/after/hZvtx");
        if (hZVtx) {
            nEvents = hZVtx->GetEntries();
            std::cout << "  Extracted nEvents (from hZVtx): " << nEvents << std::endl;
        } else {
            TH1D *hTVX = (TH1D *)topDirInput->Get("BC/hTVXCounter");
            if (hTVX) {
                nEvents = hTVX->GetBinContent(6);
                std::cout << "  Extracted nEvents (from hTVX bin 6): " << nEvents << std::endl;
            }
        }

        std::vector<TString> particles = {"Omega", "Phi"};
        std::vector<TString> categories = {"All", "Acc", "Reco"};

        for (const auto& part : particles) {
            for (const auto& cat : categories) {
                TString path;
                TString outNameBase;
                if (cat == "All") {
                    path = Form("Generated/VM/All/%s", part.Data());
                    outNameBase = Form("Generated_All/%s", part.Data());
                } else if (cat == "Acc") {
                    path = Form("Generated/VM/Acc/%s", part.Data());
                    outNameBase = Form("Generated_Acc/%s", part.Data());
                } else if (cat == "Reco") {
                    path = Form("Pair/sm/%s2ll", part.Data());
                    outNameBase = Form("Reconstructed/%s", part.Data());
                }

                TString outName = topDirName + "/" + outNameBase;
                TDirectory *dirCheck = (TDirectory*)topDirInput->Get(path);
                if (!dirCheck) continue;

                std::cout << "  Processing " << cat << " for " << part << " at " << path << std::endl;
                
                // Create recursive directories in output
                TString currentPath = "";
                TObjArray *parts = outName.Tokenize("/");
                TDirectory *current2DDir = outfile_2D;
                TDirectory *current1DDir = outfile_1D;
                for (int i = 0; i < parts->GetEntries(); i++) {
                    TString sub = ((TObjString*)parts->At(i))->GetString();
                    TDirectory *next2D = current2DDir->GetDirectory(sub);
                    if (!next2D) next2D = current2DDir->mkdir(sub);
                    current2DDir = next2D;

                    TDirectory *next1D = current1DDir->GetDirectory(sub);
                    if (!next1D) next1D = current1DDir->mkdir(sub);
                    current1DDir = next1D;
                }
                delete parts;

                THnSparse *hn = nullptr;
                TString hsPath = (cat == "Reco") ? path + "/uls/hs" : path + "/hs";
                TObject *obj = topDirInput->Get(hsPath);
                
                if (!obj || !obj->InheritsFrom(THnSparse::Class())) {
                    TDirectory *searchDir = (TDirectory*)topDirInput->Get(path);
                    if (searchDir) {
                        if (cat == "Reco") {
                            TDirectory *ulsDir = (TDirectory*)searchDir->Get("uls");
                            if (ulsDir) searchDir = ulsDir;
                        }
                        TIter next(searchDir->GetListOfKeys());
                        TKey *k;
                        while ((k = (TKey*)next())) {
                            TObject *tmp = searchDir->Get(k->GetName());
                            if (tmp && tmp->InheritsFrom(THnSparse::Class())) {
                                hn = (THnSparse*)tmp;
                                break;
                            }
                        }
                    }
                } else {
                    hn = (THnSparse*)obj;
                }

                Making2DmassfromHnSparse_mc(topDirInput, current2DDir, path, current1DDir);
                if (hn) {
                    make1Dfrom2Dhist_mc(current2DDir, current1DDir, hn, nEvents);
                }
            }
        }
    }
    
    outfile_2D->Close();
    outfile_1D->Close();
    delete outfile_2D;
    delete outfile_1D;
}
