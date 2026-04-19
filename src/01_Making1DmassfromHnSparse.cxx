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
            std::cout << "Warning: Could not find HnSparse at " << dirName << " for sign " << sign << std::endl;
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

void make1Dfrom2Dhist(TDirectory *inputDir, TDirectory *outfile_projection){
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

    // All pt
    TDirectory *allpt_Dir = outfile_projection->mkdir("Pt_0to10");
    allpt_Dir->cd();
    (new TParameter<double>("ptrange", 10))->Write();
    (new TParameter<double>("ptmin", 0))->Write();
    (new TParameter<double>("ptmax", 10))->Write();
    Projection(0, 10, SEPM2D, SEPP2D, SEMM2D, MEPM2D, MEPP2D, MEMM2D);    

    std::vector<std::pair<double, double>> pt_bins;
    for (double pt = 0.0; pt < 10.0; pt += 0.5) {
        pt_bins.push_back({pt, pt + 0.5});
    }

    for (const auto& bin : pt_bins)
    {
        double pt_min = bin.first;
        double pt_max = bin.second;
        double ptrange = pt_max - pt_min;
        TString dirName = Form("Pt_%.1fto%.1f", pt_min, pt_max);
        TDirectory *newDir = outfile_projection->mkdir(dirName);
        newDir->cd();
        Projection(pt_min, pt_max, SEPM2D, SEPP2D, SEMM2D, MEPM2D, MEPP2D, MEMM2D);
        (new TParameter<double>("ptrange", ptrange))->Write();
        (new TParameter<double>("ptmin", pt_min))->Write();
        (new TParameter<double>("ptmax", pt_max))->Write();
    }

    delete SEPM2D;
    delete SEPP2D;
    delete SEMM2D;
    delete MEPM2D;
    delete MEPP2D;
    delete MEMM2D;
}

void Making2DmassfromHnSparse_mc(TFile *inputROOT, TDirectory *outfileDir, TString dirName, TDirectory *outfile1DDir = nullptr) {
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

void Making1DfromHnSparse(TFile *inputAnalysisResult_EM){
    TFile *outfile_2D = new TFile("Made2Dhist.root","RECREATE");
    TFile *outfile_1D = new TFile("Made1Dhist.root","RECREATE");
    
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
                make1Dfrom2Dhist(dir2D, dir1D);
            }
        }
    }
    
    outfile_2D->Close();
    outfile_1D->Close();
    delete outfile_2D;
    delete outfile_1D;
}

void make1Dfrom2Dhist_mc(TDirectory *inputDir, TDirectory *outfile_projection, THnSparse *hn) {
    if (!hn) return;

    // カスタム pT ビン: 0~0.5, 0.5~1.0, 1.0~2.0, 2.0~3.0, 3.0~6.0, 6.0~10.0
    std::vector<std::pair<double, double>> pt_bins = {
        {0.0, 0.5}, {0.5, 1.0}, {1.0, 2.0}, {2.0, 3.0}, {3.0, 6.0}, {6.0, 10.0}
    };

    for (const auto& bin : pt_bins) {
        double pt_min = bin.first;
        double pt_max = bin.second;
        TString dirName = Form("Pt_%.1fto%.1f", pt_min, pt_max);
        TDirectory *newDir = outfile_projection->mkdir(dirName);
        newDir->cd();

        // pT 軸 (dim 1) で範囲を選択
        int ptDim = 1;
        int massDim = 0;
        int yDim = 2;

        TAxis *ptAxis = hn->GetAxis(ptDim);
        int binMin = ptAxis->FindBin(pt_min + 1e-6);
        int binMax = ptAxis->FindBin(pt_max - 1e-6);
        ptAxis->SetRange(binMin, binMax);

        // Mass 射影
        TH1D *hMass = (TH1D*)hn->Projection(massDim);
        hMass->SetName("SEPM");
        hMass->SetTitle(Form("Mass distribution (%.1f < p_{T} < %.1f)", pt_min, pt_max));
        hMass->Write();
        delete hMass;

        // pT 射影
        TH1D *hPt = (TH1D*)hn->Projection(ptDim);
        hPt->SetName("Pt");
        hPt->SetTitle(Form("p_{T} distribution (%.1f < p_{T} < %.1f)", pt_min, pt_max));
        hPt->Write();
        delete hPt;

        // パラメータを保存
        (new TParameter<double>("ptmin", pt_min))->Write();
        (new TParameter<double>("ptmax", pt_max))->Write();
        (new TParameter<double>("ptrange", pt_max - pt_min))->Write();

        // y 射影 (次元が存在する場合)
        if (hn->GetNdimensions() > yDim) {
            TH1D *hY = (TH1D*)hn->Projection(yDim);
            hY->SetName("Y");
            hY->SetTitle(Form("y distribution (%.1f < p_{T} < %.1f)", pt_min, pt_max));
            hY->Write();
            delete hY;
        }

        // 範囲をリセット
        ptAxis->SetRange(0, 0);
    }
}

void Making1DfromHnSparse_mc(TFile *inputAnalysisResult_EM_mc){
    TFile *outfile_2D = new TFile("Made2Dhist_mc.root","RECREATE");
    TFile *outfile_1D = new TFile("Made1Dhist_mc.root","RECREATE");
    
    std::vector<TString> particles = {"Omega", "Phi"};
    std::vector<TString> categories = {"All", "Acc", "Reco"};

    for (const auto& part : particles) {
        for (const auto& cat : categories) {
            TString path;
            TString outName;
            if (cat == "All") {
                path = Form("dimuon-mc/Generated/VM/All/%s", part.Data());
                outName = Form("Generated_All/%s", part.Data());
            } else if (cat == "Acc") {
                path = Form("dimuon-mc/Generated/VM/Acc/%s", part.Data());
                outName = Form("Generated_Acc/%s", part.Data());
            } else if (cat == "Reco") {
                path = Form("dimuon-mc/Pair/sm/%s2ll", part.Data());
                outName = Form("Reconstructed/%s", part.Data());
            }

            TDirectory *dirCheck = (TDirectory*)inputAnalysisResult_EM_mc->Get(path);
            if (!dirCheck) {
                std::cout << "Warning: Directory " << path << " not found!" << std::endl;
                continue;
            }

            std::cout << "Processing MC " << cat << " directory: " << path << std::endl;
            TDirectory *dir2D = outfile_2D->mkdir(outName);
            TDirectory *dir1D = outfile_1D->mkdir(outName);
            
            // THnSparse を取得して直接渡す
            THnSparse *hn = nullptr;
            TString hsPath = (cat == "Reco") ? path + "/uls/hs" : path + "/hs";
            TObject *obj = inputAnalysisResult_EM_mc->Get(hsPath);
            
            if (!obj || !obj->InheritsFrom(THnSparse::Class())) {
                // 見つからない場合はディレクトリ内を再帰的に探索
                TDirectory *searchDir = (TDirectory*)inputAnalysisResult_EM_mc->Get(path);
                if (searchDir) {
                    // Recoの場合はまず "uls" を見る
                    if (cat == "Reco") {
                        TDirectory *ulsDir = (TDirectory*)searchDir->Get("uls");
                        if (ulsDir) searchDir = ulsDir;
                    }
                    
                    TIter next(searchDir->GetListOfKeys());
                    TKey *key;
                    while ((key = (TKey*)next())) {
                        TObject *tmp = searchDir->Get(key->GetName());
                        if (tmp && tmp->InheritsFrom(THnSparse::Class())) {
                            hn = (THnSparse*)tmp;
                            break;
                        }
                    }
                }
            } else {
                hn = (THnSparse*)obj;
            }

            Making2DmassfromHnSparse_mc(inputAnalysisResult_EM_mc, dir2D, path, dir1D);
            if (hn) {
                make1Dfrom2Dhist_mc(dir2D, dir1D, hn);
            } else {
                std::cout << "Warning: No THnSparse found for " << outName << ", skipping custom pT projection." << std::endl;
            }
        }
    }
    
    outfile_2D->Close();
    outfile_1D->Close();
    delete outfile_2D;
    delete outfile_1D;
}
