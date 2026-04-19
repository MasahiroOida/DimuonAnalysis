#ifndef MYANALYSIS
#define MYANALYSIS

void Making1DfromHnSparse(TFile *inputAnalysisResult_EM);
void Making1DfromHnSparse_mc(TFile *inputAnalysisResult_EM_mc);
void Check_mc_mass();
void CalAcceptance_Efficiency();
void LikeSignMethod(TFile *input_SeparetePt);
void PeakFit_Gauss(TFile *input_Sig);
void PeakFit_CrystalBall(TFile *input_Sig);
void PeakFit_CrystalBall_pol4(TFile *input_Sig);
void YieldCalcuration_Gauss(TFile *input_FitResult);
void YieldCalcuration_CrystalBall(TFile *input_FitResult);
void YieldCalcuration_CrystalBall_pol4(TFile *input_FitResult);
void Significance_CrystalBall(TFile *input_FitResult);
void Significance_CrystalBall_pol4(TFile *input_FitResult);
void Significance_GlobalChi2_pol4(TFile *input_FitResult);
void Significance_CrystalBall_pol4_minpt(TFile *input_FitResult);
void Significance_MFTnCluster(TFile *input_FitResult);
void Significance_MCHnCluster(TFile *input_FitResult);
void _06_dR_MFTnCluster_Scan_2D();
void _06_dR_MCHnCluster_Scan_2D();
void _06_dR_Chi2_Scan_2D_pol4();
void _06_dR_minpt_Scan_2D_pol4();
void correction(TFile *input_Lumi,TFile *input_Yield, TFile *input_efficiency);

#endif // MYANALYSIS
