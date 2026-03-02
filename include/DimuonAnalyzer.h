#ifndef DIMUON_ANALYZER_H
#define DIMUON_ANALYZER_H

#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <tuple>
#include <string>

class DimuonAnalyzer {
public:
    DimuonAnalyzer() = default;
    ~DimuonAnalyzer() = default;

    // Create 1D mass distribution from 2D histogram for a specific pT range
    TH1D* ProjectTo1D(TH2* hist2D, double ptMin, double ptMax, const std::string& name);

    // Calculate background and signal using the Like-Sign method
    // Returns: background and signal histograms
    std::tuple<TH1D*, TH1D*> CalculateLikeSign(
        TH1D* SEPM, TH1D* SEPP, TH1D* SEMM,
        TH1D* MEPM = nullptr, TH1D* MEPP = nullptr, TH1D* MEMM = nullptr,
        const std::string& suffix = "");

    // Perform peak fitting
    // Returns: total fit, background fit, omega fit, and phi fit
    std::tuple<TF1*, TF1*, TF1*, TF1*> PerformPeakFit(TH1D* sigHist, const std::string& suffix = "");

private:
    // Discontinuous background function
    static Double_t StepBGFunction(Double_t *x, Double_t *par);
};

#endif
