import ROOT

# ROOT drawing settings
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)

# Open processed data
inFile = ROOT.TFile("Processed_Histograms.root", "READ")

# Go to specific pT bin directory and get objects
directory = inFile.Get("Pt_0.0_1.0")
histSig = directory.Get("LikeSignSig_Pt_0.0_1.0")
totalFit = directory.Get("TotalFit_Pt_0.0_1.0")
bgFit = directory.Get("BGFit_Pt_0.0_1.0")
omegaFit = directory.Get("OmegaFit_Pt_0.0_1.0")
phiFit = directory.Get("PhiFit_Pt_0.0_1.0")

# Create canvas and adjust appearance
canvas = ROOT.TCanvas("c1", "Fit Result", 800, 600)
canvas.SetGrid()
canvas.SetLeftMargin(0.15)

# Histogram settings
histSig.SetLineColor(ROOT.kBlack)
histSig.SetMarkerStyle(8)
histSig.SetMarkerSize(0.7)
histSig.GetXaxis().SetTitle("Dimuon invariant mass GeV/c^{2}")
histSig.GetYaxis().SetTitle("dN/dm")
histSig.GetXaxis().SetRangeUser(0.5, 1.3)

# Fit function settings
totalFit.SetLineColor(ROOT.kRed)
bgFit.SetLineColor(ROOT.kBlack)
bgFit.SetLineStyle(2)
omegaFit.SetLineColor(ROOT.kGreen)
phiFit.SetLineColor(ROOT.kBlue)

# Draw
histSig.Draw("PE")
totalFit.Draw("SAME")
bgFit.Draw("SAME")
omegaFit.Draw("SAME")
phiFit.Draw("SAME")

# Add legend
leg = ROOT.TLegend(0.65, 0.65, 0.85, 0.85)
leg.SetBorderSize(0)
leg.AddEntry(histSig, "Signal", "lep")
leg.AddEntry(totalFit, "Total Fit", "l")
leg.AddEntry(omegaFit, "#omega peak", "l")
leg.AddEntry(phiFit, "#phi peak", "l")
leg.Draw()

# Save as PDF
canvas.SaveAs("FitResult_Pt_0.0_1.0.pdf")
inFile.Close()
