#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>

void analyse_chi2(const char* ifname = "fits_full", const char* ofname = "bias_reso_vs_pt"){
  TFile* infile = TFile::Open(Form("%s.root", ifname));

  TH1D* hBias = new TH1D("hBias", ";1/#it{p}_{T} (#it{c}/GeV);#delta_{1/#it{p}_{T}} (#it{c}/GeV)", 40, 0., 4.);
  TH1D* hReso = new TH1D("hReso", ";1/#it{p}_{T} (#it{c}/GeV);#sigma_{1/#it{p}_{T}} (#it{c}/GeV)", 40, 0., 4.);
  for (int iB{4}; iB < 41; ++iB) {
    TH2D* chi2_prof = static_cast<TH2D*>(infile->Get(Form("chi2_prof_%d", iB)));

    int bin_x_min = 0, bin_y_min = 0, bin_z_min = 0;
    chi2_prof->GetMinimumBin(bin_x_min, bin_y_min, bin_z_min);
    double x_min = chi2_prof->GetXaxis()->GetBinCenter(bin_x_min);
    double y_min = chi2_prof->GetYaxis()->GetBinCenter(bin_y_min);

    hBias->SetBinContent(iB, x_min);
    hBias->SetBinError(iB, chi2_prof->GetXaxis()->GetBinWidth(1));
    hReso->SetBinContent(iB, y_min);
    hReso->SetBinError(iB, chi2_prof->GetYaxis()->GetBinWidth(1));
  }

  TFile* fo = TFile::Open(Form("%s.root", ofname), "recreate");
  fo->cd();
  hBias->Write();
  hReso->Write();
  fo->Close();
}
