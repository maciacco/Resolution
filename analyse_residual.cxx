#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <iostream>

int main(const int argc, const char* argv[]) {

  gROOT->SetBatch();

  TFile *f_dat = TFile::Open("foo_dat.root");
  TFile *fo = TFile::Open("residuals.root", "recreate");

  if (!f_dat) {
    std::cout << "Files not found!" << std::endl;
    return 1;
  }

  TH1D* hBias = new TH1D("hBias_0", ";1/#it{p}_{T} (#it{c}/GeV);#delta_{1/#it{p}_{T}} (#it{c}/GeV)", 40, 0., 4.);
  TH1D* hReso = new TH1D("hReso_0", ";1/#it{p}_{T} (#it{c}/GeV);#sigma_{1/#it{p}_{T}} (#it{c}/GeV)", 40, 0., 4.);

  for (int iB{4}; iB < 41; ++iB){
    std::cout << "Processing bin n. " << iB << std::endl;

    TH2D *_reso = static_cast<TH2D*>(f_dat->Get("hPtInvReso"));
    if (!_reso) {
      std::cout << "Residual histogram not found!" << std::endl;
      return 1;
    }

    TH1D* _res = static_cast<TH1D*>(_reso->ProjectionY("_reso_1d", iB, iB));

    hBias->SetBinContent(iB, _res->GetMean());
    hBias->SetBinError(iB, 0.);
    hReso->SetBinContent(iB, _res->GetStdDev());
    hReso->SetBinError(iB, 0.);
  }

  hBias->Write();
  hReso->Write();
  fo->Close();

  return 0;
}
