#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <iostream>

bool relative[]{true, true, true};

Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  Double_t u   = (x-mu)/width;
  Double_t A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  Double_t A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  Double_t B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  Double_t B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  Double_t result{1};
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

Double_t DoubleSidedCB(double* x, double *par)
{
#ifdef _MC_
  return par[0] * TMath::Exp( -0.5 * std::pow((x[0] - par[1]) / par[2], 2.) );
#else
  return par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]);
#endif
}

int main(const int argc, const char* argv[]) {

  gROOT->SetBatch();

  TFile *f_dat = TFile::Open("foo_dat_1.root");
  TFile *fo = TFile::Open("residuals_inv.root", "recreate");

  if (!f_dat) {
    std::cout << "Files not found!" << std::endl;
    return 1;
  }

  const double binEdge[]{0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 5.};
  const int nBins = 21;

  const std::string histoName[]{"hPtReso", "hPhiReso", "hEtaReso"};
  for (int iH{0}; iH < 1; ++iH){
    TH1D* hBias = new TH1D(Form("%s_bias", histoName[iH].data()), ";1/#it{p}_{T} (#it{c}/GeV);#delta_{1/#it{p}_{T}} (#it{c}/GeV)", nBins, binEdge); //40, 0., 4.);
    TH1D* hReso = new TH1D(histoName[iH].data(), ";1/#it{p}_{T} (#it{c}/GeV);#sigma_{1/#it{p}_{T}} (#it{c}/GeV)", nBins, binEdge); //40, 0., 4.);

    TAxis refAxis(nBins, binEdge);
    TH2D *_reso = static_cast<TH2D*>(f_dat->Get(histoName[iH].data()));

    if (!_reso) {
      std::cout << "Residual histogram not found!" << std::endl;
      return 1;
    }

    for (int iB{4}; iB < nBins + 1; ++iB){
      double p_low = refAxis.GetBinLowEdge(iB);
      double p_up = refAxis.GetBinUpEdge(iB);
      int iB_low = _reso->GetXaxis()->FindBin(p_low + 0.005f);
      int iB_up = _reso->GetXaxis()->FindBin(p_up - 0.005f);

      std::cout << "Processing bin n. " << iB << std::endl;

      TH1D* _res = static_cast<TH1D*>(_reso->ProjectionY("_reso_1d", iB_low, iB_up));
      _res->Scale(1. / _res->Integral(1, _res->GetNbinsX()));

      _res->Fit("gaus", "LMQR+", "", _res->GetMean() - 5. * _res->GetStdDev(), _res->GetMean() + 5. * _res->GetStdDev());
      float mu = _res->GetFunction("gaus")->GetParameter(1);
      float sg = _res->GetFunction("gaus")->GetParameter(2);
      std::string fitfname = "fit";
      fitfname = fitfname + std::to_string(iH) + "_" + std::to_string(iB);

//      TF1 fit(fitfname.data(), "gaus", -10., 10.);
      TF1 fit(fitfname.data(), DoubleSidedCB, -10., 10., 7);
      fit.SetParLimits(0, 0., 1.);
      fit.SetParLimits(1, -0.01, 0.01);
      fit.SetParLimits(2, 0.001, 0.02);
      fit.SetParLimits(3, 0.1, 5.);
      fit.SetParLimits(4, 0., 10.);
      fit.SetParLimits(5, 0.1, 5.);
      fit.SetParLimits(6, 0., 10.);
      fit.SetParameter(0, 0.5 * (0. + 1.));
      fit.SetParameter(1, 0.5 * (-0.02 + 0.02));
      fit.SetParameter(2, 0.5 * (0. + 0.02));
      fit.SetParameter(3, 0.5 * (0. + 5.));
      fit.SetParameter(4, 0.5 * (0. + 10.));
      fit.SetParameter(5, 0.5 * (0. + 5.));
      fit.SetParameter(6, 0.5 * (0. + 10.));

      for(int i{0}; i < 1; ++i)_res->Fit(fitfname.data(), "LRMQ+", "", -0.05, 0.05); //mu - 4. * sg, mu + 4. * sg);
      _res->Write();
      float scale = relative[iH] ? 1. / hBias->GetBinCenter(iB) : 1.;

      hBias->SetBinContent(iB,  _res->GetFunction(fitfname.data())->GetParameter(1) * scale); // _res->GetFunction(fitfname.data())->GetParameter(1) * scale); //_res->GetMean());
      hBias->SetBinError(iB, /* 0.); */_res->GetFunction(fitfname.data())->GetParError(1) * scale); //0.);
      hReso->SetBinContent(iB,  _res->GetFunction(fitfname.data())->GetParameter(2) * scale); //_res->GetFunction(fitfname.data())->GetParameter(2) * scale); //_res->GetStdDev());
      hReso->SetBinError(iB, /* 0.);*/_res->GetFunction(fitfname.data())->GetParError(2) * scale); //0.);

    }

    TF1 f(Form("ffit_%d", iH), "[0]*std::pow(x, [1]) + [2]", 0., 10.);
    if (iH > 0) {
      hReso->Fit(Form("ffit_%d", iH), "MQ+");
    }


    hBias->Write();
    hReso->Write();
  }

  fo->Close();

  return 0;
}
