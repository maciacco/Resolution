#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>

#include <stdlib.h>
#include <TROOT.h>

/*
TH1D* mass_mc;
double fsig(double *x, double *par) {
  double xx = x[0];
  int bin = mass_mc->GetXaxis()->FindBin(xx);
  return par[0] * mass_mc->GetBinContent(bin) + par[1] * std::exp(xx * par[2] + par[3]);

}
*/

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
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + std::exp(par[7] + par[8] * x[0]));
}

int main(const int argc, const char* argv[]) {

  int dMin = atoi(argv[1]);
  int dMax = atoi(argv[2]);
  gROOT->SetBatch();

  TFile *f_mc = TFile::Open("templates2.root");
  TFile *f_dat = TFile::Open("foo_dat_0.root");

  TFile *fo = TFile::Open(Form("fits_%d_%d.root", dMin, dMax), "recreate");

  if (!f_mc || !f_dat) {
    std::cout << "Files not found!" << std::endl;
    return 1;
  }


  for (int iB{4}; iB < 41; ++iB){
    std::cout << "Processing bin n. " << iB << std::endl;

    TH2D *_mass_2d = static_cast<TH2D*>(f_dat->Get("hMassK0s"));
    if (!_mass_2d) {
      std::cout << "Mass histogram not found!" << std::endl;
      return 1;
    }

    TH1D* _mass = static_cast<TH1D*>(_mass_2d->ProjectionY(Form("_mass_1d_%d", iB), iB, iB));
    _mass->Scale(1. / _mass->GetBinContent(_mass->GetMaximumBin()));

    TH2D* chi2_prof = new TH2D(Form("chi2_prof_%d", iB), ";#delta;#sigma", 120, -0.010125, 0.019875, 74, 0.0005, 0.0745);

    TDirectory* dir = fo->mkdir(Form("%d_fits", iB));

    TF1* tmp_dscb = new TF1("tmp_dscb", DoubleSidedCB, 0.45, 0.5, 9);
    tmp_dscb->SetParLimits(0, 1.e-7, 1.e5);
    tmp_dscb->SetParLimits(1, 0.47, 0.5);
    tmp_dscb->SetParLimits(2, 0.0005, 0.01);
    tmp_dscb->SetParLimits(3, 10., 20.);
    tmp_dscb->SetParLimits(4, 0., 40.);
    tmp_dscb->SetParLimits(5, 10., 20.);
    tmp_dscb->SetParLimits(6, 0., 40.);

    for (int i{0}; i < 1; ++i)_mass->Fit(tmp_dscb, "QNRM+", "", 0.47, 0.53);
    tmp_dscb->SetParLimits(3, 0.5, 3.);
    tmp_dscb->SetParLimits(5, 0.5, 3.);
    tmp_dscb->SetParameter(4, 2.);
    tmp_dscb->SetParameter(6, 2.);
    for (int i{0}; i < 2; ++i)_mass->Fit(tmp_dscb, "QRM+", "", 0.47, 0.53);

    for (int iD{dMin}; iD < dMax; ++iD) {
      for (int iS{1}; iS < 75; ++iS) {
        TH1D* mass = new TH1D(*_mass);
        mass->SetName(Form("hMassFit_%d_%d", iD, iS));
        double chi2_tmp = 1.e6;

        TH2D* mass_mc_2 = static_cast<TH2D*>(f_mc->Get(Form("hMass_%d_%d", iD, iS)));

        double x_bin = chi2_prof->GetXaxis()->GetBinCenter(iD);
        double y_bin = chi2_prof->GetYaxis()->GetBinCenter(iS);

        if (mass_mc_2 != nullptr) {
          TH1D* mass_mc = static_cast<TH1D*>(mass_mc_2->ProjectionY(Form("mass_mc_%d_%d", iD, iS), iB, iB));

//          TF1* fitfun = new TF1("fun", fsig, 0., 0.6, 4);
//          fitfun->FixParameter(1, 0.);
//          fitfun->SetNpx(1000);

          TF1* fitfun_ = new TF1("fun_", "gaus", 0., 0.6);
//          fitfun_->SetNpx(1000);
          mass_mc->Fit(fitfun_, "RQ+", "", 0., 0.6);
          fitfun_->FixParameter(1, fitfun_->GetParameter(1));
          fitfun_->FixParameter(2, fitfun_->GetParameter(2));

          // gaussian smoothening
          for (int ii{1}; ii < mass_mc->GetNbinsX(); ++ii) {
            mass_mc->SetBinContent(ii, fitfun_->Eval(mass_mc->GetXaxis()->GetBinCenter(ii)));
          }

          TF1* fitfun = new TF1("fun", DoubleSidedCB, 0., 0.6, 9);
          fitfun->FixParameter(1, fitfun_->GetParameter(1));
          fitfun->FixParameter(2, fitfun_->GetParameter(2));
          fitfun->FixParameter(3, tmp_dscb->GetParameter(3));
          fitfun->FixParameter(4, tmp_dscb->GetParameter(4));
          fitfun->FixParameter(5, tmp_dscb->GetParameter(5));
          fitfun->FixParameter(6, tmp_dscb->GetParameter(6));
          fitfun->FixParameter(7, tmp_dscb->GetParameter(7));
          fitfun->FixParameter(8, tmp_dscb->GetParameter(8));

          fitfun->SetParLimits(0, 0., 1.e2);
          mass->Fit(fitfun, "QMB+", "", 0.47, 0.53);

          // get chi2
          chi2_tmp = fitfun->GetChisquare();
//          dir->cd();
//          mass_mc->Write();

          delete fitfun_;
          delete fitfun;
        }

        delete mass;

        chi2_prof->SetBinContent(iD + 1, iS, chi2_tmp);

//      dir->cd();
//      mass->Write();
      }
    }

    fo->cd();
    chi2_prof->Write();
  }

  fo->Close();

  return 0;
}
