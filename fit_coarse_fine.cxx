#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <stdlib.h>
#include <TROOT.h>

//#define _MC_

bool invpT = false;

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
  return par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + par[7] * std::exp(par[8] * x[0]);
#endif
}

const double binEdge[]{0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 5.};
const int nBins = 21;

int main(const int argc, const char* argv[]) {

  int dMin = atoi(argv[1]);
  int dMax = atoi(argv[2]);
  std::string in_fname = argv[3];
  in_fname = in_fname + ".root";
  gROOT->SetBatch();

  TFile *f_mc = TFile::Open("templates__pt_smallFlatFine.root"); // template_pt
  TFile *f_dat = TFile::Open(in_fname.data());

  TFile *fo = TFile::Open(Form("fits_%d_%d.root", dMin, dMax), "recreate");

  if (!f_mc || !f_dat) {
    std::cerr << "Files not found!" << std::endl;
    return 1;
  }

  TH3D *_mass_3d = static_cast<TH3D*>(f_dat->Get(invpT ? "hMassK0s" : "hMassK0sPt"));
    if (!_mass_3d) {
    std::cerr << "Mass histogram not found!" << std::endl;
    return 1;
  }

  TAxis refAxis(nBins, binEdge);

  for (int iB{4}; iB < nBins + 1; ++iB){
    double p_low = refAxis.GetBinLowEdge(iB);
    double p_up = refAxis.GetBinUpEdge(iB);
    int iB_low = _mass_3d->GetYaxis()->FindBin(p_low + 0.005f);
    int iB_up = _mass_3d->GetYaxis()->FindBin(p_up - 0.005f);
    std::cout << "Processing bin n. " << iB << std::endl;

    TH3D* chi2_prof = new TH3D(Form("chi2_prof_%d", iB), ";Sample;#delta;#sigma", 20, 0., 20., 160, invpT ? -0.010125 : -0.030125, invpT ? 0.02025 : 0.01025, 400, -0.00025, 0.20025);
    TDirectory* dir = fo->mkdir(Form("%d_fits", iB));


    for (int iSamp{dMin}; iSamp < dMax; ++iSamp){

      TH1D* _mass = static_cast<TH1D*>(_mass_3d->ProjectionZ(Form("_mass_1d_1_%d_%d_%d", iSamp, iB, dMin), 2 * (iSamp) + 1, 2 * (iSamp) + 2, iB_low, iB_up));
//      _mass->Scale(1. / _mass->Integral(1, _mass->GetNbinsX()));

      TF1* tmp_dscb = new TF1("tmp_dscb", DoubleSidedCB, 0.46, 0.53, 9);
      tmp_dscb->SetParLimits(0, 1., 1.e7); // 1.e5 for data // 1.e9);
      tmp_dscb->SetParLimits(1, 0.48, 0.51); // 0.48, 0.505 for data
      tmp_dscb->SetParLimits(2, 0.0005, 0.02);
      tmp_dscb->SetParameter(0, _mass->GetEntries()); // 1.e3 for data
      tmp_dscb->SetParameter(1, 0.5 * (0.48 + 0.51)); // 0.48, 0.505 for data
      tmp_dscb->SetParameter(2, 0.001 ); //(0.0003 + 0.02));

#ifndef _MC_
      tmp_dscb->SetParLimits(3, 1., 30.); // 0.5 for data
      tmp_dscb->SetParLimits(4, 1., 10.);
      tmp_dscb->SetParLimits(5, 1., 30.);
      tmp_dscb->SetParLimits(6, 1., 10.);
      tmp_dscb->SetParameter(3, 5.);
      tmp_dscb->SetParameter(4, 5.);
      tmp_dscb->SetParameter(5, 9.);
      tmp_dscb->SetParameter(6, 5.);

      tmp_dscb->SetParLimits(7, 0., 10.);
      tmp_dscb->SetParLimits(8, -5., 0.);
      tmp_dscb->SetParameter(7, 5.);
      tmp_dscb->SetParameter(8, -2.);
#endif

      TFitResultPtr res = _mass->Fit(tmp_dscb, "LQSRM+", "", 0.46, 0.53);
      int i = 0;
      while (res->CovMatrixStatus() < 3 && i < 10) {
        res = _mass->Fit(tmp_dscb, "LQSRM+", "", 0.46, 0.53);
        i += 1;
      }

//      for (int i{0}; i < 2; ++i)_mass->Fit(tmp_dscb, "LQRM+", "", 0.47, 0.53);

      dir->cd();
      _mass->Write();

      TH1D* mass = static_cast<TH1D*>(_mass_3d->ProjectionZ(Form("_mass_1d_%d_%d_%d", iSamp, iB, dMin), 2 * (iSamp) + 1, 2 * (iSamp) + 2, iB_low, iB_up));
//      mass->Scale(1. / mass->Integral(1, mass->GetNbinsX()));

     int min_coarse_x= 0, min_coarse_y = 0, min_coarse_z;
     int min_fine_x=0, min_fine_y=0, min_fine_z=0;
    for (int iI{0}; iI < 3; iI++){
      int iDmin = iI == 0 ? 0 : min_coarse_x - 12 - 1;
      int iDmax = iI == 0 ? 160 : min_coarse_x + 12 - 1;
      int iSmin = iI == 0 ? 0 : min_coarse_y - 12 - 1;
      int iSmax = iI == 0 ? 400 : min_coarse_y + 12 - 1;

      if (iDmin < 0) iDmin = 0;
      if (iSmin < 0) iSmin = 0;
      if (iDmax > 160) iDmax = 160;
      if (iSmax > 400) iSmax = 400;
      for (int iD{iDmin}; iD < iDmax; ++iD) {
        for (int iS{iSmin}; iS < iSmax; ++iS) {
          if ((iD % 10 != 0 || iS % 10 != 0) && iI == 0) {
            chi2_prof->SetBinContent(iSamp + 1, iD + 1, iS + 1, 1.e6);
            continue;
          }

          double chi2_tmp = 1.e6;
          TH2D* mass_mc_2 = static_cast<TH2D*>(f_mc->Get(Form("hMass_%d_%d%s", iD, iS, invpT ? "" : "_vs_pt")));

          if (mass_mc_2 != nullptr) {

            int iB_low_mc = mass_mc_2->GetXaxis()->FindBin(p_low + 0.005f);
            int iB_up_mc = mass_mc_2->GetXaxis()->FindBin(p_up - 0.005f);

            TH1D* mass_mc = static_cast<TH1D*>(mass_mc_2->ProjectionY(Form("mass_mc_%d_%d_%d", iSamp, iD, iS), iB_low_mc, iB_up_mc));
//            mass_mc->Scale(1. / mass_mc->Integral(1, mass_mc->GetNbinsX()));


//            TF1* fitfun_ = new TF1(Form("fun__%d_%d", iD, iS), "gaus", 0., 0.6);
//            fitfun_->SetParameter(0, 1.);
//            fitfun_->SetParameter(1, 0.5);
//            fitfun_->SetParameter(2, 0.005);
//            fitfun_->SetParLimits(0, 0., 1.e2);
//            fitfun_->SetParLimits(1, 0., 1.);
//            fitfun_->SetParLimits(2, 0., 0.1);

        //    for(int i{0}; i < 2; ++i)mass_mc->Fit(fitfun_, "LRQ+", "", 0., 0.6);


// TODO: check if this has any effect
//            mass_mc->Fit("gaus", "LQM+");

            //fitfun_->FixParameter(1, fitfun_->GetParameter(1));
            //fitfun_->FixParameter(2, fitfun_->GetParameter(2));

            // gaussian smoothening
            //for (int ii{1}; ii < mass_mc->GetNbinsX(); ++ii) {
            //  mass_mc->SetBinContent(ii, fitfun_->Eval(mass_mc->GetXaxis()->GetBinCenter(ii)));
           // }

            TF1* fitfun = new TF1(Form("fun_%d_%d", iD, iS), DoubleSidedCB, 0., 0.6, 9);
            // fitfun->FixParameter(0, tmp_dscb->GetParameter(0));
            fitfun->FixParameter(1, mass_mc->GetMean() ); //fitfun_->GetParameter(1));
            fitfun->FixParameter(2, mass_mc->GetStdDev() ); //fitfun_->GetParameter(2));

#ifndef _MC_
            for (int iP{3}; iP < 9; ++iP) {
//              if (iP == 7) continue;
              fitfun->FixParameter(iP, tmp_dscb->GetParameter(iP));
            }
#endif

            //fitfun->SetParLimits(0, 0., 1.e2);
            //fitfun->SetParameter(0, 0.5 * (0. + 1.e2));
            //fitfun->FixParameter();
            fitfun->FixParameter(0, tmp_dscb->GetParameter(0));
//            double lim_inf = - tmp_dscb->GetParameter(3) * tmp_dscb->GetParameter(2) + tmp_dscb->GetParameter(1);
//            double lim_sup = tmp_dscb->GetParameter(5) * tmp_dscb->GetParameter(2) + tmp_dscb->GetParameter(1);
            mass->Fit(fitfun, "LQMR+", "", 0.46, 0.53); // lim_inf, lim_sup);
//            fitfun->FixParameter(0, fitfun->GetParameter(0));

            if (iI == 2){
              if (iD + 1 == min_fine_x && iS + 1 == min_fine_y){
                TH1D* mass__ = static_cast<TH1D*>(_mass_3d->ProjectionZ(Form("_mass_1d_%d_%d_%d_", iSamp, iB, dMin), 2 * (iSamp) + 1, 2 * (iSamp) + 2, iB_low, iB_up));
  //              mass__->Scale(1. / mass__->Integral(1, mass__->GetNbinsX()));
                mass__->Fit(fitfun, "LRQM+", "", 0.46, 0.53);
                fitfun->SetLineColor(kBlue);
                fitfun->SetNpx(1000);
                mass__->Write();
                mass_mc->Write();
              }
            }

            // get chi2
            chi2_tmp = fitfun->GetChisquare();

//            delete fitfun_;
            delete fitfun;
          }

          chi2_prof->SetBinContent(iSamp + 1, iD + 1, iS + 1, chi2_tmp);
        }
      }
      if (iI == 0){
        chi2_prof->GetXaxis()->SetRangeUser(iSamp, iSamp + 1);
        TH2D* chi2_proj = (TH2D*)chi2_prof->Project3D("zy");
        chi2_proj->GetMinimumBin(min_coarse_x, min_coarse_y, min_coarse_z);
        chi2_prof->GetXaxis()->SetRangeUser(0, 20);
      }
      if (iI ==1){
        chi2_prof->GetXaxis()->SetRangeUser(iSamp, iSamp + 1);
        TH2D* chi2_proj = (TH2D*)chi2_prof->Project3D("zy");
        chi2_proj->GetMinimumBin(min_fine_x, min_fine_y, min_fine_z);
        chi2_prof->GetXaxis()->SetRangeUser(0, 20);
      }
      if (iI ==2)
        delete mass, _mass;

   } // iI
    }

    fo->cd();
    chi2_prof->Write();
  }

  fo->Close();

  return 0;
}
