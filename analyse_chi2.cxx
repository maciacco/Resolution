#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>

int nSmp = 10;

//#define _FITCHI2_

bool relative = true;

void analyse_chi2(const char* ifname = "fitsCoarseFineM", const char* ofname = "bias_reso_vs_pt_coarse_fine_m"){
  TFile* infile = TFile::Open(Form("%s.root", ifname));

  const double binEdge[]{0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 5.};
  const int nBins = 21;

  double x_mins[10]{0};
  double y_mins[10]{0};

  double x_pitch{0};
  double y_pitch{0};

  double x_minn, x_minn_err, y_minn, y_minn_err;

  TH1D* hBias = new TH1D("hBias", ";1/#it{p}_{T} (#it{c}/GeV);#delta_{1/#it{p}_{T}} (#it{c}/GeV)", nBins, binEdge); //40, 0., 4.);
  TH1D* hReso = new TH1D("hReso", ";1/#it{p}_{T} (#it{c}/GeV);#sigma_{1/#it{p}_{T}} (#it{c}/GeV)", nBins, binEdge); //40, 0., 4.);
  for (int iB{4}; iB < nBins + 1; ++iB) {
    TH3D* chi2_prof_3d = static_cast<TH3D*>(infile->Get(Form("chi2_prof_%d", iB)));

    double mean_d{0}, sigma_d{0};
    double mean_s{0}, sigma_s{0};

    for (int iS{0}; iS < nSmp; ++iS){
      chi2_prof_3d->GetXaxis()->SetRangeUser(iS, iS + 1);
      TH2D* chi2_prof = static_cast<TH2D*>(chi2_prof_3d->Project3D("zy"));
      int bin_x_min = 0, bin_y_min = 0, bin_z_min = 0;
      chi2_prof->GetMinimumBin(bin_x_min, bin_y_min, bin_z_min);

      if (iB == 4) {
        x_pitch = chi2_prof->GetXaxis()->GetBinWidth(1);
        y_pitch = chi2_prof->GetYaxis()->GetBinWidth(1);
      }
      double x_min = chi2_prof->GetXaxis()->GetBinCenter(bin_x_min);
      double y_min = chi2_prof->GetYaxis()->GetBinCenter(bin_y_min);

      // chi2_prof->GetXaxis()->SetRangeUser(x_min - chi2_prof->GetXaxis()->GetBinWidth(1) * 6., x_min + chi2_prof->GetXaxis()->GetBinWidth(1) * 6.);
      // chi2_prof->GetYaxis()->SetRangeUser(y_min - chi2_prof->GetYaxis()->GetBinWidth(1) * 6., y_min + chi2_prof->GetYaxis()->GetBinWidth(1) * 6.);

#ifdef _FITCHI2_
      TF2 fChi2("fChi2", "((x-[0])*(x-[0])/[1]/[1]-2.*[2]*(x-[0])*(y-[3])/[1]/[4]+(y-[3])*(y-[3])/[4]/[4])/(2.*(1.-[2]*[2]))+[5]", -0.01, 0.01, 0., 0.05);
      fChi2.SetParLimits(0, -0.03, 0.01);
      fChi2.SetParameter(0, x_min);
      fChi2.SetParLimits(1, 0., 0.1);
      fChi2.SetParLimits(3, 0., 0.2);
      fChi2.SetParameter(3, y_min);
      fChi2.SetParLimits(4, 0., .1);

      TFitResultPtr res = chi2_prof->Fit("fChi2", "QMS+"); int i = 0;
      while (res->CovMatrixStatus() < 3 && i < 5) {
        chi2_prof->Fit("fChi2", "QSM+");
        i++;
      }

      double mux = fChi2.GetParameter(0);
      double sigx = fChi2.GetParameter(1);
      double muy = fChi2.GetParameter(3);
      double sigy = fChi2.GetParameter(4);

      x_minn = x_min;
      x_minn_err = sigx;
      y_minn = y_min;
      y_minn_err = sigy;

#endif // _FITCHI2_

      mean_d += x_min;
      mean_s += y_min;
      std::cout << "sigma_est " << iB << " " << iS << " " << y_min << std::endl;
      sigma_d += (x_min * x_min);
      sigma_s += (y_min * y_min);

      x_mins[iS] = x_min;
      y_mins[iS] = y_min;
    }

    double x{0}, y{0}, x_s{0}, y_s{0};
    for (int iS{0}; iS < nSmp; ++iS){
      x += x_mins[iS];
      y += y_mins[iS];
    }
    x /= nSmp; y /= nSmp;
    for (int iS{0}; iS < nSmp; ++iS){
      x_s += std::pow(x_mins[iS] - x, 2.);
      y_s += std::pow(y_mins[iS] - y, 2.);
    }
    x_s = std::sqrt(x_s / nSmp / (nSmp - 1));
    y_s = std::sqrt(y_s / nSmp / (nSmp - 1));

    if (x_s < x_pitch / std::sqrt(12.)) x_s = x_pitch / std::sqrt(12.);
    if (y_s < y_pitch / std::sqrt(12.)) y_s = y_pitch / std::sqrt(12.);

    float scale = relative ? 1. / hBias->GetXaxis()->GetBinCenter(iB) : 1.;

    hBias->SetBinContent(iB, x * scale);
    hBias->SetBinError(iB, x_s * scale);
    hReso->SetBinContent(iB, y * scale);
    hReso->SetBinError(iB, y_s * scale);

  }

  TFile* fo = TFile::Open(Form("%s.root", ofname), "recreate");
  fo->cd();
  hBias->Write();
  hReso->Write();
  fo->Close();
}

