Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2) {
  const Double_t u = (x - mu) / width;

  const Double_t A1 = TMath::Power(p1 / TMath::Abs(a1), p1) * TMath::Exp(-a1 * a1 / 2);
  const Double_t A2 = TMath::Power(p2 / TMath::Abs(a2), p2) * TMath::Exp(-a2 * a2 / 2);
  const Double_t B1 = p1 / TMath::Abs(a1) - TMath::Abs(a1);
  const Double_t B2 = p2 / TMath::Abs(a2) - TMath::Abs(a2);

  if (u < -a1) return A1 * TMath::Power(B1 - u, -p1);
  if (u < a2)  return TMath::Exp(-u * u / 2);
  return A2 * TMath::Power(B2 + u, -p2);
}

Double_t DoubleSidedCB(double* x, double* par) {
#ifdef _MC_
  return par[0] * TMath::Exp(-0.5 * std::pow((x[0] - par[1]) / par[2], 2.));
#else
  return par[0] * DoubleSidedCB2(x[0], par[1], par[2], par[3], par[4], par[5], par[6]) + par[7] * x[0] + par[8];
#endif
}

double binEdges[] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 5.5, 6., 6.5, 7., 8., 10.};
const int nBins = 35;

void configureDSCB(TF1* f) {
  f->SetParLimits(0, 1.e-5, 1.e8);
  f->SetParLimits(1, 0.49, 0.50);
  f->SetParLimits(2, 0.001, 0.005);
  f->SetParameter(2, 0.003);
  f->SetParLimits(3, 0.5, 2.);
  f->SetParameter(3, 1.5);
  f->SetParLimits(4, 0., 30.);
  f->SetParameter(4, 5.);
  f->SetParLimits(5, 0.5, 4.);
  f->SetParameter(5, 1.5);
  f->SetParLimits(6, 0., 30.);
  f->SetParameter(6, 10.);
  f->SetParLimits(7, 0., 100.);
  f->SetParameter(7, 5.);
  f->SetParLimits(8, -10., 10.);
  f->SetParameter(8, -1.);

  f->SetParError(0, 1.e-5);
  f->SetParError(1, 1.e-5);
  f->SetParError(2, 1.e-5);
}

void compareMass() {
  TFile* finSignl = TFile::Open("foo_dat_0.root");
  TFile* finToyMC = TFile::Open("check_k0s_toy_test_1.root");

  TH1D* width[2];
  TH1D* centre[2];
  for (int i = 0; i < 2; ++i) {
    width[i] = new TH1D(Form("Width_%d", i), ";#it{p}_{T} (GeV/#it{c});Width (GeV/#it{c}^{2})", nBins, binEdges);
    centre[i] = new TH1D(Form("Centre_%d", i), ";#it{p}_{T} (GeV/#it{c});Centre (GeV/#it{c}^{2})", nBins, binEdges);
  }

  auto* hSignl = static_cast<TH2D*>(finSignl->Get("hMassK0sPtK0s"));
  auto* hToyMC = static_cast<TH2D*>(finToyMC->Get("hMass_vs_pt"));

  TFile* fout = TFile::Open("compareMass2.root", "recreate");
  TCanvas* c[50];

  for (int iB = 7; iB < nBins + 1; ++iB) {
    int binLow = hSignl->GetXaxis()->FindBin(binEdges[iB - 1]);
    int binUp = hSignl->GetXaxis()->FindBin(binEdges[iB]);

    auto* hSignlProj = hSignl->ProjectionY(Form("hSignl_%d", iB), binLow, binUp);
    auto* hToyProj = hToyMC->ProjectionY(Form("hToyMC_%d", iB), binLow, binUp);

    const double mean = hToyProj->GetMean();
    const double stddev = hToyProj->GetStdDev();

    const double integralSignl = hSignlProj->Integral(hSignlProj->FindBin(mean - 1.5 * stddev), hSignlProj->FindBin(mean + 1.5 * stddev));
    const double integralMC = hToyProj->Integral(hToyProj->FindBin(mean - 1.5 * stddev), hToyProj->FindBin(mean + 1.5 * stddev));

    hToyProj->Scale(integralSignl / integralMC);
    hToyProj->SetLineColor(kBlue);
    hSignlProj->SetLineColor(kRed);

    hToyProj->Fit("gaus", "LRMQ+");

    auto* dscb = new TF1("dscb", DoubleSidedCB, 0.465, 0.53, 9);
    configureDSCB(dscb);

    TFitResultPtr fitResult = hSignlProj->Fit(dscb, "LQNSRM+", "", 0.465, 0.53);
    int retryCount = 0;
    while (fitResult->CovMatrixStatus() < 3 && fitResult->Status() != 0 && retryCount++ < 5) {
      fitResult = hSignlProj->Fit(dscb, "NLQSRM+", "", 0.465, 0.53);
    }

    hSignlProj->Fit(dscb, "LSRM+", "", 0.465, 0.53);

    width[0]->SetBinContent(iB, hToyProj->GetFunction("gaus")->GetParameter(2));
    width[0]->SetBinError(iB, hToyProj->GetFunction("gaus")->GetParError(2));
    width[1]->SetBinContent(iB, dscb->GetParameter(2));
    width[1]->SetBinError(iB, dscb->GetParError(2));

    centre[0]->SetBinContent(iB, hToyProj->GetFunction("gaus")->GetParameter(1));
    centre[0]->SetBinError(iB, hToyProj->GetFunction("gaus")->GetParError(1));
    centre[1]->SetBinContent(iB, dscb->GetParameter(1));
    centre[1]->SetBinError(iB, dscb->GetParError(1));

    hToyProj->SetLineWidth(2);
    hSignlProj->SetLineWidth(2);
    hToyProj->SetFillStyle(3004);
    hToyProj->SetFillColor(kBlue);

    // Plot
    const std::string cname = "c_" + std::to_string(iB);
    c[iB] = new TCanvas(cname.c_str(), cname.c_str(), 600, 600);
    c[iB]->cd();
    hToyProj->Draw("histoe");
    hSignlProj->Draw("histoesame");

    const double ymax = hSignlProj->GetBinContent(hSignlProj->GetMaximumBin());
    TLine(mean - 1.5 * stddev, 0, mean - 1.5 * stddev, ymax).Draw("same");
    TLine(mean + 1.5 * stddev, 0, mean + 1.5 * stddev, ymax).Draw("same");

    fout->cd();
    c[iB]->Write();
  }

  width[0]->Write();
  width[1]->Write();
  centre[0]->Write();
  centre[1]->Write();
  fout->Close();
}
