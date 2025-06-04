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

double binEdges[] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 5.5, 6., 7., 8.};
const int nBins = 34;

void configureDSCB(TF1* f) {
  f->SetParLimits(0, 1.e-5, 1.e7);
  f->SetParameter(0, 1.e4);
  f->SetParLimits(1, 0.48, 0.50);
  f->SetParLimits(2, 0.001, 0.006);
  f->SetParameter(2, 0.004);
  f->SetParameter(1, 0.49);
  f->SetParLimits(3, 0.5, 2.);
  f->SetParameter(3, 1.5);
  f->SetParLimits(4, 0., 50.);
  f->SetParameter(4, 5.);
  f->SetParLimits(5, 0.5, 4.);
  f->SetParameter(5, 1.5);
  f->SetParLimits(6, 0., 50.);
  f->SetParameter(6, 5.);
  f->SetParLimits(7, 0., 100.);
  f->SetParameter(7, 5.);
  f->SetParLimits(8, -10., 10.);
  f->SetParameter(8, -1.);

  f->SetParError(0, 1.e-5);
  f->SetParError(1, 1.e-5);
  f->SetParError(2, 1.e-5);
}

void compareMass() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile* finSignl = TFile::Open("foo_dat_1.root");
  TFile* finToyMC = TFile::Open("check_k0s_toy_test_4.root");

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
    while (fitResult->CovMatrixStatus() < 3 && fitResult->Status() != 0 && retryCount++ < 10) {
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

    c[iB]->SetLeftMargin(0.18);
    c[iB]->SetBottomMargin(0.15);
    c[iB]->SetTopMargin(0.03);
    c[iB]->SetRightMargin(0.03);
    hToyProj->GetXaxis()->SetRangeUser(0.465, 0.53);
    hToyProj->GetXaxis()->SetTitleFont(45);
    hToyProj->GetXaxis()->SetTitleSize(30);
    hToyProj->GetYaxis()->SetTitleFont(45);
    hToyProj->GetYaxis()->SetTitleSize(30);
    hToyProj->GetYaxis()->SetTitle("Entries");
    hToyProj->GetYaxis()->SetTitleOffset(1.7);
    hToyProj->Draw("histoe");
    hSignlProj->Draw("histoesame");

    const double ymax = hSignlProj->GetBinContent(hSignlProj->GetMaximumBin());
    TLine ll(mean - 1.5 * stddev, 0, mean - 1.5 * stddev, ymax);

    ll.SetLineStyle(kDashed);
    ll.Draw("same");
    TLine lr(mean + 1.5 * stddev, 0, mean + 1.5 * stddev, ymax);
    lr.SetLineStyle(kDashed);
    lr.Draw("same");

    fout->cd();
    c[iB]->Write();
    std::string cname_file = cname + ".pdf";
    c[iB]->Print(cname_file.c_str());
  }

  for (int i{0}; i < 2; ++i){
    width[i]->SetLineWidth(2);
    centre[i]->SetLineWidth(2);
  }

  // sigma
  TCanvas cw("cw", "cw", 600, 600);
  cw.cd();
  cw.SetLeftMargin(0.18);
  cw.SetBottomMargin(0.15);
  cw.SetTopMargin(0.03);
  cw.SetRightMargin(0.03);
  width[1]->SetLineColor(kRed);
  width[0]->SetLineColor(kBlue);
  width[1]->SetMarkerColor(kRed);
  width[0]->SetMarkerColor(kBlue);
  width[1]->SetMarkerStyle(20);
  width[0]->SetMarkerStyle(54);

  TLegend lw(0.2, 0.2, 0.5, 0.4);

  width[1]->GetYaxis()->SetTitle("#sigma_{#it{M}} (GeV/#it{c}^{2})");
  width[1]->GetYaxis()->SetTitleFont(45);
  width[1]->GetXaxis()->SetTitleFont(45);
  width[1]->GetYaxis()->SetTitleSize(30);
  width[1]->GetXaxis()->SetTitleSize(30);
  width[1]->GetYaxis()->SetTitleOffset(1.8);

  width[1]->Draw("pe");
  width[0]->Draw("pesame");

  lw.SetTextFont(45);
  lw.SetTextSize(25);
  lw.SetHeader("ALICE Work in progress, pp 13.6 TeV");
  lw.AddEntry(width[1], "Data, LHC24_pass1_minBias", "pe");
  lw.AddEntry(width[0], "Toy MC, #delta(#it{p}_{T}) and #sigma(#it{p}_{T}) from data", "pe");
  lw.Draw("same");
  cw.Print("cw.pdf");

  // mu
  TCanvas cm("cm", "cm", 600, 600);
  cm.cd();
  cm.SetLeftMargin(0.18);
  cm.SetBottomMargin(0.15);
  cm.SetTopMargin(0.03);
  cm.SetRightMargin(0.03);
  centre[1]->GetYaxis()->SetRangeUser(0.495, 0.4966);
  centre[1]->SetLineColor(kRed);
  centre[0]->SetLineColor(kBlue);
  centre[1]->SetMarkerColor(kRed);
  centre[0]->SetMarkerColor(kBlue);
  centre[1]->SetMarkerStyle(20);
  centre[0]->SetMarkerStyle(54);

  TLegend lm(0.2, 0.2, 0.5, 0.4);

  centre[1]->GetYaxis()->SetTitle("#mu_{#it{M}} (GeV/#it{c}^{2})");
  centre[1]->GetYaxis()->SetTitleFont(45);
  centre[1]->GetXaxis()->SetTitleFont(45);
  centre[1]->GetYaxis()->SetTitleSize(30);
  centre[1]->GetXaxis()->SetTitleSize(30);
  centre[1]->GetYaxis()->SetTitleOffset(1.8);

  centre[1]->Draw("pe");
  centre[0]->Draw("pesame");

  lm.SetTextFont(45);
  lm.SetTextSize(25);
  lm.SetHeader("ALICE Work in progress, pp 13.6 TeV");
  lm.AddEntry(centre[1], "Data, LHC24_pass1_minBias", "pe");
  lm.AddEntry(centre[0], "Toy MC, #delta(#it{p}_{T}) and #sigma(#it{p}_{T}) from data", "pe");
  lm.Draw("same");
  cm.Print("cm.pdf");


  width[0]->Write();
  width[1]->Write();
  centre[0]->Write();
  centre[1]->Write();
  fout->Close();
}
