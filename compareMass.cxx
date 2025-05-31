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
  return par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + par[7] * x[0] + par[8];
#endif
}

double binEdge[]{0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 5.5, 6., 6.5, 7., 8., 10.};
int nBins = 35;

void compareMass(){
  TFile *finSignl = TFile::Open("foo_dat_0.root");
  TFile *finToyMC = TFile::Open("check_k0s_toy_test_2.root");

  TH1D* width[2];
  TH1D* centre[2];
  for (int i{0}; i < 2; ++i){
    width[i] = new TH1D(Form("Width_%d", i), ";#it{p}_{T} (GeV/#it{c});Width (GeV/#it{c}^{2})", nBins, binEdge);
    centre[i] = new TH1D(Form("Centre_%d", i), ";#it{p}_{T} (GeV/#it{c});Centre (GeV/#it{c}^{2})", nBins, binEdge);
  }

  TH2D* hSignl = static_cast<TH2D*>(finSignl->Get("hMassK0sPtK0s"));
  TH2D* hToyMC = static_cast<TH2D*>(finToyMC->Get("hMass_vs_pt"));

  TFile *fout = TFile::Open("compareMass2.root", "recreate");
  TCanvas *c[50];
  for (int iB{7}; iB < nBins + 1; ++iB) {
    int binLow = hSignl->GetXaxis()->FindBin(binEdge[iB - 1]);
    int binUp = hSignl->GetXaxis()->FindBin(binEdge[iB]);
    TH1D* hPySignl = static_cast<TH1D*>(hSignl->ProjectionY(Form("hSignl_%d", iB), binLow, binUp));
    TH1D* hPyToyMC = static_cast<TH1D*>(hToyMC->ProjectionY(Form("hToyMC_%d", iB), binLow, binUp));
    std::string cname{"c_"};
    cname += std::to_string(iB);
    c[iB] = new TCanvas(cname.data(), cname.data(), 600, 600);
    c[iB]->cd();
//    hPySignl->Rebin(2);
//    hPyToyMC->Rebin(2);
    double mean = hPyToyMC->GetMean();
    double stddev = hPyToyMC->GetStdDev();
    double intSig = hPySignl->Integral(hPySignl->FindBin(mean - 1.5 * stddev), hPySignl->FindBin(mean + 1.5 * stddev));
    double intMC = hPyToyMC->Integral(hPySignl->FindBin(mean - 1.5 *stddev), hPyToyMC->FindBin(mean + 1.5 * stddev));
    hPyToyMC->Scale(intSig / intMC);
    hPyToyMC->SetLineColor(kBlue);
    hPySignl->SetLineColor(kRed);

    hPyToyMC->Fit("gaus", "LRMQ+");

    TF1* tmp_dscb = new TF1("tmp_dscb", DoubleSidedCB, 0.465, 0.53, 9);
    tmp_dscb->SetParLimits(0, 1.e-5, 1.e8);
    tmp_dscb->SetParLimits(1, 0.49, 0.50); // 0.47, 0.5 for data
    tmp_dscb->SetParLimits(2, 0.001, 0.005);
    tmp_dscb->SetParameter(2, 0.003);
    tmp_dscb->SetParLimits(3, 0.5, 2.);
    tmp_dscb->SetParameter(3, 1.5);
    tmp_dscb->SetParLimits(4, 0., 30.);
    tmp_dscb->SetParameter(4, 5.);
    tmp_dscb->SetParLimits(5, 0.5, 4.);
    tmp_dscb->SetParameter(5, 1.5);
    tmp_dscb->SetParLimits(6, 0., 30.);
    tmp_dscb->SetParameter(6, 10.);
    tmp_dscb->SetParLimits(7, 0., 100.);
    tmp_dscb->SetParameter(7, 5.);
    tmp_dscb->SetParLimits(8, -10., 10.);
    tmp_dscb->SetParameter(8, -1.);

    tmp_dscb->SetParError(0, 1.e-5);
    tmp_dscb->SetParError(1, 1.e-5);
    tmp_dscb->SetParError(2, 1.e-5);

    TFitResultPtr res = hPySignl->Fit(tmp_dscb, "LQNSRM+", "", 0.465, 0.53);
    int i = 0;
    while (res->CovMatrixStatus() < 3 && res->Status() != 0) {
      res = hPySignl->Fit(tmp_dscb, "NLQSRM+", "", 0.465, 0.53);
      i += 1;
    }
    hPySignl->Fit(tmp_dscb, "LSRM+", "", 0.465, 0.53);

    width[0]->SetBinContent(iB, hPyToyMC->GetFunction("gaus")->GetParameter(2));
    width[0]->SetBinError(iB, hPyToyMC->GetFunction("gaus")->GetParError(2));
    width[1]->SetBinContent(iB, tmp_dscb->GetParameter(2));
    width[1]->SetBinError(iB, tmp_dscb->GetParError(2));

    centre[0]->SetBinContent(iB, hPyToyMC->GetFunction("gaus")->GetParameter(1));
    centre[0]->SetBinError(iB, hPyToyMC->GetFunction("gaus")->GetParError(1));
    centre[1]->SetBinContent(iB, tmp_dscb->GetParameter(1));
    centre[1]->SetBinError(iB, tmp_dscb->GetParError(1));

    hPyToyMC->SetLineWidth(2);
    hPySignl->SetLineWidth(2);
    hPyToyMC->SetFillStyle(3004);
    hPyToyMC->SetFillColor(kBlue);
    hPyToyMC->Draw("histoe");
    hPySignl->Draw("histoesame");
    TLine lleft(mean - 1.5 * stddev, 0, mean - 1.5 * stddev, hPySignl->GetBinContent(hPySignl->GetMaximumBin()));
    TLine lright(mean + 1.5 * stddev, 0, mean + 1.5 * stddev, hPySignl->GetBinContent(hPySignl->GetMaximumBin()));
    lleft.Draw("same");
    lright.Draw("same");
    fout->cd();
    c[iB]->Write();
  }

  width[0]->Write();
  width[1]->Write();
  centre[0]->Write();
  centre[1]->Write();


  fout->Close();
}
