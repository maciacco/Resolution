//const double lims[]{-0.01, 0.0025};
//const char *ax_y{"#delta(#it{p}_{T})/#it{p}_{T} from K^{0}_{S}"};

const double lims[]{0., 0.025};
const char *ax_y{"#sigma(#it{p}_{T})/#it{p}_{T} from K^{0}_{S}"};

const char *ax_x{"#it{p}_{T} (GeV/#it{c})"};

void plots(const char* hname = "hReso"){
  gStyle->SetOptStat(0);

  TFile *fda = TFile::Open("bias_reso_vs_pt_coarse_fine_d.root");
  TFile *fmc = TFile::Open("bias_reso_vs_pt_coarse_fine_m.root");

  TLegend leg(0.2, 0.2, 0.5, 0.4);

  TCanvas c("c", "c", 600, 600);
  c.SetLeftMargin(0.18);
  c.SetBottomMargin(0.15);
  c.SetTopMargin(0.03);
  c.SetRightMargin(0.03);
  TH1D* h_da = static_cast<TH1D*>(fda->Get(hname));
  TH1D* h_mc = static_cast<TH1D*>(fmc->Get(hname));

  h_da->SetLineWidth(2);
  h_da->SetLineColor(kRed);
  h_mc->SetLineWidth(2);
  h_mc->SetLineColor(kBlue);
  h_da->SetMarkerStyle(20);
  h_mc->SetMarkerStyle(21);
  h_da->SetMarkerSize(0.8);
  h_mc->SetMarkerSize(0.8);
  h_da->SetMarkerColor(kRed);
  h_mc->SetMarkerColor(kBlue);

  c.cd();
  h_da->Draw("pe");
  h_mc->Draw("samepe");

  leg.SetTextFont(45);
  leg.SetTextSize(25);
  leg.SetHeader("ALICE Work in progress, pp 13.6 TeV");
  leg.AddEntry(h_da, "Data, LHC24_pass1_minBias", "pe");
  leg.AddEntry(h_mc, "MC, LHC24j7 + LHC25a3", "pe");
  leg.Draw("same");

  h_da->GetYaxis()->SetRangeUser(lims[0], lims[1]);
  h_da->GetYaxis()->SetTitle(ax_y);
  h_da->GetXaxis()->SetTitle(ax_x);
  h_da->GetYaxis()->SetTitleFont(45);
  h_da->GetXaxis()->SetTitleFont(45);
  h_da->GetYaxis()->SetTitleSize(30);
  h_da->GetXaxis()->SetTitleSize(30);
  h_da->GetYaxis()->SetTitleOffset(1.8);

  c.Print(Form("%s_da_vs_mc.pdf", hname));
}
