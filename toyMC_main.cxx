#include <TF1.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h> // ugly, but TGenPhaseSpace uses this.
#include <TMath.h>
#include <TRandom.h>

constexpr double kK0sMass = 0.497648;
constexpr double kPiMass = 0.1395703918;

const double binEdge[]{0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 2.5, 3., 4., 5.};
const int nBins = 21;

bool relative = true;

bool shiftPt = true;

int main(int const argc, const char *argv[]){
  const char* name = argv[1];
  const int kNtrials = std::atoi(argv[2]);
  const int sample = std::atoi(argv[3]);
  gRandom->SetSeed(sample);

  TFile *myFile = TFile::Open(Form("templates_%s_%d.root", name, sample), "RECREATE");

  // get resolution parameters of phi and eta from file
  TFile *reso_file = TFile::Open("residuals_inv.root");
  // TH1D* resoEta = static_cast<TH1D*>(reso_file->Get("hEtaReso"));
  // TH1D* resoPhi = static_cast<TH1D*>(reso_file->Get("hPhiReso"));
  // TF1* fResoEta = resoEta->GetFunction("ffit_2");
  // TF1* fResoPhi = resoPhi->GetFunction("ffit_1");

  TLorentzVector mother, pip, pim;
  TGenPhaseSpace gen2Pi;
  const double massesDau[2]{kPiMass, kPiMass};

  TF1 mtExpo("mtExpo","[0]*x*std::exp(-std::hypot([2], x)/[1])", 0.1, 15.);
  mtExpo.SetParameter(0, 1.);
  mtExpo.SetParameter(1, 0.5);
  mtExpo.SetParameter(2, kK0sMass);

  double shift, sigma;

  for (int iD{0}; iD < 160; ++iD) {
    shift = -0.03/* : -0.01 */ + iD * 0.00025;
    std::cout << shift << "\n";

    for (int iS{0}; iS < 400; ++iS) {
      sigma = 0.0 + iS * 0.0005;

      // TH2D hPtPosNeg("hPrPtTmp", ";#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", 200, 0.f, 10.f, 200, 0.f, 10.f);
      // TH2D hArm("hArm", ";;", 200, -1., 1., 100., 0., 0.5);
      // TH2D hMass(Form("hMass_%d_%d", iD, iS), ";1 / #it{p}_{T} (#it{c}/GeV);#it{M} (GeV/#it{c}^{2})", 100, 0., 10., 1200, 0, 0.6);
      TH2D hMassPt(Form("hMass_%d_%d_vs_pt", iD, iS), ";#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2})", 100, 0., 10., 1200, 0, 0.6);
      // TH2D hEtaReso(Form("hEtaReso_%d_%d_vs_pt", iD, iS), ":#it{p}_{T} (GeV/#it{c});#eta^{rec} - #eta^{gen}", 100, 0., 10., 10000, -.5, .5);

      for (uint64_t i = 0; i < kNtrials; ++i) {
        float pT_k0s = gRandom->Uniform(0.1, 15.); // mtExpo.GetRandom();
        float eta = gRandom->Uniform(-0.3, 0.3);
        float phi = gRandom->Uniform(0, TMath::TwoPi());

        // generate mother and decay
        mother.SetPtEtaPhiM(pT_k0s, eta, phi, kK0sMass);

        gen2Pi.SetDecay(mother, 2, massesDau);
        gen2Pi.Generate();

        auto pos = gen2Pi.GetDecay(0);
        auto neg = gen2Pi.GetDecay(1);

        float pos_eta_in = pos->Eta();
        float pos_pt = pos->Pt() + gRandom->Gaus(shift, sigma); // * pos->Pt();
        float neg_pt = neg->Pt() + gRandom->Gaus(shift, sigma); // * neg->Pt();
        float pos_eta = pos->Eta(); //+ gRandom->Gaus(0., resoEta->GetBinContent(pos->Pt()) * resoEta->GetXaxis()->GetBinCenter(pos->Pt()));
        float neg_eta = neg->Eta(); //+ gRandom->Gaus(0., resoEta->GetBinContent(neg->Pt()) * resoEta->GetXaxis()->GetBinCenter(neg->Pt()));
        float pos_phi = pos->Phi(); //+ gRandom->Gaus(0., resoPhi->GetBinContent(pos->Pt()) * resoPhi->GetXaxis()->GetBinCenter(pos->Pt()));
        float neg_phi = neg->Phi(); //+ gRandom->Gaus(0., resoPhi->GetBinContent(neg->Pt()) * resoPhi->GetXaxis()->GetBinCenter(neg->Pt()));

        pos->SetPtEtaPhiM(pos_pt, pos_eta, pos_phi, kPiMass);
        neg->SetPtEtaPhiM(neg_pt, neg_eta, neg_phi, kPiMass);

        float pT_pos = pos->Pt();
        // float pT_neg = neg->Pt();

        float pT_inv_pos = 1. / pT_pos;

        float px_pos = pos->Px();
        float py_pos = pos->Py();
        float pz_pos = pos->Pz();

        float px_neg = neg->Px();
        float py_neg = neg->Py();
        float pz_neg = neg->Pz();

        float e_pos = std::sqrt(kPiMass * kPiMass + px_pos * px_pos + py_pos * py_pos + pz_pos * pz_pos);
        float e_neg = std::sqrt(kPiMass * kPiMass + px_neg * px_neg + py_neg * py_neg + pz_neg * pz_neg);

        float px = px_pos + px_neg;
        float py = py_pos + py_neg;
        float pz = pz_pos + pz_neg;

        float p2 = px * px + py * py + pz * pz;

        float alpha = ( (px_pos - px_neg) * px + (py_pos - py_neg) * py + (pz_pos - pz_neg) * pz ) / ( (px_pos + px_neg) * px + (py_pos + py_neg) * py + (pz_pos + pz_neg) * pz );
        float mass = std::sqrt( std::pow(e_pos + e_neg, 2.) - p2);

        if (std::abs(alpha) > 0.05) continue;

        // hPtPosNeg.Fill(pT_pos, pT_neg);
        // hArm.Fill(alpha, qt);
        // hMass.Fill(pT_inv_pos, mass);
        hMassPt.Fill(pT_pos, mass);
        // hEtaReso.Fill(pT_pos, pos->Eta() - pos_eta_in);
      }

      myFile->cd();
      // hPtPosNeg.Write();
      // hArm.Write();
      hMassPt.Write();
      // hEtaReso.Write();
    }
  }

  myFile->Close();
}
