#include <TF1.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>
#include <TLorentzVector.h> // ugly, but TGenPhaseSpace uses this.
#include <TMath.h>
#include <TRandom.h>

constexpr uint64_t kNtrials{6000000};
constexpr double kK0sMass = 0.497648;
constexpr double kPiMass = 0.1395703918;

void toyMC() {

  TFile *myFile = TFile::Open("templates2.root", "RECREATE");
  gRandom->SetSeed(1234);

  // lorentz vectors to save particles
  TLorentzVector mother, pip, pim;
  // GenPhaseSpace to generate decay
  TGenPhaseSpace gen2Pi;
  //masses of resulting particles
  const double massesDau[2]{kPiMass, kPiMass};

  // simulate decay
  // loop over number of trials
  TF1 mtExpo("mtExpo","[0]*x*std::exp(-std::hypot([2], x)/[1])", 0.1, 30.);
  mtExpo.SetParameter(0, 1.);
  mtExpo.SetParameter(1, 0.5);
  mtExpo.SetParameter(2, kK0sMass);

  // double delta[100], sigma[100];

  double shift, sigma;

  for (int iD{0}; iD < 80; ++iD) {
    shift = -0.01 + iD * 0.00025;
    for (int iS{0}; iS < 50; ++iS) {
      sigma = 0.0 + iS * 0.001;

      std::cout << shift << "\t" << sigma << std::endl;

      // TH2D hPtPosNeg("hPrPtTmp", ";#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", 200, 0.f, 10.f, 200, 0.f, 10.f);
      // TH2D hArm("hArm", ";;", 200, -1., 1., 100., 0., 0.5);
      TH2D hMass(Form("hMass_%d_%d", iD, iS), ";1 / #it{p}_{T} (#it{c}/GeV);#it{M} (GeV/#it{c}^{2})", 100, 0., 10., 1200, 0, 0.6);
      TH2D hMassPt(Form("hMass_%d_%d_vs_pt", iD, iS), ";#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2})", 100, 0., 10., 1200, 0, 0.6);
      for (uint64_t i = 0; i < kNtrials; ++i) {
        // initialise variables
        float pT_k0s = mtExpo.GetRandom();
        float eta = gRandom->Uniform(-0.3, 0.3);
        float phi = gRandom->Uniform(0, TMath::TwoPi());

        // generate mother and decay
        mother.SetPtEtaPhiM(pT_k0s, eta, phi, kK0sMass);

        gen2Pi.SetDecay(mother, 2, massesDau);
        gen2Pi.Generate();

        auto pos = gen2Pi.GetDecay(0);
        auto neg = gen2Pi.GetDecay(1);

        float pos_ptinv = 1. / pos->Pt() + gRandom->Gaus(shift, sigma);
        float neg_ptinv = 1. / neg->Pt() + gRandom->Gaus(shift, sigma);

        pos->SetPtEtaPhiM(1. / pos_ptinv, pos->Eta(), pos->Phi(), kPiMass);
        neg->SetPtEtaPhiM(1. / neg_ptinv, neg->Eta(), neg->Phi(), kPiMass);

        float pT_pos = pos->Pt();
        float pT_neg = neg->Pt();

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

        float eta_rec = 0.5 * std::log( (std::sqrt(p2) + pz) / (std::sqrt(p2) - pz) );
        float qt = std::sqrt( std::pow(py_pos * pz - pz_pos * py, 2.) + std::pow(pz_pos * px - px_pos * pz, 2.) + std::pow(px_pos * py - py_pos * px, 2.)) / std::sqrt( px * px + py * py + pz * pz );
        float alpha = ( (px_pos - px_neg) * px + (py_pos - py_neg) * py + (pz_pos - pz_neg) * pz ) / ( (px_pos + px_neg) * px + (py_pos + py_neg) * py + (pz_pos + pz_neg) * pz );
        float mass = std::sqrt( std::pow(e_pos + e_neg, 2.) - p2);

        if (std::abs(alpha) > 0.05) continue;
        //if (pT_pos < 0.2 || pT_pos > 0.1) continue;
        //if (pT_neg < 0.2 || pT_neg > 0.1) continue;

        //hPtPosNeg.Fill(pT_pos, pT_neg);
        //hArm.Fill(alpha, qt);
        hMass.Fill(pT_inv_pos, mass);
        hMassPt.Fill(pT_pos, mass);
      }

      // create directon for this pT value
      myFile->cd();
      //hPtPosNeg.Write();
      //hArm.Write();
      hMass.Write();
    }
  }

  myFile->Close();
}
