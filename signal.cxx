#include <iostream>
#include <string>

#include <TFile.h>
#include <TDirectoryFile.h>
#include <ROOT/RDataFrame.hxx>

constexpr std::string _tname = "O2mcv0tableap";
const std::string _fname = "AO2D_mc_merged";
constexpr double _pimass= 0.1395703918; // GeV/c/c

void signal(){
  ROOT::EnableImplicitMT();

  std::string fname = _fname + ".root";
  TFile *inf = TFile::Open(fname.data());
  if (!inf) {
    std::cout << "No file \"" << fname.data() << "\" found!" << std::endl;
    return;
  }

  auto lKey = inf->GetListOfKeys();
  std::string kname = "";
  for (int iK{0}; iK < inf->GetNkeys(); ++iK) {
    std::string tmp_kname = lKey->At(iK)->GetName();
    auto id = tmp_kname.find("DF");
    if (id < tmp_kname.size()) {
      kname = tmp_kname;
      break;
    }
  }

  TDirectoryFile *d = static_cast<TDirectoryFile*>(inf->Get(kname.data()));
  if (!d) {
    std::cout << "Data directory \"" << kname.data() << "\" not found!" << std::endl;
    return;
  }

  kname = kname + "/" + _tname.data();
  ROOT::RDataFrame df(kname.data(), fname.data());
  auto df_rec = df.Filter("fIsReco")
                  .Define("fPx", "fPxPos + fPxNeg")
                  .Define("fPy", "fPyPos + fPyNeg")
                  .Define("fPz", "fPzPos + fPzNeg")
                  .Define("fP2", "fPx * fPx + fPy * fPy + fPz * fPz")
  		  .Define("fPtPos", "std::hypot(fPxPos, fPyPos)")
                  .Define("fPtNeg", "std::hypot(fPxNeg, fPyNeg)")
                  .Define("fP2Pos", "fPxPos * fPxPos + fPyPos * fPyPos + fPzPos * fPzPos")
                  .Define("fP2Neg", "fPxNeg * fPxNeg + fPyNeg * fPyNeg + fPzNeg * fPzNeg")
                  .Define("fEPos", "std::sqrt(_pimass * _pimass + fP2Pos)")
                  .Define("fENeg", "std::sqrt(_pimass * _pimass + fP2Neg)")
                  .Define("fK0sPt", "std::hypot(fPx, fPy)")
                  .Define("fMass", "std::sqrt(std::pow(fEPos + fENeg, 2.) - fP2Pos - fP2Neg - 2. * (fPxPos * fPxNeg + fPyPos * fPyNeg + fPzPos * fPzNeg))")
                  .Define("fAlpha", "((fPxPos - fPxNeg) * fPx + (fPyPos - fPyNeg) * fPy + (fPzPos - fPzNeg) * fPz) / ((fPxPos + fPxNeg) * fPx + (fPyPos + fPyNeg) * fPy + (fPzPos + fPzNeg) * fPz)")
                  .Define("fQt", "std::sqrt( std::pow(fPyPos * fPz - fPzPos * fPy, 2.) + std::pow(fPzPos * fPx - fPxPos * fPz, 2.) + std::pow(fPxPos * fPy - fPyPos * fPx, 2.) ) / std::sqrt(fP2)")
                  // .Filter("std::hypot(fPxPos, fPyPos) > 0.3 && std::hypot(fPxPos, fPyPos) < 0.4 && std::hypot(fPxNeg, fPyNeg) > 0.3 && std::hypot(fPxNeg, fPyNeg) < 0.4")
                  .Filter("fMass > 0.46 && fMass < 0.53 && std::abs(fAlpha) < 0.05 && std::abs(fEta) < 0.3")
                  .Define("RndmId", "fPx - fPx + gRandom->Integer(20)")
                  .Define("fPtInv", "1. / fPtPos")
                  .Define("fPtPosMC", "std::hypot(fPxPosMC, fPyPosMC)")
                  .Define("fPtReso", "fPtPos - fPtPosMC")
                  .Define("fPtInvReso", "1. / fPtPos - 1. / fPtPosMC")
                  .Define("fPhiReso", "std::atan(fPyPos / fPxPos) - std::atan(fPyPosMC / fPxPosMC)")
                  .Define("fEtaReso", "-0.5 * std::log( std::tan( 0.5 * std::atan(fPtPos / fPzPos) ) ) + 0.5 * std::log( std::tan( 0.5 * std::atan(std::hypot(fPxPosMC, fPyPosMC) / fPzPosMC) ) )");
  auto h = df_rec.Histo3D({"hMassK0s", ";Subsample;#frac{1}{#it{p}_{T}^{rec.}} (#it{c}/GeV);#it{M} (GeV/#it{c}^{2})", 20, 0., 20., 100, 0., 10., 1200, 0., 0.6}, "RndmId", "fPtInv", "fMass");
  auto h1 = df_rec.Histo3D({"hMassK0sPt", ";Subsample;#it{p}_{T} (GeV/#it{c});#it{M} (GeV/#it{c}^{2})", 20, 0., 20., 100, 0., 10., 1200, 0., 0.6}, "RndmId", "fPtPos", "fMass");
  auto h2 = df_rec.Histo2D({"hArmenteros", "hArmenteros", 100, -1, 1, 100, 0, 0.5}, "fAlpha", "fQt");
  auto h3 = df_rec.Histo2D({"hPPosNeg", "hPPosNeg", 100, 0, 4, 100, 0, 4}, "fPtPos", "fPtNeg");
  auto h4 = df_rec.Histo2D({"hPtInvReso", ";#frac{1}{#it{p}_{T}^{rec.}} (#it{c}/GeV);#frac{1}{#it{p}_{T}^{rec.}} - #frac{1}{#it{p}_{T}^{gen.}} (#it{c}/GeV)", 100, 0., 10., 500, -0.2, 0.2}, "fPtInv", "fPtInvReso");
  auto h5 = df_rec.Histo2D({"hPtReso", ";#it{p}_{T}^{rec.} (GeV/#it{c});#it{p}_{T}^{rec.} - #it{p}_{T}^{gen.} (GeV/#it{c})", 100, 0., 10., 5000, -1., 1.}, "fPtPos", "fPtReso");
  auto h6 = df_rec.Histo2D({"hPhiReso", ";#it{p}_{T}^{rec.} (GeV/#it{c});#phi^{rec.} - #phi^{gen.}", 100, 0., 10., 10000, -1., 1.}, "fPtPosMC", "fPhiReso");
  auto h7 = df_rec.Histo2D({"hEtaReso", ";#it{p}_{T}^{rec.} (GeV/#it{c});#eta^{rec.} - #eta^{gen.}", 100, 0., 10., 10000, -.5, .5}, "fPtPosMC", "fEtaReso");
  auto h8 = df_rec.Histo2D({"hEtaPtReso", ";#eta^{rec.} - #eta^{gen.};#it{p}_{T}^{rec.} - #it{p}_{T}^{gen.} (GeV/#it{c})", 10000, -.5, .5, 5000, -1., 1.}, "fEtaReso", "fPtReso");

  auto h9 = df_rec.Histo2D({"hMassK0sPtK0s", ";#it{p}_{T}^{K0s} (GeV/#it{c});#it{M} (GeV/#it{c}^{2})", 100, 0., 10., 1200, 0., 0.6}, "fK0sPt", "fMass");

  TFile *fo = TFile::Open("foo_dat_2.root", "recreate");
  fo->cd();
  h->Write();
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  h7->Write();
  h8->Write();

  h9->Write();

  fo->Close();
  inf->Close();
}
