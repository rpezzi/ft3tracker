#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "FT3Track.h"
#include "Field/MagneticField.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>

#endif

// FT3Tools

#include "ft3tools/HistosHelpers.C"
#include "ft3tools/MagField.C"

using o2::MCTrackT;
using o2::ft3::FT3Track;
using o2::itsmft::Hit;
using eventFoundTracks = std::vector<bool>;
using std::vector;

bool DEBUG_VERBOSE = false;
bool EXPORT_HISTOS_IMAGES = false;
Int_t minHitsPerTrack = 10;

//_________________________________________________________________________________________________
int FT3TrackerChecker(const Char_t *trkFile = "ft3tracks.root",
                      const Char_t *o2sim_KineFile = "o2sim_Kine.root",
                      const Char_t *HitsFT3File = "o2sim_HitsFT3.root") {

  // Histos parameters
  Double_t pMin = 0.0;
  Double_t pMax = 100.0;
  Double_t deltaetaMin = -.1;
  Double_t deltaetaMax = +.1;
  Double_t etaMin = -5.0;
  Double_t etaMax = +5.0;
  Double_t deltaphiMin = -.2;
  Double_t deltaphiMax = .2;
  Double_t deltatanlMin = -2.0;
  Double_t deltatanlMax = 2.0;

  // Seed configuration
  std::string seed_cfg{trkFile};
  std::string trk_start{"ft3tracks_"};
  std::string trk_ext{".root"};
  std::string trk_trk{"ft3tracks"};
  if (seed_cfg.find(trk_start) < seed_cfg.length())
    seed_cfg.replace(seed_cfg.find(trk_start), trk_start.length(), "");
  if (seed_cfg.find(trk_ext) < seed_cfg.length())
    seed_cfg.replace(seed_cfg.find(trk_ext), trk_ext.length(), "");
  if (seed_cfg.find(trk_trk) < seed_cfg.length())
    seed_cfg.replace(seed_cfg.find(trk_trk), trk_trk.length(), "");
  std::cout << seed_cfg << std::endl;

  // histos
  gStyle->SetOptStat("emr");
  gStyle->SetStatW(.28);
  gStyle->SetStatH(.26);
  gStyle->SetPalette(1, 0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(4);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(3);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(3);
  gStyle->SetLabelSize(0.06, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  gStyle->SetTitleSize(0.06, "xyz");
  gStyle->SetTitleSize(0.08, "o");
  gStyle->SetTitleOffset(.95, "Y");
  gStyle->SetTitleFillColor(10);
  gStyle->SetStatColor(10);

  enum TH2HistosCodes {
    kFT3TrackDeltaXYVertex,
    kFT3TrackDeltaXYVertex0_1,
    kFT3TrackDeltaXYVertex1_4,
    kFT3TrackDeltaXYVertex4plus,
    kFT3rackQPRec_MC,
    kFT3rackPtResolution,
    kFT3rackInvPtResolution,
    kMCTracksEtaZ
  };

  std::map<int, const char *> TH2Names{
      {kFT3TrackDeltaXYVertex, "FT3 Tracks at vertex"},
      {kFT3TrackDeltaXYVertex0_1, "FT3 Tracks Vertex at Z = 0 Pt0_1"},
      {kFT3TrackDeltaXYVertex1_4, "FT3 Tracks Vertex at Z = 0 Pt1_4"},
      {kFT3TrackDeltaXYVertex4plus, "FT3 Tracks Vertex at Z = 0 Pt4plus"},
      {kFT3rackQPRec_MC, "FT3 Track QP FITxMC"},
      {kFT3rackPtResolution, "FT3 Track Pt Resolution"},
      {kFT3rackInvPtResolution, "FT3 Track InvPt Resolution"},
      {kMCTracksEtaZ, "MCTracks_eta_z"}};

  std::map<int, const char *> TH2Titles{
      {kFT3TrackDeltaXYVertex, "FT3 Tracks at Z_vertex"},
      {kFT3TrackDeltaXYVertex0_1, "FT3 Tracks at Z_vertex (pt < 1)"},
      {kFT3TrackDeltaXYVertex1_4, "FT3 Tracks at Z_vertex (1 < pt < 4)"},
      {kFT3TrackDeltaXYVertex4plus, "FT3 Tracks at Z_vertex (pt > 4)"},
      {kFT3rackQPRec_MC, "Charged Momentum: Reconstructed vs MC"},
      {kFT3rackPtResolution, "Pt Resolution"},
      {kFT3rackInvPtResolution, "InvPt Resolution"},
      {kMCTracksEtaZ, "MC Tracks: Pseudorapidity vs zVertex"}};

  std::map<int, std::array<double, 6>> TH2Binning{
      {kFT3TrackDeltaXYVertex, {100, -.5, .5, 100, -.5, .5}},
      {kFT3TrackDeltaXYVertex0_1, {100, -.5, .5, 100, -.5, .5}},
      {kFT3TrackDeltaXYVertex1_4, {100, -.5, .5, 100, -.5, .5}},
      {kFT3TrackDeltaXYVertex4plus, {100, -.5, .5, 100, -.5, .5}},
      {kFT3rackQPRec_MC, {100, -10, 10, 100, -10, 10}},
      {kFT3rackPtResolution, {14, 0, 7, 250, 0, 25}},
      {kFT3rackInvPtResolution, {14, 0, 7, 300, -2, 2}},
      {kMCTracksEtaZ, {31, -15, 16, 25, etaMin, etaMax}}};

  std::map<int, const char *> TH2XaxisTitles{
      {kFT3TrackDeltaXYVertex, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaXYVertex0_1, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaXYVertex1_4, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaXYVertex4plus, "\\Delta x ~[mm]"},
      {kFT3rackQPRec_MC, "(q.p)_{MC} [GeV]"},
      {kFT3rackPtResolution, "pt_{MC} [GeV]"},
      {kFT3rackInvPtResolution, "pt_{MC} [GeV]"},
      {kMCTracksEtaZ, "Vertex PosZ [cm]"}};

  std::map<int, const char *> TH2YaxisTitles{
      {kFT3TrackDeltaXYVertex, "\\Delta y ~[mm]"},
      {kFT3TrackDeltaXYVertex0_1, "\\Delta y ~[mm]"},
      {kFT3TrackDeltaXYVertex1_4, "\\Delta y ~[mm]"},
      {kFT3TrackDeltaXYVertex4plus, "\\Delta y ~[mm]"},
      {kFT3rackQPRec_MC, "(q.p)_{fit} [GeV]"},
      {kFT3rackPtResolution, "pt_{fit} / pt_{MC}"},
      {kFT3rackInvPtResolution, "(1/(p_t)_{fit} - 1/(p_t)_{MC})*(p_t)_{MC}"},
      {kMCTracksEtaZ, "\\eta"}};

  enum TH1HistosCodes {
    kFT3TrackDeltaXErr,
    kFT3TrackDeltaYErr,
    kFT3TrackDeltaPhiErr,
    kFT3TrackDeltaTanLErr,
    kFT3TrackDeltainvQPtErr,
    kFT3TracksP,
    kFT3TrackDeltaTanl,
    kFT3TrackDeltaTanl0_1,
    kFT3TrackDeltaTanl1_4,
    kFT3TrackDeltaTanl4plus,
    kFT3TrackDeltaPhi,
    kFT3TrackDeltaPhi0_1,
    kFT3TrackDeltaPhi1_4,
    kFT3TrackDeltaPhi4plus,
    kFT3TrackDeltaPhiDeg,
    kFT3TrackDeltaPhiDeg0_1,
    kFT3TrackDeltaPhiDeg1_4,
    kFT3TrackDeltaPhiDeg4plus,
    kFT3TrackDeltaInvQPt,
    kFT3TrackDeltaInvQPtSeed,
    kFT3TrackDeltaX,
    kFT3TrackDeltaX0_1,
    kFT3TrackDeltaX1_4,
    kFT3TrackDeltaX4plus,
    kFT3TrackDeltaY,
    kFT3TrackR,
    kFT3TrackQ,
    kFT3TrackQ0_1,
    kFT3TrackQ1_4,
    kFT3TrackQ4plus,
    kFT3TrackChi2,
    kMCTrackspT,
    kMCTracksp,
    kMCTrackEta
  };

  std::map<int, const char *> TH1Names{
      {kFT3TracksP, "FT3 Tracks Fitted p"},
      {kFT3TrackDeltaXErr, "Delta X / SigmaX"},
      {kFT3TrackDeltaYErr, "Delta Y / SigmaY"},
      {kFT3TrackDeltaPhiErr, "Delta Phi at Vertex / SigmaPhi"},
      {kFT3TrackDeltaTanLErr, "Delta_Tanl / SigmaTanl"},
      {kFT3TrackDeltainvQPtErr, "Delta_InvQPt / Sigma_{q/pt}"},
      {kFT3TrackDeltaTanl, "FT3 Tracks Fitted Delta_tanl"},
      {kFT3TrackDeltaTanl0_1, "FT3 Tracks tanl (pt < 1)"},
      {kFT3TrackDeltaTanl1_4, "FT3 Tracks tanl (1 < pt < 4)"},
      {kFT3TrackDeltaTanl4plus, "FT3 Tracks tanl (pt > 4)"},
      {kFT3TrackDeltaPhi, "FT3 Tracks Fitted Phi at Vertex"},
      {kFT3TrackDeltaPhi0_1, "FT3 Tracks Fitted Phi at Vertex [rad] (pt < 1)"},
      {kFT3TrackDeltaPhi1_4,
       "FT3 Tracks Fitted Phi at Vertex [rad] (1 < pt < 4)"},
      {kFT3TrackDeltaPhi4plus,
       "FT3 Tracks Fitted Phi at Vertex [rad] (pt > 4)"},
      {kFT3TrackDeltaPhiDeg, "FT3 Tracks Fitted Phi at Vertex [deg]"},
      {kFT3TrackDeltaPhiDeg0_1,
       "FT3 Tracks Fitted Phi at Vertex [deg] (pt < 1)"},
      {kFT3TrackDeltaPhiDeg1_4,
       "FT3 Tracks Fitted Phi at Vertex [deg] (1 < pt < 4)"},
      {kFT3TrackDeltaPhiDeg4plus,
       "FT3 Tracks Fitted Phi at Vertex [deg] (pt > 4)"},
      {kFT3TrackDeltaInvQPt, "FT3 Tracks invQPt"},
      {kFT3TrackDeltaInvQPtSeed, "FT3 Tracks invQPt Seed"},
      {kFT3TrackDeltaX, "FT3 Tracks Delta X"},
      {kFT3TrackDeltaX0_1, "FT3 Tracks Delta X (pt < 1)"},
      {kFT3TrackDeltaX1_4, "FT3 Tracks Delta X (1 < pt < 4)"},
      {kFT3TrackDeltaX4plus, "FT3 Tracks Delta X (pt > 4)"},
      {kFT3TrackDeltaY, "FT3 Tracks Delta Y"},
      {kFT3TrackR, "FT3 Tracks Delta R"},
      {kFT3TrackQ, "Charge Match"},
      {kFT3TrackQ0_1, "Charge Match (pt < 1)"},
      {kFT3TrackQ1_4, "Charge Match (1 < pt < 4)"},
      {kFT3TrackQ4plus, "Charge Match (pt > 4)"},
      {kFT3TrackChi2, "Tracks Chi2"},
      {kMCTrackspT, "MC Tracks p_T"},
      {kMCTracksp, "MC Tracks p"},
      {kMCTrackEta, "MC Tracks eta"}};

  std::map<int, const char *> TH1Titles{
      {kFT3TracksP, "Standalone FT3 Tracks P"},
      {kFT3TrackDeltaXErr, "\\Delta X / \\sigma_X"},
      {kFT3TrackDeltaYErr, "\\Delta Y / \\sigma_Y"},
      {kFT3TrackDeltaPhiErr, "\\Delta \\phi / \\sigma_\\phi"},
      {kFT3TrackDeltaTanLErr, "\\Delta TanL / \\sigma_{TanL} "},
      {kFT3TrackDeltainvQPtErr, "\\Delta(q/Pt) / \\sigma_{q/pt}"},
      {kFT3TrackDeltaTanl, "tanl_{Fit} - tanl_{MC} "},
      {kFT3TrackDeltaTanl0_1, "tanl_{Fit} - tanl_{MC} (pt < 1)"},
      {kFT3TrackDeltaTanl1_4, "tanl_{Fit} - tanl_{MC} (1 < p_t < 4)"},
      {kFT3TrackDeltaTanl4plus, "tanl_{Fit} - tanl_{MC} (p_t > 4)"},
      {kFT3TrackDeltaPhi, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhi0_1, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhi1_4, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhi4plus, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhiDeg, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhiDeg0_1, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhiDeg1_4, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaPhiDeg4plus, "\\phi _{Fit} - \\phi_{MC}"},
      {kFT3TrackDeltaInvQPt, "FT3 Tracks \\Delta invQPt"},
      {kFT3TrackDeltaInvQPtSeed, "FT3 Tracks \\Delta invQPt Seed"},
      {kFT3TrackDeltaX, "FT3 Tracks Delta X at Z_vertex"},
      {kFT3TrackDeltaX0_1, "FT3 Tracks Delta X at Z_vertex"},
      {kFT3TrackDeltaX1_4, "FT3 Tracks Delta X at Z_vertex"},
      {kFT3TrackDeltaX4plus, "FT3 Tracks Delta X at Z_vertex"},
      {kFT3TrackDeltaY, "FT3 Tracks Delta Y at Z_vertex"},
      {kFT3TrackR, "FT3 Tracks Delta R at Z_vertex"},
      {kFT3TrackQ, "FT3 Tracks Charge Match"},
      {kFT3TrackQ0_1, "FT3 Tracks Charge Match (pt < 1)"},
      {kFT3TrackQ1_4, "FT3 Tracks Charge Match (1 < pt < 4)"},
      {kFT3TrackQ4plus, "FT3 Tracks Charge Match (pt > 4)"},
      {kFT3TrackChi2, "FT3 Tracks ~ \\chi^2"},
      {kMCTrackspT, "MC Tracks p_T"},
      {kMCTracksp, "MC Tracks p"},
      {kMCTrackEta, "MC Tracks Pseudorapidity"}};

  std::map<int, std::array<double, 3>> TH1Binning{
      {kFT3TracksP, {500, pMin, pMax}},
      {kFT3TrackDeltaXErr, {500, -10, 10}},
      {kFT3TrackDeltaYErr, {500, -10, 10}},
      {kFT3TrackDeltaPhiErr, {500, -10, +10}},
      {kFT3TrackDeltaTanLErr, {500, -10, +10}},
      {kFT3TrackDeltainvQPtErr, {500, -50, +50}},
      {kFT3TrackDeltaTanl, {1000, deltatanlMin, deltatanlMax}},
      {kFT3TrackDeltaTanl0_1, {1000, deltatanlMin, deltatanlMax}},
      {kFT3TrackDeltaTanl1_4, {1000, deltatanlMin, deltatanlMax}},
      {kFT3TrackDeltaTanl4plus, {1000, deltatanlMin, deltatanlMax}},
      {kFT3TrackDeltaPhi, {1000, deltaphiMin, deltaphiMax}},
      {kFT3TrackDeltaPhi0_1, {1000, deltaphiMin, deltaphiMax}},
      {kFT3TrackDeltaPhi1_4, {1000, deltaphiMin, deltaphiMax}},
      {kFT3TrackDeltaPhi4plus, {1000, deltaphiMin, deltaphiMax}},
      {kFT3TrackDeltaPhiDeg,
       {1000, TMath::RadToDeg() * deltaphiMin,
        TMath::RadToDeg() * deltaphiMax}},
      {kFT3TrackDeltaPhiDeg0_1,
       {1000, TMath::RadToDeg() * deltaphiMin,
        TMath::RadToDeg() * deltaphiMax}},
      {kFT3TrackDeltaPhiDeg1_4,
       {1000, TMath::RadToDeg() * deltaphiMin,
        TMath::RadToDeg() * deltaphiMax}},
      {kFT3TrackDeltaPhiDeg4plus,
       {1000, TMath::RadToDeg() * deltaphiMin,
        TMath::RadToDeg() * deltaphiMax}},
      {kFT3TrackDeltaInvQPt, {1000, -10., 10.}},
      {kFT3TrackDeltaInvQPtSeed, {1000, -0.5, 0.5}},
      {kFT3TrackDeltaX, {1000, -5, 5}},
      {kFT3TrackDeltaX0_1, {1000, -5, 5}},
      {kFT3TrackDeltaX1_4, {1000, -5, 5}},
      {kFT3TrackDeltaX4plus, {1000, -5, 5}},
      {kFT3TrackDeltaY, {1000, -5, 5}},
      {kFT3TrackR, {250, 0, 5}},
      {kFT3TrackQ, {5, -2.1, 2.1}},
      {kFT3TrackQ0_1, {5, -2.1, 2.1}},
      {kFT3TrackQ1_4, {5, -2.1, 2.1}},
      {kFT3TrackQ4plus, {5, -2.1, 2.1}},
      {kFT3TrackChi2, {10000, 0, 1000}},
      {kMCTrackspT, {5000, 0, 50}},
      {kMCTracksp, {1000, pMin, pMax}},
      {kMCTrackEta, {1000, etaMin, etaMax}}};

  std::map<int, const char *> TH1XaxisTitles{
      {kFT3TracksP, "p [GeV]"},
      {kFT3TrackDeltaXErr, "\\Delta x  /\\sigma_{x}"},
      {kFT3TrackDeltaYErr, "\\Delta y  /\\sigma_{y}"},
      {kFT3TrackDeltaPhiErr, "\\Delta \\phi  /\\sigma_{\\phi}"},
      {kFT3TrackDeltaTanLErr, "\\Delta tanl /\\sigma_{tanl}"},
      {kFT3TrackDeltainvQPtErr, "\\Delta (q/p_t)/\\sigma_{q/Pt}"},
      {kFT3TrackDeltaTanl, "\\Delta tanl"},
      {kFT3TrackDeltaTanl0_1, "\\Delta tanl"},
      {kFT3TrackDeltaTanl1_4, "\\Delta tanl"},
      {kFT3TrackDeltaTanl4plus, "\\Delta tanl"},
      {kFT3TrackDeltaPhi, "\\Delta \\phi ~[rad]"},
      {kFT3TrackDeltaPhi0_1, "\\Delta \\phi ~[rad]"},
      {kFT3TrackDeltaPhi1_4, "\\Delta \\phi ~[rad]"},
      {kFT3TrackDeltaPhi4plus, "\\Delta \\phi ~[rad]"},
      {kFT3TrackDeltaPhiDeg, "\\Delta \\phi ~[deg]"},
      {kFT3TrackDeltaPhiDeg0_1, "\\Delta \\phi ~[deg]"},
      {kFT3TrackDeltaPhiDeg1_4, "\\Delta \\phi ~[deg]"},
      {kFT3TrackDeltaPhiDeg4plus, "\\Delta \\phi ~[deg]"},
      {kFT3TrackDeltaInvQPt, "\\Delta invQPt"},
      {kFT3TrackDeltaInvQPtSeed, "\\Delta invQPt Seed"},
      {kFT3TrackDeltaX, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaX0_1, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaX1_4, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaX4plus, "\\Delta x ~[mm]"},
      {kFT3TrackDeltaY, "\\Delta y ~[mm]"},
      {kFT3TrackR, "\\Delta r ~[mm]"},
      {kFT3TrackQ, "q_{fit}-q_{MC}"},
      {kFT3TrackQ0_1, "q_{fit}-q_{MC}"},
      {kFT3TrackQ1_4, "q_{fit}-q_{MC}"},
      {kFT3TrackQ4plus, "q_{fit}-q_{MC}"},
      {kFT3TrackChi2, "\\chi^2"},
      {kMCTrackspT, "p_t [GeV]"},
      {kMCTracksp, "p [GeV]"},
      {kMCTrackEta, " \\eta"}};

  // Create histograms
  const int nTH1Histos = TH1Names.size();
  std::vector<std::unique_ptr<TH1F>> TH1Histos(nTH1Histos);
  auto nHisto = 0;
  for (auto &h : TH1Histos) {
    h = std::make_unique<TH1F>(TH1Names[nHisto], TH1Titles[nHisto],
                               (int)TH1Binning[nHisto][0],
                               TH1Binning[nHisto][1], TH1Binning[nHisto][2]);
    h->GetXaxis()->SetTitle(TH1XaxisTitles[nHisto]);
    ++nHisto;
  }

  const int nTH2Histos = TH2Names.size();
  std::vector<std::unique_ptr<TH2F>> TH2Histos(nTH2Histos);
  auto n2Histo = 0;
  for (auto &h : TH2Histos) {
    h = std::make_unique<TH2F>(TH2Names[n2Histo], TH2Titles[n2Histo],
                               (int)TH2Binning[n2Histo][0],
                               TH2Binning[n2Histo][1], TH2Binning[n2Histo][2],
                               (int)TH2Binning[n2Histo][3],
                               TH2Binning[n2Histo][4], TH2Binning[n2Histo][5]);
    h->GetXaxis()->SetTitle(TH2XaxisTitles[n2Histo]);
    h->GetYaxis()->SetTitle(TH2YaxisTitles[n2Histo]);

    h->SetOption("COLZ");
    ++n2Histo;
  }

  // Profiles histograms
  auto PtRes_Profile = new TProfile("Pt_res_prof", "Profile of pt{fit}/pt{MC}",
                                    14, 0, 7, 0, 20, "s");
  PtRes_Profile->GetXaxis()->SetTitle("pt_{MC}");
  PtRes_Profile->GetYaxis()->SetTitle("mean(Pt_{Fit}/Pt_{MC})");

  auto DeltaX_Profile = new TProfile("DeltaX_prof", "Vertexing resolution", 14,
                                     0, 7, -10000., 10000., "s");
  DeltaX_Profile->GetXaxis()->SetTitle("pt_{MC} [GeV]");
  DeltaX_Profile->GetYaxis()->SetTitle("\\sigma_x ~[\\mu m]");

  // TEfficiency histogram
  TEfficiency *qMatchEff = new TEfficiency(
      "QMatchEff", "Charge Match;p_t [GeV];#epsilon", 10, 0, 10);

  // Counters
  Int_t nChargeMatch = 0;
  Int_t nChargeMiss = 0;
  Int_t nChargeMatch0_1 = 0;
  Int_t nChargeMiss0_1 = 0;
  Int_t nChargeMatch1_4 = 0;
  Int_t nChargeMiss1_4 = 0;
  Int_t nChargeMatch4plus = 0;
  Int_t nChargeMiss4plus = 0;

  // Files & Trees
  // MC
  TFile *o2sim_KineFileIn = new TFile(o2sim_KineFile);
  TTree *o2SimKineTree = (TTree *)o2sim_KineFileIn->Get("o2sim");

  vector<MCTrackT<float>> *mcTr = nullptr;
  o2SimKineTree->SetBranchAddress("MCTrack", &mcTr);
  o2::dataformats::MCEventHeader *eventHeader = nullptr;
  o2SimKineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  // FT3 Hits
  TFile *HitsFT3FileIn = new TFile(HitsFT3File);
  TTree *o2FT3HitsTree = (TTree *)HitsFT3FileIn->Get("o2sim");
  vector<Hit> *ft3hit = nullptr;
  o2FT3HitsTree->SetBranchAddress("FT3Hit", &ft3hit);

  Int_t numberOfEvents = o2SimKineTree->GetEntries();
  Int_t numberOfFT3Events = o2FT3HitsTree->GetEntries();
  if (numberOfEvents == numberOfFT3Events)
    std::cout << "numberOfEvents = " << numberOfEvents << std::endl;
  else {
    std::cout << "ERROR: Inconsistent number of entries on " << o2sim_KineFile
              << " and " << HitsFT3File << std::endl;
    return -1;
  }

  // FT3 Tracks
  TFile *trkFileIn = new TFile(trkFile);
  TTree *ft3TrackTree = (TTree *)trkFileIn->Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackFT3Vec, *trackFT3VecP = &trackFT3Vec;
  ft3TrackTree->SetBranchAddress("FT3Track", &trackFT3VecP);

  vector<Int_t> *recoTrackIDs = nullptr;
  ft3TrackTree->SetBranchAddress("FT3TrackID", &recoTrackIDs);

  ft3TrackTree->GetEntry(0);
  o2SimKineTree->GetEntry(0);

  auto field_z = getZField(0, 0, 0); // Get field at Center of ALICE

  std::string outfilename = "Fittercheck_" + std::string(trkFile);

  TFile outFile(outfilename.c_str(), "RECREATE");

  // Reconstructed FT3 Tracks
  std::cout << "Loop over events and reconstructed FT3 Tracks!" << std::endl;
  // TracksFT3 - Identify reconstructed tracks
  auto totalTracks = 0;
  for (int iEvent = 0; iEvent < numberOfEvents; iEvent++) {
    ft3TrackTree->GetEntry(iEvent);
    o2SimKineTree->GetEntry(iEvent);
    auto iTrack = 0;
    if (DEBUG_VERBOSE)
      std::cout << "Event = " << iEvent << std::endl;
    o2SimKineTree->GetEntry(iEvent);
    for (auto &trackFT3 : trackFT3Vec) {
      auto trackID = recoTrackIDs->at(iTrack);
      if (trackFT3.getNumberOfPoints() < minHitsPerTrack) {
        iTrack++;
        continue;
      }
      if (1) {
        if (DEBUG_VERBOSE) {

          std::cout << "  Track #" << iTrack << ": TrackID = " << trackID
                    << std::endl;
        }

        MCTrackT<float> *thisTrack = &(*mcTr).at(trackID);
        auto vx_MC = thisTrack->GetStartVertexCoordinatesX();
        auto vy_MC = thisTrack->GetStartVertexCoordinatesY();
        auto vz_MC = thisTrack->GetStartVertexCoordinatesZ();
        auto Pt_MC = thisTrack->GetPt();
        auto P_MC = thisTrack->GetP();
        auto phi_MC = TMath::ATan2(thisTrack->Py(), thisTrack->Px());
        auto eta_MC = atanh(thisTrack->GetStartVertexMomentumZ() / P_MC);
        auto tanl_MC = thisTrack->Pz() / thisTrack->GetPt();
        auto pdgcode_MC = thisTrack->GetPdgCode();

        int Q_MC;
        if (TDatabasePDG::Instance()->GetParticle(pdgcode_MC)) {
          Q_MC =
              TDatabasePDG::Instance()->GetParticle(pdgcode_MC)->Charge() / 3;
        }

        else {
          iTrack++;
          continue;
          Q_MC = 0;
          // std::cout << " => pdgcode ERROR " << Q_MC <<  "\n";
        }

        trackFT3.propagateToZhelix(vz_MC, field_z);

        auto Q_fit = trackFT3.getCharge();
        auto dx = trackFT3.getX() - vx_MC;
        auto dy = trackFT3.getY() - vy_MC;
        auto d_eta = trackFT3.getEta() - eta_MC;
        auto d_tanl = trackFT3.getTanl() - tanl_MC;
        auto Pt_fit = trackFT3.getPt();
        auto invQPt_seed = trackFT3.getInvQPtSeed();
        auto d_invQPt = Q_fit / Pt_fit - Q_MC / Pt_MC;
        auto d_invQPtSeed = invQPt_seed - Q_fit / Pt_fit;
        auto P_fit = trackFT3.getP();
        auto P_res = P_fit / P_MC;
        auto Pt_res = Pt_fit / Pt_MC;
        auto d_Phi = trackFT3.getPhi() - phi_MC;
        auto d_Charge = Q_fit - Q_MC;
        auto trackChi2 = trackFT3.getTrackChi2();
        TH1Histos[kFT3TracksP]->Fill(trackFT3.getP());
        TH1Histos[kFT3TrackDeltaTanl]->Fill(d_tanl);
        TH1Histos[kFT3TrackDeltaPhi]->Fill(d_Phi);
        TH1Histos[kFT3TrackDeltaInvQPt]->Fill(d_invQPt);
        TH1Histos[kFT3TrackDeltaInvQPtSeed]->Fill(d_invQPtSeed);
        TH1Histos[kFT3TrackDeltaPhiDeg]->Fill(TMath::RadToDeg() * d_Phi);
        TH1Histos[kFT3TrackDeltaX]->Fill(10. * dx);

        TH1Histos[kFT3TrackDeltaXErr]->Fill(
            dx / sqrt(trackFT3.getCovariances()(0, 0)));
        TH1Histos[kFT3TrackDeltaYErr]->Fill(
            dy / sqrt(trackFT3.getCovariances()(1, 1)));
        TH1Histos[kFT3TrackDeltaPhiErr]->Fill(
            d_Phi / sqrt(trackFT3.getCovariances()(2, 2)));
        TH1Histos[kFT3TrackDeltaTanLErr]->Fill(
            d_tanl / sqrt(trackFT3.getCovariances()(3, 3)));
        TH1Histos[kFT3TrackDeltainvQPtErr]->Fill(
            d_invQPt / sqrt(trackFT3.getCovariances()(4, 4)));

        DeltaX_Profile->Fill(Pt_MC, dx * 1e4);
        TH1Histos[kFT3TrackDeltaY]->Fill(10 * dy);
        TH1Histos[kFT3TrackR]->Fill(10.0 * sqrt(dx * dx + dy * dy));
        TH1Histos[kFT3TrackQ]->Fill(d_Charge);
        TH1Histos[kFT3TrackChi2]->Fill(trackChi2);
        TH2Histos[kFT3TrackDeltaXYVertex]->Fill(10.0 * dx, 10.0 * dy);
        TH2Histos[kFT3rackQPRec_MC]->Fill(P_MC * Q_MC, P_fit * Q_fit);
        TH2Histos[kFT3rackPtResolution]->Fill(Pt_MC, Pt_fit / Pt_MC);
        PtRes_Profile->Fill(Pt_MC, Pt_fit / Pt_MC);
        TH2Histos[kFT3rackInvPtResolution]->Fill(
            Pt_MC, (1.0 / Pt_fit - 1.0 / Pt_MC) * Pt_MC);

        // MC histos
        TH1Histos[kMCTrackspT]->Fill(Pt_MC);
        TH1Histos[kMCTracksp]->Fill(P_MC);
        TH1Histos[kMCTrackEta]->Fill(eta_MC);
        TH2Histos[kMCTracksEtaZ]->Fill(vz_MC, eta_MC);

        // Differential histos
        if (Pt_MC <= 1.0) {
          TH2Histos[kFT3TrackDeltaXYVertex0_1]->Fill(10.0 * dx, 10.0 * dy);
          TH1Histos[kFT3TrackDeltaTanl0_1]->Fill(d_tanl);
          TH1Histos[kFT3TrackDeltaPhi0_1]->Fill(d_Phi);
          TH1Histos[kFT3TrackDeltaPhiDeg0_1]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kFT3TrackDeltaX0_1]->Fill(10.0 * dx);
          TH1Histos[kFT3TrackQ0_1]->Fill(d_Charge);
          d_Charge ? nChargeMiss0_1++ : nChargeMatch0_1++;
        }
        if (Pt_MC > 1.0 and Pt_MC <= 4) {
          TH2Histos[kFT3TrackDeltaXYVertex1_4]->Fill(10.0 * dx, 10.0 * dy);
          TH1Histos[kFT3TrackDeltaTanl1_4]->Fill(d_tanl);
          TH1Histos[kFT3TrackDeltaPhi1_4]->Fill(d_Phi);
          TH1Histos[kFT3TrackDeltaPhiDeg1_4]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kFT3TrackDeltaX1_4]->Fill(10.0 * dx);
          TH1Histos[kFT3TrackQ1_4]->Fill(d_Charge);
          d_Charge ? nChargeMiss1_4++ : nChargeMatch1_4++;
        }
        if (Pt_MC > 4.0) {
          TH2Histos[kFT3TrackDeltaXYVertex4plus]->Fill(10.0 * dx, 10.0 * dy);
          TH1Histos[kFT3TrackDeltaTanl4plus]->Fill(d_tanl);
          TH1Histos[kFT3TrackDeltaPhi4plus]->Fill(d_Phi);
          TH1Histos[kFT3TrackDeltaPhiDeg4plus]->Fill(TMath::RadToDeg() * d_Phi);
          TH1Histos[kFT3TrackDeltaX4plus]->Fill(10.0 * dx);
          TH1Histos[kFT3TrackQ4plus]->Fill(d_Charge);
          d_Charge ? nChargeMiss4plus++ : nChargeMatch4plus++;
        }

        d_Charge ? nChargeMiss++ : nChargeMatch++;
        qMatchEff->Fill(!d_Charge, Pt_MC);
      }
      iTrack++;
    } // Loop on TracksFT3
    totalTracks += iTrack;
  } // Loop over events

  // Customize histograms
  TH1Histos[kFT3TrackQ]->SetTitle(
      Form("nChargeMatch = %d (%.2f%%)", nChargeMatch,
           100. * nChargeMatch / (nChargeMiss + nChargeMatch)));
  TH1Histos[kFT3TrackQ0_1]->SetTitle(
      Form("nChargeMatch = %d (%.2f%%)", nChargeMatch0_1,
           100. * nChargeMatch0_1 / (nChargeMiss0_1 + nChargeMatch0_1)));
  TH1Histos[kFT3TrackQ1_4]->SetTitle(
      Form("nChargeMatch = %d (%.2f%%)", nChargeMatch1_4,
           100. * nChargeMatch1_4 / (nChargeMiss1_4 + nChargeMatch1_4)));
  TH1Histos[kFT3TrackQ4plus]->SetTitle(
      Form("nChargeMatch = %d (%.2f%%)", nChargeMatch4plus,
           100. * nChargeMatch4plus / (nChargeMiss4plus + nChargeMatch4plus)));

  qMatchEff->SetTitle(Form("Charge match = %.2f%%",
                           100. * nChargeMatch / (nChargeMiss + nChargeMatch)));

  // Remove stat boxes
  TH2Histos[kFT3rackQPRec_MC]->SetStats(0);
  TH2Histos[kFT3rackPtResolution]->SetStats(0);
  TH2Histos[kFT3rackInvPtResolution]->SetStats(0);
  TH2Histos[kMCTracksEtaZ]->SetStats(0);
  PtRes_Profile->SetStats(0);
  DeltaX_Profile->SetStats(0);
  TH1Histos[kFT3TrackQ]->SetStats(0);

  // Fit Slices: Pt resolution
  FitSlicesy(*TH2Histos[kFT3rackInvPtResolution], *TH2Histos[kFT3rackQPRec_MC]);
  FitSlicesy(*TH2Histos[kFT3rackPtResolution], *TH2Histos[kFT3rackQPRec_MC]);

  // sigmaX resultion Profile
  TH1D *DeltaX_Error = new TH1D();
  DeltaX_Error = DeltaX_Profile->ProjectionX("DeltaX_Error", "C=E");

  // Summary Canvases

  auto param_resolution = summary_report_3x2(
      *TH2Histos[kFT3TrackDeltaXYVertex], *TH2Histos[kFT3rackPtResolution],
      *PtRes_Profile, *DeltaX_Error, *TH2Histos[kFT3rackQPRec_MC], *qMatchEff,
      "Param Summary", seed_cfg, 0, 0, 0, 0, 0, 0,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3TrackDeltaXYVertex]->Integral() /
                         TH2Histos[kFT3TrackDeltaXYVertex]->GetEntries()),
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackPtResolution]->Integral() /
                         TH2Histos[kFT3rackPtResolution]->GetEntries()),
      "-", "-",
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackQPRec_MC]->Integral() /
                         TH2Histos[kFT3rackQPRec_MC]->GetEntries()),
      "-");

  auto covariances_summary = summary_report_3x2(
      *TH1Histos[kFT3TrackDeltaXErr], *TH1Histos[kFT3TrackDeltaPhiErr],
      *TH1Histos[kFT3TrackDeltainvQPtErr], *TH1Histos[kFT3TrackDeltaYErr],
      *TH1Histos[kFT3TrackDeltaTanLErr], *TH2Histos[kFT3rackQPRec_MC],
      "Covariances Summary", seed_cfg, 1, 1, 1, 1, 1, 0,
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaXErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaXErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltainvQPtErr]->Integral() /
                         TH1Histos[kFT3TrackDeltainvQPtErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaYErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaYErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanLErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanLErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackQPRec_MC]->Integral() /
                         TH2Histos[kFT3rackQPRec_MC]->GetEntries()));

  auto long_summary = summary_report_3x3(
      *TH2Histos[kFT3TrackDeltaXYVertex], *TH1Histos[kFT3TrackDeltaXErr],
      *TH1Histos[kFT3TrackDeltaYErr], *DeltaX_Error,
      *TH2Histos[kFT3rackQPRec_MC], *TH1Histos[kFT3TrackDeltaPhiErr],
      *qMatchEff, *TH1Histos[kFT3TrackDeltainvQPtErr],
      *TH1Histos[kFT3TrackDeltaTanLErr], "Summary3x3", seed_cfg, 0, 1, 1, 0, 0,
      1, 0, 1, 1,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3TrackDeltaXYVertex]->Integral() /
                         TH2Histos[kFT3TrackDeltaXYVertex]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaXErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaXErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaYErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaYErr]->GetEntries()),
      "-",
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackQPRec_MC]->Integral() /
                         TH2Histos[kFT3rackQPRec_MC]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiErr]->GetEntries()),
      "-",
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltainvQPtErr]->Integral() /
                         TH1Histos[kFT3TrackDeltainvQPtErr]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanLErr]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanLErr]->GetEntries()));

  auto param_summary_diff_pt = summary_report_3x3(
      *TH1Histos[kFT3TrackDeltaX0_1], *TH1Histos[kFT3TrackDeltaTanl0_1],
      *TH1Histos[kFT3TrackDeltaPhiDeg0_1], *TH1Histos[kFT3TrackDeltaX1_4],
      *TH1Histos[kFT3TrackDeltaTanl1_4], *TH1Histos[kFT3TrackDeltaPhiDeg1_4],
      *TH1Histos[kFT3TrackDeltaX4plus], *TH1Histos[kFT3TrackDeltaTanl4plus],
      *TH1Histos[kFT3TrackDeltaPhiDeg4plus], "ParamSummaryVsPt", seed_cfg, 1, 1,
      1, 1, 1, 1, 1, 1, 1,
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX0_1]->Integral() /
                         TH1Histos[kFT3TrackDeltaX0_1]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanl0_1]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanl0_1]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg0_1]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg0_1]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX1_4]->Integral() /
                         TH1Histos[kFT3TrackDeltaX1_4]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanl1_4]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanl1_4]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg1_4]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg1_4]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX4plus]->Integral() /
                         TH1Histos[kFT3TrackDeltaX4plus]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanl4plus]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanl4plus]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg4plus]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg4plus]->GetEntries()));

  auto pt_resolution = summary_report(
      *TH2Histos[kFT3rackPtResolution], *TH2Histos[kFT3rackQPRec_MC],
      *PtRes_Profile, *qMatchEff, "Pt Summary", seed_cfg, 0, 0, 0, 0,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackPtResolution]->Integral() /
                         TH2Histos[kFT3rackPtResolution]->GetEntries()),
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackQPRec_MC]->Integral() /
                         TH2Histos[kFT3rackQPRec_MC]->GetEntries()));

  auto invpt_resolution = summary_report(
      *TH2Histos[kFT3rackInvPtResolution], *TH2Histos[kFT3rackQPRec_MC],
      *(TH1F *)gDirectory->Get(
          (std::string(TH2Histos[kFT3rackInvPtResolution]->GetName()) +
           std::string("_1"))
              .c_str()),
      *(TH1F *)gDirectory->Get(
          (std::string(TH2Histos[kFT3rackInvPtResolution]->GetName()) +
           std::string("_2"))
              .c_str()),
      "InvPt Summary", seed_cfg, 0, 0, 0, 0,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackInvPtResolution]->Integral() /
                         TH2Histos[kFT3rackInvPtResolution]->GetEntries()),
      Form("%.2f%%", 100.0 * TH2Histos[kFT3rackQPRec_MC]->Integral() /
                         TH2Histos[kFT3rackQPRec_MC]->GetEntries()));

  auto vertexing_resolution = summary_report(
      *TH2Histos[kFT3TrackDeltaXYVertex], *TH1Histos[kFT3TrackDeltaX],
      *DeltaX_Error, *TH1Histos[kFT3TrackDeltaPhiDeg], "Vertexing Summary",
      seed_cfg, 0, 1, 0, 1,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3TrackDeltaXYVertex]->Integral() /
                         TH2Histos[kFT3TrackDeltaXYVertex]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX]->Integral() /
                         TH1Histos[kFT3TrackDeltaX]->GetEntries()),
      Form("-"),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg]->GetEntries()));

  auto vertexing_resolution0_1 = summary_report(
      *TH2Histos[kFT3TrackDeltaXYVertex0_1], *TH1Histos[kFT3TrackDeltaX0_1],
      *TH1Histos[kFT3TrackDeltaTanl0_1], *TH1Histos[kFT3TrackDeltaPhiDeg0_1],
      "Vertexing Summary pt < 1", seed_cfg, 0, 1, 1, 1,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3TrackDeltaXYVertex0_1]->Integral() /
                         TH2Histos[kFT3TrackDeltaXYVertex0_1]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX0_1]->Integral() /
                         TH1Histos[kFT3TrackDeltaX0_1]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanl0_1]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanl0_1]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg0_1]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg0_1]->GetEntries()));

  auto vertexing_resolution1_4 = summary_report(
      *TH2Histos[kFT3TrackDeltaXYVertex1_4], *TH1Histos[kFT3TrackDeltaX1_4],
      *TH1Histos[kFT3TrackDeltaTanl1_4], *TH1Histos[kFT3TrackDeltaPhiDeg1_4],
      "Vertexing Summary 1 < p_t < 4", seed_cfg, 0, 1, 1, 1,
      Form("%.2f%%", 100.0 * TH2Histos[kFT3TrackDeltaXYVertex1_4]->Integral() /
                         TH2Histos[kFT3TrackDeltaXYVertex1_4]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX1_4]->Integral() /
                         TH1Histos[kFT3TrackDeltaX1_4]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanl1_4]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanl1_4]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg1_4]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg1_4]->GetEntries()));

  auto vertexing_resolution4plus = summary_report(
      *TH2Histos[kFT3TrackDeltaXYVertex4plus], *TH1Histos[kFT3TrackDeltaX4plus],
      *TH1Histos[kFT3TrackDeltaTanl4plus],
      *TH1Histos[kFT3TrackDeltaPhiDeg4plus], "Vertexing Summary p_t > 4",
      seed_cfg, 0, 1, 1, 1,
      Form("%.2f%%", 100.0 *
                         TH2Histos[kFT3TrackDeltaXYVertex4plus]->Integral() /
                         TH2Histos[kFT3TrackDeltaXYVertex4plus]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaX4plus]->Integral() /
                         TH1Histos[kFT3TrackDeltaX4plus]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaTanl4plus]->Integral() /
                         TH1Histos[kFT3TrackDeltaTanl4plus]->GetEntries()),
      Form("%.2f%%", 100.0 * TH1Histos[kFT3TrackDeltaPhiDeg4plus]->Integral() /
                         TH1Histos[kFT3TrackDeltaPhiDeg4plus]->GetEntries()));

  // Write histograms to file and export images

  outFile.mkdir("MoreHistos");
  outFile.cd("MoreHistos");

  for (auto &h : TH2Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  for (auto &h : TH1Histos) {
    h->Write();
    if (EXPORT_HISTOS_IMAGES)
      exportHisto(*h);
  }

  PtRes_Profile->Write();
  DeltaX_Profile->Write();
  DeltaX_Error->Write();
  qMatchEff->Write();
  outFile.Close();

  std::cout << "---------------------------------------------------"
            << std::endl;
  std::cout << "-------------   Fitting Summary   -----------------"
            << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;
  std::cout << " P_mean = " << TH1Histos[kFT3TracksP]->GetMean() << std::endl;
  std::cout << " P_StdDev = " << TH1Histos[kFT3TracksP]->GetStdDev()
            << std::endl;
  std::cout << " Tanl_mean = " << TH1Histos[kFT3TrackDeltaTanl]->GetMean()
            << std::endl;
  std::cout << " Tanl_StdDev = " << TH1Histos[kFT3TrackDeltaTanl]->GetStdDev()
            << std::endl;
  std::cout << " Tanl_StdDev(pt<1) = "
            << TH1Histos[kFT3TrackDeltaTanl0_1]->GetStdDev() << std::endl;
  std::cout << " Tanl_StdDev(1<pt<4) = "
            << TH1Histos[kFT3TrackDeltaTanl1_4]->GetStdDev() << std::endl;
  std::cout << " Tanl_StdDev(pt>4) = "
            << TH1Histos[kFT3TrackDeltaTanl4plus]->GetStdDev() << std::endl;
  std::cout << " Phi_mean = " << TH1Histos[kFT3TrackDeltaPhi]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDev = " << TH1Histos[kFT3TrackDeltaPhi]->GetStdDev()
            << std::endl;
  std::cout << " Phi_StdDev(pt<1) = "
            << TH1Histos[kFT3TrackDeltaPhi0_1]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDev(1<pt<4) = "
            << TH1Histos[kFT3TrackDeltaPhi1_4]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDev(pt>4) = "
            << TH1Histos[kFT3TrackDeltaPhi4plus]->GetStdDev() << std::endl;
  std::cout << " Phi_meanDeg = " << TH1Histos[kFT3TrackDeltaPhiDeg]->GetMean()
            << std::endl;
  std::cout << " Phi_StdDevDeg = "
            << TH1Histos[kFT3TrackDeltaPhiDeg]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(pt<1) = "
            << TH1Histos[kFT3TrackDeltaPhiDeg0_1]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(1<pt<4) = "
            << TH1Histos[kFT3TrackDeltaPhiDeg1_4]->GetStdDev() << std::endl;
  std::cout << " Phi_StdDevDeg(pt>4) = "
            << TH1Histos[kFT3TrackDeltaPhiDeg4plus]->GetStdDev() << std::endl;
  std::cout << " DeltaX_mean = " << TH1Histos[kFT3TrackDeltaX]->GetMean()
            << std::endl;
  std::cout << " DeltaX_StdDev = " << TH1Histos[kFT3TrackDeltaX]->GetStdDev()
            << std::endl;
  std::cout << " DeltaX_StdDev(pt<1) = "
            << TH1Histos[kFT3TrackDeltaX0_1]->GetStdDev() << std::endl;
  std::cout << " DeltaX_StdDev(1<pt<4) = "
            << TH1Histos[kFT3TrackDeltaX1_4]->GetStdDev() << std::endl;
  std::cout << " DeltaX_StdDev(pt>4) = "
            << TH1Histos[kFT3TrackDeltaX4plus]->GetStdDev() << std::endl;
  std::cout << " DeltaY_mean = " << TH1Histos[kFT3TrackDeltaY]->GetMean()
            << std::endl;
  std::cout << " DeltaY_StdDev = " << TH1Histos[kFT3TrackDeltaY]->GetStdDev()
            << std::endl;
  std::cout << " R_mean = " << TH1Histos[kFT3TrackR]->GetMean() << std::endl;
  std::cout << " R_StdDev = " << TH1Histos[kFT3TrackR]->GetStdDev()
            << std::endl;
  std::cout << " Charge_mean = " << TH1Histos[kFT3TrackDeltaY]->GetMean()
            << std::endl;
  std::cout << " nChargeMatch = " << nChargeMatch << " ("
            << 100. * nChargeMatch / (nChargeMiss + nChargeMatch) << "%)"
            << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;

  return 0;
}
