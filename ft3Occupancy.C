#include "TFile.h"
#include "TTree.h"
#include <TH2F.h>
#include <TH3F.h>

#include "ITSMFTSimulation/Hit.h"
#include "ft3tools/HistosHelpers.C"

#include <map>
#include <set>

using o2::itsmft::Hit;

bool DEBUG_VERBOSE = false;

// Hit density scaled according to QED cross section on pureQED events
// And number of directions (Backward and/or Forward layers)
Float_t xSectionQED = 29533.4;
Float_t xSectionHad = 8;
Bool_t pureQED = true;
Float_t QEDScale = pureQED ? xSectionQED / xSectionHad : 1.0;
std::set<int> computeDirections = {1};                 // Directions on which to compute hits
                                                       // direction == 1 => Forward layers
                                                       // direction == 0 => Backward layers
Float_t nDirectionsScaling = computeDirections.size(); // Histogram scales according to number of directions

void ft3Occupancy() {

  std::string hitfile = "o2sim_HitsFT3.root";
  std::string tr3Occupancyfile = "ft3Occupancy.root";
  std::string annotation = "FT3 occupancy";

  enum TH2HistosCodes {
    kLay0Occupancy = 0,
    kLay1Occupancy = 1,
    kLay2Occupancy = 2,
    kLay3Occupancy = 3,
    kLay4Occupancy = 4,
    kLay5Occupancy = 5,
    kLay6Occupancy = 6,
    kLay7Occupancy = 7,
    kLay8Occupancy = 8,
    kLay9Occupancy = 9
  };

  std::map<int, const char *> TH2Names {
    { kLay0Occupancy, "Layer0Occupancy" },
    { kLay1Occupancy, "Layer1Occupancy" },
    { kLay2Occupancy, "Layer2Occupancy" },
    { kLay3Occupancy, "Layer3Occupancy" },
    { kLay4Occupancy, "Layer4Occupancy" },
    { kLay5Occupancy, "Layer5Occupancy" },
    { kLay6Occupancy, "Layer6Occupancy" },
    { kLay7Occupancy, "Layer7Occupancy" },
    { kLay8Occupancy, "Layer8Occupancy" },
    { kLay9Occupancy, "Layer9Occupancy" }
    };

    std::map<int, std::array<double, 6>> TH2Binning {
    { kLay0Occupancy, { 60, -3, 3, 60, -3, 3 } },
    { kLay1Occupancy, { 60, -3, 3, 60, -3, 3 } },
    { kLay2Occupancy, { 60, -3, 3, 60, -3, 3 } },
    { kLay3Occupancy, { 700/10, -35, 35, 700/10, -35, 35 } },
    { kLay4Occupancy, { 700/10, -35, 35, 700/10, -35, 35 } },
    { kLay5Occupancy, { 700/10, -35, 35, 700/10, -35, 35 } },
    { kLay6Occupancy, { 2000/10, -100, 100, 2000/10, -100, 100 } },
    { kLay7Occupancy, { 2000/10, -100, 100, 2000/10, -100, 100 } },
    { kLay8Occupancy, { 2000/10, -100, 100, 2000/10, -100, 100 } },
    { kLay9Occupancy, { 2000/10, -100, 100, 2000/10, -100, 100 } }
  };

  const int nTH2Histos = TH2Names.size();
    std::vector<std::unique_ptr<TH2F>> TH2Histos(nTH2Histos);
    auto n2Histo = 0;
    for (auto &h : TH2Histos) {
      h = std::make_unique<TH2F>(TH2Names[n2Histo], TH2Names[n2Histo],
                                 (int)TH2Binning[n2Histo][0],
                                 TH2Binning[n2Histo][1], TH2Binning[n2Histo][2],
                                 (int)TH2Binning[n2Histo][3],
                                 TH2Binning[n2Histo][4], TH2Binning[n2Histo][5]);
     h->SetOption("COLZ");
     ++n2Histo;
    }

  std::unique_ptr<TH3F> HitMap3D;
  HitMap3D = std::make_unique<TH3F> ("HitMap3D", "HitMap3D",
                                     200, -100, 100, 200, -100, 100, 2*560, -280, 280);

  TFile *HitFileIn = new TFile(hitfile.c_str());
  TTree *hitTree = (TTree *)HitFileIn->Get("o2sim");

  vector<Hit> *hit = nullptr;
  hitTree->SetBranchAddress("FT3Hit", &hit);

  Int_t nEvents = hitTree->GetEntries();
  std::cout << " Computing hits from " << nEvents << " events" << std::endl;
  for (Int_t event = 0; event < nEvents; event++) {
    hitTree->GetEntry(event);

    Int_t nHits = hit->size();
    if (DEBUG_VERBOSE) {
      std::cout << " Processing event " << event << std::endl;
    }

    for (Int_t iHit = 0; iHit < nHits; iHit++) { // Fill occupancy histograms
      Hit *thisHit = &(*hit)[iHit];
      Int_t myDir = thisHit->GetDetectorID() % 2;
      Int_t nLayer = thisHit->GetDetectorID() / 2;
      if ( computeDirections.find(myDir) != computeDirections.end()) { // Compute hit only if direction is enabled
        if (DEBUG_VERBOSE) {
          std::cout << " Direction/Layer: " << myDir << "/" << nLayer << " z = " << thisHit->GetStartZ() << std::endl;
        }
        TH2Histos[nLayer]->Fill(thisHit->GetStartX(),thisHit->GetStartY());
        HitMap3D->Fill(thisHit->GetStartX(),thisHit->GetStartY(),thisHit->GetStartZ());
      }
    }
  }

  TFile outFile(tr3Occupancyfile.c_str(), "RECREATE");

  Int_t nLayer = 0;
  for (auto &h : TH2Histos) {
    h->SetStats(0);
    h->Scale(QEDScale*1.0/nEvents/nDirectionsScaling); // Histos scaled by (xSectionQED / xSectionHad) if QEDpure is set
                                        // Divided by 2 if counting hits from Backward and Forward directions
    std::cout << " Layer " << nLayer << " Max " <<  h->GetMaximum() << " hits/event/bin " << std::endl;
    h->Write();
    nLayer++;
  }

  // Covariances summary 3x3
  auto occup_summary3x3 = summary_report_3x3(
      *TH2Histos[0], *TH2Histos[1], *TH2Histos[2], *TH2Histos[3], *TH2Histos[4],
      *TH2Histos[5], *TH2Histos[6], *TH2Histos[7], *TH2Histos[8],
      "occup_summary3x3", annotation, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      Form("Max %.2f hits/cm^2/event", 100 * TH2Histos[0]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", 100 * TH2Histos[1]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", 100 * TH2Histos[2]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", TH2Histos[3]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", TH2Histos[4]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", TH2Histos[5]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", TH2Histos[6]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", TH2Histos[7]->GetMaximum()),
      Form("Max %.2f hits/cm^2/event", TH2Histos[8]->GetMaximum())
      );
HitMap3D->Write();
outFile.Close();
}
