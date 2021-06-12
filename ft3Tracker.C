#include "FT3Track.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackFitter.h"
#include "ft3tools/MagField.C"

#include <map>

using o2::ft3::FT3Track;
using o2::itsmft::Hit;

bool DEBUG_VERBOSE = false;
std::vector<float> layersX2X0 = {};

std::vector<float_t> loadx2X0fromFile(std::string configFileName);

void ft3Tracker() {

  o2::ft3::TrackFitter fitter;
  fitter.mVerbose = DEBUG_VERBOSE;
  auto field_z = getZField(0, 0, 0); // Get field at Center of ALICE
  fitter.setBz(field_z);
  fitter.setLayersx2X0(loadx2X0fromFile("FT3_layout.cfg"));

  std::string hitfile = "o2sim_HitsFT3.root";
  std::string tr3Tracksfile = "ft3tracks.root";

  vector<map<int, FT3Track>> allFwdTracks;
  vector<vector<o2::mft::TrackMFT>> recoFwdTracks;
  vector<vector<Int_t>> recoFwdTrackIDs;

  TFile *HitFileIn = new TFile(hitfile.c_str());

  TTree *hitTree = (TTree *)HitFileIn->Get("o2sim");

  vector<Hit> *hit = nullptr;
  hitTree->SetBranchAddress("FT3Hit", &hit);

  Int_t nEvents = hitTree->GetEntries();

  allFwdTracks.resize(nEvents);
  recoFwdTracks.resize(nEvents);
  recoFwdTrackIDs.resize(nEvents);


  auto nPrintTracks = 0;

  for (Int_t event = 0; event < nEvents; event++) {

    hitTree->GetEntry(event);

    Int_t nHits = hit->size();
    if (fitter.mVerbose) {
      std::cout << " Building tracks at event " << event << std::endl;
    }

    auto &FwdTracksMap = allFwdTracks.at(event);
    for (Int_t iHit = 0; iHit < nHits; iHit++) { // Attach hits to tracks
      Hit *thisHit = &(*hit)[iHit];
      Int_t myID = thisHit->GetTrackID();
      FwdTracksMap[myID] = FwdTracksMap[myID];
      FwdTracksMap[myID].addHit(*thisHit, iHit);
    }
  }

  TFile *FT3TracksFileOut = new TFile(tr3Tracksfile.c_str(), "RECREATE");
  TTree *FT3Tree = new TTree("o2sim", "Tree with FT3 Tracks");
  vector<o2::mft::TrackMFT> *recoTracks;
  vector<Int_t> *recoTrackIDs;
  FT3Tree->Branch("FT3Track", &recoTracks);
  FT3Tree->Branch("FT3TrackID", &recoTrackIDs);

  for (Int_t event = 0; event < nEvents; event++) { // Fit and save tracks
    recoTracks = &recoFwdTracks[event];
    recoTrackIDs = &recoFwdTrackIDs[event];
    auto &FwdTracksMap = allFwdTracks.at(event);

    auto iter = FwdTracksMap.begin();
    while (iter != FwdTracksMap.end()) {
      auto &trackID = iter->first;
      auto &track = iter->second;
      auto nHits = track.getNumberOfPoints();
      if (nHits > 3) {
        fitter.initTrack(track);
        fitter.fit(track);
        //if (std::abs(track.getTanl()) < 40 and track.getP() < 500) { // filter
          recoTracks->emplace_back(track);
          recoTrackIDs->emplace_back(trackID);
        //}
      }
      if (fitter.mVerbose) {
        cout << "\n[TrackID = " << trackID << ", nHits = " << nHits
             << "] LayerIDs, zCoord => ";
        for (auto i = 0; i < nHits; i++) {
          std::cout << " " << track.getLayers()[i] << ", "
                    << track.getZCoordinates()[i] << " ";
        if(nPrintTracks > 10) fitter.mVerbose = false;
        }
        std::cout << std::endl;
      }
      ++iter;
      nPrintTracks++;
    }

    FT3Tree->Fill();
  }
  FT3TracksFileOut->Write();
}


std::vector<float_t> loadx2X0fromFile(std::string configFileName = "FT3_layout.cfg")
{
  std::vector<float_t> Layersx2X0;
  std::ifstream ifs(configFileName.c_str());
  if (!ifs.good()) {
    LOG(FATAL) << " Invalid FT3Base.configFile!";
  }
  std::string tempstr;
  float z_layer, r_in, r_out, Layerx2X0;
  char delimiter;
  int layerNumber = 0;
  while (std::getline(ifs, tempstr)) {
    if (tempstr[0] == '#') {
      LOG(INFO) << " Comment: " << tempstr;
      continue;
    }
    LOG(INFO) << " Line: " << tempstr;
    std::istringstream iss(tempstr);
    iss >> z_layer;
    iss >> r_in;
    iss >> r_out;
    iss >> Layerx2X0;

    Layersx2X0.push_back(Layerx2X0);
    LOG(INFO) << " loadx2X0fromFile z =  " << z_layer << " ; x/X0 = " << Layerx2X0 << std::endl;
  }
  return Layersx2X0;
}
