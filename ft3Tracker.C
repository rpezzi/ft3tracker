#include "FT3Track.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackFitter.h"

#include <map>

using o2::ft3::FT3Track;
using o2::itsmft::Hit;

void ft3Tracker() {

  o2::ft3::TrackFitter fitter;
  fitter.setBz(-5.0);

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

  for (Int_t event = 0; event < nEvents; event++) {

    hitTree->GetEntry(event);

    Int_t nHits = hit->size();
    std::cout << " Event " << event << " has " << nHits << " hits" << std::endl;

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
        if (std::abs(track.getTanl()) < 40 and track.getP() < 500) { // filter
          recoTracks->emplace_back(track);
          recoTrackIDs->emplace_back(trackID);
        }
      }
      cout << "\n[TrackID = " << trackID << ", nHits = " << nHits
           << "] LayerIDs, zCoord => ";
      for (auto i = 0; i < nHits; i++) {
        std::cout << " " << track.getLayers()[i] << ", "
                  << track.getZCoordinates()[i] << " ";
      }
      std::cout << std::endl;
      ++iter;
    }

    FT3Tree->Fill();
  }
  FT3TracksFileOut->Write();
}
