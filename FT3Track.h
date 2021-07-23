// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///

#ifndef O2_FT3_TRACK_H_
#define O2_FT3_TRACK_H_

#include "TRandom3.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include <vector>

namespace o2
{
namespace ft3
{

using o2::itsmft::Hit;

class FT3Track : public o2::mft::TrackMFT
{
 public:
  FT3Track() = default;
  FT3Track(const FT3Track& t) = default;
  ~FT3Track() = default;
  const std::vector<Float_t>& getXCoordinates() const { return mX; }
  const std::vector<Float_t>& getYCoordinates() const { return mY; }
  const std::vector<Float_t>& getZCoordinates() const { return mZ; }
  const std::vector<Float_t>& getSigmasX2() const { return mSigmaX2; }
  const std::vector<Float_t>& getSigmasY2() const { return mSigmaY2; }
  const std::vector<Int_t>& getLayers() const { return mLayer; }
  const std::vector<Int_t>& getHitsId() const { return mHitId; }
  void addHit(const Hit& ht, const Int_t hitId, const Float_t sigma, Bool_t smear);
  void sort();

 private:
  std::vector<Float_t> mX;
  std::vector<Float_t> mY;
  std::vector<Float_t> mZ;
  std::vector<Float_t> mSigmaX2;
  std::vector<Float_t> mSigmaY2;
  std::vector<Int_t> mLayer;
  std::vector<Int_t> mHitId;

  ClassDefNV(FT3Track, 1);
};

//_________________________________________________________________________________________________
inline void FT3Track::addHit(const Hit& ht, const Int_t hitId, const Float_t sigma = 8.44e-4, Bool_t smear = true)
{
  TRandom3 rnd(0);
  auto x = ht.GetStartX() + (smear ? rnd.Gaus(0, sigma) : 0);
  auto y = ht.GetStartY() + (smear ? rnd.Gaus(0, sigma) : 0);
  auto sigma2 = sigma * sigma;
  mX.emplace_back(x);
  mY.emplace_back(y);
  mZ.emplace_back(ht.GetStartZ());
  mSigmaX2.emplace_back(sigma2);
  mSigmaY2.emplace_back(sigma2);
  mLayer.emplace_back(ht.GetDetectorID());
  mHitId.emplace_back(hitId);
  setNumberOfPoints(mX.size());
}

//_________________________________________________________________________________________________
inline void FT3Track::sort()
{
  // Orders elements along z position
  struct HitData {
    Float_t x;
    Float_t y;
    Float_t z;
    Int_t layer;
    Int_t hitId;
    MCCompLabel label;
  };
  std::vector<HitData> points;

  // Loading hit data
  for (auto point = 0; point < (int)getNumberOfPoints(); ++point) {
    auto& somepoint = points.emplace_back();
    somepoint.x = mX[point];
    somepoint.y = mY[point];
    somepoint.z = mZ[point];
    somepoint.layer = mLayer[point];
    somepoint.hitId = mHitId[point];
  }

  // Sorting hit data
  std::sort(points.begin(), points.end(),
            [](HitData a, HitData b) { return std::abs(a.z) < std::abs(b.z); });

  // Storing sorted hit data
  for (auto point = 0; point < (int)getNumberOfPoints(); ++point) {
    mX[point] = points[point].x;
    mY[point] = points[point].y;
    mZ[point] = points[point].z;
    mLayer[point] = points[point].layer;
    mHitId[point] = points[point].hitId;
  }
}

} // namespace ft3
} // namespace o2
#endif /* O2_FT3_TRACK_H_ */
