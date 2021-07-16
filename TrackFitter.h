// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackFitter.h
/// \brief Definition of a class to fit a track to a set of clusters
///
/// \author Philippe Pillot, Subatech; adapted by Rafael Pezzi, UFRGS

#ifndef ALICEO2_FT3_TrackFitter_H_
#define ALICEO2_FT3_TrackFitter_H_

#include "DataFormatsMFT/TrackMFT.h"
#include "FT3Track.h"
#include "MFTTracking/Cluster.h"
#include <TLinearFitter.h>
#include <list>

namespace o2
{
namespace ft3
{

class TrackFitter
{

  using SMatrix55Sym =
    ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

 public:
  TrackFitter() = default;
  ~TrackFitter() = default;

  TrackFitter(const TrackFitter&) = delete;
  TrackFitter& operator=(const TrackFitter&) = delete;
  TrackFitter(TrackFitter&&) = delete;
  TrackFitter& operator=(TrackFitter&&) = delete;

  void setBz(Double_t bZ);

  bool initTrack(FT3Track& track, bool outward = false);
  bool fit(FT3Track& track, bool outward = false);
  void MinuitFit(FT3Track& track);
  void setLayersx2X0(std::vector<Double_t> v) { mLayersx2X0 = v; }

  /// Return the maximum chi2 above which the track can be considered as
  /// abnormal
  static constexpr double getMaxChi2() { return SMaxChi2; }
  bool mVerbose = false;

  //Vectors for Minuit Fitter
  static std::vector<float, std::allocator<float>> PosX;
  static std::vector<float, std::allocator<float>> PosY;
  static std::vector<float, std::allocator<float>> PosZ;
  static std::vector<float, std::allocator<float>> ErrorsX;
  static std::vector<float, std::allocator<float>> ErrorsY;

 private:
  bool computeCluster(FT3Track& track, int cluster);

  Double_t mBZField; // kiloGauss.
  static constexpr double SMaxChi2 =
    2.e10; ///< maximum chi2 above which the track can be considered as
           ///< abnormal

  bool mFieldON = true;
  std::vector<Double_t> mLayersx2X0;
};

// Functions to estimate momentum and charge from track curvature
Double_t invQPtFromFCF(const FT3Track& track, Double_t bFieldZ, Double_t& chi2);
Bool_t LinearRegression(Int_t nVal, Double_t* xVal, Double_t* yVal,
                        Double_t* yErr, Double_t& a, Double_t& ae, Double_t& b,
                        Double_t& be);

// Function used by Minuit Fitter and vector declarations
void myFitFcn(Int_t&, Double_t*, Double_t& fval, Double_t* p, Int_t);
std::vector<float, std::allocator<float>> TrackFitter::PosX;
std::vector<float, std::allocator<float>> TrackFitter::PosY;
std::vector<float, std::allocator<float>> TrackFitter::PosZ;
std::vector<float, std::allocator<float>> TrackFitter::ErrorsX;
std::vector<float, std::allocator<float>> TrackFitter::ErrorsY;

} // namespace ft3
} // namespace o2

#include "TrackFitter.cxx"

#endif // ALICEO2_FT3_TrackFitter_H_
