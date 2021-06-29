
#ifndef O2_FWD_FAT_PROBE_H_
#define O2_FWD_FAT_PROBE_H_

#include <TMath.h>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"

using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;
using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using o2::track::TrackParCovFwd;

class FwdFATProbe : public o2::track::TrackParCovFwd
{
 public:
  FwdFATProbe() = default;
  FwdFATProbe(const FwdFATProbe& t) = default;
  ~FwdFATProbe() = default;
  void init(const Double_t phi0, const Double_t tanl0, const Double_t invqpt0, const Double_t lastZ, const Double_t zField)
  {
    mZField = zField;
    SMatrix5 parameters = {0, 0, phi0, tanl0, invqpt0};
    mStartingParameters.setParameters(parameters);
    mStartingParameters.setZ(0.0);
    //mStartingParameters.propagateParamToZhelix(lastZ,mZField);

    setParameters(parameters);
    setZ(0.0);
    propagateParamToZhelix(lastZ, mZField);

    SMatrix55Sym covariances;
    Double_t qptsigma = TMath::Max((Double_t)std::abs(getInvQPt()), .5);
    Double_t tanlsigma = TMath::Max((Double_t)std::abs(getTanl()), .5);

    covariances(0, 0) = 1;                              // <X,X>
    covariances(1, 1) = 1;                              // <Y,Y>
    covariances(2, 2) = TMath::Pi() * TMath::Pi() / 16; // <PHI,PHI>
    covariances(3, 3) = 10. * tanlsigma * tanlsigma;    // <TANL,TANL>
    covariances(4, 4) = 10. * qptsigma * qptsigma;      // <INVQPT,INVQPT>
    setCovariances(covariances);
    if (mVerbose) {
      std::cout << std::endl
                << "  *** Init fwd FAT Probe: " << std::endl;
      print();
    }
  }

  void updateFAT(const Float_t nextZ, Float_t sigma2, Double_t Layerx2X0)
  {
    if (mVerbose) {
      std::cout << std::endl
                << "  *** UpdateFat Starting: currentZ = " << getZ() << " ; nextZ = " << nextZ << " ; Layerx2X0 = " << Layerx2X0 << " ; sigma = " << std::sqrt(sigma2) << std::endl;
      print();
    };
    propagateToZhelix(nextZ, mZField);
    if (mVerbose) {
      std::cout << "  UpdateFat After Propagation: " << std::endl;
      print();
    }
    const std::array<Float_t, 2>& pos = {(const Float_t)getX(), (const Float_t)getY()};
    const std::array<Float_t, 2>& cov = {sigma2, sigma2};

    update(pos, cov);
    if (mVerbose) {
      std::cout << "  UpdateFat After Kalman: " << std::endl;
      print();
    }

    if (Layerx2X0 != 0) {
      addMCSEffect(Layerx2X0);
      if (mVerbose) {
        std::cout << "  UpdateFat After MCS: " << std::endl;
        print();
      }
    }
  }

  Double_t getInvQPtResolution()
  {
    return (std::sqrt(getCovariances()(4, 4)) * getPt());
  }

  Double_t getVertexSigmaXResolution()
  {
    FwdFATProbe tempprobe(*this);
    tempprobe.propagateToZ(mStartingParameters.getZ(), mZField);
    return (1.0e4 * std::sqrt(tempprobe.getCovariances()(0, 0))); // microns
  }

  void print()
  {
    std::cout << "  ProbeParams: X = " << getX() << " Y = " << getY() << " Z = " << getZ() << " Tgl = " << getTanl() << "  Phi = " << getPhi() << " q/pt = " << getInvQPt() << std::endl;
    std::cout << "  Covariance: " << std::endl
              << getCovariances() << std::endl
              << std::endl;
  }

  bool mVerbose = false;
  o2::track::TrackParFwd mStartingParameters;
  Double_t mZField;
};

#endif /* O2_FWD_FAT_PROBE_H_ */
