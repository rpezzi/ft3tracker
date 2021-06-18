#include <vector>
#include "FwdFATProbe.h"

std::vector<Float_t> zPositionsMFT{-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76., -77.5};
std::vector<Float_t> zPositionsFT3{-16., -20., -24., -77., -100., -122., -150., -180., -220., -279.};
std::vector<Float_t> x2X0FT3{0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

Float_t EtaToTanl(Float_t eta)
{
  return tan(TMath::Pi() / 2 - 2 * atan(exp(-eta)));
}

void FWDFat(float eta = 2.6)
{
  Float_t tanl0 = EtaToTanl(eta);
  Float_t phi0 = 0;

  vector<Float_t> ptResolutions;
  vector<Float_t> pts;
  auto ptmin = 0.1;
  auto ptmax = 10.;
  auto nPtPoints = 10;

  for (int j = 0; j <= nPtPoints; j++) {
    auto pt = ptmin + j * (ptmax - ptmin) / nPtPoints;
    FwdFATProbe probe;
    probe.mVerbose = !true;
    //auto phi0 = 2.03355;
    //auto tanl0 = -6;
    //auto invqpt = 1;
    Float_t sigma = 8.44e-4;
    Float_t sigma2 = sigma * sigma;
    Float_t zField = 5.0;

    Float_t invqpt = 1.0 / pt;

    probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

    for (int i = zPositionsFT3.size(); i--;) {
      auto z = zPositionsFT3[i];
      auto x2X0Layer = x2X0FT3[i];
      //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

      probe.updateFAT(z, sigma2, x2X0Layer);
      //std::cout << "   partial Pt Resolution = " << probe.getPtResolution() << std::endl
      //          << std::endl;
    }
    probe.propagateToZ(0, zField);
    //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
    //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

    //std::cout << "   Pt Resolution = " << probe.getPtResolution() * 100. << " %" << std::endl
    //          << std::endl;
    ptResolutions.push_back(probe.getPtResolution());
    pts.push_back(pt);
  }
  std::cout << " ====> FT3 FAT: eta = " << eta << std::endl;
  for (int i = pts.size(); i--;) {
    std::cout << " pt = " << pts[i] << " ==> resolution = " << ptResolutions[i] * 100.0 << " % \n";
  }
}

float FT3FATPtRes(float pt, float eta, float sigma = 8.44e-4, float zField = 5.0)
{
  Float_t tanl0 = EtaToTanl(eta);
  Float_t phi0 = 0;
  float sigma2 = sigma * sigma;

  FwdFATProbe probe;
  probe.mVerbose = !true;

  Float_t invqpt = 1.0 / pt;

  probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

  for (int i = zPositionsFT3.size(); i--;) {
    auto z = zPositionsFT3[i];
    auto x2X0Layer = x2X0FT3[i];
    //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

    probe.updateFAT(z, sigma2, x2X0Layer);
    //std::cout << "   partial Pt Resolution = " << probe.getPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZ(0, zField);
  //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
  //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

  //std::cout << "   Pt Resolution = " << probe.getPtResolution() * 100. << " %" << std::endl
  //          << std::endl;
  return probe.getPtResolution();
}
