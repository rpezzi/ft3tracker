#include "FwdFATProbe.h"
#include <vector>

std::vector<Float_t> zPositionsMFT{-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76., -77.5};
std::vector<Float_t> zPositionsFT3{-16., -20., -24., -77., -100., -122., -150., -180., -220., -279.};
std::vector<Float_t> x2X0FT3{0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
//std::vector<Float_t> x2X0FT3(10,1e-7);




Float_t EtaToTanl(Float_t eta)
{
  return tan(TMath::Pi() / 2 - 2 * atan(exp(-eta)));
}

void FWDFat_(float eta = 2.6, Float_t sigma = 8.44e-4)
{

  Float_t tanl0 = EtaToTanl(eta);
  Float_t phi0 = 0;

  std::vector<Float_t> ptResolutions;
  std::vector<Float_t> vtxXResolutions;
  std::vector<Float_t> pts;
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
    Float_t sigma2 = sigma * sigma;
    Float_t zField = 5.0;

    Float_t invqpt = 1.0 / pt;

    probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

    for (int i = zPositionsFT3.size(); i--;) {
      auto z = zPositionsFT3[i];
      auto x2X0Layer = x2X0FT3[i];
      //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

      probe.updateFAT(z, sigma2, x2X0Layer);
      //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
      //          << std::endl;
    }
    probe.propagateToZ(0, zField);
    //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
    //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

    //std::cout << "   Pt Resolution = " << probe.getInvQPtResolution() * 100. << " %" << std::endl
    //          << std::endl;
    ptResolutions.push_back(probe.getInvQPtResolution());
    vtxXResolutions.push_back(probe.getVertexSigmaXResolution());
    pts.push_back(pt);
  }
  std::cout << " ====> FT3 FAT: eta = " << eta << std::endl;
  for (int i = pts.size(); i--;) {
    std::cout << " pt = " << pts[i] << " \n  ==> Pt Resolution = " << ptResolutions[i] * 100.0 << " % \n"
              << "  ==> Vertex sigma_x = " << vtxXResolutions[i] << " um \n";
  }
}

void FWDFat(float eta = 2.6, float pt = 1.0, Float_t sigma = 8.44e-4)
{

  Float_t tanl0 = EtaToTanl(eta);
  Float_t phi0 = 0;

  std::vector<Float_t> ptResolutions;
  std::vector<Float_t> vtxXResolutions;
  std::vector<Float_t> pts;
  auto ptmin = 0.1;
  auto ptmax = 10.;
  auto nPtPoints = 10;

  FwdFATProbe probe;

  probe.mVerbose = !true;
  //auto phi0 = 2.03355;
  //auto tanl0 = -6;
  //auto invqpt = 1;
  Float_t sigma2 = sigma * sigma;
  Float_t zField = 5.0;

  Float_t invqpt = 1.0 / pt;

  probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

  for (int i = zPositionsFT3.size(); i--;) {
    auto z = zPositionsFT3[i];
    auto x2X0Layer = x2X0FT3[i];
    //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

    probe.updateFAT(z, sigma2, x2X0Layer);
    //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZ(0, zField);
  //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
  //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

  //std::cout << "   Pt Resolution = " << probe.getInvQPtResolution() * 100. << " %" << std::endl
  //          << std::endl;
  ptResolutions.push_back(probe.getInvQPtResolution());
  vtxXResolutions.push_back(probe.getVertexSigmaXResolution());
  pts.push_back(pt);

  std::cout << " ====> FT3 FAT: eta = " << eta << std::endl;
  for (int i = pts.size(); i--;) {
    std::cout << " pt = " << pts[i] << " \n  ==> Pt Resolution = " << ptResolutions[i] * 100.0 << " % \n"
              << "  ==> Vertex sigma_x = " << vtxXResolutions[i] << " um \n";
  }
}

float FT3FATPtRes(float pt, float eta, float sigma = 8.44e-4, bool enableMCS = true, float zField = 5.0)
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
    float x2X0Layer = enableMCS ? x2X0FT3[i] : 0.0;
    //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

    probe.updateFAT(z, sigma2, x2X0Layer);
    //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZ(0, zField);
  //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
  //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

  //std::cout << "   Pt Resolution = " << probe.getInvQPtResolution() * 100. << " %" << std::endl
  //          << std::endl;
  return probe.getInvQPtResolution();
}

float FT3FATvtxXRes(float pt, float eta, float sigma = 8.44e-4, float zField = 5.0)
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
    //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZ(0, zField);
  //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
  //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

  //std::cout << "   Pt Resolution = " << probe.getInvQPtResolution() * 100. << " %" << std::endl
  //          << std::endl;
  return probe.getVertexSigmaXResolution();
}

void checkFT3FATPtRes(float pt, float eta, float sigma = 8.44e-4, float zField = 5.0)
{
  auto ptres_noMC = FT3FATPtRes(pt, eta, sigma, false);
  auto ptres_onlyMC = FT3FATPtRes(pt, eta, 1e-7, true);
  auto ptres_combined = sqrt(ptres_noMC * ptres_noMC + ptres_onlyMC * ptres_onlyMC);
  auto ptres_FAT = FT3FATPtRes(pt, eta, sigma);
  std::cout << " ptres_noMC = " << ptres_noMC << std::endl;
  std::cout << " ptres_onlyMC = " << ptres_onlyMC << std::endl;
  std::cout << " ptres_combined = " << ptres_combined << std::endl;
  std::cout << " ptres_FAT = " << ptres_FAT << std::endl;
  return;
}

//_________________________________________________________________________________________
std::vector<float> getFATPtRes_pts_at_eta(std::vector<float> pts, float eta, float sigma, bool enableMCS = true)
{
  std::vector<float> ptResolutions;
  for (auto pt : pts) {
    auto ptres = FT3FATPtRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " ptRes = " << ptres << std::endl;
    ptResolutions.push_back(ptres);
  }
  return ptResolutions;
}

//_________________________________________________________________________________________
std::vector<float> getFATPtRes_etas_at_pt(std::vector<float> etas, float pt, float sigma, bool enableMCS = true)
{
  std::vector<float> ptResolutions;
  for (auto eta : etas) {
    auto ptres = FT3FATPtRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " ptRes = " << ptres << std::endl;
    ptResolutions.push_back(ptres);
  }
  return ptResolutions;
}

//_________________________________________________________________________________________
std::vector<float> getFATvtxXRes_etas_at_pt(std::vector<float> etas, float pt, float sigma, bool enableMCS = true)
{
  std::vector<float> vtxResolutions;
  for (auto eta : etas) {
    auto vtxRes = FT3FATvtxXRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " vtxRes = " << vtxRes << std::endl;
    vtxResolutions.push_back(vtxRes);
  }
  return vtxResolutions;
}

//_________________________________________________________________________________________
std::vector<float> getFATvtxXRes_pts_at_eta(std::vector<float> pts, float eta, float sigma, bool enableMCS = true)
{
  std::vector<float> vtxResolutions;
  for (auto pt : pts) {
    auto vtxRes = FT3FATvtxXRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " vtxRes = " << vtxRes << std::endl;
    vtxResolutions.push_back(vtxRes);
  }
  return vtxResolutions;
}