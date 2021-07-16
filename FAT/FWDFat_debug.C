#include <vector>
#include <fstream>


#ifndef O2_FWD_FAT_PROBE_H_
#define O2_FWD_FAT_PROBE_H_

#include <TMath.h>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"

using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;
using SMatrix55Sym = ROOT::Math::SMatrix<Double_t, 5, 5, ROOT::Math::MatRepSym<Double_t, 5>>;
using o2::track::TrackParCovFwd;

class FwdFATProbe : public o2::track::TrackParCovFwd
{
 public:
  FwdFATProbe() = default;
  FwdFATProbe(const FwdFATProbe& t) = default;
  ~FwdFATProbe() = default;
  void init(const Double_t phi0, const Double_t tanl0_, const Double_t invqpt0, const Double_t lastZ, const Double_t zField)
  {
    Double_t tanl0 = std::copysign(tanl0_,lastZ); // Ensure tanl0 and lastZ are compatible
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

  void updateFAT(const Double_t nextZ, Double_t sigma2, Double_t Layerx2X0)
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

   //     if (Layerx2X0 != 0) {
   //   addMCSEffect(.0001 * Layerx2X0);
   //   if (mVerbose) {
   //     std::cout << "  UpdateFat After MCS: " << std::endl;
   //     print();
   //   }
   // }

    const std::array<Double_t, 2>& pos = {(const Double_t)getX(), (const Double_t)getY()};
    const std::array<Double_t, 2>& cov = {sigma2, sigma2};

    update(pos, cov);
    if (mVerbose) {
      std::cout << "  UpdateFat After Kalman: " << std::endl;
      print();
    }

    if (Layerx2X0 != 0) {
      addMCSEffect(1. * Layerx2X0);
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
    tempprobe.propagateToZhelix(mStartingParameters.getZ(), mZField);
    return (1.0e4 * std::sqrt(tempprobe.getCovariances()(0, 0))); // microns
  }

  void print()
  {
    std::cout << "  ProbeParams: X = " << getX() << " Y = " << getY() << " Z = " << getZ() << " Tgl = " << getTanl() << "  Phi = " << getPhi() << " q/pt = " << getInvQPt() << std::endl;
    std::cout << "  Variances: " 
              << getCovariances()(0,0) << " "  << getCovariances()(1,1) << " "  << getCovariances()(2,2) << " "  << getCovariances()(3,3) << " "  << getCovariances()(4,4) << " " << std::endl <<
              "  Correlations: \n" 
              << "    0,1: " << getCorrelation(0,1)
              << " ; 0,2: " << getCorrelation(0,2) << " "
              << " ; 0,3: " << getCorrelation(0,3) << " "
              << " ; 0,4: " << getCorrelation(0,4) << " "
              << " ; 1,2: " << getCorrelation(1,2) << " "
              << " ; 1,3: " << getCorrelation(1,3) << " "
              << " ; 1,4: " << getCorrelation(1,4) << " "
              << " ; 2,3: " << getCorrelation(2,3) << " "
              << " ; 2,4: " << getCorrelation(2,4) << " "
              << " ; 3,4: " << getCorrelation(3,4) << std::endl
              << std::endl;
  }

  Double_t getCorrelation(int i, int j) {
   return  getCovariances()(i,j)/TMath::Sqrt(getCovariances()(i,i) * getCovariances()(j,j));
  }

  bool mVerbose = false;
  o2::track::TrackParFwd mStartingParameters;
  Double_t mZField;
};

#endif /* O2_FWD_FAT_PROBE_H_ */



std::vector<Double_t> zPositionsMFT{-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76., -77.5};
std::vector<Double_t> zPositionsFT3{-16., -20., -24., -77., -100., -122., -150., -180., -220., -279.};
//std::vector<Double_t> x2X0FT3{0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

std::vector<Double_t> x2X0FT3;



Double_t EtaToTanl(Double_t eta)
{
  return tan(TMath::Pi() / 2 - 2 * atan(exp(-eta)));
}

void FWDFat_(Double_t eta = 2.6, Double_t sigma = 8.44e-4)
{

  Double_t tanl0 = EtaToTanl(eta);
  Double_t phi0 = 0;

  std::vector<Double_t> ptResolutions;
  std::vector<Double_t> vtxXResolutions;
  std::vector<Double_t> pts;
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
    Double_t sigma2 = sigma * sigma;
    Double_t zField = 5.0;

    Double_t invqpt = 1.0 / pt;

    probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

    for (int i = zPositionsFT3.size(); i--;) {
      auto z = zPositionsFT3[i];
      auto x2X0Layer = x2X0FT3[i];
      //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

      probe.updateFAT(z, sigma2, x2X0Layer);
      //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
      //          << std::endl;
    }
    probe.propagateToZhelix(0, zField);
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



void FWDFat2(Double_t eta = 2.6, Double_t pt = 1.0, Double_t sigma = 8.44e-4)
{

  Double_t tanl0 = EtaToTanl(eta);
  Double_t phi0 = 0;

  std::vector<Double_t> ptResolutions;
  std::vector<Double_t> vtxXResolutions;
  std::vector<Double_t> pts;
  auto ptmin = 0.1;
  auto ptmax = 10.;
  auto nPtPoints = 10;

  FwdFATProbe probe;

  probe.mVerbose = !true;
  //auto phi0 = 2.03355;
  //auto tanl0 = -6;
  //auto invqpt = 1;
  Double_t sigma2 = sigma * sigma;
  Double_t zField = 5.0;

  Double_t invqpt = 1.0 / pt;

  probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

  for (int i = zPositionsFT3.size(); i--;) {
    auto z = zPositionsFT3[i];
    auto x2X0Layer = x2X0FT3[i];
    //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

    probe.updateFAT(z, sigma2, x2X0Layer);
    //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZhelix(0, zField);
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

Double_t FT3FATPtRes(Double_t pt, Double_t eta, Double_t sigma = 8.44e-4, bool enableMCS = true, Double_t zField = 5.0, bool verbose = false)
{
  Double_t tanl0 = EtaToTanl(eta);
  Double_t phi0 = 1;
  Double_t sigma2 = sigma * sigma;
  if (!x2X0FT3.size()) {
      x2X0FT3 = {0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  }

  FwdFATProbe probe;
  probe.mVerbose = verbose;

  Double_t invqpt = 1.0 / pt;
  if (verbose) {
    std::cout << "FT3FATPtRes. Probe starting parameters: phi0 = " << phi0 << " tanl0 = " << tanl0 << " invqpt = " << invqpt << std::endl;
  }
  probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

  for (int i = zPositionsFT3.size(); i--;) {
    auto z = zPositionsFT3[i];
    Double_t x2X0Layer = enableMCS ? x2X0FT3[i] : 0.0;
    if (verbose) std::cout << "===> layer " << i << " ; z = " << z << std::endl;

    probe.updateFAT(z, sigma2, x2X0Layer);
    //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZhelix(0, zField);
  //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
  //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

  //std::cout << "   Pt Resolution = " << probe.getInvQPtResolution() * 100. << " %" << std::endl
  //          << std::endl;
  return probe.getInvQPtResolution();
}

Double_t FT3FATvtxXRes(Double_t pt, Double_t eta, Double_t sigma = 8.44e-4, Double_t zField = 5.0)
{ 
  
  Double_t tanl0 = EtaToTanl(eta);
  Double_t phi0 = 0;
  Double_t sigma2 = sigma * sigma;
    if (!x2X0FT3.size()) {
      x2X0FT3 = {0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  }

  FwdFATProbe probe;
  probe.mVerbose = !true;

  Double_t invqpt = 1.0 / pt;

  probe.init(phi0, tanl0, invqpt, zPositionsFT3.back(), zField);

  for (int i = zPositionsFT3.size(); i--;) {
    auto z = zPositionsFT3[i];
    auto x2X0Layer = x2X0FT3[i];
    //std::cout << "===> layer " << i << " ; z = " << z << std::endl;

    probe.updateFAT(z, sigma2, x2X0Layer);
    //std::cout << "   partial Pt Resolution = " << probe.getInvQPtResolution() << std::endl
    //          << std::endl;
  }
  probe.propagateToZhelix(0, zField);
  //std::cout << "  ProbeInit:   X = " << probe.mStartingParameters.getX() << " Y = " << probe.mStartingParameters.getY() << " Z = " << probe.mStartingParameters.getZ() << " Tgl = " << probe.mStartingParameters.getTanl() << "  Phi = " << probe.mStartingParameters.getPhi() << " q/pt = " << probe.mStartingParameters.getInvQPt() << std::endl;
  //std::cout << "  ProbeParams: X = " << probe.getX() << " Y = " << probe.getY() << " Z = " << probe.getZ() << " Tgl = " << probe.getTanl() << "  Phi = " << probe.getPhi() << " q/pt = " << probe.getInvQPt() << std::endl;

  //std::cout << "   Pt Resolution = " << probe.getInvQPtResolution() * 100. << " %" << std::endl
  //          << std::endl;
  return probe.getVertexSigmaXResolution();
}

void checkFT3FATPtRes(Double_t pt, Double_t eta, Double_t sigma = 8.44e-4, Double_t zField = 5.0, bool verbose = false)
{
  auto ptres_noMCS = FT3FATPtRes(pt, eta, sigma, false, zField, verbose);
  auto ptres_onlyMCS = FT3FATPtRes(pt, eta, 1e-5, true, zField, verbose);
  auto ptres_combined = sqrt(ptres_noMCS * ptres_noMCS + ptres_onlyMCS * ptres_onlyMCS);
  auto ptres_FAT = FT3FATPtRes(pt, eta, sigma, true, zField, verbose);
  std::cout << " ptres_noMCS = " << ptres_noMCS << std::endl;
  std::cout << " ptres_onlyMCS = " << ptres_onlyMCS << std::endl;
  std::cout << " ptres_combined = " << ptres_combined << std::endl;
  std::cout << " ptres_FAT = " << ptres_FAT << std::endl;
  return;
}


void checkFT3FATPtResCSV(Double_t pt, Double_t eta, Double_t sigma = 8.44e-4, Double_t zField = 5.0, bool verbose = false)
{ 
  static bool firstRun = true;
  if(firstRun) {
    std::cout << " pt" << " ; " << "eta" << " ; " << "ptres_noMCS" << " ; " 
  << "ptres_onlyMCS" << " ; " 
  << "ptres_combined" << " ; " 
  << "ptres_FAT" <<
   " ; " <<  "ptres_combined / ptres_FAT " << std::endl;
   firstRun = false;

  }
  auto ptres_noMCS = FT3FATPtRes(pt, eta, sigma, false, zField, verbose);
  auto ptres_onlyMCS = FT3FATPtRes(pt, eta, 1e-6, true, zField, verbose);
  auto ptres_combined = sqrt(ptres_noMCS * ptres_noMCS + ptres_onlyMCS * ptres_onlyMCS);
  auto ptres_FAT = FT3FATPtRes(pt, eta, sigma, true, zField, verbose);
  std::cout <<  pt  << " ; " <<  eta  << " ; " <<  ptres_noMCS  << " ; " 
  <<  ptres_onlyMCS  << " ; " 
  <<  ptres_combined  << " ; " 
  <<  ptres_FAT  <<
   " ; " << ptres_combined / ptres_FAT << std::endl;

  return;
}


//_________________________________________________________________________________________
std::vector<Double_t> getFATPtRes_pts_at_eta(std::vector<Double_t> pts, Double_t eta, Double_t sigma, bool enableMCS = true)
{
  std::vector<Double_t> ptResolutions;
  for (auto pt : pts) {
    auto ptres = FT3FATPtRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " ptRes = " << ptres << std::endl;
    ptResolutions.push_back(ptres);
  }
  return ptResolutions;
}

//_________________________________________________________________________________________
std::vector<Double_t> getFATPtRes_etas_at_pt(std::vector<Double_t> etas, Double_t pt, Double_t sigma, bool enableMCS = true)
{
  std::vector<Double_t> ptResolutions;
  for (auto eta : etas) {
    auto ptres = FT3FATPtRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " ptRes = " << ptres << std::endl;
    ptResolutions.push_back(ptres);
  }
  return ptResolutions;
}

//_________________________________________________________________________________________
std::vector<Double_t> getFATvtxXRes_etas_at_pt(std::vector<Double_t> etas, Double_t pt, Double_t sigma, bool enableMCS = true)
{
  std::vector<Double_t> vtxResolutions;
  for (auto eta : etas) {
    auto vtxRes = FT3FATvtxXRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " vtxRes = " << vtxRes << std::endl;
    vtxResolutions.push_back(vtxRes);
  }
  return vtxResolutions;
}

//_________________________________________________________________________________________
std::vector<Double_t> getFATvtxXRes_pts_at_eta(std::vector<Double_t> pts, Double_t eta, Double_t sigma, bool enableMCS = true)
{
  std::vector<Double_t> vtxResolutions;
  for (auto pt : pts) {
    auto vtxRes = FT3FATvtxXRes(pt, eta, sigma, enableMCS);
    //std::cout << "pt =  " << pt << " eta = " << eta << " vtxRes = " << vtxRes << std::endl;
    vtxResolutions.push_back(vtxRes);
  }
  return vtxResolutions;
}




void FWDFat_debug() {
Double_t pt = 1.007 ; 
Double_t eta = 0.6;
Double_t sigma = 8.44e-4; 
bool enableMCS = false; 
Double_t field = 5.0; 
bool verbose = true;
/*
Double_t ptres;
ptres = FT3FATPtRes(pt, eta, sigma, enableMCS, field, verbose);
std::cout << "FT3 FAT Pt Resolution without MCS for pt = " << pt << " ;  eta = " << eta << " ; sigma = " << sigma << " ==> " << ptres << std::endl;
eta = -eta;
ptres = FT3FATPtRes(pt, eta, sigma, enableMCS, field, verbose);
std::cout << "FT3 FAT Pt Resolution without MCS for pt = " << pt << " ;  eta = " << eta << " ; sigma = " << sigma << " ==> " << ptres << std::endl;
*/

for (auto eta : {0.2, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.4, 2.8, 3.0, 3.4, 3.6}) { 
  checkFT3FATPtResCSV(pt,eta);
}
}