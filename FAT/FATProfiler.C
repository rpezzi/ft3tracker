#include <TFile.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TH2.h>

#include "TEfficiency.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TTree.h"
#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TProfile2D.h>
#include <TMultiGraph.h>

#include "FWDFat.C"

double etaWindow = 0.05;
double ptWindow = 0.2;

//FAT configuration
float sigma = 8.44e-4; // 1e-6;// 8.44e-4
float etamin = 3.8;
float etamax = 2.6;
int nEtaPoints = 10;
float ptmin = 0.01;
float ptmax = 10;
int nPtPoints = 200;

// Estimages pt and vertexing resolution from TH3 histograms produced by FT3TrackerChecker
//
//_________________________________________________________________________________________
void FATProfiler()
{

  gStyle->SetHistLineWidth(3);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetLineWidth(3);
  gStyle->SetPalette(kRainBow);

  std::vector<double> etaList({3.7, 3.4, 3.1, 2.8});
  std::vector<double> etaListFAT({3.7 + 0.5 * etaWindow, 3.4 + 0.5 * etaWindow, 3.1 + 0.5 * etaWindow, 2.8 + 0.5 * etaWindow});
  std::vector<double> ptList({0.2, 0.4, 1., 2., 4., 9., 20.});
  std::vector<double> ptListFAT({0.2 + 0.5 * ptWindow, 0.4 + 0.5 * ptWindow, 1. + 0.5 * ptWindow, 2. + 0.5 * ptWindow, 4. + 0.5 * ptWindow, 9. + 0.5 * ptWindow, 20 + 0.5 * ptWindow});
  // Tracks with 7 hits minium
  bool Min7Hits = false;
  if (Min7Hits) {
    etaList = {4., 3.5, 2.75, 2.};
    etaListFAT = {4. + 0.5 * etaWindow, 3.5 + 0.5 * etaWindow, 2.75 + 0.5 * etaWindow, 2. + 0.5 * etaWindow};
    etamin = 4.;
    etamax = 2.0;
  }

  std::vector<float> etas;
  for (auto i = 0; i < nEtaPoints; i++) {
    etas.emplace_back(i * (etamax - etamin) / nEtaPoints + etamin);
  }

  TFile* chkFileIn = new TFile("Fittercheck_ft3tracks.root");

  // FAT pt resolution vs. pt
  auto CPtResInvPt = new TCanvas();
  TMultiGraph* ptResGraphs = new TMultiGraph("PtRes", "PtRes");

  for (auto eta : etaListFAT) {
    std::vector<float> pts;
    for (auto i = 0; i < nPtPoints; i++) {
      pts.emplace_back(i * (ptmax - ptmin) / nPtPoints + ptmin);
    }

    auto ptres = getFATPtRes_pts_at_eta(pts, eta, sigma);
    TGraph* ptResolution = new TGraph(pts.size(), &pts[0], &ptres[0]);
    ptResolution->SetTitle(Form("\\text{FAT: }\\eta = %1.3f", eta));
    ptResolution->SetDrawOption("A SAME PLC PMC");
    ptResGraphs->Add(ptResolution);
  }
  ptResGraphs->GetYaxis()->SetTitle("(q/pt) resolution");
  ptResGraphs->GetXaxis()->SetTitle("\\text{p}_t (GeV/c)");
  ptResGraphs->Draw("A PLC");

  CPtResInvPt->BuildLegend();

  // FAT pt resolution vs. eta
  auto CPtResInvPtEta = new TCanvas();
  TMultiGraph* ptResGraphsEta = new TMultiGraph("PtRes", "PtRes");

  for (auto pt : ptListFAT) {

    auto ptres = getFATPtRes_etas_at_pt(etas, pt, sigma);
    TGraph* ptResolution = new TGraph(etas.size(), &etas[0], &ptres[0]);
    ptResolution->SetTitle(Form("\\text{FAT: } p_t = %1.2f", pt));
    ptResolution->SetDrawOption("A SAME PLC PMC");
    ptResGraphsEta->Add(ptResolution);
  }

  { // FAT: Line at infinite sensor position resolution
    auto pt = 1.0;

    auto ptres = getFATPtRes_etas_at_pt(etas, pt, 1e-10);
    TGraph* ptResolution = new TGraph(etas.size(), &etas[0], &ptres[0]);
    ptResolution->SetTitle(Form("FAT: MCS only"));
    ptResolution->SetDrawOption("A SAME PLC PMC");
    ptResolution->SetLineStyle(kDashed);
    ptResGraphsEta->Add(ptResolution);
  }

  // FAT: no MCS
  for (auto pt : ptListFAT) {

    auto ptres = getFATPtRes_etas_at_pt(etas, pt, sigma, false);
    TGraph* ptResolution = new TGraph(etas.size(), &etas[0], &ptres[0]);
    ptResolution->SetTitle(Form("\\text{FAT: } p_t = %1.2f (w.o. MCS)", pt));
    ptResolution->SetDrawOption("A SAME PMC");
    ptResolution->SetLineStyle(kDashed);
    ptResGraphsEta->Add(ptResolution);
  }

  ptResGraphsEta->GetYaxis()->SetTitle("(q/pt) resolution");
  ptResGraphsEta->GetXaxis()->SetTitle("\\eta");
  ptResGraphsEta->Draw("A PLC");
  CPtResInvPtEta->SetTicky();
  CPtResInvPtEta->BuildLegend();

  /*


  // FAT X vertexing resolution vs. eta
  auto CPtResVertEta = new TCanvas();
  TMultiGraph* vtxXResGraphsEta = new TMultiGraph("PtRes", "PtRes");

  for (auto pt : ptListFAT) {

    auto ptres = getFATvtxXRes_etas_at_pt(etas, pt, sigma);
    TGraph* ptResolution = new TGraph(etas.size(), &etas[0], &ptres[0]);

    ptResolution->SetTitle(Form("\\text{FAT: } pt = %1.2f", pt));
    ptResolution->SetDrawOption("A SAME PLC PMC");

    vtxXResGraphsEta->Add(ptResolution);
  }
  vtxXResGraphsEta->GetYaxis()->SetTitle("\\sigma_x \\text{ @ vertex} (\\mu m) ");
  vtxXResGraphsEta->GetXaxis()->SetTitle("\\eta");
  vtxXResGraphsEta->Draw("A PLC");
  CPtResVertEta->BuildLegend();

  auto CVertexResEta = new TCanvas();

  // FAT vertexing resolution vs. pt
  TMultiGraph* vtxXResGraphsPt = new TMultiGraph("VtxXRes", "VtxXRes");

  for (auto eta : etaListFAT) {
    std::vector<float> pts;

    for (auto i = 0; i < nPtPoints; i++) {
      pts.emplace_back(i * (ptmax - ptmin) / nPtPoints + ptmin);
    }

    auto ptres = getFATvtxXRes_pts_at_eta(pts, eta, sigma);
    TGraph* ptResolution = new TGraph(pts.size(), &pts[0], &ptres[0]);
    ptResolution->SetTitle(Form("\\text{FAT: } \\eta = %1.3f", eta));
    ptResolution->SetDrawOption("A SAME PLC PMC");
    vtxXResGraphsPt->Add(ptResolution);
  }
  vtxXResGraphsPt->GetYaxis()->SetTitle("\\sigma_x \\text{ @ Vertex } (\\mu m)");
  vtxXResGraphsPt->GetXaxis()->SetTitle("p_t (GeV/c)");
  vtxXResGraphsPt->Draw("A PLC");

  CVertexResEta->BuildLegend();
  */
}
