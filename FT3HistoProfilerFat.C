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

std::vector<float> getFATPtRes_pts_at_eta(std::vector<float> pts, float eta);
std::vector<float> getFATPtRes_etas_at_pt(std::vector<float> etas, float pt);

// Estimages pt and vertexing resolution from TH3 histograms produced by FT3TrackerChecker
//
//_________________________________________________________________________________________
void FT3HistoProfilerFat()
{

  gStyle->SetHistLineWidth(3);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetLineWidth(3);
  gStyle->SetPalette(107);

  TFile* chkFileIn = new TFile("Fittercheck_ft3tracks.root");
  //auto C3D = new TCanvas();
  //FT3TrackPtResolutionPtEta->RebinX(2);
  //FT3TrackPtResolutionPtEta->Draw();
  bool first = true;
  int marker = kFullCircle;

  // InvPtResolution vs pt
  auto FT3TrackInvPtResolutionPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackInvPtResolutionPtEta");
  //FT3TrackPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

  //auto C3DInvPt = new TCanvas();
  //FT3TrackInvPtResolutionPtEta->RebinX(2);
  //FT3TrackInvPtResolutionPtEta->Draw();
  first = true;
  auto CPtResInvPt = new TCanvas();

  // FAT pt resolution vs. pt
  TMultiGraph* ptResGraphs = new TMultiGraph("PtRes", "PtRes");

  for (auto eta : {2.85, 2.95, 3.05, 3.15, 3.25, 3.35, 3.45, 3.55}) {
    //for (auto eta : {2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6}) {

    std::vector<float> pts;
    auto ptmin = 0.1;
    auto ptmax = 10;
    auto npoints = 100;
    for (auto i = 0; i < npoints; i++) {
      pts.emplace_back(i * (ptmax - ptmin) / npoints + ptmin);
    }

    auto ptres = getFATPtRes_pts_at_eta(pts, eta);
    TGraph* ptResolution = new TGraph(pts.size(), &pts[0], &ptres[0]);

    //ptResolution->SetMarkerStyle(21);
    if (first) {
      //ptResolution->SetDrawOption("A PLC PMC");
      //ptResolution->GetYaxis()->SetTitle("(1/pt) resolution");
      //th2InvPtResolution->SetStats(0);
      //th2InvPtResolution->GetYaxis()->SetTitle("(1/pt) resolution");
      //ptResolution->Draw("PLC PMC");
      first = false;
    }
    ptResolution->SetTitle(Form("\\text{FAT: }\\eta = %1.2f", eta));
    ptResolution->SetDrawOption("A SAME PLC PMC");

    //ptResolution->SetLineColor(2);
    //ptResolution->SetLineWidth(4);
    //ptResolution->SetFillStyle(0);
    //ptResolution->SetMarkerColor(2);
    //ptResolution->SetMarkerSize(2);

    ptResGraphs->Add(ptResolution);
  }
  ptResGraphs->GetYaxis()->SetTitle("(1/pt) resolution");
  ptResGraphs->GetXaxis()->SetTitle("p_t (GeV/c)");
  ptResGraphs->Draw("A PLC");

  marker = kFullCircle;
  for (auto etamin : {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5}) { // 10 Hits
    //for (auto etamin : {2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5}) { // 7 Hits

    auto etamax = etamin + 0.1;

    FT3TrackInvPtResolutionPtEta->GetYaxis()->SetRangeUser(etamin, etamax);
    //FT3TrackInvPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

    auto title = Form("InvQPtRes_%1.1f_%1.1f_xz", etamin, etamax);
    auto a = (TH2F*)FT3TrackInvPtResolutionPtEta->Project3D(title);
    a->GetXaxis()->SetRangeUser(-1, 20);

    //a->SetTitle(Form("%f < \\eta < %f",etamin,etamax));
    //new TCanvas();
    //a->Draw("colz");

    a->FitSlicesX(0, 0, -1, 1);
    auto th2InvPtResolution = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2InvPtResolution->SetTitle(Form("%1.1f < \\eta < %1.1f", etamin, etamax));
    th2InvPtResolution->SetMarkerStyle(marker++);

    th2InvPtResolution->SetStats(0);
    th2InvPtResolution->Draw("PLC PMC same");
  }
  CPtResInvPt->BuildLegend();

  // FAT pt resolution vs. eta

  auto CPtResInvPtEta = new TCanvas();
  TMultiGraph* ptResGraphsEta = new TMultiGraph("PtRes", "PtRes");

  for (auto pt : {1.5, 2.5, 4.5, 6.5, 8.5, 9.5}) {
    //for (auto eta : {2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6}) {

    std::vector<float> etas;
    auto etamin = 3.6;
    auto etamax = 2.8;
    auto npoints = 10;
    for (auto i = 0; i < npoints; i++) {
      etas.emplace_back(i * (etamax - etamin) / npoints + etamin);
    }

    auto ptres = getFATPtRes_etas_at_pt(etas, pt);
    TGraph* ptResolution = new TGraph(etas.size(), &etas[0], &ptres[0]);

    //ptResolution->SetMarkerStyle(21);
    if (first) {
      //ptResolution->SetDrawOption("A PLC PMC");
      //ptResolution->GetYaxis()->SetTitle("(1/pt) resolution");
      //th2InvPtResolution->SetStats(0);
      //th2InvPtResolution->GetYaxis()->SetTitle("(1/pt) resolution");
      //ptResolution->Draw("PLC PMC");
      first = false;
    }
    ptResolution->SetTitle(Form("\\text{FAT: } pt = %1.1f", pt));
    ptResolution->SetDrawOption("A SAME PLC PMC");

    //ptResolution->SetLineColor(2);
    //ptResolution->SetLineWidth(4);
    //ptResolution->SetFillStyle(0);
    //ptResolution->SetMarkerColor(2);
    //ptResolution->SetMarkerSize(2);

    ptResGraphsEta->Add(ptResolution);
  }
  ptResGraphsEta->GetYaxis()->SetTitle("(1/pt) resolution");
  ptResGraphsEta->GetXaxis()->SetTitle("\\eta");
  ptResGraphsEta->Draw("A PLC");

  // InvPtResolution vs eta
  //FT3TrackPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);
  FT3TrackInvPtResolutionPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackInvPtResolutionPtEta");
  FT3TrackInvPtResolutionPtEta->GetYaxis()->SetRange(0, 0);
  //FT3TrackInvPtResolutionPtEta->GetXaxis()->SetRange(2.8,3.6);

  //auto C3DInvPtEta = new TCanvas();
  //FT3TrackInvPtResolutionPtEta->RebinX(2);
  //FT3TrackInvPtResolutionPtEta->Draw();
  first = true;
  marker = kFullCircle;
  for (auto ptmin : {1., 2., 4., 6., 8., 9.}) {
    auto ptmax = ptmin + 1.;

    FT3TrackInvPtResolutionPtEta->GetXaxis()->SetRangeUser(ptmin, ptmax);
    //FT3TrackInvPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

    auto title = Form("InvQPtResEta_%1.1f_%1.1f_yz", ptmin, ptmax);
    auto a = (TH2F*)FT3TrackInvPtResolutionPtEta->Project3D(title);
    a->GetXaxis()->SetRangeUser(-1, 10);

    //a->SetTitle(Form("%f < p_t < %f",ptmin,ptmax));
    //new TCanvas();
    //a->Draw("colz");

    a->FitSlicesX(0, 0, -1, 1);
    auto th2InvPtResolutionEta = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2InvPtResolutionEta->SetTitle(Form("%1.1f < p_t < %1.1f", ptmin, ptmax));
    th2InvPtResolutionEta->SetMarkerStyle(marker++);

    th2InvPtResolutionEta->SetStats(0);
    th2InvPtResolutionEta->Draw("PLC PMC same");
  }
  CPtResInvPtEta->BuildLegend();
}

std::vector<float> getFATPtRes_pts_at_eta(std::vector<float> pts, float eta)
{
  std::vector<float> ptResolutions;
  for (auto pt : pts) {
    auto ptres = FT3FATPtRes(pt, eta);
    std::cout << "pt =  " << pt << " eta = " << eta << " ptRes = " << ptres << std::endl;
    ptResolutions.push_back(ptres);
  }
  return ptResolutions;
}

std::vector<float> getFATPtRes_etas_at_pt(std::vector<float> etas, float pt)
{
  std::vector<float> ptResolutions;
  for (auto eta : etas) {
    auto ptres = FT3FATPtRes(pt, eta);
    std::cout << "pt =  " << pt << " eta = " << eta << " ptRes = " << ptres << std::endl;
    ptResolutions.push_back(ptres);
  }
  return ptResolutions;
}
