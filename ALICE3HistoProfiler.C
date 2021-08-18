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

// Estimages pt and vertexing resolution from TH3 histograms produced by FT3TrackerChecker
//
//_________________________________________________________________________________________
void ALICE3HistoProfiler(const Char_t* FittercheckFile = "Fittercheck_ALICE3tracks.root")
{
  gStyle->SetHistLineWidth(3);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetLineWidth(3);
  gStyle->SetPalette(107);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridY(1);
  TFile* chkFileIn = new TFile(FittercheckFile);
  auto FT3TrackPtResolutionPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackPtResolutionPtEta");
  //auto C3D = new TCanvas();
  //FT3TrackPtResolutionPtEta->RebinX(2);
  //FT3TrackPtResolutionPtEta->Draw();
  bool first = true;
  auto CPtRes = new TCanvas();
  int marker = kFullCircle;
  for (auto etamin : {0.5, 0.9, 1.9, 2.9, 3.8}) {
    auto etamax = etamin + 0.1;

    FT3TrackPtResolutionPtEta->GetYaxis()->SetRangeUser(etamin, etamax);
    auto title = Form("PtRes_%1.1f_%1.1f_xz", etamin, etamax);
    auto a = (TH2F*)FT3TrackPtResolutionPtEta->Project3D(title);

    //a->SetTitle(Form("%f < \\eta < %f",etamin,etamax));
    // new TCanvas();
    //a->Draw("colz");

    a->FitSlicesX(0, 0, -1, 1);
    auto th2PtResolution = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2PtResolution->SetTitle(Form("%1.1f < \\eta < %1.1f", etamin, etamax));
    th2PtResolution->SetMarkerStyle(marker++);
    if (first) {
      th2PtResolution->SetStats(0);
      th2PtResolution->GetYaxis()->SetTitle("pt resolution");

      th2PtResolution->Draw("PLC PMC");
      first = false;
    } else {
      th2PtResolution->SetStats(0);
      th2PtResolution->Draw("PLC PMC same");
    }
  }
  CPtRes->BuildLegend();

  // InvPtResolution vs pt
  auto FT3TrackInvQPtResolutionPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackInvQPtResolutionPtEta");

  first = true;
  auto CPtResInvPy = new TCanvas();
  marker = kFullCircle;
  for (auto etamin : {0.5, 0.9, 1.9, 2.9, 3.8}) {
    auto etamax = etamin + 0.1;

    FT3TrackInvQPtResolutionPtEta->GetYaxis()->SetRangeUser(etamin, etamax);
    //FT3TrackInvQPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

    auto title = Form("InvQPtRes_%1.1f_%1.1f_xz", etamin, etamax);
    auto a = (TH2F*)FT3TrackInvQPtResolutionPtEta->Project3D(title);
    a->GetXaxis()->SetRangeUser(-20, 20);

    a->FitSlicesX(0, 0, -1, 1);
    auto th2InvPtResolution = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2InvPtResolution->SetTitle(Form("%1.1f < \\eta < %1.1f", etamin, etamax));
    th2InvPtResolution->SetMarkerStyle(marker++);

    if (first) {
      th2InvPtResolution->SetStats(0);
      th2InvPtResolution->GetYaxis()->SetTitle("(q/pt) resolution");
      th2InvPtResolution->Draw("PLC PMC");
      first = false;
    } else {
      th2InvPtResolution->SetStats(0);
      th2InvPtResolution->Draw("PLC PMC same");
    }
  }
  CPtResInvPy->BuildLegend();

  // InvPtResolution vs eta
  //FT3TrackPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);
  FT3TrackInvQPtResolutionPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackInvQPtResolutionPtEta");
  FT3TrackInvQPtResolutionPtEta->GetYaxis()->SetRange(0, 0);
  //FT3TrackInvQPtResolutionPtEta->GetXaxis()->SetRange(2.8,3.6);

  //auto C3DInvPyEta = new TCanvas();
  //FT3TrackInvQPtResolutionPtEta->RebinX(2);
  //FT3TrackInvQPtResolutionPtEta->Draw();
  first = true;
  auto CPtResInvPyEta = new TCanvas();
  marker = kFullCircle;
  for (auto ptmin : {1., 2., 4., 6., 8., 9.}) {
    auto ptmax = ptmin + 1.;

    FT3TrackInvQPtResolutionPtEta->GetXaxis()->SetRangeUser(ptmin, ptmax);
    //FT3TrackInvQPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

    auto title = Form("InvQPtResEta_%1.1f_%1.1f_yz", ptmin, ptmax);
    auto a = (TH2F*)FT3TrackInvQPtResolutionPtEta->Project3D(title);
    //a->GetXaxis()->SetRangeUser(-10, 10);

    //a->SetTitle(Form("%f < p_t < %f",ptmin,ptmax));
    //new TCanvas();
    //a->Draw("colz");

    a->FitSlicesX(0, 0, -1, 1);
    auto th2InvPtResolutionEta = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2InvPtResolutionEta->SetTitle(Form("%1.1f < p_t < %1.1f", ptmin, ptmax));
    th2InvPtResolutionEta->SetMarkerStyle(marker++);

    if (first) {
      th2InvPtResolutionEta->SetStats(0);
      th2InvPtResolutionEta->GetYaxis()->SetTitle("(q/pt) resolution");
      th2InvPtResolutionEta->Draw("PLC PMC");
      first = false;
    } else {
      th2InvPtResolutionEta->SetStats(0);
      th2InvPtResolutionEta->Draw("PLC PMC same");
    }
  }
  CPtResInvPyEta->BuildLegend();

  //
  // Vertexing resolution vs eta
  //FT3TrackPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);
  auto FT3TrackDeltaXVertexPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackDeltaXVertexPtEta");
  FT3TrackDeltaXVertexPtEta->GetYaxis()->SetRange(0, 0);
  //FT3TrackDeltaXVertexPtEta->GetXaxis()->SetRange(2.8,3.6);

  //auto C3DInvPyEta = new TCanvas();
  FT3TrackDeltaXVertexPtEta->RebinY(2);
  //FT3TrackDeltaXVertexPtEta->Draw();
  first = true;
  auto CPtResVertEta = new TCanvas();
  marker = kFullCircle;
  //for (auto ptmin : {1., 2., 4., 6., 8., 9.}) {
  for (auto ptmin : {1., 2., 4., 6., 9., 15., 20.}) {

    auto ptmax = ptmin + 0.5;

    FT3TrackDeltaXVertexPtEta->GetXaxis()->SetRangeUser(ptmin, ptmax);
    //FT3TrackDeltaXVertexPtEta->GetXaxis()->SetRangeUser(-1,20);

    auto title = Form("VertXResEta_%1.1f_%1.1f_yz", ptmin, ptmax);
    auto a = (TH2F*)FT3TrackDeltaXVertexPtEta->Project3D(title);
    //a->GetXaxis()->SetRangeUser(-1, 10);

    //a->SetTitle(Form("%f < p_t < %f",ptmin,ptmax));
    //new TCanvas();
    //a->Draw("colz");

    a->FitSlicesX(0, 0, -1, 1);
    auto th2VertEta = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2VertEta->SetTitle(Form("%1.1f < p_t < %1.1f", ptmin, ptmax));
    th2VertEta->SetMarkerStyle(marker++);

    if (first) {
      th2VertEta->SetStats(0);
      th2VertEta->GetYaxis()->SetTitle("\\sigma_x \\text{ @ Vertex } (\\mu m)");
      th2VertEta->Draw("PLC PMC");
      first = false;
    } else {
      th2VertEta->SetStats(0);
      th2VertEta->Draw("PLC PMC same");
    }
  }
  CPtResVertEta->BuildLegend();

  // Vertexing resolution vs pt
  FT3TrackDeltaXVertexPtEta = (TH3F*)chkFileIn->Get("MoreHistos/FT3TrackDeltaXVertexPtEta");
  FT3TrackDeltaXVertexPtEta->GetYaxis()->SetRange(0, 0);
  FT3TrackDeltaXVertexPtEta->GetXaxis()->SetRange(0, 0);

  //auto C3DInvPy = new TCanvas();
  //FT3TrackInvQPtResolutionPtEta->RebinX(2);
  //FT3TrackInvQPtResolutionPtEta->Draw();
  first = true;
  auto CVertexResEta = new TCanvas();
  marker = kFullCircle;
  for (auto etamin : {0.5, 0.9, 1.9, 2.9, 3.8}) {
    auto etamax = etamin + 0.1;

    FT3TrackDeltaXVertexPtEta->GetYaxis()->SetRangeUser(etamin, etamax);
    //FT3TrackInvQPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

    auto title = Form("VertResPt_%1.1f_%1.1f_xz", etamin, etamax);
    auto a = (TH2F*)FT3TrackDeltaXVertexPtEta->Project3D(title);
    //a->GetXaxis()->SetRangeUser(-1,20);

    a->SetTitle(title);
    //new TCanvas();
    //a->Draw("colz");

    a->FitSlicesX(0, 0, -1, 1);
    auto th2VertResolutionPt = (TH2F*)gDirectory->Get((std::string(a->GetName()) + std::string("_2")).c_str());
    th2VertResolutionPt->SetTitle(Form("%1.1f < \\eta < %1.1f", etamin, etamax));
    th2VertResolutionPt->SetMarkerStyle(marker++);

    if (first) {
      th2VertResolutionPt->SetStats(0);
      th2VertResolutionPt->GetYaxis()->SetTitle("\\sigma_x \\text{ @ Vertex } (\\mu m)");
      th2VertResolutionPt->Draw("PLC PMC");
      first = false;
    } else {
      th2VertResolutionPt->SetStats(0);
      th2VertResolutionPt->Draw("PLC PMC same");
    }
  }
  CVertexResEta->BuildLegend();
}