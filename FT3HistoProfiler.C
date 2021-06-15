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



void profiler() {

TFile *chkFileIn = new TFile("Fittercheck_ft3tracks.root");
auto FT3TrackInvPtResolutionPtEta = (TH3F *)chkFileIn->Get("MoreHistos/FT3TrackInvPtResolutionPtEta");
auto C3D = new TCanvas();
//FT3TrackInvPtResolutionPtEta->RebinX(2);
FT3TrackInvPtResolutionPtEta->Draw();
bool first = true;
auto CPtRes = new TCanvas();

for (auto etamin : {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4}) {
auto etamax = etamin + 0.1;

FT3TrackInvPtResolutionPtEta->GetYaxis()->SetRangeUser(etamin,etamax);
auto title = Form("InvQPtRes_%1.1f_%1.1f_xz",etamin, etamax);
auto a = (TH2F *)FT3TrackInvPtResolutionPtEta->Project3D(title);


//a->SetTitle(Form("%f < \\eta < %f",etamin,etamax));
 //new TCanvas();
//a->Draw("colz");



a->FitSlicesX(0,0,-1,1);
auto th2InvPtResolution = (TH2F *)gDirectory->Get((std::string(a->GetName())+std::string("_2")).c_str());
th2InvPtResolution->SetTitle(Form("%1.1f < \\eta < %1.1f",etamin,etamax));
if(first) {
    th2InvPtResolution->SetStats(0);
    th2InvPtResolution->Draw("PLC PMC");
    first=false;
} else {
    th2InvPtResolution->SetStats(0);
    th2InvPtResolution->Draw("PLC PMC same");
}
}
CPtRes->BuildLegend();

}

void func3() {
gStyle->SetHistLineWidth(3);
gStyle->SetFrameLineWidth(3);
gStyle->SetLineWidth(3);
gStyle->SetPalette(107);
TFile *chkFileIn = new TFile("Fittercheck_ft3tracks.root");
auto FT3TrackPtResolutionPtEta = (TH3F *)chkFileIn->Get("MoreHistos/FT3TrackPtResolutionPtEta");
auto C3D = new TCanvas();
//FT3TrackPtResolutionPtEta->RebinX(2);
FT3TrackPtResolutionPtEta->Draw();
bool first = true;
auto CPtRes = new TCanvas();
int marker = kFullCircle;
for (auto etamin : {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5}) {
auto etamax = etamin + 0.1;

FT3TrackPtResolutionPtEta->GetYaxis()->SetRangeUser(etamin,etamax);
auto title = Form("PtRes_%1.1f_%1.1f_xz",etamin, etamax);
auto a = (TH2F *)FT3TrackPtResolutionPtEta->Project3D(title);


//a->SetTitle(Form("%f < \\eta < %f",etamin,etamax));
// new TCanvas();
//a->Draw("colz");


a->FitSlicesX(0,0,-1,1);
auto th2PtResolution = (TH2F *)gDirectory->Get((std::string(a->GetName())+std::string("_2")).c_str());
th2PtResolution->SetTitle(Form("%1.1f < \\eta < %1.1f",etamin,etamax));
th2PtResolution->SetMarkerStyle(marker++);
if(first) {
    th2PtResolution->SetStats(0);
    th2PtResolution->GetYaxis()->SetTitle("pt resolution");

    th2PtResolution->Draw("PLC PMC");
    first=false;
} else {
    th2PtResolution->SetStats(0);
    th2PtResolution->Draw("PLC PMC same");
}
}
CPtRes->BuildLegend();

// InvPtResolution
auto FT3TrackInvPtResolutionPtEta = (TH3F *)chkFileIn->Get("MoreHistos/FT3TrackInvPtResolutionPtEta");
//FT3TrackPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

auto C3DInvPy = new TCanvas();
//FT3TrackInvPtResolutionPtEta->RebinX(2);
FT3TrackInvPtResolutionPtEta->Draw();
first = true;
auto CPtResInvPy = new TCanvas();
marker = kFullCircle;
for (auto etamin : {2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5}) {
auto etamax = etamin + 0.1;

FT3TrackInvPtResolutionPtEta->GetYaxis()->SetRangeUser(etamin,etamax);
//FT3TrackInvPtResolutionPtEta->GetXaxis()->SetRangeUser(-1,20);

auto title = Form("InvQPtRes_%1.1f_%1.1f_xz",etamin, etamax);
auto a = (TH2F *)FT3TrackInvPtResolutionPtEta->Project3D(title);
a->GetXaxis()->SetRangeUser(-1,20);

//a->SetTitle(Form("%f < \\eta < %f",etamin,etamax));
//new TCanvas();
//a->Draw("colz");

a->FitSlicesX(0,0,-1,1);
auto th2InvPtResolution = (TH2F *)gDirectory->Get((std::string(a->GetName())+std::string("_2")).c_str());
th2InvPtResolution->SetTitle(Form("%1.1f < \\eta < %1.1f",etamin,etamax));
th2InvPtResolution->SetMarkerStyle(marker++);

if(first) {
    th2InvPtResolution->SetStats(0);
    th2InvPtResolution->GetYaxis()->SetTitle("(1/pt) resolution");
    th2InvPtResolution->Draw("PLC PMC");
    first=false;
} else {
    th2InvPtResolution->SetStats(0);
    th2InvPtResolution->Draw("PLC PMC same");
}
}
CPtResInvPy->BuildLegend();

}