#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TFile.h>
#include <TTree.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TCut.h>
#include <TText.h>
#include <TLine.h>
#include <TPaveText.h>
#include "ERTools.h"
#include "EdbLog.h"

using namespace std;

struct
{
  bool process = false;     // process raw.root file
  bool interactive = false; // run interactive session
} DO;

struct
{
  string filename;
  TDatime file_creation_date;
  TDatime file_modification_date;
  Long_t file_size;
  float scan_time;    // total scanning time in seconds
  float scanned_area; // total scanned area in cm2 (binned)
  float empty_frac_top;
  float empty_frac_bot;
  float density_top; // density as total/area excluding lost (zero) bins [/mm2]
  float density_bot;
  float requested_area; // requested area in cm2
  float covered_area;   // covered area in cm2 - really scanned from the first to the last scanned point
} RES;                  // Result to be exported as report

struct
{
  float fraction;                   // fraction of the area to be covered
  float xminR, xmaxR, yminR, ymaxR; // requested area
  float areaR;                      // requested area
  float fractionS;                  // fraction of the scanned area
  float xminS, xmaxS, yminS, ymaxS; // scanned area
  float areaS;                      // scanned area
} AREA;                             // requested and scanned area

struct
{
  int nx;
  float xmin, xmax, xbin;
  int ny;
  float ymin, ymax, ybin;
} B; // view step binning

struct
{
  int prof = 1;
  float X0 = 100000;
  float Y0 = 100000;
  float width = 1000;
} P; // profile definition

struct
{
  float dzx1;
  float dzx2;
  float dzy1;
  float dzy2;
  float szx1;
  float szx2;
  float szy1;
  float szy2;
} Z; // mean dz calculated from profiles

struct
{
  TProfile2D *nseg_top;
  TProfile2D *nseg_bot;
  TProfile2D *thick_top;
  TProfile2D *thick_bot;
  TH2D *thick_base;
  TProfile2D *glass;
} H; // oputput histos

struct
{
  float bin = 20;
} FFT;

using namespace ERTools;

TEnv cenv("qualityenv");

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout << "\nUsage: \n\t  emthickness -input=file.raw.root [ -out=report -I -v=DEBUG] \n";
  cout << "\t  Generate quality report for raw data\n";
  cout << "\t\t -input      -  input raw root file with tracks\n";
  cout << "\t\t -out        -  output report filename (default: raw_quality_report[.json],[.png] )\n";
  cout << "\t\t -I          -  open interactive root session\n";
  cout << "\t\t -v          -  verbosity level\n";
  cout << endl;
}

//----------------------------------------------------------------------------------------
void set_default_env()
{
  cenv.SetValue("quality.input", "");
  cenv.SetValue("quality.output", "raw_quality_report");

  cenv.SetValue("quality.request.Xmin", 5000.);
  cenv.SetValue("quality.request.Ymin", 5000.);
  cenv.SetValue("quality.request.Xmax", 185000.);
  cenv.SetValue("quality.request.Ymax", 185000.);
  cenv.SetValue("quality.request.fraction", 0.9);

  cenv.SetValue("quality.prof", 1);
  cenv.SetValue("quality.prof.X0", 100000);
  cenv.SetValue("quality.prof.Y0", 100000);
  cenv.SetValue("quality.prof.width", 1000);

  cenv.SetValue("quality.frameal", 1);

  cenv.SetValue("quality.EdbDebugLevel", 1);
}

//----------------------------------------------------------------------------------------
void generate_json_report(TH2 *h, std::ofstream &report, bool first_plot)
{
  if (!first_plot)
    report << ",\n";

  // Calculate statistics
  double sum = 0, sum2 = 0;
  int total_bins = h->GetNbinsX() * h->GetNbinsY();
  int empty_bins = 0;
  int count = 0;

  for (int binx = 1; binx <= h->GetNbinsX(); ++binx)
  {
    for (int biny = 1; biny <= h->GetNbinsY(); ++biny)
    {
      if (h->GetBinContent(binx, biny) == 0)
      {
        empty_bins++;
        continue;
      }
      double val = h->GetBinContent(binx, biny);
      sum += val;
      sum2 += val * val;
      count++;
    }
  }

  double mean = (count > 0) ? sum / count : 0;
  double rms = (count > 0) ? TMath::Sqrt(sum2 / count - mean * mean) : 0;
  double empty_frac = (double)empty_bins / total_bins;
  double density = mean / B.xbin / B.ybin * 1000 * 1000; // density in mm2

  if (!strcmp(h->GetName(), "nseg_top"))
  {
    RES.empty_frac_top = empty_frac;
    RES.density_top = density;
  }
  if (!strcmp(h->GetName(), "nseg_bot"))
  {
    RES.empty_frac_bot = empty_frac;
    RES.density_bot = density;
  }

  int nsig = 6;
//  h->GetZaxis()->SetRangeUser(TMath::Max(0., TMath::Max(mean - nsig/2 * rms, h->GetMinimum())),
//                              TMath::Min(3 * mean, TMath::Min(mean + nsig * rms, h->GetMaximum())));
  h->GetZaxis()->SetRangeUser(TMath::Max(mean - nsig/2 * rms, h->GetMinimum()),
                              TMath::Min(mean + nsig * rms, h->GetMaximum()));

  h->SetTitle(Form("%s: mean = %.2f  RMS = %.2f", h->GetTitle(), mean, rms));

  // JSON Output
  report << "    {\n";
  report << "      \"name\": \"" << h->GetName() << "\",\n";
  report << "      \"title\": \"" << h->GetTitle() << "\",\n";
  report << "      \"statistics\": {\n";
  report << "        \"mean\": " << mean << ",\n";
  report << "        \"rms\": " << rms << ",\n";
  report << "        \"empty_bins\": " << empty_bins << ",\n";
  report << "        \"total_bins\": " << total_bins << ",\n";
  report << "        \"empty_fraction\": " << empty_frac << "\n";
  report << "      }\n";
  report << "    }";
}

//----------------------------------------------------------------------------------------
TH1D *get_step_fft(TH1D *h)
{
  // Perform the Fourier Transform
  TH1 *h_fft = h->FFT(nullptr, Form("MAG_%s", h->GetName())); // Compute magnitude of FFT
  h_fft->SetTitle("Fourier Transform;Frequency;Magnitude");

  auto peaks = FindPeaksWithIntegral(h_fft, 5, 0.1, 2, h_fft->GetNbinsX() / 2); // 5-bin window

  // Convert bin to frequency
  double frequency = peaks[0].peak_position;
  float xmin = h->GetXaxis()->GetXmin();
  float xmax = h->GetXaxis()->GetXmax();
  double step_size = (xmax - xmin) / (frequency - 1); // Step size
  if (gEDBDEBUGLEVEL > 1)
  {
    for (int i = 0; i < 3; i++)
      printf("%d %f %f %f\n",
             i, peaks[i].peak_position, peaks[i].peak_height, peaks[i].window_integral);
    std::cout << "Estimated step size: " << step_size << std::endl;
  }

  int n = TMath::Ceil((xmax - xmin) / step_size);
  xmin -= step_size / 2;
  xmax += step_size / 2;
  n += 1;

  TH1D *hstep = new TH1D(Form("hstep_%s", h->GetName()), "hstep", n, xmin, xmax);

  if (0)
  {
    // Visualize the results
    TCanvas *c = new TCanvas(h->GetName(), "Step Size Estimation with FFT", 1200, 600);
    c->Divide(2, 1);
    c->cd(1);
    h->Draw(); // Original binned data
    c->cd(2);
    h_fft->Draw(); // Fourier Transform
  }
  return hstep;
}

//----------------------------------------------------------------------------------------
void define_steps(TTree *tree)
{
  TH1D *hx = get_h_var(tree, "headers.eXview", "hx", FFT.bin, "!(eXview==0&&eYview==0)");
  TH1D *hy = get_h_var(tree, "headers.eYview", "hy", FFT.bin, "!(eXview==0&&eYview==0)");
  // hx->Draw();

  TH1D *hstepx = get_step_fft(hx);
  TH1D *hstepy = get_step_fft(hy);
  B.nx = hstepx->GetXaxis()->GetNbins();
  B.xmin = hstepx->GetXaxis()->GetXmin();
  B.xmax = hstepx->GetXaxis()->GetXmax();
  B.xbin = (B.xmax - B.xmin) / B.nx;
  B.xmin -= 5 * B.xbin;
  B.nx += 5;
  B.xmax += 5 * B.xbin;
  B.nx += 5;
  B.ny = hstepy->GetXaxis()->GetNbins();
  B.ymin = hstepy->GetXaxis()->GetXmin();
  B.ymax = hstepy->GetXaxis()->GetXmax();
  B.ybin = (B.ymax - B.ymin) / B.ny;
  B.ymin -= 5 * B.ybin;
  B.ny += 5;
  B.ymax += 5 * B.ybin;
  B.ny += 5;

  AREA.xminS = hx->GetXaxis()->GetXmin() - B.xbin;
  AREA.xmaxS = hx->GetXaxis()->GetXmax() + B.xbin;
  AREA.yminS = hy->GetXaxis()->GetXmin() - B.ybin;
  AREA.ymaxS = hy->GetXaxis()->GetXmax() + B.ybin;
  AREA.areaS = (AREA.xmaxS - AREA.xminS) * (AREA.ymaxS - AREA.yminS);
  AREA.fractionS = AREA.areaS / AREA.areaR;
  Log(1, "define_steps",
      "Determined step sizes:  B.xbin= %.1f  B.ybin= %.1f\n"
      "Scanned area: X: %.1f , %.1f  Y: %.1f , %.1f   Area= %.1f cm2  Fraction= %.3f",
      B.xbin, B.ybin,
      AREA.xminS, AREA.xmaxS, AREA.yminS, AREA.ymaxS, AREA.areaS / 10000. / 10000., AREA.fractionS);
}

//----------------------------------------------------------------------------------------
void make_histos(TTree *tree)
{
  TCut sideTop = "eNframesTop>0";
  TCut sideBot = "eNframesBot>0";
  TProfile2D *hview = new TProfile2D("hview", "hview", B.nx, B.xmin, B.xmax, B.ny, B.ymin, B.ymax);

  H.nseg_top = (TProfile2D *)(hview->Clone("nseg_top"));
  H.nseg_top->SetTitle("nseg_top");
  H.nseg_bot = (TProfile2D *)(hview->Clone("nseg_bot"));
  H.nseg_bot->SetTitle("nseg_bot");
  H.thick_top = (TProfile2D *)(hview->Clone("thick_top"));
  H.thick_top->SetTitle("thick_top");
  H.thick_bot = (TProfile2D *)(hview->Clone("thick_bot"));
  H.thick_bot->SetTitle("thick_bot");
  H.glass = (TProfile2D *)(hview->Clone("glass"));
  H.glass->SetTitle("glass");
  TProfile2D *hz2 = (TProfile2D *)(hview->Clone("hz2"));
  TProfile2D *hz3 = (TProfile2D *)(hview->Clone("hz3"));
  H.thick_base = new TH2D("thick_base", "thick_base", B.nx, B.xmin, B.xmax, B.ny, B.ymin, B.ymax);

  tree->Draw("eNsegments:eYview:eXview>>nseg_top", sideTop, "goff prof colz");
  tree->Draw("eNsegments:eYview:eXview>>nseg_bot", sideBot, "goff prof colz");
  tree->Draw("eZ1-eZ2:eYview:eXview>>thick_top", sideTop, "goff prof colz");
  tree->Draw("eZ3-eZ4:eYview:eXview>>thick_bot", sideBot, "goff prof colz");
  tree->Draw("eZ2:eYview:eXview>>hz2", sideTop, "goff prof colz");
  tree->Draw("eZ3:eYview:eXview>>hz3", sideBot, "goff prof colz");
  tree->Draw("eZ4:eYview:eXview>>glass", sideBot, "goff prof colz");

  DiffProfile2D(hz2, hz3, H.thick_base);

  RES.scanned_area = AREA.areaS / 10000. / 10000.; // cm2

  tree->Draw("eEvent>>htime(10000)", "eViewID>0", "goff");
  TH1 *htime = (TH1 *)(gDirectory->Get("htime"));
  if (htime)
    RES.scan_time = htime->Integral() * htime->GetMean() / 1000.; // scan time in seconds

  SafeDelete(hview);
  SafeDelete(hz2);
  SafeDelete(hz3);
}

void make_canvas(const char *nameo = "ccc")
{
  gStyle->SetNumberContours(256);
  gStyle->SetPalette(107);
  gStyle->SetOptStat("");
  gStyle->SetPadRightMargin(0.12);

  TCanvas *cc = new TCanvas("thickness", Form("thickness at %s", nameo), 1920, 1080);

  TPad *header = new TPad("header", "Header", 0.0, 0.91, 1, 1); // xlow, ylow, xup, yup
  header->Draw();
  header->cd();
  TPaveText *tp = new TPaveText(0.01, 0.01, 0.99, 0.99, "NDC");
  tp->AddText(Form("%s  of  %s   %.1f GB",
                   RES.filename.c_str(),
                   RES.file_creation_date.AsString(),
                   RES.file_size / 1024. / 1024. / 1024.));
  tp->AddText(Form("Xrange:  %.1f   %.1f     Yrange:  %.1f   %.1f     Step:  %.1f   %.1f   ScanTime: %.2f h",
                   B.xmin, B.xmax, B.ymin, B.ymax, B.xbin, B.ybin, RES.scan_time / 60. / 60.));
  tp->AddText(Form(
      "Requested/Scanned:  %.1f / %.1f cm2    EmptyTop: %.2f %%  EmptyBot: %.2f %%   Density/mm2: Top = %.1f  Bot = %.1f",
      RES.requested_area, RES.scanned_area, RES.empty_frac_top * 100, RES.empty_frac_bot * 100, RES.density_top, RES.density_bot));
  if (AREA.fractionS * (1. - RES.empty_frac_top) < AREA.fraction || AREA.fractionS * (1. - RES.empty_frac_bot) < AREA.fraction)
  {
    tp->SetTextColor(kRed);
    // tp->AddText("WARNING: Scanned area or empty fraction below requested!");
  }
  tp->Draw();

  cc->cd(0);
  TPad *pad = new TPad("pad", "", 0., 0., 1., 0.9);
  pad->Draw();
  pad->cd();
  pad->Divide(3, 2, 0.01, 0.01);

  pad->cd(1);
  H.nseg_top->Draw("prof colz");
  pad->cd(4);
  H.nseg_bot->Draw("prof colz");
  pad->cd(2);
  H.thick_top->Draw("prof colz");
  pad->cd(5);
  H.thick_bot->Draw("prof colz");
  pad->cd(3);
  H.thick_base->Draw("prof colz");
  pad->cd(6);
  H.glass->Draw("prof colz");
  pad->cd(0);
  TDatime time;
  TText *t = new TText();
  t->SetTextSize(0.015);
  t->DrawText(0.25, 0.0001, Form("%s/%s    %s", gSystem->WorkingDirectory(), RES.filename.c_str(), time.AsString()));
  if (gROOT->IsBatch())
    cc->SaveAs(Form("%s.png", nameo));
}

int make_report(const char *output_file)
{
  // Open JSON report
  std::ofstream report(Form("%s.json", output_file));
  if (!report.is_open())
  {
    Log(1, "quality_report", "Error opening output file: %s", output_file);
    return 3;
  }

  RES.requested_area = AREA.areaR / 10000. / 10000.; // cm2
  RES.covered_area = AREA.areaS / 10000. / 10000.;   // cm2

  report << "{\n";
  report << "  \"metadata\": {\n";
  report << "    \"input_file\": \"" << RES.filename << "\",\n";
  report << "    \"file_creation_date\": \"" << RES.file_creation_date.AsString() << "\",\n";
  report << "    \"file_modification_date\": \"" << RES.file_modification_date.AsString() << "\",\n";
  report << "    \"file_size\": \"" << RES.file_size << "\",\n";
  report << "    \"timestamp\": \"" << TDatime().AsString() << "\"\n";
  report << "  },\n";
  report << "  \"plots\": [\n";

  // Process all 2D histograms
  bool first_plot = true;
  TIter next(gDirectory->GetList());
  while (TObject *obj = next())
  {
    if (obj->InheritsFrom("TH2") || obj->InheritsFrom("TProfile2D"))
    {
      generate_json_report((TH2 *)obj, report, first_plot);
      first_plot = false;
    }
  }
  report << "\n  ],\n";
  report << "  \"result\": {\n";
  report << "    \"scan_time\": \"" << RES.scan_time << "\",\n";
  report << "    \"scanned_area\": \"" << RES.scanned_area << "\",\n";
  report << "    \"empty_frac_top\": \"" << RES.empty_frac_top << "\",\n";
  report << "    \"empty_frac_bot\": \"" << RES.empty_frac_bot << "\",\n";
  report << "    \"density_top\": \"" << RES.density_top << "\",\n";
  report << "    \"density_bot\": \"" << RES.density_bot << "\"\n";
  report << "  }\n";
  report << "}\n";
  report.close();

  Log(2, "quality_report", "Report generated: %s", output_file);
  return 0;
}

//----------------------------------------------------------------------------------------
void draw_frame_align()
{
  TCut cut("cut", "nsg>2000");
  TTree *tfa = (TTree *)gFile->Get("FrameAlign");
  if (tfa)
  {
    TH1D *z1s1 = get_h_var(tfa, "z1", "z1s1", 1., "side==0");
    TH1D *z1s2 = get_h_var(tfa, "z1", "z1s2", 1., "side==1");

    TH2D *dxs1 = new TH2D("dxs1", "dx:z  side1", 1000, z1s1->GetXaxis()->GetXmin(), z1s1->GetXaxis()->GetXmax(), 500, -2., 2.);
    TH2D *dxs2 = new TH2D("dxs2", "dx:z  side2", 1000, z1s2->GetXaxis()->GetXmin(), z1s2->GetXaxis()->GetXmax(), 500, -2., 2.);
    TH2D *dys1 = new TH2D("dys1", "dy:z  side1", 1000, z1s1->GetXaxis()->GetXmin(), z1s1->GetXaxis()->GetXmax(), 500, -2., 2.);
    TH2D *dys2 = new TH2D("dys2", "dy:z  side2", 1000, z1s2->GetXaxis()->GetXmin(), z1s2->GetXaxis()->GetXmax(), 500, -2., 2.);

    TCanvas *c = new TCanvas("fral", "frames align", 1000, 800);
    c->Divide(2, 2);
    c->cd(1)->SetGrid();
    tfa->Draw("dx:z1>>dxs1", cut && "side==0", "colz");
    c->cd(2)->SetGrid();
    tfa->Draw("dy:z1>>dys1", cut && "side==0", "colz");
    c->cd(3)->SetGrid();
    tfa->Draw("dx:z1>>dxs2", cut && "side==1", "colz");
    c->cd(4)->SetGrid();
    tfa->Draw("dy:z1>>dys2", cut && "side==1", "colz");
  }
}

void SetLargeLabelsStyle()
{
  gStyle->SetOptStat("e");
  // Axis titles and labels
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetLabelSize(0.06, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  // Statistics box (if enabled)
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatFont(42);
  // Legend
  gStyle->SetLegendTextSize(0.05);
  gStyle->SetLegendFont(42);
}

//----------------------------------------------------------------------------------------
void make_snd_profiles(const char *nameo = "snd_profiles")
{
  gStyle->SetNumberContours(256);
  gStyle->SetPalette(1);
  gStyle->SetOptStat("");
  gStyle->SetPadRightMargin(0.12);
  // gStyle->SetPadTopMargin(0.12);

  P.X0 = cenv.GetValue("quality.prof.X0", 100000);
  P.Y0 = cenv.GetValue("quality.prof.Y0", 100000);
  P.width = cenv.GetValue("quality.prof.width", 1000);

  TTree *Views = (TTree *)gFile->Get("Views");
  TH1D *hz4 = get_h_var(Views, "headers.eZ4", "hz4", 1, "eNframesBot>0");
  float z4mean = hz4->GetMean();
  float glassmin = z4mean - 100;
  float glassmax = z4mean + 500;

  cout << "Check scanned data surface and profiles" << endl;

  TCut cTop = "eNframesTop!=0";
  TCut cBot = "eNframesBot!=0";
  TCut cutXline = Form("abs(eYview-%.1f)<%.1f", P.Y0, P.width);
  TCut cutYline = Form("abs(eXview-%.1f)<%.1f", P.X0, P.width);

  TCanvas *cc = new TCanvas("check_snd", gFile->GetName(), 1920, 1080);
  TPad *header = new TPad("header", "Header", 0.005, 0.94, 1, 1); // xlow, ylow, xup, yup
  header->Draw();
  header->cd();
  TPaveText *tp = new TPaveText(0.01, 0.01, 0.99, 0.99, "NDC");
  tp->AddText(Form("%s  of  %s", RES.filename.c_str(), RES.file_creation_date.AsString()));
  tp->AddText(Form("profiles at %.1f %.1f +- %.1f ", P.X0, P.Y0, P.width));
  tp->Draw();

  cc->cd(0);
  TPad *pad = new TPad("pad", "", 0., 0., 1., 0.94);
  pad->Draw();
  pad->cd();
  pad->Divide(3, 2, 0.01, 0.01);

  TLine *lineX1 = new TLine(B.xmin, P.Y0 - P.width, B.xmax, P.Y0 - P.width);
  lineX1->SetLineStyle(3);
  TLine *lineX2 = new TLine(B.xmin, P.Y0 + P.width, B.xmax, P.Y0 + P.width);
  lineX2->SetLineStyle(3);
  TLine *lineY1 = new TLine(P.X0 - P.width, B.ymin, P.X0 - P.width, B.ymax);
  lineY1->SetLineStyle(3);
  TLine *lineY2 = new TLine(P.X0 + P.width, B.ymin, P.X0 + P.width, B.ymax);
  lineY2->SetLineStyle(3);
  pad->cd(1);
  Views->Draw(Form("eNsegments:eYview:eXview>>hxy1(%d,%f,%f,%d,%f,%f)", B.nx, B.xmin, B.xmax, B.ny, B.ymin, B.ymax), cTop, "prof colz");
  lineX1->Draw();
  lineX2->Draw();
  lineY1->Draw();
  lineY2->Draw();
  pad->cd(4);
  Views->Draw(Form("eNsegments:eYview:eXview>>hxy2(%d,%f,%f,%d,%f,%f)", B.nx, B.xmin, B.xmax, B.ny, B.ymin, B.ymax), cBot, "prof colz");
  lineX1->Draw();
  lineX2->Draw();
  lineY1->Draw();
  lineY2->Draw();

  gStyle->SetOptStat("n");
  pad->cd(2)->SetGrid();
  Views->Draw(Form("eNcl:eZframe:eXview>>hpx(%d,%f,%f,%d,%f,%f)", B.nx, B.xmin, B.xmax, 1200, glassmin, glassmax), cutXline, "prof colz");
  Views->SetMarkerColor(6);
  Views->SetMarkerStyle(7);
  Views->Draw("eZ1:eXview", cTop && cutXline, "same");
  Views->Draw("eZ2:eXview", cTop && cutXline, "same");
  Views->Draw("eZ3:eXview", cBot && cutXline, "same");
  Views->Draw("eZ4:eXview", cBot && cutXline, "same");
  Views->SetMarkerColor(kBlack);
  pad->cd(3)->SetGrid();
  Views->Draw(Form("eNcl:eZframe:eYview>>hpy(%d,%f,%f,%d,%f,%f)", B.ny, B.ymin, B.ymax, 1200, glassmin, glassmax), cutYline, "prof colz");
  Views->SetMarkerColor(6);
  Views->SetMarkerStyle(7);
  Views->Draw("eZ1:eYview", cTop && cutYline, "same");
  Views->Draw("eZ2:eYview", cTop && cutYline, "same");
  Views->Draw("eZ3:eYview", cBot && cutYline, "same");
  Views->Draw("eZ4:eYview", cBot && cutYline, "same");
  // Views->SetMarkerColor(kBlack);

  SetLargeLabelsStyle();

  TVirtualPad *px = pad->cd(5);
  px->Divide(1, 2);
  px->cd(1);
  Views->Draw(Form("segments.eDz:eXview>>hdz51(%d,%f,%f,90,0,90)", B.nx, B.xmin, B.xmax), cTop && cutXline, "colz");
  Views->Draw(Form("eZ1-eZ2:eXview>>hdz510(%d,%f,%f,90,0,90)", B.nx, B.xmin, B.xmax), cTop && cutXline, "same");
  Views->Draw("eZ1-eZ2>>hdz511(90,0,90)", cTop && cutYline, "goff");
  TH2F *hdz51 = (TH2F *)gDirectory->Get("hdz51");
  Z.dzx1 = hdz51->GetMean(2);
  TH1F *hdz511 = (TH1F *)gDirectory->Get("hdz511");
  Z.szx1 = hdz511->GetMean();
  TText *t51 = new TText(0.15, 0.2, Form("surf_dZ - seg_dZ = %.1f - %.1f = %.1f um", Z.szx1, Z.dzx1, Z.szx1 - Z.dzx1));
  t51->SetNDC();
  t51->SetTextSize(0.07);
  t51->Draw();

  px->cd(2);
  Views->Draw(Form("segments.eDz:eXview>>hdz52(%d,%f,%f,90,0,90)", B.nx, B.xmin, B.xmax), cBot && cutXline, "colz");
  Views->Draw(Form("eZ3-eZ4:eXview>>hdz520(%d,%f,%f,90,0,90)", B.nx, B.xmin, B.xmax), cBot && cutXline, "same");
  Views->Draw("eZ1-eZ2>>hdz521(90,0,90)", cTop && cutYline, "goff");
  TH2F *hdz52 = (TH2F *)gDirectory->Get("hdz52");
  Z.dzx2 = hdz52->GetMean(2);
  TH1F *hdz521 = (TH1F *)gDirectory->Get("hdz521");
  Z.szx2 = hdz521->GetMean();
  TText *t52 = new TText(0.15, 0.2, Form("surf_dZ - seg_dZ = %.1f - %.1f = %.1f um", Z.szx2, Z.dzx2, Z.szx2 - Z.dzx2));
  t52->SetNDC();
  t52->SetTextSize(0.07);
  t52->Draw();

  TVirtualPad *py = pad->cd(6);
  py->Divide(1, 2);
  py->cd(1);
  Views->Draw(Form("segments.eDz:eYview>>hdz61(%d,%f,%f,90,0,90)", B.ny, B.ymin, B.ymax), cTop && cutYline, "colz");
  Views->Draw(Form("eZ1-eZ2:eYview>>hdz610(%d,%f,%f,90,0,90)", B.ny, B.ymin, B.ymax), cTop && cutYline, "same");
  Views->Draw("eZ1-eZ2>>hdz611(90,0,90)", cTop && cutYline, "goff");
  TH2F *hdz61 = (TH2F *)gDirectory->Get("hdz61");
  Z.dzy1 = hdz61->GetMean(2);
  TH1F *hdz611 = (TH1F *)gDirectory->Get("hdz611");
  Z.szy1 = hdz611->GetMean();
  TText *t61 = new TText(0.15, 0.2, Form("surf_dZ - seg_dZ = %.1f - %.1f = %.1f um", Z.szy1, Z.dzy1, Z.szy1 - Z.dzy1));
  t61->SetNDC();
  t61->SetTextSize(0.07);
  t61->Draw();

  py->cd(2);
  Views->Draw(Form("segments.eDz:eYview>>hdz62(%d,%f,%f,90,0,90)", B.ny, B.ymin, B.ymax), cBot && cutYline, "colz");
  Views->Draw(Form("eZ3-eZ4:eYview>>hdz620(%d,%f,%f,90,0,90)", B.ny, B.ymin, B.ymax), cBot && cutYline, "same");
  Views->Draw("eZ3-eZ4>>hdz621(90,0,90)", cBot && cutYline, "goff");
  TH2F *hdz62 = (TH2F *)gDirectory->Get("hdz62");
  Z.dzy2 = hdz62->GetMean(2);
  TH1F *hdz621 = (TH1F *)gDirectory->Get("hdz621");
  Z.szy2 = hdz621->GetMean();
  TText *t62 = new TText(0.15, 0.2, Form("surf_dZ - seg_dZ = %.1f - %.1f = %.1f um", Z.szy2, Z.dzy2, Z.szy2 - Z.dzy2));
  t62->SetNDC();
  t62->SetTextSize(0.07);
  t62->Draw();

  pad->cd(0);
  TDatime time;
  TText *t = new TText();
  t->SetTextSize(0.015);
  t->DrawText(0.25, 0.0001, Form("%s/%s    %s", gSystem->WorkingDirectory(), RES.filename.c_str(), time.AsString()));

  if (gROOT->IsBatch())
    cc->SaveAs(Form("%s.prof.png", nameo));
}

//----------------------------------------------------------------------------------------
int process_file(const char *input_file, const char *output_file)
{
  TFile *file = TFile::Open(input_file);
  if (!file || file->IsZombie())
  {
    Log(1, "process_file", "Error opening input file: %s", input_file);
    return 1;
  }
  RES.filename = input_file;

  Long_t id, flags, time;
  gSystem->GetPathInfo(input_file, &id, &(RES.file_size), &flags, &time);
  RES.file_modification_date.Set(time);
  RES.file_creation_date = file->GetCreationDate();
  std::cout << "File: " << input_file << "\n";

  TTree *tree = (TTree *)file->Get("Views");
  if (!tree)
  {
    Log(1, "process_file", "Error: No 'Views' tree found in file");
    return 2;
  }

  AREA.xmaxR = cenv.GetValue("quality.request.Xmax", 185000.);
  AREA.xminR = cenv.GetValue("quality.request.Xmin", 5000.);
  AREA.ymaxR = cenv.GetValue("quality.request.Ymax", 185000.);
  AREA.yminR = cenv.GetValue("quality.request.Ymin", 5000.);
  AREA.fraction = cenv.GetValue("quality.request.fraction", 0.9);
  AREA.areaR = (AREA.xmaxR - AREA.xminR) * (AREA.ymaxR - AREA.yminR);

  define_steps(tree);

  P.prof = cenv.GetValue("quality.prof", 1);

  make_histos(tree);
  if (P.prof)
    make_snd_profiles(output_file);
  make_report(output_file);
  make_canvas(output_file);
  if (DO.interactive)
    draw_frame_align();

  return 0;
}

//----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    print_help_message();
    return 0;
  }

  set_default_env();
  gEDBDEBUGLEVEL = cenv.GetValue("quality.EdbDebugLevel", 1);

  const char *input_file = "";
  const char *output_file = cenv.GetValue("quality.output", "quality_report.json");

  for (int i = 1; i < argc; i++)
  {
    char *key = argv[i];
    if (!strncmp(key, "-input=", 7))
    {
      input_file = key + 7;
      DO.process = true;
    }
    if (!strncmp(key, "-out=", 5))
    {
      output_file = key + 5;
    }
    if (!strncmp(key, "-I", 2))
    {
      DO.interactive = true;
    }
    if (!strncmp(key, "-v=", 3))
    {
      gEDBDEBUGLEVEL = atoi(key + 3);
    }
  }

  cenv.ReadFile("quality.rootrc", kEnvLocal);
  cenv.WriteFile("quality.save.rootrc");

  if (DO.interactive)
  {
    gROOT->SetBatch(false);
    if (!gApplication)
    {
      int dummy_argc = 1;
      char *dummy_argv[] = {"-l"};
      new TRint("APP", &dummy_argc, dummy_argv);
    }
  }

  int iproc = 0;
  if (DO.process)
  {
    iproc = process_file(input_file, output_file);
  }
  else
  {
    print_help_message();
    return 0;
  }

  if (DO.interactive)
  {
    gApplication->Run();
  }
  cenv.WriteFile("quality.save.rootrc");
  return iproc;
}
