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
#include <TPaveText.h>
#include "ERTools.h"
#include "EdbLog.h"

using namespace std;

struct
{
  bool process=false;        // process raw.root file
  bool interactive=false;    // run interactive session
} DO;
struct
{
  string  filename;
  TDatime file_creation_date;
  TDatime file_modification_date;
  Long_t  file_size;
  float   scan_time;           // total scanning time in seconds
  float   requested_area;
  float   scanned_area_top;
  float   scanned_area_bot;
  float   empty_frac_top;
  float   empty_frac_bot;
  float   density_top;         // density as total/area excluding lost (zero) bins [/mm2]
  float   density_bot;
} RES;                         // Result to be exported as report

struct 
{
  int nx;
  float xmin, xmax, xbin;
  int ny;
  float ymin, ymax, ybin;
  double rxmin, rxmax;
  double rymin, rymax;
}B;                         // view step binning

struct
{
  TProfile2D *nseg_top;
  TProfile2D *nseg_bot;
  TProfile2D *thick_top;
  TProfile2D *thick_bot;
  TH2D       *thick_base;
  TProfile2D *glass;
}H;                          //oputput histos

struct
{
  float bin=20;
} FFT;

using namespace ERTools;

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout<< "\nUsage: \n\t  emthickness -input=file.raw.root [ -out=report -I -v=DEBUG] \n";
  cout<< "\t  Generate quality report for raw data\n";
  cout<< "\t\t -input      -  input raw root file with tracks\n";
  cout<< "\t\t -out        -  output report filename (default: raw_quality_report[.json],[.png] )\n";
  cout<< "\t\t -I          -  open interactive root session\n";
  cout<< "\t\t -v          -  verbosity level\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default(TEnv &cenv)
{
  cenv.SetValue("quality.input"          , "");
  cenv.SetValue("quality.output"         , "raw_quality_report");
  cenv.SetValue("quality.EdbDebugLevel"  , 1);
  cenv.SetValue("quality.requestXmin"    , 5000);
  cenv.SetValue("quality.requestYmin"    , 5000);
  cenv.SetValue("quality.requestXmax"    , 185000);
  cenv.SetValue("quality.requestYmax"    , 185000);
}

//----------------------------------------------------------------------------------------
void generate_json_report(TH2* h, std::ofstream &report, bool first_plot)
{
  if (!first_plot) report << ",\n";
  
  // Calculate statistics
  double sum = 0, sum2 = 0;
  int total_bins = h->GetNbinsX() * h->GetNbinsY();
  int empty_bins = 0;
  int count = 0;
  
  for (int binx = 1; binx <= h->GetNbinsX(); ++binx) {
    for (int biny = 1; biny <= h->GetNbinsY(); ++biny) {
      if (h->GetBinContent(binx, biny) == 0) {
	empty_bins++;
	continue;
      }
      double val = h->GetBinContent(binx, biny);
      sum += val;
      sum2 += val * val;
      count++;
    }
  }
  int rnbinsx = int((B.rxmax-B.rxmin)/B.xbin) + 1;
  int rnbinsy = int((B.rymax-B.rymin)/B.ybin) + 1;
  RES.requested_area = B.xbin*B.ybin*rnbinsx*rnbinsy;
  double scanned_area = count*B.xbin*B.ybin;
  double mean = (count > 0) ? sum / count : 0;
  double rms = (count > 0) ? TMath::Sqrt(sum2/count - mean*mean) : 0;
  //double empty_frac = (double)empty_bins / total_bins;
  double empty_frac = (RES.requested_area - scanned_area)/(RES.requested_area);
  double density = mean/B.xbin/B.ybin*1000*1000;  // density in mm2
  
  if( !strcmp( h->GetName(),"nseg_top") ) 
  {
    RES.scanned_area_top = scanned_area;
    RES.empty_frac_top = empty_frac;
    RES.density_top = density;
  }
  if( !strcmp( h->GetName(),"nseg_bot") )
  {
    RES.scanned_area_bot =scanned_area;
    RES.empty_frac_bot = empty_frac;
    RES.density_bot = density;
  }
  
  int nsig=3;
  h->GetZaxis()->SetRangeUser( TMath::Max(0.,     TMath::Max( mean-3*rms, h->GetMinimum() )), 
			       TMath::Min(3*mean, TMath::Min( mean+5*rms, h->GetMaximum() )) );
  h->SetTitle( Form("%s: mean = %.2f  RMS = %.2f",h->GetTitle(), mean, rms) );
  
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
TH1D *get_step_fft( TH1D *h )
{
  // Perform the Fourier Transform
  TH1 *h_fft = h->FFT(nullptr, Form("MAG_%s",h->GetName()) ); // Compute magnitude of FFT
  h_fft->SetTitle("Fourier Transform;Frequency;Magnitude");

  auto peaks = FindPeaksWithIntegral(h_fft, 5, 0.1, 2, h_fft->GetNbinsX()/2); // 5-bin window

  // Convert bin to frequency
  double frequency = peaks[0].peak_position;
  float xmin=h->GetXaxis()->GetXmin();
  float xmax=h->GetXaxis()->GetXmax();
  double step_size = (xmax-xmin) / (frequency-1); // Step size 
  if(gEDBDEBUGLEVEL>1) {
    for(int i=0; i<3; i++) printf("%d %f %f %f\n",
      i,peaks[i].peak_position, peaks[i].peak_height, peaks[i].window_integral );
    std::cout << "Estimated step size: " << step_size << std::endl;
  }
  
  int n = TMath::Ceil((xmax-xmin)/step_size);
  xmin -= step_size/2;   xmax += step_size/2; n+=1;
  
  TH1D *hstep = new TH1D(Form("hstep_%s",h->GetName()),"hstep",n,xmin,xmax);
  
  if(0) {
    // Visualize the results
    TCanvas *c = new TCanvas(h->GetName(), "Step Size Estimation with FFT", 1200, 600);
    c->Divide(2, 1);
    c->cd(1); h->Draw(); // Original binned data
    c->cd(2); h_fft->Draw(); // Fourier Transform
  }
  return hstep;
}


//----------------------------------------------------------------------------------------
void define_steps(TTree *tree)
{
  TH1D *hx = get_h_var(tree,"headers.eXview","hx",FFT.bin,"1");
  TH1D *hy = get_h_var(tree,"headers.eYview","hy",FFT.bin,"1");
  //hx->Draw();
  
  TH1D *hstepx = get_step_fft(hx);
  TH1D *hstepy = get_step_fft(hy);
  B.nx   = hstepx->GetXaxis()->GetNbins();
  B.xmin = hstepx->GetXaxis()->GetXmin();
  B.xmax = hstepx->GetXaxis()->GetXmax();
  B.xbin = (B.xmax-B.xmin)/B.nx;
  B.ny   = hstepy->GetXaxis()->GetNbins();
  B.ymin = hstepy->GetXaxis()->GetXmin();
  B.ymax = hstepy->GetXaxis()->GetXmax();
  B.ybin = (B.ymax-B.ymin)/B.ny;

}

//----------------------------------------------------------------------------------------
void  make_histos(TTree *tree)
{
  TCut sideTop="eNframesTop>0";
  TCut sideBot="eNframesBot>0";
  TProfile2D *hview = new TProfile2D("hview", "hview", B.nx,B.xmin,B.xmax, B.ny, B.ymin,B.ymax);
  
  H.nseg_top   = (TProfile2D*)(hview->Clone("nseg_top"));   H.nseg_top->SetTitle("nseg_top");
  H.nseg_bot   = (TProfile2D*)(hview->Clone("nseg_bot"));   H.nseg_bot->SetTitle("nseg_bot");
  H.thick_top  = (TProfile2D*)(hview->Clone("thick_top"));  H.thick_top->SetTitle("thick_top");
  H.thick_bot  = (TProfile2D*)(hview->Clone("thick_bot"));  H.thick_bot->SetTitle("thick_bot");
  H.glass      = (TProfile2D*)(hview->Clone("glass"));      H.glass->SetTitle("glass");
  TProfile2D *hz2 = (TProfile2D*)(hview->Clone("hz2"));
  TProfile2D *hz3 = (TProfile2D*)(hview->Clone("hz3"));
  H.thick_base    = new TH2D("thick_base", "thick_base", B.nx,B.xmin,B.xmax, B.ny, B.ymin,B.ymax);
  
  tree->Draw("eNsegments:eYview:eXview>>nseg_top",sideTop,"goff prof colz");
  tree->Draw("eNsegments:eYview:eXview>>nseg_bot",sideBot,"goff prof colz");
  tree->Draw("eZ1-eZ2:eYview:eXview>>thick_top",sideTop,"goff prof colz");
  tree->Draw("eZ3-eZ4:eYview:eXview>>thick_bot",sideBot,"goff prof colz");
  tree->Draw("eZ2:eYview:eXview>>hz2",sideTop,"goff prof colz");
  tree->Draw("eZ3:eYview:eXview>>hz3",sideBot,"goff prof colz");
  tree->Draw("eZ4:eYview:eXview>>glass",sideBot,"goff prof colz");
  
  DiffProfile2D( hz2, hz3, H.thick_base );
  
  float dxcm = (B.xmax-B.xmin)/10000;
  float dycm = (B.ymax-B.ymin)/10000;
  //RES.scanned_area = dxcm*dycm;        // cm2
 
  tree->Draw("eEvent>>htime(10000)","","goff");
  TH1 *htime = (TH1*)(gDirectory->Get("htime"));
  if(htime)  RES.scan_time = htime->Integral()*htime->GetMean()/1000.;  // scan time in seconds
  
  SafeDelete(hview);
  SafeDelete(hz2);
  SafeDelete(hz3);
}

void make_canvas(const char *nameo="ccc")
{
  gStyle->SetNumberContours(256);
  gStyle->SetPalette(107);
  gStyle->SetOptStat("");
  gStyle->SetPadRightMargin(0.12);
  
  TCanvas *cc = new TCanvas("thickness",Form("thickness at %s",nameo),1920,1080);
  
  TPad *header = new TPad("header", "Header", 0.0, 0.91, 1, 1);  // xlow, ylow, xup, yup
  header->Draw();
  header->cd();  
  TPaveText *tp = new TPaveText(0.01,0.01,0.99,0.99, "NDC");
  tp->AddText( Form("%s  of  %s   %.1f GB",
		    RES.filename.c_str(), 
		    RES.file_creation_date.AsString(), 
		    RES.file_size/1024./1024./1024. ) );
  tp->AddText( Form("Xrange:  %.1f   %.1f     Yrange:  %.1f   %.1f     Requested area  (%.1f, %.1f) x (%.1f, %.1f)     Step:  %.1f   %.1f   ScanTime: %.2f h",
		    B.xmin, B.xmax, B.ymin, B.ymax, B.rxmin, B.rxmax, B.rymin, B.rymax, B.xbin, B.ybin, RES.scan_time/60./60.) );
  tp->AddText( Form(
    "Scanned area top:  %.1f cm2  EmptyTop: %.2f %%  Scanned area bot: %.1f cm2  EmptyBot: %.2f %%   Density/mm2: Top = %.1f  Bot = %.1f",
    RES.scanned_area_top/1E8, RES.empty_frac_top*100, RES.scanned_area_bot/1E8, RES.empty_frac_bot*100, RES.density_top, RES.density_bot) );
  if ( RES.scanned_area_top < 0.96*RES.requested_area || RES.scanned_area_bot < 0.96*RES.requested_area) {
    TText *tt = tp->AddText("Warning! Scanned area is lower than 95 % of the requested area");
    tt->SetTextColor(kRed);
  }
  tp->Draw();
  TProfile2D* dummy= new TProfile2D("dummy", "", 1, 0, 200000, 1, 0, 200000);
  dummy->SetStats(0);
  
  cc->cd(0);
  TPad    *pad = new TPad("pad","",0.,0.,1.,0.9);
  pad->Draw();
  pad->cd();
  pad->Divide(3,2, 0.01, 0.01);
  
  //pad->cd(1);  dummy->SetTitle(H.nseg_top->GetTitle());dummy->SetMaximum(H.nseg_top->GetMaximum());dummy->SetMinimum(H.nseg_top->GetMinimum());dummy->DrawCopy("AXIS COLZ");gPad->Update();H.nseg_top->Draw("same prof colz");gPad->Update();
  pad->cd(1);  dummy->SetTitle(H.nseg_top->GetTitle());dummy->SetMaximum(H.nseg_top->GetMaximum()*1.2);dummy->SetMinimum(0.);dummy->DrawCopy("AXIS COLZ");gPad->Update();H.nseg_top->Draw("same prof colz");gPad->Update();
  //pad->cd(4);  dummy->SetTitle(H.nseg_bot->GetTitle());dummy->SetMaximum(H.nseg_bot->GetMaximum());dummy->SetMinimum(H.nseg_bot->GetMinimum());dummy->DrawCopy("AXIS COLZ");gPad->Update();H.nseg_bot->Draw("same prof colz");gPad->Update();
  pad->cd(4);  dummy->SetTitle(H.nseg_bot->GetTitle());dummy->SetMaximum(H.nseg_bot->GetMaximum()*1.2);dummy->SetMinimum(0.);dummy->DrawCopy("AXIS COLZ");gPad->Update();H.nseg_bot->Draw("same prof colz");gPad->Update();
  pad->cd(2);  dummy->SetTitle(H.thick_top->GetTitle());dummy->SetMaximum(H.thick_top->GetMaximum());dummy->SetMinimum(H.thick_top->GetMinimum());dummy->DrawCopy("AXIS COLZ");gPad->Update();H.thick_top->Draw("same prof colz");gPad->Update();
  pad->cd(5);  dummy->SetTitle(H.thick_bot->GetTitle());dummy->SetMaximum(H.thick_bot->GetMaximum());dummy->SetMinimum(H.thick_bot->GetMinimum());dummy->DrawCopy("AXIS COLZ");gPad->Update();H.thick_bot->Draw("same prof colz");gPad->Update();
  pad->cd(3);  dummy->SetTitle(H.thick_base->GetTitle());dummy->SetMaximum(H.thick_base->GetMaximum());dummy->SetMinimum(H.thick_base->GetMinimum());dummy->DrawCopy("AXIS COLZ");gPad->Update();H.thick_base->Draw("same prof colz");gPad->Update();
  pad->cd(6);  dummy->SetTitle(H.glass->GetTitle());dummy->SetMaximum(H.glass->GetMaximum());dummy->SetMinimum(H.glass->GetMinimum());dummy->DrawCopy("AXIS COLZ");gPad->Update();H.glass->Draw("same prof colz");gPad->Update();
  pad->cd(0);
  TDatime time;
  TText *t = new TText();
  t->SetTextSize(0.015);
  t->DrawText(0.25,0.0001, Form("%s/%s    %s",gSystem->WorkingDirectory(),
	      RES.filename.c_str(),time.AsString()) );
  if(gROOT->IsBatch()) cc->SaveAs(Form("%s.png",nameo));
}

int make_report( const char *output_file )
{
  // Open JSON report
  std::ofstream report( Form("%s.json",output_file) );
  if (!report.is_open()) {
    Log(1, "quality_report", "Error opening output file: %s", output_file);
    return 3;
  }
  
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
      while (TObject* obj = next()) {
	if (obj->InheritsFrom("TH2") || obj->InheritsFrom("TProfile2D")) {
	  generate_json_report((TH2*)obj, report, first_plot);
	  first_plot = false;
      }
    }
  report << "\n  ],\n";
  report << "  \"result\": {\n";
  report << "    \"scan_time\": \"" << RES.scan_time << "\",\n";
  report << "    \"scanned_area_top\": \"" << RES.scanned_area_top << "\",\n";
  report << "    \"scanned_area_bot\": \"" << RES.scanned_area_bot << "\",\n";
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
  TCut cut("cut","nsg>2000");
  TTree *tfa = (TTree*)gFile->Get("FrameAlign");

  TH1D *z1s1 = get_h_var( tfa, "z1", "z1s1", 1., "side==0" );
  TH1D *z1s2 = get_h_var( tfa, "z1", "z1s2", 1., "side==1" );

  TH2D *dxs1 = new TH2D("dxs1","dx:z  side1", 1000,z1s1->GetXaxis()->GetXmin(),z1s1->GetXaxis()->GetXmax(), 500,-2.,2.);
  TH2D *dxs2 = new TH2D("dxs2","dx:z  side2", 1000,z1s2->GetXaxis()->GetXmin(),z1s2->GetXaxis()->GetXmax(), 500,-2.,2.);
  TH2D *dys1 = new TH2D("dys1","dy:z  side1", 1000,z1s1->GetXaxis()->GetXmin(),z1s1->GetXaxis()->GetXmax(), 500,-2.,2.);
  TH2D *dys2 = new TH2D("dys2","dy:z  side2", 1000,z1s2->GetXaxis()->GetXmin(),z1s2->GetXaxis()->GetXmax(), 500,-2.,2.);
   
  TCanvas *c = new TCanvas("fral","frames align",1000,800);
  c->Divide(2,2);
  c->cd(1)->SetGrid();
  tfa->Draw("dx:z1>>dxs1", cut && "side==0","colz");
  c->cd(2)->SetGrid();
  tfa->Draw("dy:z1>>dys1", cut && "side==0","colz");
  c->cd(3)->SetGrid();
  tfa->Draw("dx:z1>>dxs2", cut &&"side==1","colz");
  c->cd(4)->SetGrid();
  tfa->Draw("dy:z1>>dys2", cut && "side==1","colz");
}
  
//----------------------------------------------------------------------------------------
int process_file(const char* input_file, const char* output_file)
{
  TFile *file = TFile::Open(input_file);
  if (!file || file->IsZombie()) {
    Log(1, "process_file", "Error opening input file: %s", input_file);
    return 1;
  }
  RES.filename=input_file;
  
  Long_t  id, flags, time;
  gSystem->GetPathInfo(input_file, &id, &(RES.file_size), &flags, &time);
  RES.file_modification_date.Set(time);
  RES.file_creation_date = file->GetCreationDate();
  std::cout << "File: " << input_file << "\n";
  
  TTree *tree = (TTree*)file->Get("Views");
  if (!tree) {
    Log(1, "process_file", "Error: No 'Views' tree found in file");
    return 2;
  }
  
  define_steps(tree);
  make_histos(tree);
  make_report(output_file);
  make_canvas(output_file);
  if(DO.interactive) draw_frame_align();
  
  return 0;
}

//----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argc < 2) { print_help_message(); return 0; }
  
  TEnv cenv("qualityenv");
  set_default(cenv);
  gEDBDEBUGLEVEL = cenv.GetValue("quality.EdbDebugLevel", 1);
  
  const char *input_file = "";
  const char *output_file = cenv.GetValue("quality.output", "quality_report.json");
  
  for(int i=1; i<argc; i++) {
    char *key = argv[i];
    if(!strncmp(key,"-input=",7)) {
      input_file = key+7;
      DO.process = true;
    }
    if(!strncmp(key,"-out=",5)) {
      output_file = key+5;
    }
    if(!strncmp(key,"-I",2)) {
      DO.interactive = true;
    }
    if(!strncmp(key,"-v=",3)) {
      gEDBDEBUGLEVEL = atoi(key+3);
    }
  }
  
  cenv.ReadFile("quality.rootrc", kEnvLocal);
  B.rxmin = cenv.GetValue("quality.requestXmin", 5000);
  B.rymin = cenv.GetValue("quality.requestYmin", 5000);
  B.rxmax = cenv.GetValue("quality.requestXmax", 185000);
  B.rymax = cenv.GetValue("quality.requestYmax", 185000);
  if(DO.interactive) {
    gROOT->SetBatch(false);
    if(!gApplication) {
        int dummy_argc = 1;
        char* dummy_argv[] = {"-l"};
        new TRint("APP", &dummy_argc, dummy_argv);
    }
 }
  
  int iproc=0;
  if(DO.process) {
    iproc = process_file(input_file, output_file);
  } else {
    print_help_message();
    return 0;
  }
  
  if(DO.interactive) {
    gApplication->Run();
  }

  return iproc;
}
