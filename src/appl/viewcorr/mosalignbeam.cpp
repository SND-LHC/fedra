//-- Author :  Valeri Tioukov   20/04/2023

#include <string.h>
#include <iostream>
#include <TRint.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TChain.h>
#include <TList.h>
#include <TFile.h>
#include <TH3.h>
#include <TSystem.h>
#include "EdbLog.h"
#include "EdbRunAccess.h"
#include "EdbLinking.h"
#include "EdbScanProc.h"
#include "EdbPlateAlignment.h"
#include "EdbMosaic.h"
#include "EdbMosaicIO.h"
#include "EdbAttachPath.h"

using namespace std;
using namespace TMath;
void AlignToBeam(EdbID id, TEnv &env);
bool AlignFragmentToBeam0(EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2, float offMax, int flag = 0);
void TuneShrinkage(EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2, TEnv &env);
TH2D *ProfileAndCleanTH3(TH3 *h3, double min_entries = 20);
void AlignMicrotracksAngles(EdbPattern &p, TEnv &env);

void print_help_message()
{
  cout << "\nUsage: \n";
  cout << "\t  mosalignbeam  -id=ID   [-from=frag0 -nfrag=N  -merge -v=DEBUG] \n";
  cout << "\t  mosalignbeam  -set=ID  [-from=frag0 -nfrag=N  -merge -v=DEBUG] \n";

  cout << "\t\t  ID    - id of the data piece or data set formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout << "\t\t  frag0 - the first fragment (default: 0) \n";
  cout << "\t\t  N     - number of fragments to be processed (default: upto 1000000, stop at first empty) \n";
  cout << "\t\t  merge - merge all fragments into one cp file \n";

  cout << "\n If the data location directory if not explicitly defined\n";
  cout << " the current directory will be assumed to be the brick directory \n";
  cout << "\n If the parameters file (mosalignbeam.rootrc) is not presented - the default \n";
  cout << " parameters will be used. After the execution them are saved into mosalignbeam.save.rootrc file\n";
  cout << endl;
}

//----------------------------------------------------------------------------------------
void set_default_link(TEnv &cenv)
{
  // default parameters for the new linking

  cenv.SetValue("fedra.mosalignbeam.make_ab0", 0); // produce debug output (beam)
  cenv.SetValue("fedra.mosalignbeam.make_ab1", 0); // produce debig output (shrinkage)
  cenv.SetValue("fedra.mosalignbeam.read.ICUT", "-1     0. 120000.   0.   120000.    -1.   1.      -1.   1.       0.  500.");
  cenv.SetValue("fedra.mosalignbeam.DoAlignBeam", 1);
  cenv.SetValue("fedra.mosalignbeam.DoTuneShrinkage", 1);
  cenv.SetValue("fedra.mosalignbeam.DoFullLinking", 0);
  cenv.SetValue("fedra.mosalignbeam.DoAlignMicrotracksAngles", 0);
  cenv.SetValue("fedra.mosalignbeam.DoAlignMicrotracksAngles.bin", 100);
  cenv.SetValue("fedra.mosalignbeam.DoAlignMicrotracksAngles.minbin", 20);

  cenv.SetValue("fedra.link.AFID", 1); // 1 is usually fine for scanned data; for the db-read data use 0!
  cenv.SetValue("fedra.link.DoImageCorr", 0);
  cenv.SetValue("fedra.link.ImageCorrSide1", "1. 1. 0.");
  cenv.SetValue("fedra.link.ImageCorrSide2", "1. 1. 0.");
  cenv.SetValue("fedra.link.DoImageMatrixCorr", 0);
  cenv.SetValue("fedra.link.ImageMatrixCorrSide1", "");
  cenv.SetValue("fedra.link.ImageMatrixCorrSide2", "");

  cenv.SetValue("fedra.link.CheckUpDownOffset", 1); // check dXdY offsets between up and correspondent down views
  cenv.SetValue("fedra.link.BinOK", 6.);
  cenv.SetValue("fedra.link.NcorrMin", 100);
  cenv.SetValue("fedra.link.DoCorrectShrinkage", true);
  cenv.SetValue("fedra.link.read.UseDensityAsW", false);
  cenv.SetValue("fedra.link.RemoveDoublets", "1    2. .01   1"); // yes/no   dr  dt  checkview(0,1,2)
  cenv.SetValue("fedra.link.DumpDoubletsTree", true);
  cenv.SetValue("fedra.link.shr.NsigmaEQ", 7.5);
  cenv.SetValue("fedra.link.shr.Shr0", .85);
  cenv.SetValue("fedra.link.shr.DShr", .3);
  cenv.SetValue("fedra.link.shr.ThetaLimits", "0.05  1.");
  cenv.SetValue("fedra.link.DoCorrectAngles", true);
  cenv.SetValue("fedra.link.ang.Chi2max", 1.5);
  cenv.SetValue("fedra.link.DoFullLinking", false);
  cenv.SetValue("fedra.link.full.NsigmaEQ", 5.5);
  cenv.SetValue("fedra.link.full.DR", 20.);
  cenv.SetValue("fedra.link.full.DT", 0.1);
  cenv.SetValue("fedra.link.full.CHI2Pmax", 3.);
  cenv.SetValue("fedra.link.DoSaveCouples", false);
  cenv.SetValue("fedra.link.Sigma0", "1 1 0.007 0.007");
  cenv.SetValue("fedra.link.PulsRamp0", "6 9");
  cenv.SetValue("fedra.link.PulsRamp04", "6 9");
  cenv.SetValue("fedra.link.Degrad", 5);

  cenv.SetValue("fedra.link.LLfunction", "0.256336-0.16489*x+2.11098*x*x");
  cenv.SetValue("fedra.link.CPRankingAlg", 0);

  cenv.SetValue("emlink.reportfileformat", "pdf");
  cenv.SetValue("emlink.outdir", "..");
  cenv.SetValue("emlink.env", "link.rootrc");
  cenv.SetValue("emlink.EdbDebugLevel", 1);
}

bool do_make_ab0;
bool do_make_ab1;
int from_fragment = 0;
int n_fragments = 0;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    print_help_message();
    return 0;
  }

  TEnv cenv("mosalignbeamenv");
  gEDBDEBUGLEVEL = cenv.GetValue("mosalignbeam.EdbDebugLevel", 1);
  const char *env = cenv.GetValue("mosalignbeam.env", "mosalignbeam.rootrc");
  const char *outdir = cenv.GetValue("mosalignbeam.outdir", "..");

  bool do_single = false;
  bool do_set = false;
  bool do_merge = false;
  EdbID id;

  for (int i = 1; i < argc; i++)
  {
    char *key = argv[i];

    if (!strncmp(key, "-set=", 5))
    {
      if (strlen(key) > 5)
        if (id.Set(key + 5))
          do_set = true;
    }
    else if (!strncmp(key, "-id=", 4))
    {
      if (strlen(key) > 4)
        if (id.Set(key + 4))
          do_single = true;
    }
    else if (!strncmp(key, "-from=", 6))
    {
      if (strlen(key) > 6)
        from_fragment = atoi(key + 6);
    }
    else if (!strncmp(key, "-nfrag=", 7))
    {
      if (strlen(key) > 7)
        n_fragments = atoi(key + 7);
    }
    else if (!strncmp(key, "-merge", 6))
    {
      do_merge = true;
    }
    else if (!strncmp(key, "-v=", 3))
    {
      if (strlen(key) > 3)
        gEDBDEBUGLEVEL = atoi(key + 3);
    }
  }

  if (!(do_single || do_set))
  {
    print_help_message();
    return 0;
  }
  if (do_single && do_set)
  {
    print_help_message();
    return 0;
  }

  set_default_link(cenv);
  cenv.SetValue("mosalignbeam.env", env);
  cenv.ReadFile(cenv.GetValue("mosalignbeam.env", "mosalignbeam.rootrc"), kEnvLocal);
  cenv.SetValue("mosalignbeam.outdir", outdir);

  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("mosalignbeam.outdir", "..");
  cenv.WriteFile("mosalignbeam.save.rootrc");

  printf("\n----------------------------------------------------------------------------\n");
  printf("mosalignbeam  %s\n", id.AsString());
  printf("----------------------------------------------------------------------------\n\n");

  do_make_ab0 = cenv.GetValue("fedra.mosalignbeam.make_ab0", 0);
  do_make_ab1 = cenv.GetValue("fedra.mosalignbeam.make_ab1", 0);

  if (do_single)
  {
    AlignToBeam(id, cenv);
  }
  else if (do_set)
  {
    EdbScanSet *ss = sproc.ReadScanSet(id);
    if (ss)
    {
      int n = ss->eIDS.GetSize();
      for (int i = 0; i < n; i++)
      {
        EdbID *id_pl = ss->GetID(i);
        if (id_pl)
          AlignToBeam(*id_pl, cenv);
      }
    }
  }

  cenv.WriteFile("mosalignbeam.save.rootrc");
  return 1;
}

//-----------------------------------------------------------------------------------
void AlignToBeam(EdbID id, TEnv &cenv)
{
  EdbMosaicIO mio;
  TString mosfile;
  mosfile.Form("p%3.3d/%d.%d.%d.%d.mos.root",
               id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
  mio.Init(mosfile.Data());

  bool use_saved_alignment = true;
  int do_align_beam = cenv.GetValue("fedra.mosalignbeam.DoAlignBeam", 1);
  int do_align_microtracks_angles = cenv.GetValue("fedra.mosalignbeam.DoAlignMicrotracksAngles", 1);
  int do_correct_shrinkage = cenv.GetValue("fedra.mosalignbeam.DoTuneShrinkage", 1);
  int do_full_linking = cenv.GetValue("fedra.mosalignbeam.DoFullLinking", 0);
  mio.AddSegmentCut(1, cenv.GetValue("fedra.mosalignbeam.read.ICUT", "-1     0. 120000.   0.   120000.    -1.   1.      -1.   1.       0.  500."));

  EdbLayer *mapside1 = mio.GetCorrMap(id.ePlate, 1);
  EdbLayer *mapside2 = mio.GetCorrMap(id.ePlate, 2); // align side 2 to side 1
  mapside1->SetZ(97.5);                              // TODO take it from set.root
  mapside2->SetZ(-97.5);

  int first, last;
  int nc = mapside2->Map().Ncell();
  if (from_fragment == 0 && n_fragments == 0)
  {
    first = 0;
    last = first + nc;
  }
  else
  {
    first = from_fragment;
    last = first + n_fragments;
  }
  Log(1, "mosalignbeam::AlignToBeam", "with %d fragments [%d:%d] out of %d",
      last - first, first, last - 1, nc);

  for (int i = first; i < last; i++)
  {
    EdbLayer *l1 = mapside1->Map().GetLayer(i);
    EdbLayer *l2 = mapside2->Map().GetLayer(i);
    if (l1 && l2)
    {
      if (!use_saved_alignment)
      {
        l1->GetAffineXY()->Reset();
        l2->GetAffineXY()->Reset();
        l1->GetAffineTXTY()->Reset();
        l2->GetAffineTXTY()->Reset();
      }
      l1->SetZ(mapside1->Z()); // base thickness is considered fixed...
      l2->SetZ(mapside2->Z());
      EdbPattern *p1 = mio.GetFragment(id.ePlate, 1, i, use_saved_alignment); // get side 1
      EdbPattern *p2 = mio.GetFragment(id.ePlate, 2, i, use_saved_alignment); // get side 2
      if (p1 && p2)
      {
        p1->SetScanID(id);
        p2->SetScanID(id);
        p1->SetID(i);
        p2->SetID(i);
        p1->SetSide(1);
        p2->SetSide(2);
        p1->SetSegmentsFlag(0);
        p2->SetSegmentsFlag(0);
        Log(1, "mosalignbeam::AlignFragmentToBeam", "fragment %d: %d & %d", p1->ID(), p1->N(), p2->N());

        if (do_align_microtracks_angles)
        {
          AlignMicrotracksAngles(*p1, cenv); // locally align microtracks angles to beam
          AlignMicrotracksAngles(*p2, cenv); // locally align microtracks angles to beam
        }
        if (do_align_beam)
        {
          AlignFragmentToBeam0(*p2, *p1, *l2, *l1, 10);     // align 2 to 1 using parallel beam tracks
          AlignFragmentToBeam0(*p2, *p1, *l2, *l1, 5);      // align 2 to 1 using parallel beam tracks
          AlignFragmentToBeam0(*p2, *p1, *l2, *l1, 3, -10); // align 2 to 1 using parallel beam tracks, assign flag -10 to used segments
        }
        if (do_correct_shrinkage)
        {
          TuneShrinkage(*p2, *p1, *l2, *l1, cenv); // shrinkage correction using non-beam tracks
        }
        if (do_full_linking)
        {
          cenv.SetValue("fedra.link.DoFullLinking", true);
          cenv.SetValue("fedra.link.DoSaveCouples", true);
          cenv.SetValue("fedra.link.DoCorrectShrinkage", false);
          cenv.SetValue("fedra.link.DoCorrectAngles", false);
          cenv.SetValue("fedra.link.shr.ThetaLimits", "0.  1.");
          EdbLinking link;
          TString cpfile(Form("p%3.3d/%d.%d.%d.%d.%d.cp.root",
                              id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor, p1->ID()));
          link.InitOutputFile(cpfile.Data());
          Log(1, "mosalignbeam::AlignFragmentToBeam", "full linking -> %s", cpfile.Data());
          p1->SetSegmentsFlag(0);
          p2->SetSegmentsFlag(0);
          l1->ResetCorr();
          l2->ResetCorr();
          link.Link(*p2, *p1, *l2, *l1, cenv);
          TH3F *htx1 = (TH3F *)gROOT->Get(Form("htx_%d_%d", p1->Side(), p1->ID()));
          TH3F *hty1 = (TH3F *)gROOT->Get(Form("hty_%d_%d", p1->Side(), p1->ID()));
          TH3F *htx2 = (TH3F *)gROOT->Get(Form("htx_%d_%d", p2->Side(), p2->ID()));
          TH3F *hty2 = (TH3F *)gROOT->Get(Form("hty_%d_%d", p2->Side(), p2->ID()));
          if (htx1)
            htx1->Write();
          if (hty1)
            hty1->Write();
          if (htx2)
            htx2->Write();
          if (hty2)
            hty2->Write();
          link.CloseOutputFile();
        }
      }
      SafeDelete(p1);
      SafeDelete(p2);
    }
  }
  mio.Close();
  if (!gSystem->AccessPathName(mosfile.Data()))
  {
    if (!gSystem->AccessPathName(mosfile.Data(), kWritePermission))
    {
      mio.Init(mosfile.Data(), "UPDATE");
      mio.SaveCorrMap(id.ePlate, 1, *mapside1);
      mio.SaveCorrMap(id.ePlate, 2, *mapside2);
      mio.Close();
      Log(1, "mosalignbeam", "%s maps saved into %s", id.AsString(), mosfile.Data());
    }
    else
      Log(1, "mosalignbeam", "Error: file %s is not writable!", mosfile.Data());
  }
  else
    Log(1, "mosalignbeam", "Error: file %s is not accessible!", mosfile.Data());
}

//-----------------------------------------------------------------------
bool AlignFragmentToBeam0(EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2, float offMax, int flag)
{
  // Assume 0 angle beam here
  //
  bool success = false;
  int eMinPeak = 100;
  bool do_transform = true;

  EdbPlateAlignment av;
  av.eNoScale = 1;    // calculate shift and rotation
  av.eNoScaleRot = 0; // calculate shift only
  av.eOffsetMax = offMax;
  av.eDZ = 0.;
  av.eDPHI = 0.0;
  av.eDoFine = 1;
  if (do_make_ab0)
    av.eSaveCouples = 1;
  else
    av.eSaveCouples = 0;
  av.SetSigma(0.3, 0.025);
  av.eDoublets[0] = av.eDoublets[1] = 0.01;
  av.eDoublets[2] = av.eDoublets[3] = 0.0001;
  av.eDoCorrectAngle = false;
  av.eSaveCouples = 0;

  if (do_make_ab0)
    av.InitOutputFile(Form("p%.3d/%d_%d.ab0.root", p1.ScanID().ePlate, p1.ID(), p2.ID()));
  av.Align(p1, p2, 0, flag); //-190
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();

  float dtx1 = av.CalcMeanDiff2Const(2, 0, 0);
  float dty1 = av.CalcMeanDiff2Const(3, 0, 0);
  float dtx2 = av.CalcMeanDiff2Const(2, 1, 0);
  float dty2 = av.CalcMeanDiff2Const(3, 1, 0);

  EdbAffine2D aa1;
  aa1.ShiftX(-dtx1);
  aa1.ShiftY(-dty1);
  EdbAffine2D aa2;
  aa2.ShiftX(-dtx2);
  aa2.ShiftY(-dty2);

  // printf("\n angular offsets found: %f %f %f %f\n\n", dtx1,dty1,dtx2,dty2);

  if (av.eNcoins > eMinPeak)
  {
    if (do_transform)
    {
      p1.Transform(affXY);
      p1.TransformA(&aa1);
      p2.TransformA(&aa2);
    }
    l1.GetAffineXY()->Transform(affXY);
    l1.GetAffineTXTY()->Transform(aa1);
    l2.GetAffineTXTY()->Transform(aa2);
    success = true;
  }

  if (do_make_ab0)
    av.CloseOutputFile();
  return success;
}

//-----------------------------------------------------------------------
void TuneShrinkage(EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2, TEnv &env)
{
  EdbLinking link;
  if (do_make_ab1)
  {
    link.InitOutputFile(Form("p%.3d/%d_%d.ab1.root", p1.ScanID().ePlate, p1.ID(), p2.ID()));
  }
  else
    env.SetValue("fedra.link.DumpDoubletsTree", false);
  // calculate valid segments
  int n1_valid = 0;
  int n2_valid = 0;
  for (int i = 0; i < p1.N(); i++)
  {
    EdbSegP *s = p1.GetSegment(i);
    if (s->Flag() >= 0)
      n1_valid++;
  }
  for (int i = 0; i < p2.N(); i++)
  {
    EdbSegP *s = p2.GetSegment(i);
    if (s->Flag() >= 0)
      n2_valid++;
  }
  link.Link(p1, p2, l1, l2, env);
  l1.SetShrinkage(l1.Shr() * link.eL1.Shr());
  l2.SetShrinkage(l2.Shr() * link.eL2.Shr());

  Log(1, "mosalignbeam:TuneShrinkage", "with %d %d valid segments, shrinkages: %f %f", n1_valid, n2_valid, l1.Shr(), l2.Shr());
  if (do_make_ab1)
    link.CloseOutputFile();
}

//-----------------------------------------------------------------------
void AlignMicrotracksAngles(EdbPattern &p, TEnv &env)
{
  float bin = env.GetValue("fedra.mosalignbeam.DoAlignMicrotracksAngles.bin", 100.);
  int minbin = env.GetValue("fedra.mosalignbeam.DoAlignMicrotracksAngles.minbin", 20);
  float xmin = p.Xmin();
  float xmax = p.Xmax();
  float ymin = p.Ymin();
  float ymax = p.Ymax();
  int ibinX = (int)((xmax - xmin) / bin) + 2;
  int ibinY = (int)((ymax - ymin) / bin) + 2;
  xmin = xmin - 1. * bin;
  ymin = ymin - 1. * bin;
  xmax = xmin + ibinX * bin;
  ymax = ymin + ibinY * bin;
  if (ibinX < 1)
    ibinX = 1;
  if (ibinY < 1)
    ibinY = 1;
  Log(1, "mosalignbeam::AlignMicrotracksAngles", "fragment %d: %d microtracks, binning %d x %d", p.ID(), p.N(), ibinX, ibinY);
  if (p.N() < 10)
  {
    Log(1, "mosalignbeam::AlignMicrotracksAngles", "fragment %d: too few microtracks, skip", p.ID());
    return;
  }
  TH3F *htx = new TH3F(Form("htx_%d_%d", p.Side(), p.ID()), "htx", ibinX, xmin, xmax, ibinY, ymin, ymax, 100, -0.1, 0.1);
  TH3F *hty = new TH3F(Form("hty_%d_%d", p.Side(), p.ID()), "hty", ibinX, xmin, xmax, ibinY, ymin, ymax, 100, -0.1, 0.1);
  htx->SetDirectory(gROOT);
  hty->SetDirectory(gROOT);
  for (int i = 0; i < p.N(); i++)
  {
    EdbSegP *s = p.GetSegment(i);
    if (s->TX() < -0.1 || s->TX() > 0.1 || s->TY() < -0.1 || s->TY() > 0.1)
      continue;
    htx->Fill(s->X(), s->Y(), s->TX());
    hty->Fill(s->X(), s->Y(), s->TY());
  }

  TH2D *h2_zmean = ProfileAndCleanTH3(htx, minbin);
  TH2D *h2_zmean_ty = ProfileAndCleanTH3(hty, minbin);

  for (int i = 0; i < p.N(); i++)
  {
    EdbSegP *s = p.GetSegment(i);
    float txcorr = h2_zmean->Interpolate(s->X(), s->Y());
    float tycorr = h2_zmean_ty->Interpolate(s->X(), s->Y());
    s->SetTX(s->TX() - txcorr);
    s->SetTY(s->TY() - tycorr);
  }
}

//-----------------------------------------------------------------------
// Converts TH3(x,y,z) -> cleaned TH2D(x,y) where bin content = mean z
TH2D *ProfileAndCleanTH3(TH3 *h3, double min_entries)
{
  int nBinsX = h3->GetNbinsX();
  int nBinsY = h3->GetNbinsY();
  int nBinsZ = h3->GetNbinsZ();

  // 1. Create a 2D histogram matching the X and Y axes of the TH3
  TH2D *h2_zmean = new TH2D(
      Form("%s_zmean", h3->GetName()),
      "Mean Z Map;x;y",
      nBinsX, h3->GetXaxis()->GetXmin(), h3->GetXaxis()->GetXmax(),
      nBinsY, h3->GetYaxis()->GetXmin(), h3->GetYaxis()->GetXmax());
  h2_zmean->SetDirectory(0);

  // Track total entries per (x,y) cell to identify sparse bins
  std::vector<std::vector<double>> cell_entries(nBinsX + 1, std::vector<double>(nBinsY + 1, 0.0));

  // 2. Step 1: Calculate mean z for each (x, y) cell
  for (int ix = 1; ix <= nBinsX; ++ix)
  {
    for (int iy = 1; iy <= nBinsY; ++iy)
    {
      double sum_wz = 0.0;
      double sum_w = 0.0;

      // find maximum z bin content for this (x, y) cell
      double max_z_content = 0.0;
      for (int iz = 1; iz <= nBinsZ; ++iz)
      {
        double content = h3->GetBinContent(ix, iy, iz);
        if (content > max_z_content)
          max_z_content = content;
      }

      // Integrate over Z bins gt max_z_content/2 for this (x, y) cell
      for (int iz = 1; iz <= nBinsZ; ++iz)
      {
        double content = h3->GetBinContent(ix, iy, iz);
        double z_center = h3->GetZaxis()->GetBinCenter(iz);
        if (content < max_z_content / 2)
          continue; // skip low content bins
        sum_wz += content * z_center;
        sum_w += content;
      }

      cell_entries[ix][iy] = sum_w;

      // If the cell meets the minimum entry threshold, store mean z
      if (sum_w >= min_entries)
      {
        h2_zmean->SetBinContent(ix, iy, sum_wz / sum_w);
      }
    }
  }

  // 3. Step 2: Fill sparse/empty (x, y) cells using nearest non-empty 2D neighbors
  int max_radius = std::max(nBinsX, nBinsY);

  for (int ix = 1; ix <= nBinsX; ++ix)
  {
    for (int iy = 1; iy <= nBinsY; ++iy)
    {

      // Check if cell is below threshold
      if (cell_entries[ix][iy] < min_entries)
      {
        double sum_weighted_z = 0.0;
        double sum_weight = 0.0;
        bool found_neighbors = false;

        // Search expanding 2D rings around (ix, iy)
        for (int r = 1; r <= max_radius; ++r)
        {
          for (int dx = -r; dx <= r; ++dx)
          {
            for (int dy = -r; dy <= r; ++dy)
            {
              if (std::abs(dx) != r && std::abs(dy) != r)
                continue;

              int nx = ix + dx;
              int ny = iy + dy;
              if (nx < 1 || nx > nBinsX || ny < 1 || ny > nBinsY)
                continue;

              if (cell_entries[nx][ny] >= min_entries)
              {
                double dist = std::sqrt(dx * dx + dy * dy);
                double weight = 1.0 / dist; // Inverse distance weight

                sum_weighted_z += h2_zmean->GetBinContent(nx, ny) * weight;
                sum_weight += weight;
                found_neighbors = true;
              }
            }
          }
          if (found_neighbors)
            break; // Stop at closest distance ring
        }

        if (found_neighbors)
        {
          h2_zmean->SetBinContent(ix, iy, sum_weighted_z / sum_weight);
        }
      }
    }
  }

  return h2_zmean;
}