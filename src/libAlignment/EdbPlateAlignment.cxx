//-- Author :  Valeri Tioukov   19-06-2009
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// EdbPlateAlignment                                                      //
//                                                                        //
// Plate-to-plate alignment, normally with the basetracks                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "EdbLog.h"
#include "EdbPlateAlignment.h"
#include "EdbLayer.h"
#include "EdbCouplesTree.h"
#include "EdbTrackFitter.h"
#include "EdbScanCond.h"

ClassImp(EdbPlateAlignment)

using namespace TMath;

//---------------------------------------------------------------------
EdbPlateAlignment::EdbPlateAlignment()
{
  SetSigma( 10, 0.005 );
  eOffsetMax = 500.;
  eDZ        = 120;
  eDPHI      = 0.009;
  eDoublets[0]=eDoublets[1]=5;
  eDoublets[2]=eDoublets[3]=0.003;

  eStatus   =false;
  eDoTestAl =false; eTestAlOK =true;
  eDoCoarse =true;  eCoarseOK =false;
  eDoFine   =true;  eFineOK   =false;
  eSaveCouples=false;
  eRankCouples=false;
  eDoCorrectBeforeSaving=false;
  eDoCorrectAngle=true;               // keep true for compatibility 
  eNoScale=0;
  
  eNcoins    = 0;
  eCoarseMin = 5;
  eFineMin   = 5;
  eNiter     = 3;
}

//---------------------------------------------------------------------
EdbPlateAlignment::~EdbPlateAlignment()
{
  CloseOutputFile();
}

//---------------------------------------------------------------------
void EdbPlateAlignment::Align(EdbPattern &p1, EdbPattern &p2, float dz, int flag)
{
  // if flag is different from 0 assign it to the selected segments
  //
  eStatus = false;
  eHxy.Delete();
  eH_zphi_coarse.Delete();
  eH_xy_coarse.Delete();
  eH_xy_final.Delete();

  Log(2,"EdbPlateAlignment::Align","pattern 1 with %d in limits x:(%12.2f %12.2f) y:(%12.2f %12.2f)",
        p1.N(), p1.Xmin(),p1.Xmax(),p1.Ymin(),p1.Ymax() );
  Log(2,"EdbPlateAlignment::Align","pattern 2 with %d in limits x:(%12.2f %12.2f) y:(%12.2f %12.2f)",
        p2.N(), p2.Xmin(),p2.Xmax(),p2.Ymin(),p2.Ymax() );
  
  if(p1.N()<1||p2.N()<1) return;
  DoubletsFilterOut(p1,p2);

  SetParTestAl(dz);
  if(eDoTestAl)     {
    TestAl(p1,p2);
    if(!eTestAlOK)    goto END;
  }

  CorrToCoG(1,p2);

  if(eDoCoarse) {
    SetParCoarseAl( eCorr[0].V(2) , eOffsetMax, 3*eSigma[1], eDZ, eDPHI );
    CoarseAl(p1,p2);
    if(!eCoarseOK)    goto END;
  }

 if(eDoFine) {
    SetParFineAl();
    for(int i=0; i<eNiter; i++)   FineAl(p1,p2);

    eUseAffCorr=true;
    Corr2Aff( eCorrL[0] );

    for(int i=0; i<eNiter; i++)   FineAlAff(p1,p2,eCorrL[0] );
    eH_xy_final.Init(eHxy);
    if(!eFineOK)    goto END;
  }
  
 END:
 
  if(eRankCouples) {
//    eCorr[0].SetV(2,0);  eCorr[1].SetV(2, 0); ??
    RankCouples( eS[0], eS[1] );
  }
  ProduceReport();
  if(flag) AssignFlagToS(flag);           // assign particular flag to the selected segments
  if(eSaveCouples) SaveCouplesTree();
}

//---------------------------------------------------------------------
void EdbPlateAlignment::SaveCouplesTree()
{
  EdbCouplesTree ect;
  ect.InitCouplesTree("couples",0,"NEW");
  int nseg = CheckEqualArr(eS[0],eS[1]);
  for(int i=0; i<nseg; i++)
    if(eDoCorrectBeforeSaving) eCorrL[0].CorrectSeg(*(EdbSegP*)eS[0].UncheckedAt(i));
  for(int i=0; i<nseg; i++) {
    if( eRankCouples ) {
      EdbSegCouple *sc = (EdbSegCouple *)eSegCouples.At(i);
      ect.Fill( sc->eS1, sc->eS2, sc->eS, sc );
    }
    else ect.Fill( (EdbSegP*)eS[0].UncheckedAt(i), (EdbSegP*)eS[1].UncheckedAt(i) );
  }
  ect.eTree->AutoSave();
}

//---------------------------------------------------------------------
void EdbPlateAlignment::ProduceReport()
{
  TH2F *zphiCoarse=0, *xyCoarse=0;
  TH2F *xyFine=0;
  gStyle->SetPalette(1);
  if( eDoTestAl&&eTestAlOK ) {
    Log(1,"Alignment Report","TestAl OK");
  }
  if( eDoCoarse&&eCoarseOK ) {
    zphiCoarse = eH_zphi_coarse.DrawH2("zphi_coarse", "Z vs Phi after coarse align");
    xyCoarse   = eH_xy_coarse.DrawH2("xy_coarse", "dY vs dX after coarse align");
    float dx,dy;
    int npk = (int)eH_xy_coarse.FindPeak9( dx, dy );
    eH_xy_coarse.ProbPeak();

    if     (npk<eCoarseMin)                       eStatus = false;   // too few coinsidences
    else if(npk < 3*Sqrt(eH_xy_coarse.eMean[0]) ) eStatus = false;   // too small signal/noise
    else                                          eStatus = true;    // good alignment

    if(eStatus) 
      Log(2,"Alignment Report","Coarse %s: %d coins with peak of %d/%6.3f","OK",
	  npk, (int)eH_xy_coarse.ePeak[0], eH_xy_coarse.eMean[0]);    
    else {
      Log(1,"Alignment Report","Coarse %s: %d coins with peak of %d/%6.3f","FAILED!",
	  npk, (int)eH_xy_coarse.ePeak[0], eH_xy_coarse.eMean[0]);
      return;
    }
  }


  if( eDoFine&&eFineOK ) {
    xyFine = eH_xy_final.DrawH2("xy_final","dY vs dX after final alignment");
    eNcoins = eH_xy_final.Integral();
    if( eNcoins < eFineMin )               eStatus = false;   // too few coinsidences
    else                                   eStatus = true;    // good alignment

     if(eStatus) 
       Log(2,"Alignment Report","Fine   %s: %d coins inside (%6.1f %6.1f) with angular acc: (%6.3f %6.3f)", "OK",
	(int)(eH_xy_final.Integral()), 
	eH_xy_final.Xmax()-eH_xy_final.Xmin(),eH_xy_final.Ymax()-eH_xy_final.Ymin(),
	eDVsame[2], eDVsame[3] );
     else {
       Log(2,"Alignment Report","Fine   %s: %d coins inside (%6.1f %6.1f) with angular acc: (%6.3f %6.3f)", "FAILED!",
	(int)(eH_xy_final.Integral()), 
	eH_xy_final.Xmax()-eH_xy_final.Xmin(),eH_xy_final.Ymax()-eH_xy_final.Ymin(),
	eDVsame[2], eDVsame[3] );
	return;
     }
  }

    EdbAffine2D *aXY = eCorrL[0].GetAffineXY();
    float xcenter1 = (ePC[0].Xmin()+ePC[0].Xmax())/2.;
    float ycenter1 = (ePC[0].Ymin()+ePC[0].Ymax())/2.;
    float xoffset= aXY->A11()*xcenter1 + aXY->A12()*ycenter1 + aXY->B1() - xcenter1;
    float yoffset= aXY->A21()*xcenter1 + aXY->A22()*ycenter1 + aXY->B2() - ycenter1;
    
    eCorrL[0].SetXY( xcenter1, ycenter1 );
    
    const char *str = Form( "Nfinal= %5d    Peak: %5d/%6.3f    dx,dy,dz = %7.3f  %7.3f %7.3f",
			    eS[0].GetEntries(),(int)eH_xy_coarse.ePeak[0], eH_xy_coarse.eMean[0], 
			    xoffset, yoffset, eCorrL[0].Zcorr() );
    Log(1,"Alignment Report","%s", str);

  if(eOutputFile) {
    Log(2,"Alignment Report","Save to file %s", eOutputFile->GetName());
    gStyle->SetPalette(1);
    bool batch = gROOT->IsBatch();
    gROOT->SetBatch();
    
    TCanvas *cc = new TCanvas("cc","Alignment report",800,800);
    
    
    TPaveText *ctit = new TPaveText(0.01,0.943,0.99,0.998);
    ctit->AddText( Form("Alignment of  %s",eOutputFile->GetName()) );
    ctit->AddText( str );
    ctit->Draw();
    
    TPad *c = new TPad("c","plots",0.01,0.05,0.99,0.94);
    c->Divide(3,3);
    c->Draw();
    TH2F *hp1 = ePC[0].DrawH2("hp1","Pattern1"); c->cd(1); hp1->SetStats(0); hp1->Draw("colz");
    TH2F *hp2 = ePC[1].DrawH2("hp2","Pattern2"); c->cd(2); hp2->SetStats(0); hp2->Draw("colz");
    if(zphiCoarse)  { c->cd(3); zphiCoarse->SetStats(0); zphiCoarse->Draw("colz"); }
    if(xyCoarse)    { c->cd(4); xyCoarse->SetStats(0);   xyCoarse->Draw("colz");   }
    if(xyFine)      { c->cd(5); xyFine->SetStats(0);     xyFine->Draw("colz");     }


    float xmin = Min( ePC[0].Xmin(), ePC[1].Xmin() );
    float xmax = Max( ePC[0].Xmax(), ePC[1].Xmax() );
    float ymin = Min( ePC[0].Ymin(), ePC[1].Ymin() );
    float ymax = Max( ePC[0].Ymax(), ePC[1].Ymax() );

    Log(2,"Alignment Report","xmax-xmin: %f (%f-%f) ymax-ymin: %f (%f-%f)", 
	(xmax-xmin),xmax,xmin, (ymax-ymin),ymax,ymin);


    int n = CheckEqualArr(eS[0],eS[1]);

    c->cd(6);

    TH2I *hdtxy  = new TH2I("dTXdTY","dTY vs dTX after final alignment",25,-3*eSigma[1], 3*eSigma[1], 25, -3*eSigma[1], 3*eSigma[1]);
    for(int i=0; i<n; i++) {
      EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
      EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);
      float dtx  =   TX(1,*s2) - TX(0,*s1);
      float dty  =   TY(1,*s2) - TY(0,*s1);
      hdtxy->Fill(dtx,dty);
    }
    hdtxy->SetStats(0); hdtxy->Draw("colz");

    c->cd(8);

    TH2I *h2  = new TH2I("plateXY","plateXY",100,xmin-10, xmax+10, 100, ymin-10,ymax+10);
    h2->SetStats(0);
    h2->Draw();

    if(n<1000) {
      for(int i=0; i<n; i++) {
        EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
        EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);
        float x1  =  X(0,*s1),  y1 =  Y(0,*s1);
        float x2  =  X(1,*s2),  y2 =  Y(1,*s2);
        TMarker *m1 = new TMarker(x1,y1,22);
        m1->SetMarkerColor(kRed);
        m1->SetMarkerSize(0.6);
        m1->Draw();
        TMarker *m2 = new TMarker(x2,y2,23);
        m2->SetMarkerColor(kBlue);
        m2->SetMarkerSize(0.6);
        m2->Draw();
      }
    }

    c->cd(9);
    TH2I *h2t = new TH2I("plateTXTY","plateTXTY",100,-1, 1, 100, -1,1 );
    h2t->SetStats(0);
    h2t->Draw();
    if(n<1000) {
      for(int i=0; i<n; i++) {
        EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
        EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);
        float tx1 =  TX(0, *s1), ty1 =  TY(0, *s1);
        float tx2 =  TX(1, *s2), ty2 =  TY(1, *s2);
        TMarker *m1 = new TMarker(tx1,ty1,22);
        m1->SetMarkerColor(kRed);
        m1->SetMarkerSize(0.6);
        m1->Draw();
        TMarker *m2 = new TMarker(tx2,ty2,23);
        m2->SetMarkerColor(kBlue);
        m2->SetMarkerSize(0.6);
        m2->Draw();
      }
    }
    

    c->cd(7);
    int nbin=25;
    TArrayF xbins(nbin+1);
    xbins[0]= 0.15;
    float s0 = xbins[0]*xbins[0];
    for(int i=1; i<nbin+1; i++) xbins[i] = Sqrt( xbins[i-1]*xbins[i-1] + s0 );

    TH1F *thetadens = new TH1F("thetadens","theta density",nbin,xbins.GetArray() );
    for(int i=0; i<n; i++) {
      //EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
      EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);
      float tx  =   TX(1,*s2), ty  =   TY(1,*s2);
      float theta = Sqrt(tx*tx+ty*ty);
      thetadens->Fill(theta);
    }
    thetadens->Draw();

    eH_xy_coarse.Write("peak2c");
    cc->Write("report_al");
    eCorrL[0].Write("corr_layer1");
    eCorrL[1].Write("corr_layer2");
    SafeDelete(c);
    
    SafeDelete(cc);
    
    gROOT->SetBatch(batch);
  }

  SafeDelete(zphiCoarse);
  SafeDelete(xyCoarse);
  SafeDelete(xyFine);
}

//---------------------------------------------------------------------
void EdbPlateAlignment::DoubletsFilterOut(EdbPattern &p1, EdbPattern &p2)
{
  EdbAlignmentV adup;
  for(int i=0; i<4; i++) adup.eDVsame[i]=eDoublets[i];
  
  adup.FillGuessCell(p1,p1,1.,eOffsetMax);
  adup.FillCombinations();
  adup.DoubletsFilterOut(0);   // assign flag -10 to the duplicated segments

  adup.FillGuessCell(p2,p2,1.,eOffsetMax);
  adup.FillCombinations();
  adup.DoubletsFilterOut(0);   // assign flag -10 to the duplicated segments
}

//---------------------------------------------------------------------
void EdbPlateAlignment::FineAlAff(EdbPattern &p1, EdbPattern &p2, EdbLayer &la1)
{
  // assuming that the patterns are already aligned with the accuracy of a few 
  // microns or the corrections are setted.
  // find affine transfornations
  Log(3, "FineAlAff","with patterns %d %d", p1.N(), p2.N() );
  FillGuessCell(p1,p2,1.,eOffsetMax);
  float dxlim = eH[0][0].Xmax()-eH[0][0].Xmin();
  float dylim = eH[0][1].Xmax()-eH[0][1].Xmin();

  eDVsame[0] = dxlim/2.;  eDVsame[1] = dylim/2.;
  //eDVsame[0] = dxlim;  eDVsame[1] = dylim;
  eDVsame[2]=eDVsame[3]=3.*eSigma[1];
  FillCombinations( eDVsame, dxlim, dylim, 1);

  int      n[2] = { 9,9 };
  float  min[2] = { eH[0][0].Xmin(), eH[0][1].Xmin() };
  float  max[2] = { eH[0][0].Xmax(), eH[0][1].Xmax() };
  eHxy.InitH2(n, min, max);

  TObjArray sel1, sel2;
  int npk= Ncoins(eDVsame, &eHxy, 0, &sel1, &sel2);
  EdbAffine2D aff;
  if(eNoScale)         CalculateAffXYTurn(sel1,sel2,aff);
  else if(eNoScaleRot) CalculateAffXYShift(sel1,sel2,aff);
  else                 CalculateAffXY(sel1,sel2,aff);
  la1.GetAffineXY()->Transform(&aff);
  
  if(eDoCorrectAngle) {
    EdbAffine2D afftxty;
    CalculateAffTXTYTurn(sel1,sel2,afftxty);
    la1.GetAffineTXTY()->Transform(&afftxty);
  }
  
  Log(2,"FineAlAff","peak of %d  with patterns of: %d %d",  npk, p1.N(),p2.N() );
  if(gEDBDEBUGLEVEL>2)  la1.GetAffineXY()->Print();

  if(npk>3)  eFineOK = true;
  else       eFineOK = false;
}

//---------------------------------------------------------------------
void EdbPlateAlignment::FineAl(EdbPattern &p1, EdbPattern &p2)
{
  Log(3, "FineAl","with patterns %d %d", p1.N(), p2.N() );
  FillGuessCell(p1,p2,1.,eOffsetMax);
  float dxlim = eH[0][0].Xmax()-eH[0][0].Xmin();
  float dylim = eH[0][1].Xmax()-eH[0][1].Xmin();

  eDVsame[0] = dxlim/2.;  eDVsame[1] = dylim/2.;
  //eDVsame[0] = dxlim;  eDVsame[1] = dylim;

  //eDVsame[2]=eDVsame[3]=0.01;
  eDVsame[2]=eDVsame[3]=3.*eSigma[1];
  FillCombinations( eDVsame, dxlim, dylim, 1);

  int      n[2] = { 9,9 };
  float  min[2] = { eH[0][0].Xmin(), eH[0][1].Xmin() };
  float  max[2] = { eH[0][0].Xmax(), eH[0][1].Xmax() };
  eHxy.InitH2(n, min, max);

  for(int i=0; i<3; i++) FindCorrDiff(eDVsame);

  TObjArray sel1,sel2;
  Ncoins(eDVsame, &eHxy, 0, &sel1, &sel2);

  if( eDZ>0.0001 ) {
    for(int i=0; i<3; i++) {
      float dz = FineCorrZ(sel1,sel2);
      eCorr[0].AddV(2,dz);
    }
  }
  
//   for(int i=0; i<3; i++) {
//     float dphi = FineCorrPhi(sel1,sel2);
//     eCorr[0].AddV(6,dphi);
//   }

  int npk=Ncoins(eDVsame, &eHxy);

  Log(2,"FineAl","peak of %d with offsets (dx:dy:dz:dphi): %7.2f %7.2f %8.2f %7.3f  with patterns of: %d %d", 
      npk, 
      eCorr[0].V(0)-eCorr[1].V(0), eCorr[0].V(1)-eCorr[1].V(1), eCorr[0].V(2), eCorr[0].V(6), 
      p1.N(),p2.N() );


  if(npk>3)  eFineOK = true;
  else       eFineOK = false;
}

//---------------------------------------------------------------------
void EdbPlateAlignment::CoarseAl(EdbPattern &p1, EdbPattern &p2)
{
  // all pars should be already defined
  
  if( Log(3, "CoarseAl","with patterns %d %d", p1.N(), p2.N() ) ) {
    eCorr[0].Print();    eCorr[1].Print();
  }
  FillGuessCell(p1,p2,1.,eOffsetMax);
  FillCombinations( eDVsame, eDVsame[0], eDVsame[1], 1);
  int      n[2] = { eH[0][0].N()   , eH[0][1].N()    };
  float  min[2] = { eH[0][0].Xmin(), eH[0][1].Xmin() };
  float  max[2] = { eH[0][0].Xmax(), eH[0][1].Xmax() };
  eHxy.InitH2(n, min, max);

  //OptimiseVar1(0, 6, &eHxy);
  //OptimiseVar1(0, 2, &eHxy);

  EdbH2  h_z_phi;
  OptimiseVar2( 0, 6,  0, 2, h_z_phi, &eHxy );  // variate rotation and z of the p1

  eH_zphi_coarse.Init(h_z_phi);

  EdbPeak2 pk2z(h_z_phi);
  float  dphi, dz;
  pk2z.FindPeak( dphi, dz );
  pk2z.ProbPeak();
  eCorr[0].SetV(2, dz  );
  eCorr[0].SetV(6, dphi);

  Ncoins(eDVsame, &eHxy);
  Log(2,"CoarseAl","DVsame: %s",StrDVsame()); 
  eH_xy_coarse.Init(eHxy);

  EdbPeak2 pk2(eHxy);
  float  dx, dy;
  int npk = (int)pk2.FindPeak9( dx, dy );
  eCorr[0].AddV(0,dx);   eCorr[0].AddV(1,dy);
  pk2.ProbPeak();
  //pk2.Print();

  Log(2,"CoarseAl","%d coins with peak of %d/%6.3f with offsets (dx:dy:dz:dphi): %7.2f %7.2f %8.2f %7.3f  with patterns of: %d %d", 
      npk, (int)pk2.ePeak[0], pk2.eMean[0],
      eCorr[0].V(0)-eCorr[1].V(0), eCorr[0].V(1)-eCorr[1].V(1), eCorr[0].V(2), eCorr[0].V(6), 
      p1.N(),p2.N() );

  if(npk>3)  eCoarseOK = true;
  else       eCoarseOK = false;
}

//---------------------------------------------------------------------
void EdbPlateAlignment::TestAl(EdbPattern &p1, EdbPattern &p2)
{
  // alignment of 2 plates with testal-like preselection. 
  // Application: basetracks patterns with: well defined angles, big possible offsets, 
  // big possible rotations; z-variation
  // Result: find best dz and rotation and position offset and set them as: 
  // eCorr[0].SetV(2, dz),  eCorr[0].SetV(6, dphi); eCorr[0].AddV(0, dx); eCorr[0].AddV(1, dy);

  float dxlim = eH[0][0].Xmax()-eH[0][0].Xmin();
  float dylim = eH[0][1].Xmax()-eH[0][1].Xmin();
  HDistance(p1,p2, dxlim, dylim );

  float x0 = p2.Xmean(), y0 = p2.Ymean();           // centralize patterns
  eCorr[0].SetV(0,-x0);  eCorr[0].SetV(1,-y0);
  eCorr[1].SetV(0,-x0);  eCorr[1].SetV(1,-y0);
    
  int      n[2] = { eH[0][0].N()   , eH[0][1].N()        };
  float  min[2] = { eH[0][0].Xmin(), eH[0][1].Xmin() };
  float  max[2] = { eH[0][0].Xmax(), eH[0][1].Xmax() };
  eHxy.InitH2(n, min, max);

  eDVsame[0] = dxlim/2.;
  eDVsame[1] = dylim/2.;
  eDVsame[2]=eDVsame[3]=eSigma[1];

  EdbH2  h_z_phi;
  OptimiseVar2( 0, 6,  0, 2, h_z_phi, &eHxy );  // variate rotation and z of the p1

  h_z_phi.DrawH2()->Draw("colz");

  EdbPeak2 pk2z(h_z_phi);
  float  dphi, dz;
  pk2z.FindPeak( dphi, dz );
  pk2z.ProbPeak();

  eCorr[0].SetV(2, dz  );
  eCorr[0].SetV(6, dphi);

  Ncoins(eDVsame, &eHxy);
  EdbPeak2 pk2(eHxy);
  float  dx, dy;
  pk2.FindPeak( dx, dy );
  pk2.ProbPeak();
  eCorr[0].AddV(0,dx);      // set x-y corrections
  eCorr[0].AddV(1,dy);
 
  Log(2,"TestAl","peak of %d/%6.3f with offsets (dx:dy:dz:dphi): %7.2f %7.2f %8.2f %7.3f  with patterns of: %d %d", 
      (int)pk2.ePeak[0], pk2.eMean[0], 
      eCorr[0].V(0)-eCorr[1].V(0), eCorr[0].V(1)-eCorr[1].V(1), eCorr[0].V(2), eCorr[0].V(6), 
      p1.N(),p2.N() );
}

//---------------------------------------------------------------------
void EdbPlateAlignment::SetParFineAl()
{
  // narrow parameters setting for the fine alignment
  InitHx(   1,  -4.*eSigma[0],   4.*eSigma[0]   );
  InitHy(   1,  -4.*eSigma[0],   4.*eSigma[0]   );
}


//---------------------------------------------------------------------
void EdbPlateAlignment::SetParCoarseAl(float zcorr, float dpos, float dang, float dz, float dphi)
{
  // medium-wide parameters setting for the coarse alignment
  Log(3,"EdbPlateAlignment::SetParCoarseAl","zcorr = %f  dpos = %f dang=%f  dz=%f  dphi=%f",zcorr, dpos, dang, dz, dphi);
  eDVsame[0] = eDVsame[1] = dpos;
  //eDVsame[0] = eDVsame[1] = 3*eSigma[0];
  eDVsame[2] = eDVsame[3] = dang;
  int n = (int)( 2*dpos/ (4*eSigma[0]) );
  n = n%2 ? n : n+1;              // make n - odd
  Log(3, "SetParCoarseAl"," n = %d dpos = %f ", n, dpos);
  InitHx(   n,  -dpos,   dpos   );
  InitHy(   n,  -dpos,   dpos   );

  float zstepmin=0.5;   // 0.5 microns is the minimum step size
  int nz = Min( (int)( 2*dz/ zstepmin ), 51);
  nz = nz%2 ? nz : nz+1;              // make n - odd
  InitHz(   nz,  zcorr-dz, zcorr+dz );
  eCorr[0].SetV(2,zcorr);

  float phimin=0.00001;   // minimum step size
  int nphi = Min( (int)( 2*dphi/ phimin ), 51);
  nphi = nphi%2 ? nphi : nphi+1;             // make n - odd
  InitHphi( nphi,  eCorr[0].V(6)-dphi,   eCorr[0].V(6)+dphi );
}

//---------------------------------------------------------------------
void EdbPlateAlignment::SetParTestAl(float zcorr, float dz, float dphi)
{
  // wide parameters setting for the test pre-alignment
  InitHx(   81,  -2000,   2000   );
  InitHy(   81,  -2000,   2000   );
  int nbinz = Min(101, Max(21, (int)(dz/50.) ));
  InitHz(   nbinz,   zcorr-dz, zcorr+dz );
  int nbinphi = Min(101, Max(21, (int)(dphi/0.0001) ));
  InitHphi( nbinphi,  -dphi,   dphi   );
  eCorr[0].SetV(2,zcorr);
}

//---------------------------------------------------------------------
void EdbPlateAlignment::SlowAlignXY(EdbPattern &pat1, EdbPattern &pat2, EdbH2 &hxy, EdbH1 &hphi, const char *name )
{
  int n1= pat1.N();
  int n2= pat2.N();
  float x0 = (Min(pat1.Xmin(),pat2.Xmin())+Max(pat1.Xmax(),pat2.Xmax()))/2;
  float y0 = (Min(pat1.Ymin(),pat2.Ymin())+Max(pat1.Ymax(),pat2.Ymax()))/2;
  struct point{ int i; float x; float y; };
  point **p1 = new point*[n1];
  point **p2 = new point*[n2];
 

  TFile f( Form("%s.root",name),"RECREATE");
  TNtuple *ntpeak = new TNtuple("ntpeak","","phi:dx:dy:peak");
  
  for(int i=0; i<n1; i++)     p1[i] = new point;
  for(int i=0; i<n2; i++)     p2[i] = new point;
  
  for(int iphi=0; iphi<hphi.N(); iphi++) {
    float phi=hphi.X(iphi);
    for(int i=0; i<n1; i++) {
      EdbSegP *s= pat1.GetSegment(i);
      p1[i]->x = (s->X()-x0)*Cos(phi) - (s->Y()-y0)*Sin(phi);
      p1[i]->y = (s->X()-x0)*Sin(phi) + (s->Y()-y0)*Cos(phi);
      p1[i]->i = i;
    }
    for(int i=0; i<n2; i++) {
      EdbSegP *s= pat2.GetSegment(i);
      p2[i]->x = (s->X()-x0); 
      p2[i]->y = (s->Y()-y0); 
      p2[i]->i = i;
    }

    for(int i1=0; i1<n1; i1++) {
      for(int i2=0; i2<n2; i2++) {
        float dx= p2[i2]->x-p1[i1]->x;
        float dy= p2[i2]->y-p1[i1]->y;
        hxy.Fill(dx,dy);
      }
    }
    EdbPeak2 pk2(hxy);
    float  dx, dy;
    int npk = pk2.FindPeak( dx, dy );
    printf("phi = %f  dx = %f dy = %f  npk = %d\n",phi,dx,dy,npk);
    hphi.Fill(phi,npk);
    ntpeak->Fill(phi,dx,dy,npk);
    hxy.Write(Form("hxy%d",iphi));
    hxy.CleanCells();
  }
  hphi.Write("hphi");
  ntpeak->Write();
  f.Close();
}

//---------------------------------------------------------------------
void EdbPlateAlignment::RankCouples( TObjArray &arr1,TObjArray &arr2 )
{
  int n = arr1.GetEntries();
  Log(3,"RankCouples","%d" ,n);

  EdbTrackFitter tf;
  EdbScanCond cond;
  cond.SetSigma0(eSigma[0],eSigma[0],eSigma[1],eSigma[1]);
  eSegCouples.Delete();
  EdbSegP seg, seg1, seg2;
  for(int i=0; i<n; i++) {
    EdbSegP *s1 = ((EdbSegP*)arr1.UncheckedAt(i));
    EdbSegP *s2 = ((EdbSegP*)arr2.UncheckedAt(i));
    
    seg.Copy(*s1);   // to set correctly vid, aid, etc
    seg1.Copy(*s1);
    seg2.Copy(*s2);
    
    eCorr[0].ApplyCorrections(seg1);
    eCorr[1].ApplyCorrections(seg2);
   
    tf.Chi2PSeg( seg1, seg2, seg, cond, cond);
    
    s1->SetFlag(0);
    s2->SetFlag(0);
    seg.SetFlag(0);
    seg.SetSide(0);
    seg.SetVolume(seg1.Volume()+seg2.Volume());
    seg.SetDZ( Abs(seg1.eZ - seg2.eZ) );
    
    seg.SetDZem( seg1.Chi2() + seg2.Chi2() );    // HACK: use eDZem variable to keep microtracking Likelihood 
    
    EdbSegCouple *sc=new EdbSegCouple();
    sc->eS1=s1;
    sc->eS2=s2;
    sc->eS = new EdbSegP(seg);
    sc->SetCHI2P( seg.Chi2() );
    eSegCouples.Add(sc);
  }

  EdbSegCouple::SetSortFlag(0);    // sort by CHI2P
  eSegCouples.UnSort();
  eSegCouples.Sort();

  int ncp = eSegCouples.GetEntries();

  for(int i=0; i<ncp; i++) {
    EdbSegCouple *sc = (EdbSegCouple*)(eSegCouples.UncheckedAt(i));
    sc->eS1->SetFlag( sc->eS1->Flag()+1 );
    sc->eS2->SetFlag( sc->eS2->Flag()+1 );
    sc->SetN1(sc->eS1->Flag());
    sc->SetN2(sc->eS2->Flag());
  }

  for(int i=0; i<ncp; i++) {
    EdbSegCouple *sc = (EdbSegCouple*)(eSegCouples.UncheckedAt(i));
    sc->SetN1tot(sc->eS1->Flag());
    sc->SetN2tot(sc->eS2->Flag());
  }

  Log(2,"RankCouples","%d couples ok", ncp );
}
