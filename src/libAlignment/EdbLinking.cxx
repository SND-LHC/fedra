//-- Author :  Valeri Tioukov   4-03-2011
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// EdbLinking                                                             //
//                                                                        //
// microtracks linking inside one plate                                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "EdbLog.h"
#include "EdbLinking.h"
#include "EdbCouplesTree.h"
#include "EdbTrackFitter.h"
#include "EdbLayer.h"
#include "EdbScanCond.h"
#include "EdbSEQ.h"

ClassImp(EdbLinking)

using namespace TMath;

//---------------------------------------------------------------------
EdbLinking::EdbLinking()
{
}

//---------------------------------------------------------------------
void EdbLinking::GetPar(TEnv &env)
{
  eBinOK               = env.GetValue("fedra.link.BinOK"              , 6. );

  eNcorrMin            = env.GetValue("fedra.link.NcorrMin"           , 100 );

  eDoCorrectShrinkage  = env.GetValue("fedra.link.DoCorrectShrinkage" , true );
  eNsigmaEQshr         = env.GetValue("fedra.link.shr.NsigmaEQ"        , 7. );
  eShr0                = env.GetValue("fedra.link.shr.Shr0"            , .9   );
  eDShr                = env.GetValue("fedra.link.shr.DShr"            , .25  );

  eDoCorrectAngles     = env.GetValue("fedra.link.DoCorrectAngles"    , true );
  eChi2Acorr           = env.GetValue("fedra.link.ang.Chi2max"        , 1.5   );

  eDoFullLinking       = env.GetValue("fedra.link.DoFullLinking"      , true );
  eNsigmaEQlnk         = env.GetValue("fedra.link.full.NsigmaEQ"        , 5. );
  eDRfull              = env.GetValue("fedra.link.full.DR"             , 30.  );
  eDTfull              = env.GetValue("fedra.link.full.DT"             , 0.1  );
  eCHI2Pmax            = env.GetValue("fedra.link.full.CHI2Pmax"       , 3.   );

  eDoSaveCouples       = env.GetValue("fedra.link.DoSaveCouples"       , true );
}

//---------------------------------------------------------------------
void EdbLinking::GetPreselectionPar(EdbSEQ &seq, TEnv &env)
{
  const char *str=0;
  float x1,x2;
  str=env.GetValue("fedra.link.shr.XLimits","");     if(str) if(sscanf(str,"%f %f",&x1,&x2) == 2)  seq.SetXLimits(x1,x2);
  str=env.GetValue("fedra.link.shr.YLimits","");     if(str) if(sscanf(str,"%f %f",&x1,&x2) == 2)  seq.SetYLimits(x1,x2);
  str=env.GetValue("fedra.link.shr.WLimits","");     if(str) if(sscanf(str,"%f %f",&x1,&x2) == 2)  seq.SetWLimits(x1,x2);
  str=env.GetValue("fedra.link.shr.ThetaLimits",""); if(str) if(sscanf(str,"%f %f",&x1,&x2) == 2)  seq.SetThetaLimits(x1,x2);
  str=env.GetValue("fedra.link.shr.ChiLimits","");   if(str) if(sscanf(str,"%f %f",&x1,&x2) == 2)  seq.SetChiLimits(x1,x2);
}

//---------------------------------------------------------------------
Double_t EdbLinking::EstimatePatternArea(EdbPattern &p)
{
  return (p.Xmax()-p.Xmin())*(p.Ymax()-p.Ymin());
}

//---------------------------------------------------------------------
void EdbLinking::Link(EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2, TEnv &env, Double_t area1, Double_t area2 )
{
   // main linking function
  // Input: p1,p2 - patterns, l1,l2 layers, env - linking parameters

  Log(2,"EdbLinking::Link","patterns with %d and %d segments and z1 = %f  z2 = %f",p1.N(),p2.N(),l1.Z(),l2.Z());
  Log(2,"EdbLinking::Link","segments z1 = %f  z2 = %f",p1.GetSegment(0)->Z(),p2.GetSegment(0)->Z());

  GetPar(env);
  eL1 = l1;
  eL2 = l2;

  if(area1<=0) area1 = EstimatePatternArea(p1);
  if(area2<=0) area2 = EstimatePatternArea(p2);
  EdbSEQ seq;
  GetPreselectionPar(seq,env);
  seq.eNsigma=eNsigmaEQshr;
  seq.PrintLimits();

  TObjArray p1pre(p1.N()), p1shr;
  seq.PreSelection(p1,p1pre);
  seq.EqualizeMT(p1pre, p1shr, area1);
  TH1F *hTshr1 = seq.ThetaPlot(p1shr, "hTshr1","Theta plot shr, side 1 "); hTshr1->Write();
  TH1F *hTall1 = seq.ThetaPlot(p1, "hTall1","Theta plot all, side 1 ");    hTall1->Write();
  
  TObjArray p2pre(p2.N()), p2shr;
  seq.PreSelection(p2,p2pre);
  seq.EqualizeMT(p2pre, p2shr, area2);
  TH1F *hTshr2 = seq.ThetaPlot(p2shr, "hTshr2","Theta plot shr, side 2 "); hTshr2->Write();
  TH1F *hTall2 = seq.ThetaPlot(p2, "hTall2","Theta plot all, side 2 ");    hTall2->Write();

  DoubletsFilterOut(p1shr,p2shr);
  

  FillCombinationsAtMeanZ(p1shr,p2shr);

  if(eDoCorrectShrinkage) {
    eCorr[0].SetV(5,eShr0);
    eCorr[1].SetV(5,eShr0);
    CorrectShrinkage( eDShr );
    CorrectShrinkage( eDShr*0.8 );
    if( eCorr[0].V(5) > eShr0*(1+eDShr) ||  eCorr[0].V(5)< eShr0/(1+eDShr) ) {
      Log(1,"EdbLinking::Link","Shrinkage correction side 1 out of range: %f [%f %f] reset to default: %f", 
             eCorr[0].V(5), eShr0/(1+eDShr), eShr0*(1+eDShr), eShr0 );
      eCorr[0].SetV(5,eShr0);
      eDoCorrectShrinkage = false;
    }
    if( eCorr[1].V(5) > eShr0*(1+eDShr) ||  eCorr[1].V(5)< eShr0/(1+eDShr) ) {
      Log(1,"EdbLinking::Link","Shrinkage correction side 2 out of range: %f [%f %f] reset to default: %f", 
             eCorr[1].V(5), eShr0/(1+eDShr), eShr0*(1+eDShr), eShr0 );
      eCorr[1].SetV(5,eShr0);
      eDoCorrectShrinkage = false;
    }
  }
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectShrinkage) {
    CorrectShrinkage(  eDShr*0.5 );
    CorrectShrinkage(  eDShr*0.5 );
    CorrectShrinkage(  eDShr*0.5 );
    WriteShrinkagePlots();
  }
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );

  eCorr[0].Write("corr1");
  eCorr[1].Write("corr2");
  ePC[0].DrawH2("hxy_shr1","xy for the shrinkage corr sample side 1")->Write();
  ePC[1].DrawH2("hxy_shr2","xy for the shrinkage corr sample side 2")->Write();
  EdbH2 htxy1(50,-1,1,50,-1,1);
  FillThetaHist(0,htxy1);
  htxy1.DrawH2("htxy_shr1","txy plot for shr corr sample side 1");
  EdbH2 htxy2(50,-1,1,50,-1,1);
  FillThetaHist(0,htxy2);
  htxy2.DrawH2("htxy_shr2","txy plot for shr corr sample side 2");

  seq.ResetExcludeThetaRange();
  if(eDoFullLinking)   {
    seq.eNsigma = eNsigmaEQlnk;
    TObjArray p1lnk(p1pre.GetEntriesFast());
    seq.EqualizeMT(p1pre, p1lnk, area1);
    TObjArray p2lnk(p2pre.GetEntriesFast());
    seq.EqualizeMT(p2pre, p2lnk, area2);
    FullLinking(p1lnk,p2lnk);
  }

  ProduceReport();
  Corr2Aff(eCorr[0],eL1);
  Corr2Aff(eCorr[1],eL2);
  eL1.Print();                     // layers with the corrections
  eL2.Print();
}

//---------------------------------------------------------------------
void EdbLinking::FillCombinationsAtMeanZ(TObjArray &p1, TObjArray &p2)
{
  float dz = eL2.Z()-eL1.Z();
  int n1 = p1.GetEntries(), n2=p2.GetEntries();
  if(n1<1||n2<1) return;
  eCorr[0].SetV(2,dz/2.);  eCorr[1].SetV(2, -dz/2.);
  FillGuessCell(p1,p2,eBinOK,0);
  float dvcomb[4] = {50,50, 0.1, 0.1};
  FillCombinations(dvcomb, 50,50, 1);
}

//---------------------------------------------------------------------
void EdbLinking::CorrectAngles(TObjArray &p1, TObjArray &p2)
{
    Log(2,"EdbLinking::CorrectAngles","arrays with %d and %d segments",p1.GetEntries(),p2.GetEntries());

    eCorr[0].SetV(2,0);  eCorr[1].SetV(2,0);
    RankCouples( eS[0], eS[1] );

    int nc=0;
    double dtx1=0, dty1=0,dtx2=0, dty2=0;
    int ncp = eSegCouples.GetEntries();
    for(int i=0; i<ncp; i++) {
      EdbSegCouple *sc = (EdbSegCouple*)(eSegCouples.At(i));
      if(sc->CHI2P()>eChi2Acorr) continue;
      dtx1+= TX( 0, *(sc->eS1) ) - sc->eS->TX();
      dty1+= TY( 0, *(sc->eS1) ) - sc->eS->TY();
      dtx2+= TX( 1, *(sc->eS2) ) - sc->eS->TX();
      dty2+= TY( 1, *(sc->eS2) ) - sc->eS->TY();
      nc++;
    }
    if(nc<eNcorrMin) {Log(1,"EdbLinking::CorrectAngles","Warning: number of the selected segments too small: %d < %d do nothing",nc,eNcorrMin); return;}
    
    float cc=1.8;
    dtx1 /= nc;  dty1 /= nc;  dtx2 /= nc;  dty2 /= nc;
    eCorr[0].AddV(3,-cc*dtx1);     eCorr[0].AddV(4,-cc*dty1);
    eCorr[1].AddV(3,-cc*dtx2);     eCorr[1].AddV(4,-cc*dty2);
    
    Log(2,"EdbLinking::CorrectAngles","using %d segments dx1,dy1:( %6.3f %6.3f )  dx2,dy2:( %6.3f %6.3f )",
    nc, dtx1*cc, dty1*cc, dtx2*cc, dty2*cc );
}

//---------------------------------------------------------------------
void EdbLinking::FullLinking(TObjArray &p1, TObjArray &p2)
{
    float z1 = eL1.Z(),  z2 = eL2.Z(),  dz = z2-z1;
    Log(2,"EdbLinking::FullLinking","arrays with %d and %d segments and z1 = %f  z2 = %f",p1.GetEntries(),p2.GetEntries(),z1,z2);
    Log(2,"EdbLinking::FullLinking","segments z1 = %f  z2 = %f",((EdbSegP*)(p1.At(0)))->Z(),((EdbSegP*)(p2.At(0)))->Z());
     
    float dvcomb[4] = {eDRfull,eDRfull, eDTfull, eDTfull};
    eCorr[0].SetV(2,dz/2.);  eCorr[1].SetV(2,   -dz/2.);
    FillGuessCell(p1,p2,eBinOK,0);
    FillCombinations(dvcomb, eDRfull,eDRfull, 1);
    eCorr[0].SetV(2,0);  eCorr[1].SetV(2, 0);
    RankCouples( eS[0], eS[1] );
    SaveCouplesTree();
}

//---------------------------------------------------------------------
void EdbLinking::FullLinking(EdbPattern &p1, EdbPattern &p2)
{
    float z1 = eL1.Z(),  z2 = eL2.Z(),  dz = z2-z1;
    Log(2,"EdbLinking::FullLinking","patterns with %d and %d segments and z1 = %f  z2 = %f",p1.N(),p2.N(),z1,z2);
    
    float dvcomb[4] = {eDRfull,eDRfull, eDTfull, eDTfull};
    eCorr[0].SetV(2,dz/2.);  eCorr[1].SetV(2,   -dz/2.);
    FillGuessCell(p1,p2,eBinOK,0);
    FillCombinations(dvcomb, eDRfull,eDRfull, 1);
    eCorr[0].SetV(2,0);  eCorr[1].SetV(2, 0);
    RankCouples( eS[0], eS[1] );
    SaveCouplesTree();
}

//---------------------------------------------------------------------
void EdbLinking::SaveCouplesTree()
{
  EdbCouplesTree ect;
  ect.InitCouplesTree("couples",0,"NEW");
  
  ect.eTree->SetAlias("dz1"  ,"107.");
  ect.eTree->SetAlias("dz2"  ,"-107.");
  ect.eTree->SetAlias("x1"   ,"s1.eX+dz1*s1.eTX");
  ect.eTree->SetAlias("x2"   ,"s2.eX+dz2*s2.eTX");
  ect.eTree->SetAlias("y1"   ,"s1.eY+dz1*s1.eTY");
  ect.eTree->SetAlias("y2"   ,"s2.eY+dz2*s2.eTY");
  ect.eTree->SetAlias("tx1c" , Form("s1.eTX/(%f)+(%f)",eCorr[0].V(5),eCorr[0].V(3)) );
  ect.eTree->SetAlias("ty1c" , Form("s1.eTY/(%f)+(%f)",eCorr[0].V(5),eCorr[0].V(4)) );
  ect.eTree->SetAlias("tx2c" , Form("s2.eTX/(%f)+(%f)",eCorr[1].V(5),eCorr[1].V(3)) );
  ect.eTree->SetAlias("ty2c" , Form("s2.eTY/(%f)+(%f)",eCorr[1].V(5),eCorr[1].V(4)) );
  ect.eTree->SetAlias("x1c"  ,  "s1.eX+dz1*tx1c" );
  ect.eTree->SetAlias("y1c"  ,  "s1.eY+dz1*ty1c" );
  ect.eTree->SetAlias("x2c"  ,  "s2.eX+dz2*tx2c" );
  ect.eTree->SetAlias("y2c"  ,  "s2.eY+dz2*ty2c" );
  ect.eTree->SetAlias("dx"   , "x2c-x1c" );
  ect.eTree->SetAlias("dy"   , "y2c-y1c" );
  ect.eTree->SetAlias("dtx"  , "tx2c-tx1c" );
  ect.eTree->SetAlias("dty"  , "ty2c-ty1c" );
  ect.eTree->SetAlias("dr"   , "sqrt(dx*dx+dy*dy)" );
  ect.eTree->SetAlias("dt"   , "sqrt(dtx*dtx+dty*dty)" );

  EdbH2 hxy; hxy.InitH2(ePC[0]);
  TH2F htxy("htxy_cp","s.eTY vs s.eTX for all couples",50,-1,1,50,-1,1);
  
  TH1F  hchi("hchi"  ,"chi of the couples",80,0,4);
  TH1F  hchi20("hchi20","chi of the couples [s.eW>20]",80,0,4);
  TH1F  hdtx1("hdtx1","s1.eTX-s.eTX [s.eW>20]",50,-0.2,0.2);
  TH1F  hdty1("hdty1","s1.eTY-s.eTY [s.eW>20]",50,-0.2,0.2);
  TH1F  hdtx2("hdtx2","s2.eTX-s.eTX [s.eW>20]",50,-0.2,0.2);
  TH1F  hdty2("hdty2","s2.eTY-s.eTY [s.eW>20]",50,-0.2,0.2);

  int ntr = eSegCouples.GetEntries();
  for(int i=0; i<ntr; i++) {
    EdbSegCouple *sc = (EdbSegCouple *)eSegCouples.At(i);
    eCorr[0].ApplyCorrections( *(sc->eS1) );
    eCorr[1].ApplyCorrections( *(sc->eS2) );
    
    hxy.Fill(sc->eS->X(), sc->eS->Y());
    htxy.Fill(sc->eS->TX(), sc->eS->TY());
    hchi.Fill(sc->CHI2P());
    if(sc->eS->W()>20) {
      hchi20.Fill(sc->CHI2P());
      hdtx1.Fill(sc->eS1->TX()-sc->eS->TX());
      hdty1.Fill(sc->eS1->TY()-sc->eS->TY());
      hdtx2.Fill(sc->eS2->TX()-sc->eS->TX());
      hdty2.Fill(sc->eS2->TY()-sc->eS->TY());
    }
    ect.Fill( sc->eS1, sc->eS2, sc->eS, sc );
  }
  ect.eTree->Write();
  hxy.DrawH2("hxy_cp","xy of the output couples")->Write();
  htxy.Write();
  hchi.Write();
  hchi20.Write();
  hdtx1.Write();
  hdtx2.Write();
  hdty1.Write();
  hdty2.Write();
}


//---------------------------------------------------------------------
void EdbLinking::CorrectShrinkage(TObjArray &p1, TObjArray &p2, float dshr)
{
  float dz = eL2.Z()-eL1.Z();
  int n1 = p1.GetEntries(), n2=p2.GetEntries();
  Log(2,"EdbLinking::CorrectShrinkage","arrays with %d and %d segments and z1 = %f  z2 = %f",n1,n2,eL1.Z(),eL2.Z());
  if(n1<1||n2<1) return;
  //DoubletsFilterOut(p1,p2);
  eCorr[0].SetV(2,dz/2.);  eCorr[1].SetV(2, -dz/2.);
  FillGuessCell(p1,p2,eBinOK,0);
  
  CorrectShrinkage(dshr);
}

//---------------------------------------------------------------------
void EdbLinking::CorrectShrinkage(EdbPattern &p1, EdbPattern &p2, float dshr)
{
  float dz = eL2.Z()-eL1.Z();
  Log(2,"EdbLinking::CorrectShrinkage","patterns with %d and %d segments and z1 = %f  z2 = %f",p1.N(),p2.N(),eL1.Z(),eL2.Z());
  if(p1.N()<1||p2.N()<1) return;
  //DoubletsFilterOut(p1,p2);
  eCorr[0].SetV(2,dz/2.);  eCorr[1].SetV(2, -dz/2.);
  FillGuessCell(p1,p2,eBinOK,0);

  CorrectShrinkage(dshr);
}

//---------------------------------------------------------------------
void EdbLinking::CorrectShrinkage(float dshr)
{
  EdbSegCorr c0 = eCorr[0], c1 = eCorr[1];  //save initial values
  float dz = eL2.Z()-eL1.Z();

  eDVsame[0] = eDVsame[1] = 30;  eDVsame[2] = eDVsame[3] = 0.03;
  InitHshr0(   25,  eCorr[0].V(5)/(1+dshr),   eCorr[0].V(5)*(1+dshr) );
  InitHshr1(   25,  eCorr[1].V(5)/(1+dshr),   eCorr[1].V(5)*(1+dshr) );
  
  eHdxyShr[0].InitH2(7, -35, 35, 7, -35, 35);
  eCorr[0].SetV(2,dz);  eCorr[1].SetV(2,   0);
  OptimiseVar1( 0, 5, &eHdxyShr[0] );              // variate shr1
  EdbSegCorr c0_opt = eCorr[0];

  eHdxyShr[1].InitH2(7, -35, 35, 7, -35, 35);
  eCorr[0]=c0;                                     // to remove the dependency on the correction order
  eCorr[0].SetV(2, 0);  eCorr[1].SetV(2, -dz);
  OptimiseVar1( 1, 5, &eHdxyShr[1]);               // variate shr2
  eCorr[0]=c0_opt;                                 // return to the optimized value

  Log(2,"EdbLinking::CorrectShrinkage","side1: %f (%f)   side2: %f (%f)", eCorr[0].V(5), c0.V(5), eCorr[1].V(5), c1.V(5) );
}

//---------------------------------------------------------------------
void EdbLinking::WriteShrinkagePlots()
{
  eH[0][5].DrawH1("shr1")->Write();
  eH[1][5].DrawH1("shr2")->Write();

  float dz = eL2.Z()-eL1.Z();
  float dvcomb[4] = {50,50, 0.03, 0.03};
  EdbH2 hdxy; hdxy.InitH2(25, -50, 50, 25, -50, 50);
  eCorr[0].SetV(2, dz/2.);  eCorr[1].SetV(2, -dz/2.);
  Int_t nc = Ncoins(dvcomb, &hdxy );
  hdxy.DrawH2("hdxy_shr","DY vs DX at the same Z for shrinkage correction sample")->Write();

  EdbH2 htxy; htxy.InitH2(35, -0.175, 0.175, 35, -0.175, 0.175);
  dvcomb[0]=dvcomb[1]=9;
  dvcomb[2]=dvcomb[3]=0.15;
  nc = Ncoins(dvcomb, 0, &htxy );
  htxy.DrawH2("htxy_shr","s2.eTY-s1.eTY vs s2.eTX-s1.eTX for shrinkage correction sample")->Write();
}

//---------------------------------------------------------------------
void EdbLinking::RankCouples( TObjArray &arr1,TObjArray &arr2 )
{
  int n = arr1.GetEntries();
  Log(3,"RankCouples","%d" ,n);

  EdbTrackFitter tf;
  EdbScanCond cond1,cond2;
  cond1.SetSigma0(1.,1., 0.013,0.013);  cond1.SetDegrad(4); cond1.SetPulsRamp0(6,9); cond1.SetPulsRamp04(6,9);
  cond2.SetSigma0(1.,1., 0.013,0.013);  cond2.SetDegrad(4); cond2.SetPulsRamp0(6,9); cond2.SetPulsRamp04(6,9);
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
   
    //tf.Chi2SegM( *s1, *s2, seg, cond1, cond2);
    tf.Chi2ASeg( seg1, seg2, seg, cond1, cond2);
    //seg.SetChi2( tf.Chi2ACP( *s1, *s2, cond1) );   //TODO test!!
    
    if(seg.Chi2() > eCHI2Pmax)  continue;
    
    s1->SetFlag(0);
    s2->SetFlag(0);
    
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

//---------------------------------------------------------------------
void EdbLinking::ProduceReport()
{
  if(eOutputFile) {
    Log(2,"Linking Report","Save to file %s", eOutputFile->GetName());
    gStyle->SetPalette(1);
    gStyle->SetOptStat(1);
    bool batch = gROOT->IsBatch();
    gROOT->SetBatch();
    
    TH1F *h1=0;    TH2F *h2=0;
    
    gStyle->SetOptDate(1);
    TCanvas *crep1 = new TCanvas("crep1","Linking report2",900,1000);

    TPaveText *ctit = new TPaveText(0.01,0.943,0.99,0.998);
    ctit->AddText( Form("Linking of  %s",eOutputFile->GetName()) );
 
    float dx,dy;
    EdbPeak2 pk_shr1(eHdxyShr[0]);
    pk_shr1.ProbPeak( dx, dy );
    EdbPeak2 pk_shr2(eHdxyShr[1]);
    pk_shr2.ProbPeak( dx, dy );
    const char *str = Form( "Shrinkage side1 is  %5.3f with peak of  %5.0f/%6.1f      side2 is  %5.3f with peak of %5.0f/%6.1f ",
                              eCorr[0].V(5),pk_shr1.Peak(0), pk_shr1.Mean(0),
                              eCorr[1].V(5),pk_shr2.Peak(0), pk_shr2.Mean(0));
    Log(1,"Shrinkage corr:","%s", str);
    ctit->AddText(str);
    ctit->AddText( Form("angular offsets: side1  %6.3f %6.3f   side2 %6.3f %6.3f",
                         eCorr[0].V(3),eCorr[0].V(4), eCorr[1].V(3), eCorr[1].V(4) ));
    ctit->Draw();

    TPad *c = new TPad("c","plots",0.01,0.03,0.99,0.94);
    c->Divide(4,4);    c->Draw();

    float densAll=0,densShr=0;
    
    h2 = (TH2F*)eOutputFile->Get("hxy_shr1");    if(h2) {c->cd(1); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h2 = (TH2F*)eOutputFile->Get("hxy_shr2");    if(h2) {c->cd(2); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    
    h1 = (TH1F*)eOutputFile->Get("hTall1");      if(h1) {c->cd(3);  h1->Draw(); densAll=h1->Integral();h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hTshr1");      if(h1) {c->cd(3);  h1->Draw("same"); densShr=h1->Integral(); h1=0;}
    c->cd(3);
    TPaveText *lable1 = new TPaveText(0.2,0.6,0.6,0.8,"NDC");
    lable1->AddText(Form("all: %7.0f",densAll));
    lable1->AddText(Form("shr: %7.0f",densShr));
    lable1->Draw();
    h1 = (TH1F*)eOutputFile->Get("hTall2");      if(h1) {c->cd(4);  h1->Draw(); densAll=h1->Integral(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hTshr2");      if(h1) {c->cd(4);  h1->Draw("same"); densShr=h1->Integral(); h1=0;}
    c->cd(4);
    TPaveText *lable2 = new TPaveText(0.2,0.6,0.6,0.8,"NDC");
    lable2->AddText(Form("all: %7.0f",densAll));
    lable2->AddText(Form("shr: %7.0f",densShr));
    lable2->Draw();
    
    h2 = (TH2F*)eOutputFile->Get("hdxy_shr");  if(h2) {c->cd(5); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h2 = (TH2F*)eOutputFile->Get("htxy_shr");  if(h2) {c->cd(6); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h1 = (TH1F*)eOutputFile->Get("shr1");      if(h1) {c->cd(7);  h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("shr2");      if(h1) {c->cd(8);  h1->Draw(); h1=0;}
    
    h1 = (TH1F*)eOutputFile->Get("hdtx1");     if(h1) {c->cd(9);  h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hdty1");     if(h1) {c->cd(10);  h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hdtx2");     if(h1) {c->cd(11);  h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hdty2");     if(h1) {c->cd(12);  h1->Draw(); h1=0;}
    
    h2 = (TH2F*)eOutputFile->Get("hxy_cp");    if(h2) {c->cd(13);  h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h2 = (TH2F*)eOutputFile->Get("htxy_cp");   if(h2) {c->cd(14); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h1 = (TH1F*)eOutputFile->Get("hchi");      if(h1) {c->cd(15); h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hchi20");    if(h1) {c->cd(16); h1->Draw(); h1=0;}

    crep1->Write("report");
    SafeDelete(crep1);
    gROOT->SetBatch(batch); 
  }
}

//---------------------------------------------------------------------
void EdbLinking::DoubletsFilterOut(TObjArray &p1, TObjArray &p2)
{
  EdbAlignmentV adup;
  adup.eDVsame[0]=adup.eDVsame[1]=5;
  adup.eDVsame[2]=adup.eDVsame[3]=0.003;
  
  adup.FillGuessCell(p1,p1,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(0);   // assign flag -10 to the duplicated segments

  adup.FillGuessCell(p2,p2,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(0);   // assign flag -10 to the duplicated segments
}