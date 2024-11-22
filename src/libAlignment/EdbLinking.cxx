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
#include "EdbSEQ.h"
#include "EdbTraceBack.h"

ClassImp(EdbLinking)

using namespace TMath;

//---------------------------------------------------------------------
EdbLinking::EdbLinking()
{
  eSegCouples.SetOwner();
}

//---------------------------------------------------------------------
void EdbLinking::GetPar(TEnv &env)
{
  eBinOK               = env.GetValue("fedra.link.BinOK"              , 6. );

  eNcorrMin            = env.GetValue("fedra.link.NcorrMin"           , 50 );

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
  eCPRankingAlg        = env.GetValue("fedra.link.CPRankingAlg"        , 0   );

  eDoSaveCouples       = env.GetValue("fedra.link.DoSaveCouples"       , true );
  
  eDoDumpDoubletsTree  = env.GetValue("fedra.link.DumpDoubletsTree"    , false );
  
  eCond.SetSigma0(       env.GetValue("fedra.link.Sigma0"         , "1 1 0.013 0.013") );
  eCond.SetPulsRamp0(    env.GetValue("fedra.link.PulsRamp0"      , "6 9") );
  eCond.SetPulsRamp04(   env.GetValue("fedra.link.PulsRamp04"     , "6 9") );
  eCond.SetDegrad(       env.GetValue("fedra.link.Degrad"         , 5) );
  eCond.DefineLLFunction(  env.GetValue("fedra.link.LLfunction"         , "0.256336-0.16489*x+2.11098*x*x") );

  GetDoubletsPar(env);
}

//---------------------------------------------------------------------
void EdbLinking::GetDoubletsPar(TEnv &env)
{
  const char *str=env.GetValue("fedra.link.RemoveDoublets"      , "");
  if(str) 
    if(sscanf(str,"%d %f %f %d",&(eRemoveDoublets.remove),&eRemoveDoublets.dr,&eRemoveDoublets.dt,&eRemoveDoublets.checkview) == 4) return;
  eRemoveDoublets.remove=0;
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
  if(p1.N()<=0) return;
  if(p2.N()<=0) return;
  Log(2,"EdbLinking::Link","segments z1 = %f  z2 = %f",
      p1.GetSegment(0)->Z(),p2.GetSegment(0)->Z());

  GetPar(env);
  eL1.Copy(l1,true);
  eL2.Copy(l2,true);

  if(area1<=0) area1 = EstimatePatternArea(p1);
  if(area2<=0) area2 = EstimatePatternArea(p2);
  if(area1<0.0000001||area2<0.0000001) return;
  EdbSEQ seq1;  seq1.eArea = area1;
  EdbSEQ seq2;  seq2.eArea = area2;
  GetPreselectionPar(seq1,env);
  GetPreselectionPar(seq2,env);
  seq1.eNsigma=eNsigmaEQshr;
  seq2.eNsigma=eNsigmaEQshr;
  if(gEDBDEBUGLEVEL>1) seq1.PrintLimits();
  
  if(eOutputFile) 
  {
    TH1F *hTall1 = seq1.ThetaPlot(p1   , "hTall1","Theta plot all, side 1 "); 
    hTall1->Write();
    TH1F *hTall2 = seq2.ThetaPlot(p2   , "hTall2","Theta plot all, side 2 "); 
    hTall2->Write();
  }
  TObjArray  p1pre(p1.N());  seq1.PreSelection(p1,p1pre);
  TObjArray  p2pre(p2.N());  seq2.PreSelection(p2,p2pre);
  Log(2,"EdbLinking::Link","after preselection: n1pre = %d  n2pre = %d", p1pre.GetEntries(),p2pre.GetEntries() );
  if(p1pre.GetEntries()<=0) return;
  if(p2pre.GetEntries()<=0) return;
  DoubletsFilterOut(p1pre,p2pre);
  
  TObjArray  p1shr,p2shr;
  
  if(eDoCorrectShrinkage || eDoCorrectAngles) {
    seq1.EqualizeMT(p1pre, p1shr, area1);
    seq2.EqualizeMT(p2pre, p2shr, area2);
    Log(2,"EdbLinking::Link","A n1shr = %d  n2shr = %d", p1shr.GetEntries(),p2shr.GetEntries() );
    if(eOutputFile)
    {
      TH1F *hTshr1 = seq1.ThetaPlot(p1shr, "hTshr1","Theta plot shr, side 1 "); hTshr1->Write();
      seq1.eHEq.DrawH1("eHEq1","eHEq1")->Write();
      TH1F *hTshr2 = seq2.ThetaPlot(p2shr, "hTshr2","Theta plot shr, side 2 "); hTshr2->Write();
      seq2.eHEq.DrawH1("eHEq2","eHEq2")->Write();
    }
    FillCombinationsAtMeanZ(p1shr,p2shr);
  }
  if(eDoCorrectShrinkage) {
    eCorr[0].SetV(5,eShr0);
    eCorr[1].SetV(5,eShr0);
    CorrectShrinkage( eDShr );
    CorrectShrinkage( eDShr*0.8 );
    eDoCorrectShrinkage = VerifyShrinkageCorr(0) * VerifyShrinkageCorr(1);
  }
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectShrinkage) {
    CorrectShrinkage(  eDShr*0.5 );
    CorrectShrinkage(  eDShr*0.5 );
    CorrectShrinkage(  eDShr*0.5 );
    if(eOutputFile) WriteShrinkagePlots();
  }
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectAngles)    CorrectAngles( p1shr,p2shr );
  if(eDoCorrectShrinkage || eDoCorrectAngles) {
    if(eOutputFile)
    {
      ePC[0].DrawH2("hxy_shr1","xy for the shrinkage corr sample side 1")->Write();
      ePC[1].DrawH2("hxy_shr2","xy for the shrinkage corr sample side 2")->Write();
      EdbH2 htxy1(50,-1,1,50,-1,1);
      FillThetaHist(0,htxy1);
      htxy1.DrawH2("htxy_shr1","txy plot for shr corr sample side 1");
      EdbH2 htxy2(50,-1,1,50,-1,1);
      FillThetaHist(0,htxy2);
      htxy2.DrawH2("htxy_shr2","txy plot for shr corr sample side 2");
    }
  }
  
  if(eOutputFile) 
  {
    eCorr[0].Write("corr1");
    eCorr[1].Write("corr2");
  }
  
  if(eDoFullLinking)   {
    TObjArray p1lnk(p1pre.GetEntriesFast());
    TObjArray p2lnk(p2pre.GetEntriesFast());
    if(eNsigmaEQlnk>0.1) {
      seq1.eNsigma = eNsigmaEQlnk;
      seq2.eNsigma = eNsigmaEQlnk;
      seq1.EqualizeMT(p1pre, p1lnk, area1);
      seq2.EqualizeMT(p2pre, p2lnk, area2);
    } else {
      p1lnk=p1pre;
      p2lnk=p2pre;
    }

    FullLinking(p1lnk,p2lnk);
    if(eOutputFile)
    {
      TH1F *hTlnk1 = seq1.ThetaPlot(p1lnk, "hTlnk1","Theta plot lnk, side 1 ");
      hTlnk1->Write();
      TH1F *hTlnk2 = seq2.ThetaPlot(p2lnk, "hTlnk2","Theta plot lnk, side 2 ");
      hTlnk2->Write();
      ePC[0].DrawH2("hxy_full1","xy full 1")->Write();
      ePC[1].DrawH2("hxy_full2","xy full 2")->Write();
    }
  }

  if(eOutputFile) ProduceReport();
  Corr2Aff(eCorr[0],eL1);
  Corr2Aff(eCorr[1],eL2);
  //eL1.Print();                     // layers with the corrections
  //eL2.Print();
}

//---------------------------------------------------------------------
bool EdbLinking::VerifyShrinkageCorr( int side )
{
  if( eCorr[side].V(5) > eShr0*(1+eDShr) ||  eCorr[side].V(5)< eShr0/(1+eDShr) ) {
      Log(1,"EdbLinking::Link","Shrinkage correction side %d out of range: %f [%f %f] reset shrinkage to default: %f",
             side, eCorr[side].V(5), eShr0/(1+eDShr), eShr0*(1+eDShr), eShr0 );
      eCorr[side].SetV(5,eShr0);
      return  false;
  }
  if(eNshr[side]<eNcorrMin) { 
     Log(1,"EdbLinking::Link","Shrinkage correction side %d not enough coins: %d < %d reset shrinkage to default: %f",
            side, eNshr[side],eNcorrMin,eShr0 );
     eCorr[side].SetV(5,eShr0);
     return false;
   }
  
  return true;
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
    if(nc<eNcorrMin) {Log(1,"EdbLinking::CorrectAngles","Warning: number of the selected segments too small: %d < %d do not correct angles",
         nc,eNcorrMin); return;}
    
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
    if(eDoSaveCouples) SaveCouplesTree();
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
    if(eDoSaveCouples) SaveCouplesTree();
}

//---------------------------------------------------------------------
void EdbLinking::SaveCouplesTree(const char *file)
{
  EdbCouplesTree ect;
  ect.InitCouplesTree("couples",file,"NEW");
  
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
  ect.eTree->SetAlias("x1c"  , "s1.eX+dz1*tx1c" );
  ect.eTree->SetAlias("y1c"  , "s1.eY+dz1*ty1c" );
  ect.eTree->SetAlias("x2c"  , "s2.eX+dz2*tx2c" );
  ect.eTree->SetAlias("y2c"  , "s2.eY+dz2*ty2c" );
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
    //s->SetID(i);             // basetrack id will be the tree entry number
    int entr  = ect.eTree->GetEntries();
    sc->eS->SetID(entr);             // basetrack id will be the tree entry number
    EdbTraceBack::SetBaseTrackVid( *(sc->eS), 0, 0, entr );   //TODO: plate, piece if available
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
  Log(2,"EdbLinking::SaveCouplesTree","%d couples are stored into %s",ntr, ect.GetFileName() );
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

  eDVsame[0] = eDVsame[1] = 50;  eDVsame[2] = eDVsame[3] = 0.03;
  InitHshr0(   25,  eCorr[0].V(5)/(1+dshr),   eCorr[0].V(5)*(1+dshr) );
  InitHshr1(   25,  eCorr[1].V(5)/(1+dshr),   eCorr[1].V(5)*(1+dshr) );
  
  eHdxyShr[0].InitH2(7, -35, 35, 7, -35, 35);
  eCorr[0].SetV(2,dz);  eCorr[1].SetV(2,   0);
  eNshr[0] = OptimiseVar1( 0, 5, &eHdxyShr[0] );              // variate shr1
  EdbSegCorr c0_opt = eCorr[0];

  eHdxyShr[1].InitH2(7, -35, 35, 7, -35, 35);
  eCorr[0]=c0;                                     // to remove the dependency on the correction order
  eCorr[0].SetV(2, 0);  eCorr[1].SetV(2, -dz);
  eNshr[1] = OptimiseVar1( 1, 5, &eHdxyShr[1]);               // variate shr2
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
    //seg.SetChi2( tf.Chi2ACP( *s1, *s2, cond1) );   //TODO test!!
    if(eCPRankingAlg==1) tf.Chi2ASegLL( seg1, seg2, seg, eCond, eCond);
    else tf.Chi2ASeg( seg1, seg2, seg, eCond, eCond);
    
    if(seg.Chi2() > eCHI2Pmax)  continue;
    
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

  eSegCouples.SetOwner();
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
 
    if(eDoCorrectShrinkage) {
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
    }
    if(eDoCorrectAngles) {
     ctit->AddText( Form("angular offsets: side1  %6.3f %6.3f   side2 %6.3f %6.3f",
                         eCorr[0].V(3),eCorr[0].V(4), eCorr[1].V(3), eCorr[1].V(4) ));
    }
    ctit->Draw();

    TPad *c = new TPad("c","plots",0.01,0.03,0.99,0.94);
    c->Divide(4,4);    c->Draw();

    float densAll=0,densShr=0, densLnk=0;

    h2 = (TH2F*)eOutputFile->Get("hxy_shr1");
    if(!h2) h2 = (TH2F*)eOutputFile->Get("hxy_full1");  
    if(h2) {c->cd(1); h2->SetStats(0); h2->Draw("colz"); h2=0;}

    h2 = (TH2F*)eOutputFile->Get("hxy_shr2");
    if(!h2) h2 = (TH2F*)eOutputFile->Get("hxy_full2");
    if(h2) {c->cd(2); h2->SetStats(0); h2->Draw("colz"); h2=0;}

    h1 = (TH1F*)eOutputFile->Get("hTall1");      if(h1) {c->cd(3);  h1->Draw(); densAll=h1->Integral();h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hTshr1");      if(h1) {c->cd(3);  h1->Draw("same"); densShr=h1->Integral(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hTlnk1");      if(h1) {c->cd(3);  h1->Draw("same"); densLnk=h1->Integral(); h1=0;}
    c->cd(3);
    TPaveText *lable1 = new TPaveText(0.2,0.6,0.6,0.8,"NDC");
    lable1->AddText(Form("all: %7.0f",densAll));
    lable1->AddText(Form("shr: %7.0f",densShr));
    lable1->AddText(Form("lnk: %7.0f",densLnk));
    lable1->Draw();
    h1 = (TH1F*)eOutputFile->Get("hTall2");      if(h1) {c->cd(4);  h1->Draw(); densAll=h1->Integral(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hTshr2");      if(h1) {c->cd(4);  h1->Draw("same"); densShr=h1->Integral(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hTlnk2");      if(h1) {c->cd(4);  h1->Draw("same"); densLnk=h1->Integral(); h1=0;}
    c->cd(4);
    TPaveText *lable2 = new TPaveText(0.2,0.6,0.6,0.8,"NDC");
    lable2->AddText(Form("all: %7.0f",densAll));
    lable2->AddText(Form("shr: %7.0f",densShr));
    lable2->AddText(Form("lnk: %7.0f",densLnk));
    lable2->Draw();

    h2 = (TH2F*)eOutputFile->Get("hdxy_shr");  if(h2) {c->cd(5);  h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h2 = (TH2F*)eOutputFile->Get("htxy_shr");  if(h2) {c->cd(6);  h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h1 = (TH1F*)eOutputFile->Get("shr1");      if(h1) {c->cd(7);  h1->Draw("hist"); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("shr2");      if(h1) {c->cd(8);  h1->Draw("hist"); h1=0;}

    h1 = (TH1F*)eOutputFile->Get("hdtx1");     if(h1) {c->cd(9);  h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hdty1");     if(h1) {c->cd(10); h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hdtx2");     if(h1) {c->cd(11); h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hdty2");     if(h1) {c->cd(12); h1->Draw(); h1=0;}

    h2 = (TH2F*)eOutputFile->Get("hxy_cp");    if(h2) {c->cd(13); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h2 = (TH2F*)eOutputFile->Get("htxy_cp");   if(h2) {c->cd(14); h2->SetStats(0); h2->Draw("colz"); h2=0;}
    h1 = (TH1F*)eOutputFile->Get("hchi");      if(h1) {c->cd(15); h1->Draw(); h1=0;}
    h1 = (TH1F*)eOutputFile->Get("hchi20");    if(h1) {c->cd(16); h1->Draw(); h1=0;}

    crep1->Write("report");
    SafeDelete(crep1);
    SafeDelete(ctit);
    gROOT->SetBatch(batch);
  }
}

//---------------------------------------------------------------------
void EdbLinking::DoubletsFilterOut(TObjArray &p1, TObjArray &p2, bool fillhist)
{
  TH2F *hxy1=0, *htxty1=0, *hxy2=0, *htxty2=0;
  if(fillhist)  {
    hxy1   = new TH2F("dblXY1","Doublets DX DY",50,-eRemoveDoublets.dr,eRemoveDoublets.dr,50,-eRemoveDoublets.dr,eRemoveDoublets.dr);
    htxty1 = new TH2F("dblTXTY1","Doublets DTX DTY",50,-eRemoveDoublets.dt,eRemoveDoublets.dt,50,-eRemoveDoublets.dt,eRemoveDoublets.dt);
    hxy2   = new TH2F("dblXY2","Doublets DX DY",50,-eRemoveDoublets.dr,eRemoveDoublets.dr,50,-eRemoveDoublets.dr,eRemoveDoublets.dr);
    htxty2 = new TH2F("dblTXTY2","Doublets DTX DTY",50,-eRemoveDoublets.dt,eRemoveDoublets.dt,50,-eRemoveDoublets.dt,eRemoveDoublets.dt);
  }
  EdbAlignmentV adup;
  adup.eDVsame[0]=adup.eDVsame[1]= eRemoveDoublets.dr;
  adup.eDVsame[2]=adup.eDVsame[3]= eRemoveDoublets.dt;
  
  adup.FillGuessCell(p1,p1,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(eRemoveDoublets.checkview, hxy1, htxty1);   // assign flag -10 to the duplicated segments
  if(eDoDumpDoubletsTree) DumpDoubletsTree(adup,"doublets1");

  adup.FillGuessCell(p2,p2,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(eRemoveDoublets.checkview, hxy2, htxty2);   // assign flag -10 to the duplicated segments
  if(eDoDumpDoubletsTree) DumpDoubletsTree(adup,"doublets2");
  
  if(hxy1)   hxy1->Write();
  if(hxy2)   hxy2->Write();
  if(htxty1) htxty1->Write();
  if(htxty2) htxty2->Write();
  SafeDelete(hxy1);
  SafeDelete(hxy2);
  SafeDelete(htxty1);
  SafeDelete(htxty2);
}

//---------------------------------------------------------------------
void EdbLinking::DumpDoubletsTree(EdbAlignmentV &adup, const char *name)
{
  EdbCouplesTree dup;
  dup.InitCouplesTree(name, 0 ,"RECREATE");
  int n=adup.CheckEqualArr(adup.eS[0],adup.eS[1]);
  for(int i=0; i<n; i++)
  {
    EdbSegP *s0 = (EdbSegP*)(adup.eS[0].At(i));
    EdbSegP *s1 = (EdbSegP*)(adup.eS[1].At(i));
    if( s0->Flag() ==-10)
      dup.Fill( s0, s1);
  }
  dup.WriteTree();
}

///______________________________________________________________________________
void EdbLinking::CloneCouplesTree( const char *ifile, const char *ofile, EdbAffine2D *aff, TCut *cut )
{
  //read tree using eCut, apply affine transformation if any and write tree
  EdbPattern pat, p1, p2;

  EdbCouplesTree ict;
  if(cut) ict.eCut = (*cut);
  ict.InitCouplesTree("couples",ifile,"READ");
  int n = ict.GetCPData(eSegCouples);
  ict.Close();

  if(aff) {
    for(int i=0; i<n; i++) {
      EdbSegCouple *sc = (EdbSegCouple*)(eSegCouples.UncheckedAt(i));
      if(sc->eS) sc->eS->Transform(aff);
      if(sc->eS) sc->eS1->Transform(aff);
      if(sc->eS) sc->eS2->Transform(aff);
    }
  }

  SaveCouplesTree(ofile);
}

//---------------------------------------------------------------------
void EdbLinking::GetFillPattern(EdbPattern &p, int layer)
{
  int ncp = eSegCouples.GetEntries();
  for(int i=0; i<ncp; i++) {
    EdbSegCouple *sc = (EdbSegCouple*)(eSegCouples.At(i));
    if     (layer==0) p.AddSegment(*(sc->eS));
    else if(layer==1) p.AddSegment(*(sc->eS1));
    else if(layer==2) p.AddSegment(*(sc->eS2));
  }
  Log(2,"EdbLinking::GetPattern","for layer %d with %d segments",layer,ncp);
}
