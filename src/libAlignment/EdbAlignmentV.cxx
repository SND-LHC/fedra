//-- Author :  Valeri Tioukov   19-06-2009
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// EdbAlignmentV                                                          //
//                                                                        //
// General alignment of the preselected arrays (combinations) of segments //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector2.h"
#include "EdbLog.h"
#include "EdbAlignmentV.h"
#include "EdbLayer.h"

ClassImp(EdbAlignmentV)

using namespace TMath;

//---------------------------------------------------------------------
EdbAlignmentV::EdbAlignmentV()
{
  eXmarg = eYmarg = 0.0001;        // margins for the cell definition
  eDoubletsRate=0;
  eUseAffCorr = false;   // by default use eCorr (separated corrections)
  eOutputFile=0;
}

//---------------------------------------------------------------------
EdbAlignmentV::~EdbAlignmentV()
{
  eS[0].Clear();
  eS[1].Clear();
  SafeDelete(eDoubletsRate);
}

//---------------------------------------------------------------------
void EdbAlignmentV::InitOutputFile(const char *file, const char *option)
{
  CloseOutputFile();
  eOutputFile = TFile::Open(file,option);
  if (!eOutputFile || !eOutputFile->IsOpen() || eOutputFile->IsZombie()) 
    Log(1, "EdbAlignmentV::InitOutputFile", "Failed open file: %s", file);
}

//---------------------------------------------------------------------
void EdbAlignmentV::CloseOutputFile()
{
  if(eOutputFile) {
    eOutputFile->Close();
    SafeDelete(eOutputFile);
  }
}

//---------------------------------------------------------------------
Bool_t EdbAlignmentV::IsInsideDVsame(EdbSegP &s1, EdbSegP &s2)
{
  if( Abs(s1.X()-s2.X())   > eDVsame[0] )  return 0;
  if( Abs(s1.Y()-s2.Y())   > eDVsame[1] )  return 0;
  if( Abs(s1.TX()-s2.TX()) > eDVsame[2] )  return 0;
  if( Abs(s1.TY()-s2.TY()) > eDVsame[3] )  return 0;
  return 1;
}

//---------------------------------------------------------------------
int EdbAlignmentV::DoubletsFilterOut(int checkview, TH2F *hxy, TH2F *htxty )
{
  // assumed the same pattern's segments in eS[0] and eS[1]
  // checkview = 0 do not check if close tracks belongs to different views
  // checkview = 1 check views, select best and cancel the worst candidate
  // checkview = 2 check views but do not cancell worst candidate (debug option)
  
  int nout=0;
  int n = CheckEqualArr(eS[0], eS[1]);
  for(int i=0; i<n; i++) {
    EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);
    if( s1->Flag()==-10 && s2->Flag()==-10) continue;
    if(checkview>0) if( s2->Aid(0)==s1->Aid(0) && s2->Aid(1)==s1->Aid(1) && s2->Side()==s1->Side() )    continue;
    if( !IsInsideDVsame(*s1,*s2) )                                                                      continue;
    if(s1->Flag()>-10 && s2->Flag()>-10) {
      if(hxy)   hxy->Fill( s1->X()-s2->X() , s1->Y()-s2->Y());
      if(htxty) htxty->Fill(s1->TX()-s2->TX(),s1->TY()-s2->TY());
    }
    if(checkview==1||checkview==0) {
      if( s2->W()>s1->W() )                                        s1->SetFlag(-10);
      else if( s2->W() == s1->W() &&  s2->Chi2() <  s1->Chi2() )   s1->SetFlag(-10);
      else if( s2->W() == s1->W() &&  s2->Chi2() >= s1->Chi2() )   s2->SetFlag(-10);
      else if( s2->W()<s1->W() )                                   s2->SetFlag(-10);
      nout++;
    }
  }
  
  //check 
  int miss=0, wrong=0;
  for(int i=0; i<n; i++) {
    EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);
    if( s1->Flag()==-10 && s2->Flag()==-10 )                  continue;
    if( s1->Flag()!=-10 && s2->Flag()!=-10 )        { miss++; continue;}
    if( s1->Flag()==-10 && (s1->W() > s2->W()) )    {
      wrong++;
      s1->SetFlag(0);
      s2->SetFlag(-10);
    }
  }
  Log(2,"DubletsFilterOut","miss: %d wrong: %d",miss,wrong);
  
  Log(2,"DubletsFilterOut","%d segments discarded with DX,DY,DTX,DTY: (%8.4f %8.4f %7.5f %7.5f) checkview =%d", 
      nout,eDVsame[0],eDVsame[1],eDVsame[2],eDVsame[3], checkview );
  return nout;
}

//---------------------------------------------------------------------
int EdbAlignmentV::SelectIsolated()
{
  // assume the different pattern's segments in eS[0] and eS[1]
  // assume that the selection method in FillCombinations - "first s1 then s2"
  // remove all duplicated combs

  int nrem=0;
  int n = CheckEqualArr(eS[0], eS[1]);
  EdbSegP *s1=0,*s2=0;
  EdbSegP *s1ok=0,*s2ok=0;
  int   iok=-1;
  for(int i=0; i<n; i++) {
    s1 = (EdbSegP*)eS[0].UncheckedAt(i);
    s2 = (EdbSegP*)eS[1].UncheckedAt(i);
    if(s1==s1ok) {
      eS[0].RemoveAt(i);
      eS[1].RemoveAt(i);
      nrem++;
    } 
    else if(s1ok!=0) {
      eS[0].RemoveAt(iok);
      eS[1].RemoveAt(iok);
      s1ok=0;
      nrem++;
    }
    else {  s1ok = s1, s2ok = s2; iok=i; }
  }
  
  eS[0].Compress();
  eS[1].Compress();

  Log(2,"SelectIsolated","%d couples discarded with", nrem );
  return nrem;
}

//---------------------------------------------------------------------
int EdbAlignmentV::SelectBestCouple()
{
  // assume the different pattern's segments in eS[0] and eS[1]
  // assume that the selection method in FillCombinations - "first s1 then s2"

  int nrem=0;
  int n = CheckEqualArr(eS[0], eS[1]);
  EdbSegP *s1=0,*s2=0;
  EdbSegP *s1ok=0,*s2ok=0;
  int   iok=-1;
  for(int i=0; i<n; i++) {
    s1 = (EdbSegP*)eS[0].UncheckedAt(i);
    s2 = (EdbSegP*)eS[1].UncheckedAt(i);
    if(s1==s1ok) {
      if( CoupleQuality(*s1,*s2) >  CoupleQuality(*s1ok,*s2ok) ) {
	s1ok = s1, s2ok = s2;
	eS[0].RemoveAt(iok);
	eS[1].RemoveAt(iok);
	iok=i;
      }
      else { 
	eS[0].RemoveAt(i);
	eS[1].RemoveAt(i);
      }
      nrem++;
    } 
    else {  s1ok = s1, s2ok = s2; iok=i; }
  }
  
  eS[0].Compress();
  eS[1].Compress();

  Log(2,"SelectBestCouple","%d couples discarded with", nrem );
  return nrem;
}

//---------------------------------------------------------------------
float EdbAlignmentV::CoupleQuality( EdbSegP &s1, EdbSegP &s2 )
{
  // TODO
  return s1.W()+s2.W();
}

//---------------------------------------------------------------------
void EdbAlignmentV::FillThetaHist(int side, EdbH2 &htxy)
{
  if( side<0 || side>1 ) return;
  TObjArray arr(10000);
  int n = ePC[side].SelectObjects(arr);
  for(int i=0; i<n; i++) {
    EdbSegP *s = (EdbSegP*)arr.UncheckedAt(i);
    htxy.Fill( s->TY(), s->TX() );
  }
}

//---------------------------------------------------------------------
void EdbAlignmentV::AssignFlagToS(int flag)
{
  int n = CheckEqualArr(eS[0], eS[1]);
  for(int i=0; i<n; i++) {
    ((EdbSegP*)eS[0].UncheckedAt(i))->SetFlag(flag);
    ((EdbSegP*)eS[1].UncheckedAt(i))->SetFlag(flag);
 }
}

//---------------------------------------------------------------------
int EdbAlignmentV::FillCombinations()
{
  float dxMax = 2*Max( ePC[0].Xbin(), ePC[1].Xbin() );
  float dyMax = 2*Max( ePC[0].Ybin(), ePC[1].Ybin() );
  return FillCombinations( eDVsame, dxMax, dyMax, 1);
}

//---------------------------------------------------------------------
int EdbAlignmentV::FillCombinations(float dv[4], float dxMax, float dyMax, bool doFill)
{
  // input:  dv = coinsidence condition (dx,dy,dtx,dty)
  //         dxMax, dyMax - the area around each segment for the combinations 
  //         selection, it can be selected bigger then coinsidence condition due to 
  //         corrections displacements
  // 
  // the cells must be already filled

  eS[0].Clear();
  eS[1].Clear();
  int ir2[2] = { int((dxMax-0.0001)/ePC[0].Xbin())+1, int((dyMax-0.0001)/ePC[0].Ybin())+1 };

  int nout=0;

  TObjArray arr1(10000);
  TObjArray arr2(100);
  float v[2];
  EdbSegP *s1,*s2;

  int ncomb;                      // combinations rate counter
  int n1 = ePC[0].SelectObjects(arr1);
  for(int i=0; i<n1; i++) {
    ncomb=0;
    s1 = (EdbSegP*)arr1.UncheckedAt(i);
    if(s1->Flag()==-10)                                   continue;
    arr2.Clear();
    v[0] = X( 0, *s1);
    v[1] = Y( 0, *s1);
    int n2 = ePC[1].SelectObjectsC(v,ir2,arr2);

    if(n2<1)                                              continue;
    float tx1 = TX(0, *s1);
    float ty1 = TY(0, *s1);
    for(int i2=0; i2<n2; i2++) {
      s2 = (EdbSegP*)arr2.UncheckedAt(i2);
      if(s2->Flag()==-10)                                 continue;
      if( s2==s1 )                                        continue;
      if( Abs(X(1, *s2) - v[0]) > dv[0] )                 continue;
      if( Abs(Y(1, *s2) - v[1]) > dv[1] )                 continue;
      if( Abs(TX(1, *s2)- tx1)  > dv[2] )                 continue;
      if( Abs(TY(1, *s2)- ty1)  > dv[3] )                 continue;

      if(doFill) {
        eS[0].Add(s1);
        eS[1].Add(s2);
      }
      ncomb++;
    }
    nout+=ncomb;
    if(eDoubletsRate) eDoubletsRate->Fill(ncomb);
  }

  arr1.Clear();
  arr2.Clear();

  Log(3,"FillCombinations","%d selected with the acceptance: %7.2f %7.2f (%d,%d) and tolerance (%7.2f %7.2f %7.4f %7.4f)",
      nout, dxMax,dyMax, ir2[0], ir2[1], dv[0], dv[1], dv[2], dv[3] );
  return nout;
}

//---------------------------------------------------------------------
int EdbAlignmentV::OptimiseVar1(int side, int ivar,  EdbH2 *hdxy, EdbH2 *hdtxy)
{
  // use already preselected "couples" - useful to test small variations of parameters
  if(!SideOK(side)) return 0;

  EdbH1 &h = eH[side][ivar];
  if( h.N() < 1 ) {
    Log(1,"OptimiseVar1","ERROR: var %d of side %d have %d steps - skipped", ivar,side, h.N() );
    return 0;
  }

  h.CleanCells();

  int   npk0=0;
  float dx0=0,dy0=0;  //the position of the highest peak in xy
  float var0=0;      //the position of the highest peak on var
  for(int i=0; i<h.N(); i++) 
    {
      float var= h.X(i);
      eCorr[side].SetV( ivar, var );

      int npk = Ncoins(eDVsame, hdxy);

      float dx=0, dy=0;
      if(hdxy){
	EdbPeak2 pk2(*hdxy);
	npk = (int)(pk2.ProbPeak( dx, dy ));
        //hdxy->DrawH2(Form("hdxy%d_%2.2d",side,i),Form("OptimiseVar1: ivar=%d var=%f npk=%d",ivar, var,npk) )->Write();
      }
      h.SetBin( i, npk );
      if( npk > npk0 ) { npk0=npk; dx0 = dx; dy0=dy; var0 = var; }
      Log(4,"OptimiseVar1"," %d %f %d (%f %f)", i, h.X(i), h.Bin(i), dx,dy );
    }

  //eCorr[side].AddV(0,dx0);      // set x-y corrections if any
  //eCorr[side].AddV(1,dy0);
  eCorr[side].SetV(ivar,var0);
  if(hdxy) Ncoins(eDVsame, hdxy);  // to keep the optimal histogram in hdxy
  
  Log(3,"OptimiseVar1"," side:var (%d:%d): found %d at var=%f  xy:(%f %f)", side, ivar, npk0, var0, dx0,dy0 );
  return npk0;
}


//---------------------------------------------------------------------
void EdbAlignmentV::OptimiseVar2( int side1, int ivar1, int side2, int ivar2, 
				  EdbH2 &h12,  EdbH2 *hdxy, EdbH2 *hdtxy )
{
  // use already preselected "couples" - useful to test small variations of parameters
  // especially after Hdistance preselection
  // Input:    2 variables defoned as side1/ivar1 and side1/ivar1
  // Output:   h12 - filled selections histo
  // Optional: hdxy, hdtxy - 2d- histograms defined outside and used for the peak selection
  if(!SideOK(side1)) return;
  if(!SideOK(side2)) return;
  EdbH1 &h1 = eH[side1][ivar1];
  EdbH1 &h2 = eH[side2][ivar2];
  if( h1.N() < 1 || h2.N() < 1) {
    Log(1,"OptimiseVar2","ERROR: var %d of side %d have %d steps - skipped", ivar1,side1, h1.N() );
    Log(1,"OptimiseVar2","ERROR: var %d of side %d have %d steps - skipped", ivar2,side2, h2.N() );
    return;
  }
  Log(3,"OptimiseVar2","side/var: %d/%d and %d/%d with  %d x %d = %d attempts", 
      side1, ivar1, side2, ivar2, h1.N(), h2.N(), h1.N()*h2.N() );

    //h12.Delete();
  int    n12[2]  = { h1.N(), h2.N() };
  float  min[2]  = { h1.Xmin(), h2.Xmin() };
  float  max[2]  = { h1.Xmax(), h2.Xmax() };
  h12.InitH2(n12, min, max);

  h1.CleanCells();
  h2.CleanCells();

  int   npk0=0;
  float dx0=0,dy0=0;  //the position of the highest peak

  for(int i1=0; i1<h1.N(); i1++)
    {
      eCorr[side1].SetV( ivar1, h1.X(i1) );
      for(int i2=0; i2<h2.N(); i2++) 
	{
	  eCorr[side2].SetV( ivar2, h2.X(i2) );
	  int nc = Ncoins(eDVsame, hdxy);
	  if(!hdxy) h12.SetBin( i1, i2, nc );
	  else {
	    EdbPeak2 pk2(*hdxy);
	    float  dx, dy;
	    int npk = pk2.FindPeak( dx, dy );
	    //int npk = pk2.FindPeak9( dx, dy );
	    if( npk > npk0 ) { npk0=npk; dx0 = dx; dy0=dy; }
	    h12.SetBin( i1, i2, npk );
	  }
	}
    }

  eCorr[side1].AddV(0,dx0);      // set x-y corrections
  eCorr[side1].AddV(1,dy0);
}

//---------------------------------------------------------------------
bool EdbAlignmentV::SideOK(int side)
{
  if(side==0||side==1) return 1; 
  else { 
    Log(1,"SideOK","ERROR: invalid side index %d",side); 
    return 0;
  }
}

//---------------------------------------------------------------------
Int_t EdbAlignmentV::CheckEqualArr(TObjArray &arr1, TObjArray &arr2)
{
  int n1 = arr1.GetEntries();
  int n2 = arr2.GetEntries();
  int n = Min(n1,n2);
  if(n1!=n2) Log(1,"CheckEqualArr","WARNING: arrays are not equal: %d %d, take first %d segs", n1,n2,n );
  return n;
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::FineCorrZ(TObjArray &arr1, TObjArray &arr2)
{
  double dzz=0, dzx=0, dzy=0;
  int n = CheckEqualArr(arr1,arr2);
  int icx=0, icy=0;
  for(int i=0; i<n; i++) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    float tx1 =  TX(0, *s1), ty1 =  TY(0, *s1);
    float x1  =  X(0,*s1),  y1 =  Y(0,*s1);
    float x2  =  X(1,*s2),  y2 =  Y(1,*s2);

    if(tx1>0.1) {
      dzx += (x2-x1)/tx1;
      icx++;
    }
    if(ty1>0.1) {
      dzy += (y2-y1)/ty1;
      icy++;
    }
  }

  if(icx+icy>3)   dzz = (dzx+dzy)/(icx+icy);

  Log(3,"FineCorrZ","dzz =%f with icx=%d and icy=%d (of %d) segments",
      (float)dzz, icx,icy, n );

  return (float)dzz;
}

//---------------------------------------------------------------------
void EdbAlignmentV::CorrToCoG(int side, EdbPattern &p)
{
  int n=p.N();
  Double_t x0=0, y0=0;
  for(int i=0; i<n; i++) {
    EdbSegP *s = p.GetSegment(i);
    x0 += X(side, *s);
    y0 += Y(side, *s);
  }
  x0/=n;   y0/=n;
  eCorr[0].AddV(0, -x0);
  eCorr[0].AddV(1, -y0);
  eCorr[1].AddV(0, -x0);
  eCorr[1].AddV(1, -y0);
}

//---------------------------------------------------------------------
void EdbAlignmentV::CorrToCoG(int side, TObjArray &p)
{
  int n=p.GetEntries();
  Double_t x0=0, y0=0;
  for(int i=0; i<n; i++) {
    EdbSegP *s = (EdbSegP*)p.UncheckedAt(i);
    x0 += X(side,*s);
    y0 += Y(side,*s);
  }
  x0/=n;   y0/=n;
  eCorr[0].AddV(0, -x0);
  eCorr[0].AddV(1, -y0);
  eCorr[1].AddV(0, -x0);
  eCorr[1].AddV(1, -y0);
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::FineCorrPhi(TObjArray &arr1, TObjArray &arr2)
{
  // for the correct answer the O,O should be in the CoG of the second pattern
  int n = CheckEqualArr(arr1,arr2);
  if(n<3) return 0;
  double x1=0, y1=0, x2=0,y2=0;
  Double_t R=0, r=0, dPhi=0;
  for(int i=0; i<n; i++) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1 =  X(0,*s1);
    y1 =  Y(0,*s1);
    x2 =  X(1,*s2);
    y2 =  Y(1,*s2);
    r = Sqrt(x1*x1+y1*y1);
    dPhi += (ATan2(-y2,-x2) - ATan2(-y1,-x1)) * r;
    R += r;
  }
  dPhi /= R;
  Log(3,"FineCorrPhi","dPhi =%f with %d segments",      (float)dPhi, n );
  return (float)dPhi;
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::CalcMeanShr(float tmin, float tmax)
{
  //  Calculate mean shrinkage (to be applied to side 0 to get side 1)
  int n = CheckEqualArr(eS[0],eS[1]);
  int ic=0;
  Double_t shr=0;
  for(int i=0; i<n; i++) {
    float t1  = eCorr[0].T(*((EdbSegP*)eS[0].UncheckedAt(i)));
    //float t1  = ((EdbSegP*)eS[0].UncheckedAt(i))->Theta();
    if(t1<tmin)   continue;
    if(t1>tmax)   continue;
    float t2  = eCorr[1].T(*((EdbSegP*)eS[1].UncheckedAt(i)));
    //float t2  = ((EdbSegP*)eS[1].UncheckedAt(i))->Theta();
    shr   += (t1/t2);                // the angle is divided by shr - that's why it is inverse
    ic++;
  }
  shr /= ic;
  return shr;
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::CalcMeanDZ(float tmin, float tmax)
{
  //  Calculate mean shrinkage (to be applied to side 0 to get side 1)
  int n = CheckEqualArr(eS[0],eS[1]);
  int ic=0;
  float dz = eCorr[0].V(2); if(Abs(dz)<0.000001) return 0;
  eCorr[0].SetV(2,0.);
  Double_t dzn=0;
  for(int i=0; i<n; i++) {
    EdbSegP *s1=(EdbSegP*)eS[0].UncheckedAt(i);
    float t1  = eCorr[0].T(*s1);
    if(t1<tmin)   continue;
    if(t1>tmax)   continue;
    EdbSegP *s2=(EdbSegP*)eS[1].UncheckedAt(i);
    float tx = ( X(1,*s2) - X(0,*s1) ) / dz;
    float ty = ( Y(1,*s2) - Y(0,*s1) ) / dz;
    float t = TMath::Sqrt(tx*tx+ty*ty);
    dzn += dz*t/t1;
    ic++;
  }
  dzn /= ic;
  eCorr[0].SetV(2,dz);
  Log(3,"CalcMeanDZ","dz = %f dzn = %f n=%d",dz,dzn,ic);
  return dzn;
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::CalcMeanDiff(int ivar)
{
  //  Calculate mean difference for the given variable taking into account corrections
  int n = CheckEqualArr(eS[0],eS[1]);
  int ic=0;
  Double_t sd=0;
  for(int i=0; i<n; i++) {
    sd   += ( Var( 1, i, ivar ) - Var( 0, i, ivar) );
    ic++;
  }
  sd /= ic;
  return sd;
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::CalcMeanDiff2Const(int ivar, int side, float mean)
{
  //  Calculate mean difference to the the given constant taking into account corrections
  int n = Ncp();
  int ic=0;
  Double_t sd=0;
  for(int i=0; i<n; i++) {
    sd   += ( Var( side, i, ivar ) - mean );
    ic++;
  }
  sd /= ic;
  return sd;
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::CalcFractMeanDiff(int ivar, float fraction)
{
  //  Calculate mean difference for the highest density region
  int n = CheckEqualArr(eS[0],eS[1]);
  TArrayF arr(n);
  for(int i=0; i<n; i++) {
    arr[i] = ( Var( 1, i, ivar ) - Var( 0, i, ivar) );
  }
  return FindDensityPeak( arr, fraction);
}

//---------------------------------------------------------------------
Float_t EdbAlignmentV::FindDensityPeak(TArrayF &arr, float fraction)
{
  if(fraction<=0||fraction>1) Log(1,"EdbAlignmentV::FindDensityPeak","ERROR! fraction out of limits: %f",fraction);
  TArrayI ind(arr.fN);
  Sort( arr.fN, arr.GetArray(), ind.GetArray(), 0 );
  int nfra = (int)(fraction*arr.fN);
  float dxmin = arr[ ind[nfra] ] - arr[ ind[0] ];
  int   ibest=0;
 for(int i=1; i<arr.fN-nfra; i++)   {
    float dx = arr[ ind[nfra+i] ] - arr[ind[i]];
    if( dx < dxmin) {
      dxmin = dx;
      ibest=i;
    }
  }

  Double_t sd=0;
  int ic=0;
  for(int i=ibest; i<ibest+nfra-1; i++) {
    sd   += arr[ ind[i] ];
    ic++;
  }
  sd /= ic;
  Log(2,"EdbAlignmentV::FindDensityPeak", "nfra/n: %d/%d     dxmin = %g   sd = %g",nfra,arr.fN,dxmin,sd);
  return (Float_t)sd;
}

//---------------------------------------------------------------------
Int_t EdbAlignmentV::CalcApplyFractMeanDiff()
{
  // by default apply it to side 0
  int side=0;
  int n = CheckEqualArr(eS[0],eS[1]);

//  if(n<nlim) return 0;           //TODO
  
  eCorr[side].AddV( 3, CalcFractMeanDiff( 2, 0.4 ) );                    // calc and  apply angles
  eCorr[side].AddV( 4, CalcFractMeanDiff( 3, 0.4 ) );
  
  //eCorr[side].SetV( 2, CalcMeanDZ()  );                      // calc and apply DZ
  
  //eCorr[side].SetV( 5, eCorr[side].V(5) * CalcMeanShr() );   // calc and apply shrinkage
  
  
  eCorr[side].AddV( 0, CalcFractMeanDiff( 0, 0.4 ) );                    // then apply coord
  eCorr[side].AddV( 1, CalcFractMeanDiff( 1, 0.4 ) );
  
  return n;
}

//---------------------------------------------------------------------
Int_t EdbAlignmentV::CalcAffFull()
{
  eUseAffCorr=true;
  eCorrL[0].SetZcorr( FineCorrZ( eS[0], eS[1] ) );
  EdbAffine2D aff1to2XY;
  CalculateAffXY( eS[0], eS[1], aff1to2XY);
  eCorrL[0].GetAffineXY()->Transform(&aff1to2XY);
  EdbAffine2D aff2to1TXTY;
  CalculateAffTXTY( eS[1], eS[0], aff2to1TXTY);
  eCorrL[1].GetAffineTXTY()->Transform(&aff2to1TXTY);
}

/*
//---------------------------------------------------------------------
Int_t EdbAlignmentV::CalcApplyMeanDiffCoord()
{
  // by default apply it to side 0
  int side=0;
  int n = CheckEqualArr(eS[0],eS[1]);

//  if(n<nlim) return 0;           //TODO
  eCorr[side].SetV( 2, CalcMeanDZ()  );                      // at first calc and apply DZ
  eCorr[side].AddV( 0, CalcMeanDiff(0) );                    // then apply coord
  eCorr[side].AddV( 1, CalcMeanDiff(1) );
  return n;
}
*/
//---------------------------------------------------------------------
Int_t EdbAlignmentV::CalcApplyMeanDiff()
{
  // by default apply it to side 0
  int side=0;
  int n = CheckEqualArr(eS[0],eS[1]);

//  if(n<nlim) return 0;           //TODO
  
  eCorr[side].AddV( 3, CalcMeanDiff(2) );                    // then apply angles
  eCorr[side].AddV( 4, CalcMeanDiff(3) );
  
  eCorr[side].SetV( 2, CalcMeanDZ()  );                      // at first calc and apply DZ
  
  eCorr[side].SetV( 5, eCorr[side].V(5) * CalcMeanShr() );   // second calc and apply shrinkage
  
  
  eCorr[side].AddV( 0, CalcMeanDiff(0) );                    // then apply coord
  eCorr[side].AddV( 1, CalcMeanDiff(1) );
  
  return n;
}

//---------------------------------------------------------------------
Int_t EdbAlignmentV::FindCorrDiff( float dvsame[4], int side, int nlim )
{
  //find diffs with the default settings and correct requested side
  if(!SideOK(side)) return 0;
  float dvf[4];
  int n = FindDiff(eS[0], eS[1], dvsame, dvf );
  if(n<nlim) return n;
  eCorr[side].AddV(0,dvf[0]);
  eCorr[side].AddV(1,dvf[1]);
  //eCorr[side].AddV(3,dvf[2]);
  //eCorr[side].AddV(4,dvf[3]);
  Log(3,"FindCorrDiff","correct side %d  with %11.7f %11.7f %10.7f %10.7f using %5d coinsidences",
      side, dvf[0] ,dvf[1] ,dvf[2] ,dvf[3],n);
  return n;
}

//---------------------------------------------------------------------
Int_t EdbAlignmentV::FindDiff(TObjArray &arr1, TObjArray &arr2, float dvlim[4], float dvfound[4] )
{
  int n = CheckEqualArr(arr1,arr2);
  int ic=0;
  EdbSegP *s1,*s2;
  Double_t sdx=0, sdy=0, sdtx=0, sdty=0;
  float dvf[4];
  for(int i=0; i<n; i++) {
    s1 = (EdbSegP*)arr1.UncheckedAt(i);
    s2 = (EdbSegP*)arr2.UncheckedAt(i);
    if( ValidCoinsidence(*s1, *s2, dvlim, dvf) ) {
      sdx   += dvf[0];
      sdy   += dvf[1];
      sdtx  += dvf[2];
      sdty  += dvf[3];
      ic++;
    }
  }
  dvfound[0]  = sdx/ic;
  dvfound[1]  = sdy/ic;
  dvfound[2]  = sdtx/ic;
  dvfound[3]  = sdty/ic;
  return ic;
}


//---------------------------------------------------------------------
Bool_t EdbAlignmentV::ValidCoinsidence(EdbSegP &s1, EdbSegP &s2, float dv[4], float dfound[4] )
{
  if(&s1==&s2)                                                                return 0;
  if(s1.Flag()==-10)                                                          return 0;
  if(s2.Flag()==-10)                                                          return 0;
  dfound[0] = X(1,s2)  - X(0,s1);  if(Abs(dfound[0]) > dv[0])   return 0;
  dfound[1] = Y(1,s2)  - Y(0,s1);  if(Abs(dfound[1]) > dv[1])   return 0;
  dfound[2] = TX(1,s2) - TX(0,s1); if(Abs(dfound[2]) > dv[2])   return 0;
  dfound[3] = TY(1,s2) - TY(0,s1); if(Abs(dfound[3]) > dv[3])   return 0;
  return 1;
}

//---------------------------------------------------------------------
Int_t EdbAlignmentV::Ncoins(float dv[4], EdbH2 *hdxy, EdbH2 *hdtxy, TObjArray *sel1, TObjArray *sel2 )
{
  // calculate all coinsidences inside dv:(dx,dy, dtx,dty)
  // fill position and angular 2-d plots if requested
  // fill arrays with selected segments sel1, sel2 if requested

  int ic = 0; // counter
  int n = CheckEqualArr(eS[0],eS[1]);

  if(hdxy)   hdxy->CleanCells();
  if(hdtxy)  hdtxy->CleanCells();

  float dfound[4];
  for(int i=0; i<n; i++)
    {
      EdbSegP *s1 = (EdbSegP*)eS[0].UncheckedAt(i);
      EdbSegP *s2 = (EdbSegP*)eS[1].UncheckedAt(i);

      if( !ValidCoinsidence(*s1, *s2, dv, dfound) )    continue;

      ic++;
      if(sel1)  sel1->Add(s1);
      if(sel2)  sel2->Add(s2);
      if(hdxy)  hdxy->Fill(dfound[0],dfound[1]);
      if(hdtxy) hdtxy->Fill(dfound[2],dfound[3]);
    }
  return ic;
}

//---------------------------------------------------------------------
void EdbAlignmentV::HDistance(EdbPattern &p1, EdbPattern &p2, float dxMax, float dyMax)
{
  Log(3,"HDistance","preselection for patterns with %d and %d segments",p1.N(), p2.N());

  int eITMAX=51;  // angular normalization (def=50) (the step is ~1./eITMAX)
  int eOCMAX=100; // occupancy (def=100)
  int nx=201, ny=201, no=eOCMAX;
 
  Long_t *ind1 = new Long_t[nx*ny*no];
  Long_t *ind2 = new Long_t[nx*ny*no];
  Int_t  *oc1  = new Int_t[nx*ny];
  Int_t  *oc2  = new Int_t[nx*ny];
  memset(ind1,0,nx*ny*no*sizeof(Long_t));
  memset(ind2,0,nx*ny*no*sizeof(Long_t));
  memset( oc1,0,nx*ny*sizeof(Int_t));
  memset( oc2,0,nx*ny*sizeof(Int_t));

  int itx,ity;  

  for(int i=0;i<eITMAX;i++) for(int j=0;j<eITMAX;j++) { oc1[i*nx+j]=0; oc2[i*nx+j]=0;} 
  for(int i=0;i<p1.N();i++){
    itx=int( (p1.GetSegment(i)->TX()+1.)*eITMAX/2. );
    ity=int( (p1.GetSegment(i)->TY()+1.)*eITMAX/2. );
    if(itx<eITMAX)
      if(ity<eITMAX)
	if(itx>=0)
	  if(ity>=0)
	    if(oc1[itx*nx+ity]<eOCMAX-1) 
	      { ind1[itx*nx*ny+ity*ny+oc1[itx*nx+ity]]=i; (oc1[itx*nx+ity])++;}
  }
  for(int i=0;i<p2.N();i++){
    itx=int( (p2.GetSegment(i)->TX()+1.)*eITMAX/2 );
    ity=int( (p2.GetSegment(i)->TY()+1.)*eITMAX/2 );
    if(itx<eITMAX)
      if(ity<eITMAX)
	if(itx>=0)
	  if(ity>=0)
	    if(oc2[itx*nx+ity]<eOCMAX-1) 
	      { ind2[itx*nx*ny+ity*ny+oc2[itx*nx+ity]]=i; (oc2[itx*nx+ity])++;}
  }
  
  eS[0].Clear();
  eS[1].Clear();

  EdbSegP *s1,*s2;
  for(int itx=0;itx<eITMAX;itx++) for(int ity=0;ity<eITMAX;ity++) 
    for(int i=0;i<oc1[itx*nx+ity];i++) {
      s1 = p1.GetSegment(ind1[itx*nx*ny+ity*ny+i]);
      for(int j=0;j<oc2[itx*nx+ity];j++)
	{
	  s2 = p2.GetSegment(ind2[itx*nx*ny+ity*ny+j]);
	  if(Abs(s2->eX-s1->eX)>dxMax) continue;
	  if(Abs(s2->eY-s1->eY)>dyMax) continue;
	  eS[0].Add( s1 );
	  eS[1].Add( s2 );
	}
    }
  delete [] ind1;
  delete [] ind2;
  delete [] oc1;
  delete [] oc2;
  Log(3,"HDistance","preselection for patterns with %d and %d segments: %d comb",p1.N(), p2.N(), eS[0].GetEntries());
}

//---------------------------------------------------------------------
void  EdbAlignmentV::InitPatCellN(EdbCell2 &cell, EdbPattern &pat, int nx, int ny)
{
  // note: no correction (do we need tis function?)
  float min[2] = {pat.Xmin()-eXmarg, pat.Ymin()-eYmarg };
  float max[2] = {pat.Xmax()+eXmarg, pat.Ymax()+eYmarg };
  int   n[2]   = {nx,ny};
  int maxcell = pat.N()/n[0]/n[1]+10;
  maxcell += (int)(5*Sqrt(maxcell));
  Log(2, "InitPatCellN", "maxcell = %d\n",maxcell);
  cell.InitCell(maxcell,n,min,max);
}

//---------------------------------------------------------------------
void EdbAlignmentV::InitPatCellBin( int side, EdbPattern &pat, float binx, float biny)
{
  if(!SideOK(side)) return;

  EdbSegP s1min(0,pat.Xmin(), pat.Ymin(), 0,0 );
  EdbSegP s1max(0,pat.Xmax(), pat.Ymax(), 0,0 );
  float min[2] = { X(side,s1min)-eXmarg, X(side,s1min)-eYmarg };
  float max[2] = { X(side,s1max)+eXmarg, X(side,s1max)+eYmarg };

  int   n[2]   = { (int)((max[0]-min[0])/binx+1), (int)((max[1]-min[1])/biny+1) };
  int maxcell = pat.N()/n[0]/n[1]+10;
  maxcell += (int)(5*Sqrt(maxcell));
  Log(2, "InitPatCellBin", "maxcell = %d\n",maxcell);
  ePC[side].InitCell(maxcell,n,min,max);
}

//---------------------------------------------------------------------
void  EdbAlignmentV::ApplyLimitsOffset( float &xmin1, float &xmax1, float &xmin2, float &xmax2, float offsetMax)
{
  if(Abs(xmin2-xmin1)>offsetMax)                                                // cut useless patterns margins
    if(xmin2>xmin1) xmin1=xmin2-offsetMax; 
    else  xmin2=xmin1-offsetMax;
  if(Abs(xmax2-xmax1)>offsetMax) 
    if(xmax2>xmax1) xmax2=xmax1+offsetMax; 
    else  xmax1=xmax2+offsetMax;
}

//---------------------------------------------------------------------
bool  EdbAlignmentV::DefineGuessCell( float xmin1, float xmax1, float ymin1, float ymax1, 
                                      float xmin2, float xmax2, float ymin2, float ymax2, int np1, int np2, float binOK)
{
  float s1 = (xmax1-xmin1)*(ymax1-ymin1);
  float xbin1 = Sqrt( s1/(np1/binOK)  );
  float s2 = (xmax2-xmin2)*(ymax2-ymin2);
  float xbin2 = Sqrt( s2/(np2/binOK)  );  
  float xbin = Min(xbin1,xbin2);
  float min[2] = { Min(xmin1,xmin2)-eXmarg, Min(ymin1,ymin2)-eYmarg };
  float max[2] = { Max(xmax1,xmax2)+eXmarg, Max(ymax1,ymax2)+eYmarg };
  int   n[2]   = { (int)((max[0]-min[0])/xbin+1), (int)((max[1]-min[1])/xbin+1) };
  if(n[0]<=0||n[1]<=0) return false;
  int maxcell1 = np1/n[0]/n[1]+10;   maxcell1 += (int)(5*Sqrt(maxcell1));
  ePC[0].InitCell(maxcell1,n,min,max);
  int maxcell2 = np2/n[0]/n[1]+10;    maxcell2 += (int)(5*Sqrt(maxcell2));
  ePC[1].InitCell(maxcell2,n,min,max);
  return true;
}

//---------------------------------------------------------------------
void  EdbAlignmentV::FillGuessCell( EdbPattern &p1, EdbPattern &p2, float binOK, float offsetMax)
{
  // binOK - is the mean number of entries requested per bin
  
  float xmin1=Xmin(0,p1), ymin1=Ymin(0,p1), xmax1=Xmax(0,p1), ymax1=Ymax(0,p1);
  float xmin2=Xmin(1,p2), ymin2=Ymin(1,p2), xmax2=Xmax(1,p2), ymax2=Ymax(1,p2);
  Log(3, "FillGuessCell", "x1:(%f %f) y1:(%f %f)",xmin1,xmax1, ymin1, ymax1);
  Log(3, "FillGuessCell", "x2:(%f %f) y2:(%f %f)",xmin2,xmax2, ymin2, ymax2);
  ApplyLimitsOffset( xmin1, xmax1, xmin2, xmax2, offsetMax);
  ApplyLimitsOffset( ymin1, ymax1, ymin2, ymax2, offsetMax);
  int np1 = p1.N(), np2 = p2.N();    
  if( xmax1<xmin1||ymax1<ymin1||xmax2<xmin2||ymax2<ymin2||np1<1||np2<1 ) {
    Log(1, "FillGuessCell", "No any overlap! x1:(%f %f) y1:(%f %f) np1:%d",
	xmin1,xmax1, ymin1, ymax1,np1);
    Log(1, "FillGuessCell", "No any overlap! x2:(%f %f) y2:(%f %f) np2:%d",
	xmin2,xmax2, ymin2, ymax2,np2);
  } else
  {
    if(DefineGuessCell( 
      xmin1, xmax1, ymin1, ymax1, 
      xmin2, xmax2, ymin2, ymax2, 
      np1, np2, binOK))
    {
      FillCell( 0, p1 );
      FillCell( 1, p2 ); 
      if( Log(3, "FillGuessCell", "with patterns of %d %d  statistics:",np1,np2) ) {
	ePC[0].PrintStat();
	ePC[1].PrintStat();
      }
    }
  }
}

//---------------------------------------------------------------------
void  EdbAlignmentV::FillGuessCell( TObjArray &p1, TObjArray &p2, float binOK, float offsetMax)
{
  // binOK - is the mean number of entries requested per bin

  float xmin1=Xmin(0,p1), ymin1=Ymin(0,p1), xmax1=Xmax(0,p1), ymax1=Ymax(0,p1);
  float xmin2=Xmin(1,p2), ymin2=Ymin(1,p2), xmax2=Xmax(1,p2), ymax2=Ymax(1,p2);
  ApplyLimitsOffset( xmin1, xmax1, xmin2, xmax2, offsetMax);
  ApplyLimitsOffset( ymin1, ymax1, ymin2, ymax2, offsetMax);
  int np1 = p1.GetEntries(), np2 = p2.GetEntries();

  Log(3, "FillGuessCell", "xmin1, xmax1, ymin1, ymax1, xmin2, xmax2, ymin2, ymax2, np1, np2, binOK: %f %f %f %f   %f %f %f %f   %d %d   %f",xmin1, xmax1, ymin1, ymax1, xmin2, xmax2, ymin2, ymax2, np1, np2, binOK);
  if( xmax1<xmin1||ymax1<ymin1||xmax2<xmin2||ymax2<ymin2||np1<1||np2<1 ) {
    Log(1, "FillGuessCell", "No any overlap! x1:(%f %f) y1:(%f %f) np1:%d",
	xmin1,xmax1, ymin1, ymax1,np1);
    Log(1, "FillGuessCell", "No any overlap! x2:(%f %f) y2:(%f %f) np2:%d",
	xmin2,xmax2, ymin2, ymax2,np2);
  } else
  {
    if(DefineGuessCell( 
      xmin1, xmax1, ymin1, ymax1, 
      xmin2, xmax2, ymin2, ymax2, 
      np1, np2, binOK))
    {     
      FillCell( 0, p1 );
      FillCell( 1, p2 );      
      if( Log(3, "FillGuessCell", "with arrays of %d %d  statistics:",np1,np2) ) {
	ePC[0].PrintStat();
	ePC[1].PrintStat();
      }
    }
  }
  
}

//---------------------------------------------------------------------
void  EdbAlignmentV::FillCell( int side, EdbPattern &pat )
{
  // assume that the cell is already initialized
  if(!SideOK(side)) return;
  EdbSegP *s=0;
  int n = pat.N();
  for(int i=0; i<n; i++) {
    s = pat.GetSegment(i);
    ePC[side].AddObject( X(side,*s), Y(side,*s), (TObject*)s );
  }
}

//---------------------------------------------------------------------
void  EdbAlignmentV::FillCell( int side, TObjArray &pat)
{
  // assume that the cell is already initialized
  if(!SideOK(side)) return;
  EdbSegP *s=0;
  int n = pat.GetEntriesFast();
  for(int i=0; i<n; i++) {
    s = (EdbSegP*)pat.UncheckedAt(i);
    ePC[side].AddObject( X(side,*s), Y(side,*s), (TObject*)s );
  }
}

//---------------------------------------------------------------------
Int_t  EdbAlignmentV::CalculateAffXY(TObjArray &arr1, TObjArray &arr2, EdbAffine2D &aff)
{
  int n = CheckEqualArr(arr1,arr2);
  TArrayF x1(n), x2(n), y1(n), y2(n);
  for(int i=0; i<n; i++ ) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1[i] = X(0,*s1);    y1[i] = Y(0,*s1);
    x2[i] = X(1,*s2);    y2[i] = Y(1,*s2);
  }
  aff.Calculate(n,x1.GetArray(),y1.GetArray(),x2.GetArray(),y2.GetArray(),0);
  return n;
}

//---------------------------------------------------------------------
Int_t  EdbAlignmentV::CalculateAffXYTurn(TObjArray &arr1, TObjArray &arr2, EdbAffine2D &aff)
{
  int n = CheckEqualArr(arr1,arr2);
  TArrayF x1(n), x2(n), y1(n), y2(n);
  for(int i=0; i<n; i++ ) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1[i] = X(0,*s1);    y1[i] = Y(0,*s1);
    x2[i] = X(1,*s2);    y2[i] = Y(1,*s2);
  }
  aff.Calculate(n,x1.GetArray(),y1.GetArray(),x2.GetArray(),y2.GetArray(),2);
  return n;
}

//---------------------------------------------------------------------
Int_t  EdbAlignmentV::CalculateAffXYShift(TObjArray &arr1, TObjArray &arr2, EdbAffine2D &aff)
{
  int n = CheckEqualArr(arr1,arr2);
  TArrayF x1(n), x2(n), y1(n), y2(n);
  for(int i=0; i<n; i++ ) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1[i] = X(0,*s1);    y1[i] = Y(0,*s1);
    x2[i] = X(1,*s2);    y2[i] = Y(1,*s2);
  }
  aff.CalculateShift(n,x1.GetArray(),y1.GetArray(),x2.GetArray(),y2.GetArray());
  return n;
}

//---------------------------------------------------------------------
Int_t  EdbAlignmentV::CalculateAffTXTY(TObjArray &arr1, TObjArray &arr2, EdbAffine2D &aff)
{
  int n = CheckEqualArr(arr1,arr2);
  TArrayF x1(n), x2(n), y1(n), y2(n);
  for(int i=0; i<n; i++ ) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1[i] = TX(0, *s1);    y1[i] = TY(0, *s1);
    x2[i] = TX(1, *s2);    y2[i] = TY(1, *s2);
  }
  aff.Calculate(n,x1.GetArray(),y1.GetArray(),x2.GetArray(),y2.GetArray(),0);
  return n;
}

//---------------------------------------------------------------------
Int_t  EdbAlignmentV::CalculateAffTXTYTurn(TObjArray &arr1, TObjArray &arr2, EdbAffine2D &aff)
{
  int n = CheckEqualArr(arr1,arr2);
  TArrayF x1(n), x2(n), y1(n), y2(n);
  for(int i=0; i<n; i++ ) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1[i] = TX(0, *s1);    y1[i] = TY(0, *s1);
    x2[i] = TX(1, *s2);    y2[i] = TY(1, *s2);
  }
  aff.Calculate(n,x1.GetArray(),y1.GetArray(),x2.GetArray(),y2.GetArray(),2);
  return n;
}

//---------------------------------------------------------------------
Int_t  EdbAlignmentV::CalculateAffTXTYShift(TObjArray &arr1, TObjArray &arr2, EdbAffine2D &aff)
{
  int n = CheckEqualArr(arr1,arr2);
  TArrayF x1(n), x2(n), y1(n), y2(n);
  for(int i=0; i<n; i++ ) {
    EdbSegP *s1 = (EdbSegP*)arr1.UncheckedAt(i);
    EdbSegP *s2 = (EdbSegP*)arr2.UncheckedAt(i);
    x1[i] = TX(0, *s1);    y1[i] = TY(0, *s1);
    x2[i] = TX(1, *s2);    y2[i] = TY(1, *s2);
  }
  aff.CalculateShift(n,x1.GetArray(),y1.GetArray(),x2.GetArray(),y2.GetArray());
  return n;
}

//---------------------------------------------------------------------
void EdbAlignmentV::Corr2Aff(EdbSegCorr &corr, EdbLayer &layer)
{
  EdbAffine2D &a0 = *(layer.GetAffineXY());
  a0.Reset();
  a0.ShiftX(  corr.V(0) );  a0.ShiftY( corr.V(1) );
  a0.Rotate(  corr.V(6) );
  EdbAffine2D &atxty = *(layer.GetAffineTXTY());
  atxty.Reset();
  atxty.ShiftX(corr.V(3)); atxty.ShiftY(corr.V(4));
  layer.SetShrinkage( corr.V(5) );
  layer.SetZcorr(     corr.V(2) );
}

//---------------------------------------------------------------------
void EdbAlignmentV::Corr2Aff(EdbLayer &layer)
{
/*  EdbAffine2D &a0 = *(layer.GetAffineXY());
  a0.Reset();
  a0.ShiftX(  eCorr[0].V(0) );  a0.ShiftY( eCorr[0].V(1) );
  a0.Rotate(  eCorr[0].V(6) );
  a0.ShiftX( -eCorr[1].V(0) );  a0.ShiftY( -eCorr[1].V(1) );
  layer.SetShrinkage( eCorr[0].V(5) );
  layer.SetZcorr(     eCorr[0].V(2) );*/
  
  Corr2Aff(eCorr[0],layer);
  layer.GetAffineXY()->ShiftX( -eCorr[1].V(0) );
  layer.GetAffineXY()->ShiftY( -eCorr[1].V(1) );
}

//---------------------------------------------------------------------
float  EdbAlignmentV::Xmin(int side,EdbPattern &p)
{
  int n=p.N();  float xmin = kMaxInt;
  if(n>0)
    for(int i=0; i<n; i++) {
      EdbSegP *s = p.GetSegment(i);
      float x=X(side,*s);
      if(x>kMinInt) xmin=(x<xmin)?x:xmin;
    }
    return xmin;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Xmax(int side,EdbPattern &p)
{
  int n=p.N();  float xmax = kMinInt;
  if(n>0)
    for(int i=0; i<n; i++) {
      EdbSegP *s = p.GetSegment(i);
      float x=X(side,*s);
      if(x<kMaxInt) xmax=(x>xmax)?x:xmax;
    }
    return xmax;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Ymin(int side,EdbPattern &p)
{
  int n=p.N();  float ymin = kMaxInt;
  if(n>0)
    for(int i=0; i<n; i++) {
      EdbSegP *s = p.GetSegment(i);
      float y=Y(side,*s);
      if(y>kMinInt) ymin=(y<ymin)?y:ymin;
    }
    return ymin;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Ymax(int side,EdbPattern &p)
{
  int n=p.N();  float ymax = kMinInt;
  if(n>0)
    for(int i=0; i<n; i++) {
      EdbSegP *s = p.GetSegment(i);
      float y=Y(side,*s);
      if(y<kMaxInt) ymax=(y>ymax)?y:ymax;
    }
    return ymax;
}


/*
//---------------------------------------------------------------------
float  EdbAlignmentV::Xmax(int side,EdbPattern &p)
{
  int n=p.N();  float xmax = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = p.GetSegment(i);
     float x=X(side,*s);
     if(!i) xmax=x;     else if( x>xmax) xmax=x;
  }
  return xmax;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Ymin(int side,EdbPattern &p)
{
  int n=p.N();  float ymin = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = p.GetSegment(i);
     float y=Y(side,*s);
     if(!i) ymin=y;     else if( y<ymin) ymin=y;
  }
  return ymin;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Ymax(int side,EdbPattern &p)
{
  int n=p.N();  float ymax = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = p.GetSegment(i);
     float y=Y(side,*s);
     if(!i) ymax=y;     else if( y>ymax) ymax=y;
  }
  return ymax;
}
*/

//---------------------------------------------------------------------
float  EdbAlignmentV::Xmin(int side,TObjArray &p)
{
  int n=p.GetEntries();  float xmin = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = (EdbSegP*)(p.At(i));
     float x=X(side,*s);
     if(!i) xmin=x;     else if( x<xmin) xmin=x;
  }
  return xmin;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Xmax(int side,TObjArray &p)
{
  int n=p.GetEntries();  float xmax = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = (EdbSegP*)(p.At(i));
     float x=X(side,*s);
     if(!i) xmax=x;     else if( x>xmax) xmax=x;
  }
  return xmax;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Ymin(int side,TObjArray &p)
{
  int n=p.GetEntries();  float ymin = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = (EdbSegP*)(p.At(i));
     float y=Y(side,*s);
     if(!i) ymin=y;     else if( y<ymin) ymin=y;
  }
  return ymin;
}
//---------------------------------------------------------------------
float  EdbAlignmentV::Ymax(int side,TObjArray &p)
{
  int n=p.GetEntries();  float ymax = 0;
  for(int i=0; i<n; i++) {
     EdbSegP *s = (EdbSegP*)(p.At(i));
     float y=Y(side,*s);
     if(!i) ymax=y;     else if( y>ymax) ymax=y;
  }
  return ymax;
}
