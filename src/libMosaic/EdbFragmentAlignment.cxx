#include "EdbLog.h"
#include "EdbMosaicPath.h"
#include "EdbFragmentAlignment.h"

ClassImp(EdbFragmentAlignment);

//-----------------------------------------------------------------------
void EdbFragmentAlignment::SetAlPar( const AlPar &ap, EdbPlateAlignment &av )
{
  av.eNoScaleRot  = ap.NoScaleRot;         // calculate shift only
  av.eOffsetMax   = ap.OffsetMax;
  av.eDZ          = ap.DZ;
  av.eDPHI        = ap.DPHI;
  av.eDoFine      = ap.DoFine;
  av.eSaveCouples = ap.DoSaveCouples;
  av.SetSigma( ap.SigmaR, ap.SigmaT );
  for(int i=0; i<4; i++) av.eDoublets[i] = ap.Doublets[i];
}

//-----------------------------------------------------------------------
void EdbFragmentAlignment::SetHarr( TObjArray &ha )
{
  eN = ha.GetEntries();
  eHarr.Expand(eN);
  for( int i=0; i<eN; i++ ) eHarr.AddAt( ha.At(i),i );
  eParr.Expand(eN);
  printf("*** eN=%d\n",eN);
  FillVC();
  printf("*** eN=%d\n",eN);
}

//-----------------------------------------------------------------------
void EdbFragmentAlignment::FillVC()
{
  eVC0 = new EdbPattern(0,0,0,eN);
  eVC  = new EdbPattern(0,0,0,eN);
  for( int i=0; i<eN; i++ ) {
    EdbViewHeader *h = GetHeader(i);
    EdbSegP *s0 = eVC0->AddSegment(i, h->GetXview(), h->GetYview(), 0, 0, 0,0 );
    EdbSegP *s  = eVC->AddSegment( i, h->GetXview(), h->GetYview(), 0, 0, 0,0 );
//    int side=0;
//    if(     h->GetNframesTop()==0) side=2;         // 2- bottom
//    else if(h->GetNframesBot()==0) side=1;         // 1- top
//    s0->SetAid( h->GetAreaID(),h->GetViewID(), side);
//    s->SetAid( h->GetAreaID(),h->GetViewID(), side);
    s0->SetAid( h->GetAreaID(),h->GetViewID(), eSide);
    s->SetAid( h->GetAreaID(),h->GetViewID(), eSide);
  }
}

//-----------------------------------------------------------------------
void EdbFragmentAlignment::AlignFragment( EdbPattern &pf )
{
  EdbMosaicPath mp(eN);
  mp.eR0=1200;
  EdbViewHeader *h = mp.FindNearest( eHarr, pf.X(), pf.Y() );
  if(h) {
    Log(1,"EdbFragmentAlignment::AlignFragment","with %d views at x0,y0 was %f %f:   %f %f",
	eN,h->GetXview(), h->GetYview(), pf.X(), pf.Y() );
    mp.InitArea( eHarr, h->GetXview(), h->GetYview(), eMinPeak );
  }
  else {
    Log(1,"EdbFragmentAlignment::AlignFragment","Warning! central view close to (%f %f) was not found! abandon this fragment",pf.X(), pf.Y());
    //mp.InitArea( eHarr, pf.X(), pf.Y(), eMinPeak );
    return;
  }
    
  AlignAndShift( mp );  
  RealignAndShift( mp );
  AlignAndShift( mp );
  RealignAndShift( mp );
  
  for(int i=0; i<eN; i++)
  {
    if( mp.OK(i) )   pf.AddPattern( *GetPattern(i) );
  }
  
  eAff.Reset();
  eAff.Calculate(eVC,eVC0);                     // from found to original
  Log(1,"EdbFragmentAlignment::AlignFragment","apply %s",eAff.AsString());
  eVC->Transform(&eAff);
  pf.Transform(&eAff);
}

//-----------------------------------------------------------------------
void EdbFragmentAlignment::FillVDT( EdbCouplesTree &vdt  )
{
  int n = eVC0->N();
  for(int i=0; i<n; i++)
  {
    vdt.Fill( eVC0->GetSegment(i), eVC->GetSegment(i), 0,0, 0,0, eID, eSide );
  }
}

//-----------------------------------------------------------------------
float EdbFragmentAlignment::CheckScaleX( float y0 )
{
  EdbMosaicPath mp(eN);
  float length = mp.InitLineX(eHarr, y0, 100. );
  mp.eR0 = 800;
  EdbAffine2D aff;
  CheckScale(mp,aff);
  return (length-aff.B1())/length;
}

//-----------------------------------------------------------------------
float EdbFragmentAlignment::CheckScaleY( float x0 )
{
  EdbMosaicPath mp(eN);
  float length = mp.InitLineY(eHarr, x0, 100. );
  mp.eR0 = 800;
  EdbAffine2D aff;
  CheckScale(mp,aff);
  return (length-aff.B2())/length;
}

//-----------------------------------------------------------------------
void EdbFragmentAlignment::CheckScale( EdbMosaicPath &mp, EdbAffine2D &aff )
{
  int n = mp.N();
  mp.SetOK( (mp.I(0)) );
  
  for(int i=1; i<n; i++) 
  {
    EdbPattern *p = (EdbPattern *)(eParr.At(mp.I(i)));
    TArrayI narr(10);
    int nb = mp.GetAlignedNeighbours( mp.I(i), narr );
    if(nb>10)
    {
      Log(1,"EdbFragmentAlignment::CheckScale","Warning! too many neigbours: %d  in dr = %f  reset to 10",nb, mp.eR0);
      nb=10;
    }
    EdbPattern  alp;
    for(int ii=0; ii<nb; ii++) {
      alp.AddPattern( *(EdbPattern *)(eParr.At(narr[ii])) );
    }
    printf("align %d -> %d  at dist %.1f \n",p->N(), alp.N(), mp.Dist(mp.I(i)) );
    if( ViewSideAl(*p, alp, aff, 0) > eMinPeak ) mp.SetOK( mp.I(i) );
  }
}

//-----------------------------------------------------------------------
void EdbFragmentAlignment::ApplyAff()
{
  Log(1,"EdbFragmentAlignment::ApplyAff","%s",eAff.AsString());
  for( int i=0; i<eN; i++ ) GetPattern(i)->Transform(&eAff);
}

//-------------------------------------------------------------------
void EdbFragmentAlignment::AlignAndShift( EdbMosaicPath &mp )
{
  Log(1,"EdbFragmentAlignment::AlignAndShift","");
  mp.SetOK( mp.I(0) );
  for(int i=1; i<eN; i++)
  {
    EdbPattern *p = GetPattern(mp.I(i));
    TArrayI narr(10);
    int nb = mp.GetAlignedNeighbours( mp.I(i), narr );
    if(nb>10) 
    { 
      Log(1,"EdbFragmentAlignment::AlignAndShift","Warning! too many neigbours: %d  in dr = %f  reset to 10",nb, mp.eR0);
      nb=10;
    }
    EdbPattern  alp;
    for(int ii=0; ii<nb; ii++) {
      alp.AddPattern( *GetPattern(narr[ii]) );
    }
    printf("align %d -> %d  at dist %.1f \n",p->N(), alp.N(), mp.Dist(mp.I(i)) );
    EdbAffine2D aff;
    int peak = ViewSideAl(*p, alp, aff, 1);
    EdbSegP *s = eVC->GetSegment( mp.I(i) );
    s->SetW(peak);
    if( peak > eMinPeak )
    {
      mp.SetOK( mp.I(i) );
      s->Transform(&aff);
    }
  }
}

//-------------------------------------------------------------------
void EdbFragmentAlignment::RealignAndShift( EdbMosaicPath &mp )
{
  //assume that most of views were already aligned
  Log(1,"EdbFragmentAlignment::RealignAndShift","");
  for(int i=0; i<eN; i++)
  {
    if(mp.OK(i)) continue;
    EdbPattern *p = GetPattern(i);
    TArrayI narr(10);
    int nb = mp.GetAlignedNeighbours( i, narr );
    if(nb>10) 
    { 
      Log(1,"EdbFragmentAlignment::AlignAndShift","Warning! too many neigbours: %d  in dr = %f  reset to 10",nb, mp.eR0);
      nb=10;
    }
    EdbPattern  alp;
    for(int ii=0; ii<nb; ii++) {
      alp.AddPattern( *GetPattern(narr[ii]) );
    }
    printf("align %d -> %d  at dist %.1f \n",p->N(), alp.N(), mp.Dist(i) );
    EdbAffine2D aff;
    int peak = ViewSideAl(*p, alp, aff, 1);
    EdbSegP *s = eVC->GetSegment( i );
    s->SetW(peak);
    if( peak > eMinPeak )
    {
      mp.SetOK( i );
      s->Transform(&aff);
    }
  }
}

//-----------------------------------------------------------------------
int EdbFragmentAlignment::ViewSideAl( EdbPattern &p1, EdbPattern &p2, EdbAffine2D &aff, bool do_shift )
{
  // transform aff instead of pattern
  EdbPlateAlignment av;
  SetAlPar( eAP, av );
  av.eSaveCouples=0;
  //av.InitOutputFile( Form( "p%.3d/%d_%d.al.vsa.root", p2.Plate(), p1.ID(), p2.ID() ) ); 
  av.Align( p1, p2, 0);
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();
  if(av.eNcoins > eMinPeak ) 
  {
    if(do_shift) p1.Transform( affXY );
    aff.Transform( affXY );
  }
  //av.CloseOutputFile();
  return av.eNcoins;  
}

void EdbFragmentAlignment::AddPatternAt( EdbPattern *p, int i ) 
{ 
  if(i>=0&&i<eParr.GetSize()) eParr.AddAt( p, i );
  else Log(1,"EdbFragmentAlignment::AddPatternAt","Error! i = %d is out of bounds(0:%d)",i,eParr.GetSize());
}