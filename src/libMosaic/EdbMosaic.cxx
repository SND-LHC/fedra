//-- Author :  Valeri Tioukov   15/02/2024

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaic                                                            //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TCut.h>
#include "EdbMosaicIO.h"
#include "EdbMosaic.h"
#include "EdbPlateAlignment.h"
#include "EdbScanProc.h"
#include "EdbLog.h"

using namespace TMath;

ClassImp(EdbMosaicAl);

//-----------------------------------------------------------------------
EdbMosaicAl::EdbMosaicAl()
{
  eCorrMap[0]=eCorrMap[1]=eCorrMap[2]=0;
  eH_XY[0] = 0;
  eH_XY[1] = 0;
  eH_XY[2] = 0;
}

//-----------------------------------------------------------------------
EdbMosaicAl::~EdbMosaicAl()
{
  for(int i=0; i<3; i++) 
  {
    if(eCorrMap[i]) delete eCorrMap[i];
    if(eH_XY[i])    delete eH_XY[i];
  }
}

//-----------------------------------------------------------------------
void EdbMosaicAl::InitH2( const TEnv &env )
{
  int    nbinsx, nbinsy;
  float  xlow, xup, ylow, yup;
  const char *str = env.GetValue( "fedra.vsa.hxy", "2000 0. 200000. 2000 0. 200000." );
  if(6 != sscanf(str,"%d %f %f %d %f %f",
           &nbinsx, &xlow, &xup, &nbinsy, &ylow, &yup)) {
    Log(1,"EdbMosaicAl::InitH2","Error parsing hxy parameters: %s");
    return;
  }
  for(int i=1; i<=2; i++) {
    if(eH_XY[i]) delete eH_XY[i];
    eH_XY[i] = new TH2F( Form("hxy%d",i), Form("xy side %d",i), nbinsx,xlow,xup, nbinsy,ylow,yup);
  }
}

//-----------------------------------------------------------------------
void EdbMosaicAl::ProcRun( EdbID id, const TEnv &env )
{
  InitH2(env);
  EdbID idset =id; idset.ePlate =0;
  EdbScanProc sproc;
  sproc.eProcDirClient="..";
  TString fin;
  sproc.MakeFileName(fin,id,"raw.root");
  TCut cut("cut", env.GetValue("fedra.vsa.HeaderCut" , "1") );
  eRAW.eAFID = 1;
  eRAW.eHeaderCut = cut.GetTitle();
  eRAW.AddSegmentCut(1,env.GetValue("fedra.vsa.ICUT"      , "-1") );
  eRAW.AddSegmentCut(0,env.GetValue("fedra.vsa.XCUT"      , "-1") );
  eRAW.SetScoreCut( env.GetValue("fedra.vsa.ScoreCUT"  , "-1.5 1.5") );
  
  if(!sproc.InitRunAccessNew(eRAW,idset,id.ePlate)) return;  
  eRAW.eDoImageMatrixCorr = env.GetValue("fedra.vsa.DoImageMatrixCorr", 0  );
  if(eRAW.eDoImageMatrixCorr) {
    eRAW.ReadImageMatrixCorrection( 1, env.GetValue("fedra.vsa.ImageMatrixCorrSide1"      , "") );
    eRAW.ReadImageMatrixCorrection( 2, env.GetValue("fedra.vsa.ImageMatrixCorrSide2"      , "") );
  }

  eRAW.eDoImageCorr = env.GetValue("fedra.vsa.DoImageCorr", 0  ); 
  if(eRAW.eDoImageCorr) {
    eRAW.SetImageCorrectionAff( 1, env.GetValue("fedra.vsa.ImageCorrSide1"      , "1 0 0 1 0 0") );
    eRAW.SetImageCorrectionAff( 2, env.GetValue("fedra.vsa.ImageCorrSide2"      , "1 0 0 1 0 0") );
  }

  float fx = env.GetValue("fedra.vsa.Xfrag" , 10000);
  float fy = env.GetValue("fedra.vsa.Yfrag" , 5000);
  eMinPeak = env.GetValue("fedra.vsa.MinPeak" , 20);  
  eR0      = env.GetValue("fedra.vsa.R0"      , 1200);
  eDoSaveAlignments = env.GetValue("fedra.vsa.DoSaveAlignments", 0);
  eNoAlign = env.GetValue("fedra.vsa.NoAlign", 0);
  SetAlPar( env, eAP );

  EdbViewMap vm;
  vm.ReadViewsHeaders(fin.Data(), cut);    // read headers from runfile, fill eViewHeaders
  int nh = vm.eViewHeaders.GetEntries();
  Log(1,"EdbMosaicAl::ProcRun","%d headers selected by cut %s", nh, cut.GetTitle() );
  if(nh>0) {
    FormFragments( fx,fy, vm.eViewHeaders );
    
    eID=id;
    eMIO.Init( Form("p%3.3d/%d.%d.%d.%d.mos.root",
	       id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor), 
	       "RECREATE");
	       
    AlignFragments();
    
    if(eCorrMap[1]) eMIO.SaveCorrMap(id.ePlate, 1, *eCorrMap[1]);
    if(eCorrMap[2]) eMIO.SaveCorrMap(id.ePlate, 2, *eCorrMap[2]);
    
    if(eH_XY[1]) eMIO.SaveSideObj( eH_XY[1], id.ePlate, 1, "hxy_" );
		if(eH_XY[2]) eMIO.SaveSideObj( eH_XY[2], id.ePlate, 2, "hxy_" );
  }
}

//-----------------------------------------------------------------------
void EdbMosaicAl::SetAlPar( const TEnv &env, AlPar &ap )
{
  ap.NoScaleRot    = 1;                // calculate shift only
  ap.OffsetMax     = env.GetValue("fedra.vsa.offsetMax",15);
  ap.DZ            = env.GetValue("fedra.vsa.DZ",0);
  ap.DPHI          = env.GetValue("fedra.vsa.DPHI",0);
  ap.DoFine        = env.GetValue("fedra.vsa.DoFine",1);
  ap.DoSaveCouples = env.GetValue("fedra.vsa.SaveCouples",1);
  ap.SigmaR        = env.GetValue("fedra.vsa.SigmaR",0.5);
  ap.SigmaT        = env.GetValue("fedra.vsa.SigmaT",0.005);
  ap.Doublets[0]   = ap.Doublets[1] = 0.25;
  ap.Doublets[2]   = ap.Doublets[3] = 0.0025;
}

//-----------------------------------------------------------------------
void EdbMosaicAl::AlignFragments()
{
  EdbCouplesTree vdt;
  //vdt.InitCouplesTree("couples", Form("p%3.3d/%d.%d.%d.%d.vdt.root",
//		  eID.ePlate, eID.eBrick, eID.ePlate, eID.eMajor, eID.eMinor),"RECREATE");
  vdt.InitCouplesTree("couples",0,"NEW"); // by default create tree in .mos.root file already opened

  int nc=eCF[1].Ncell();
  for( int side=1; side<=2; side++ ) {
   for( int i=0; i<nc; i++ ) {
     TObjArray *a = (TObjArray *)(eCF[side].GetObject(i,0));
     EdbLayer *l = eCorrMap[side]->Map().GetLayer(i);
     if(a) 
      {
	EdbPattern pf( eCF[side].Xj(i), eCF[side].Yj(i), 0, 2000 );
	pf.SetSide(side);
	pf.SetX(eCF[side].Xj(i));
	pf.SetY(eCF[side].Yj(i));
	pf.SetID(i);
	pf.SetScanID(eID);
	
	if(l) 
	{
	  l->SetXY(pf.X(), pf.Y());
	  l->SetID(pf.ID());
	}
	
	EdbFragmentAlignment fa;
	fa.SetID(i);
	fa.SetSide(side);
	fa.SetHarr(*a);
	fa.SetMinPeak( eMinPeak );
  fa.eR0 = eR0;
  fa.eDoSaveAlignments = eDoSaveAlignments;
	ReadPatterns( fa );
	fa.eAP=eAP;
	if(eNoAlign) fa.AlignFragment_dummy(pf);
	else         fa.AlignFragment(pf);
	fa.FillVDT(vdt);
	FillHXY( pf, *(eH_XY[side]) );
	eMIO.SaveFragment( pf );
      }
    }
  }
  vdt.SaveTree();
}

//-----------------------------------------------------------------------
void EdbMosaicAl::FillHXY( const EdbPattern &pf, TH2F &h )
{
  int n=pf.N();
  for(int i=0; i<n; i++) 
  {
    EdbSegP *s = pf.GetSegment(i);
    h.Fill(s->eX,s->eY);
  }
}

//-----------------------------------------------------------------------
void EdbMosaicAl::ReadPatterns( EdbFragmentAlignment &fa )
{
  for(int i=0; i<fa.N(); i++)
  {
    EdbViewHeader *h=  fa.GetHeader(i);
    EdbPattern *p = new EdbPattern( h->GetXview(), h->GetYview(), 0, 2000 );
    int nrej=0;
    int nseg = eRAW.GetPatternView(  *p, fa.Side(), h->GetStatus(), nrej );
    p->SetID(h->GetStatus());
    fa.AddPatternAt(p,i);
  }
  fa.SetPatternsOwner();
}

//-----------------------------------------------------------------------
void EdbMosaicAl::AlignFragment( EdbPattern &pf, TObjArray &harr )
{
  // do not join view patterns
  int nh = harr.GetEntries();  
  EdbMosaicPath mp(nh);
  mp.eR0=eR0;
  mp.InitArea(harr, pf.X(), pf.Y(), eMinPeak );

  Log(1,"EdbMosaicAl::AlignFragment","with %d views at x0,y0:   %f %f",nh,mp.eX0,mp.eY0);

  TObjArray parr(nh);     // read views as a patterns to this array
  for(int i=0; i<nh; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    EdbPattern *p = new EdbPattern( h->GetXview(), h->GetYview(), 0, 2000 );
    int nrej=0;
    int nseg = eRAW.GetPatternView(  *p, pf.Side(), h->GetStatus(), nrej );
    p->SetID(h->GetStatus());
    parr.Add(p);
  }
  parr.SetOwner();
  
  AlignSpot(parr, mp);
  AlignSpot(parr, mp);
 
  for(int i=0; i<nh; i++)
  {
    EdbPattern *p = (EdbPattern *)(parr.At(i));
    if( mp.OK(i) )
    {
      pf.AddPattern(*p);
    }
  }

  parr.Delete();
}


//-------------------------------------------------------------------
void EdbMosaicAl::AlignSpot(TObjArray &parr, EdbMosaicPath &mp)
{
  int nh = parr.GetEntries();  
  mp.SetOK( mp.I(0) );  
  for(int i=1; i<nh; i++)
  {
    EdbPattern *p = (EdbPattern *)(parr.At(mp.I(i)));
    if(p->N() >= eMinPeak) {
      TArrayI narr(10);
      int nb = mp.GetAlignedNeighbours( mp.I(i), narr );
      if(nb>10) { Log(1,"EdbMosaicAl::AlignSpot","too many neigbours! %d ", nb ); nb=10;}
      else { Log(3,"EdbMosaicAl::AlignSpot","%d neighbours found", nb );}
      EdbPattern  alp;
      for(int ii=0; ii<nb; ii++) {
        alp.AddPattern(*(EdbPattern *)(parr.At(narr[ii])) );
      }
      printf("align %d -> %d  at dist %.1f \n",p->N(), alp.N(), mp.Dist(mp.I(i)) );
      if( ViewSideAl0(*p, alp) > eMinPeak ) mp.SetOK( mp.I(i) );
    } else Log(2,"EdbMosaicAl::AlignSpot","skip view with %d segments", p->N() );
  }  
}

//-------------------------------------------------------------------
int EdbMosaicAl::ViewSideAl0(EdbPattern &p1, EdbPattern &p2)
{
  // align AND TRANSFORM p1 to p2 RS
  
//  gEDBDEBUGLEVEL       =  1;
  EdbPlateAlignment av;
  av.eNoScaleRot=1;                // calculate shift only
  av.SetSigma(eAP.SigmaR,eAP.SigmaT);
  av.eOffsetMax = eAP.OffsetMax;
  av.eDZ        = eAP.DZ;
  av.eDPHI      = eAP.DPHI;
  av.eDoFine = eAP.DoFine;
  av.eSaveCouples = eAP.DoSaveCouples;
  for(int i=0; i<4; i++) av.eDoublets[i]=eAP.Doublets[i];

  //av.InitOutputFile( Form( "p%.3d/%d_%d.al.vsa.root", p2.Plate(), p1.ID(), p2.ID() ) ); 
  av.Align( p1, p2, 0);
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();
  if(av.eNcoins > eMinPeak ) 
  {
    p1.Transform( affXY );
  }
  //av.CloseOutputFile();
  return av.eNcoins;
}

//-----------------------------------------------------------------------
void EdbMosaicAl::FormFragments( float fx, float fy, TObjArray &harr )
{
  float xmin,xmax,ymin,ymax;  
  int nh = harr.GetEntries();
  EdbViewHeader *h=0;
  for(int i=0; i<nh; i++)
  {
    h=(EdbViewHeader *)(harr.At(i));
    //h->SetStatus(i);                      // keep tree entry number here
    if(i==0) 
    {
      xmin=xmax = h->GetXview();
      ymin=ymax = h->GetYview();
    } 
    else
    {
      if(xmin>h->GetXview() ) xmin = h->GetXview();
      if(xmax<h->GetXview() ) xmax = h->GetXview();
      if(ymin>h->GetYview() ) ymin = h->GetYview();
      if(ymax<h->GetYview() ) ymax = h->GetYview();
    }
  }

  float xstep = 800;
  float ystep = 600;
 
  int nx = Max(1,int((xmax-xmin+xstep)/fx));  //TODO: optimize or fix bin size?
  int ny = Max(1,int((ymax-ymin+ystep)/fy));
  
  for( int iside=1; iside<=2; iside++)
  {
    eCF[iside].InitCell(nx, xmin-xstep*2./3., xmax+xstep*1./3., ny, ymin-ystep*2./3., ymax+ystep*1./3., 1);
    eCorrMap[iside] = new EdbLayer();
    eCorrMap[iside]->Map().Init(eCF[iside]);
    int nc=eCF[iside].Ncell();
    for( int i=0; i<nc; i++ ) {
      TObjArray *a = new TObjArray();
      eCF[iside].AddObject(i, (TObject*)a );
      eCF[iside].SetBin(i,0);
    }
  }
  
  for(int i=0; i<nh; i++)
  {
    h=(EdbViewHeader *)(harr.At(i));
    int iside=0;
    if(     h->GetNframesTop()==0) iside=2;         // 2- bottom
    else if(h->GetNframesBot()==0) iside=1;         // 1- top
    if(iside)
    {
      TObjArray *a = (TObjArray *)(eCF[iside].GetObject( h->GetXview(), h->GetYview(),0));
      a->Add(h);
      eCF[iside].Fill(h->GetXview(), h->GetYview());
    }
    //eCF[iside].PrintStat();
  }

}
