//-- Author :  Valeri Tioukov   15/02/2024

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaic                                                          //
//                                                                      //
//                                                //
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
void EdbMosaicAl::ProcRun( EdbID id, const TEnv &env )
{
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

  if(!sproc.InitRunAccessNew(eRAW,idset,id.ePlate)) return;  
  eRAW.eDoImageMatrixCorr = env.GetValue("fedra.vsa.DoImageMatrixCorr", 0  );
  if(eRAW.eDoImageMatrixCorr) {
    eRAW.ReadImageMatrixCorrection( 1, env.GetValue("fedra.vsa.ImageMatrixCorrSide1"      , "") );
    eRAW.ReadImageMatrixCorrection( 2, env.GetValue("fedra.vsa.ImageMatrixCorrSide2"      , "") );
  }

  float fx = env.GetValue("fedra.vsa.Xfrag" , 10000);
  float fy = env.GetValue("fedra.vsa.Yfrag" , 5000);
  eMinPeak = env.GetValue("fedra.vsa.MinPeak" , 20);  

  EdbViewMap vm;
  vm.ReadViewsHeaders(fin.Data(), cut);    // read headers from runfile, fill eViewHeaders
  
  FormFragments( fx,fy, vm.eViewHeaders );
  
  eID=id;
  eMIO.Init( Form("p%3.3d/%d.%d.%d.%d.mos.root",
		  id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor), 
	     "RECREATE");
  
  SetAlPar( env, eAP );

  AlignFragments();

  if(eCorrMap[1]) eMIO.SaveCorrMap(id.ePlate, 1, *eCorrMap[1]);
  if(eCorrMap[2]) eMIO.SaveCorrMap(id.ePlate, 2, *eCorrMap[2]);
}

//-----------------------------------------------------------------------
void EdbMosaicAl::SetAlPar( const TEnv &env, AlPar &ap )
{
  ap.NoScaleRot    = 1;                // calculate shift only
  ap.OffsetMax     = 15.;
  ap.DZ            = 0;
  ap.DPHI          = 0;
  ap.DoFine        = 1;
  ap.DoSaveCouples = 1;
  ap.SigmaR        = 0.5; 
  ap.SigmaT        = 0.005;
  ap.Doublets[0]   = ap.Doublets[1] = 0.5;
  ap.Doublets[2]   = ap.Doublets[3] = 0.005;
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
	ReadPatterns( fa );
	fa.eAP=eAP;
	fa.AlignFragment(pf);
	fa.FillVDT(vdt);
	eMIO.SaveFragment( pf );
      }
    }
  }
  vdt.SaveTree();
}

//-----------------------------------------------------------------------
void EdbMosaicAl::SetAlPar( EdbFragmentAlignment &fa )
{
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
  mp.eR0=1200;
  mp.InitArea(harr, pf.X(), pf.Y() );

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
    TArrayI narr(10);
    int nb = mp.GetAlignedNeighbours( mp.I(i), narr );
    EdbPattern  alp;
    for(int ii=0; ii<nb; ii++) {
      alp.AddPattern(*(EdbPattern *)(parr.At(narr[ii])) );
    }
    printf("align %d -> %d  at dist %.1f \n",p->N(), alp.N(), mp.Dist(mp.I(i)) );
    if( ViewSideAl0(*p, alp) > eMinPeak ) mp.SetOK( mp.I(i) );
  }  
}

//-------------------------------------------------------------------
int EdbMosaicAl::ViewSideAl0(EdbPattern &p1, EdbPattern &p2)
{
  // align AND TRANSFORM p1 to p2 RS
  
  gEDBDEBUGLEVEL       =  1;
  float      sigmaR    =  0.5;
  float      sigmaT    =  0.005;
  float      offsetMax = 15. ;
  float      DZ        =  0.;
  float      DPHI      =  0.;
  
  EdbPlateAlignment av;
  av.eNoScaleRot=1;                // calculate shift only
  av.SetSigma(sigmaR,sigmaT);
  av.eOffsetMax = offsetMax;
  av.eDZ        = DZ;
  av.eDPHI      = DPHI;
  av.eDoFine = 1;
  av.eSaveCouples = 1;
  av.eDoublets[0]=av.eDoublets[1]=0.5;
  av.eDoublets[2]=av.eDoublets[3]=0.005;

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
    eCF[iside].PrintStat();
  }

}
