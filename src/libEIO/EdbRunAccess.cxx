//-- Author :  Valeri Tioukov   15.11.2004

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbRunAccess                                                         //
//                                                                      //
// OPERA Run Access helper class                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "Riostream.h"
#include "TSystem.h"
#include "TArrayL.h"
#include "TArrayI.h"
#include "TEventList.h"
#include "EdbMath.h"
#include "EdbSegment.h"
#include "EdbCluster.h"
#include "EdbRunAccess.h"
#include "EdbLog.h"
#ifndef __CINT__
#include "libDataConversion.h"
#endif

using namespace TMath;

ClassImp(EdbRunAccess)

///_________________________________________________________________________
EdbRunAccess::EdbRunAccess()
{
  Set0();
}

///_________________________________________________________________________
EdbRunAccess::EdbRunAccess(const char *fname)
{
  Set0();
  eRunFileName=fname;
}

///_________________________________________________________________________
EdbRunAccess::~EdbRunAccess()
{
  if(eRun) { delete eRun; eRun=0; }
  for(int i=0; i<3; i++) {
    if(eVP[i])     { delete eVP[i];     eVP[i]=0;     }
    if(eLayers[i]) { delete eLayers[i]; eLayers[i]=0; }
    if(eCuts[i])   { delete eCuts[i];   eCuts[i]=0;   }
    if(eCond[i])   { delete eCond[i];   eCond[i]=0;   }
  }
}

///_________________________________________________________________________
EdbRunAccess::EdbRunAccess(EdbRun *r)
{
  Set0();
  eRun=r;
}

///_________________________________________________________________________
void EdbRunAccess::Set0()
{
  eRun=0;
  for(int i=0; i<3; i++) {
    eVP[i]=0;
    eLayers[i]=0;
    eCuts[i]=0;
    eCond[i]=0;
  }
  eAFID=0;
  eCLUST=0;
  eFirstArea = 9999999;
  eLastArea  = 0;
  eNareas=0;
  eNviewsPerArea=0;
  eXmin=eXmax=eYmin=eYmax=0;
  eUseExternalSurface=false;
  eDoViewAnalysis=true;
  eHViewXY[1].InitH2(500,-250,250, 500, -250, 250);
  eHViewXY[2].InitH2(500,-250,250, 500, -250, 250);
  eDoImageCorr=0;
  eHeaderCut="1";
  eTracking=-1;
  eWeightAlg=0;
  eViewCorr = 0;
  eAffStage2Abs=0;
  eInvertSides=0;
  for(int i=0; i<4; i++) eGraphZ[i] = 0;
  for(int i=0; i<3; i++) eGraphDZ[i] = 0;
}

///_________________________________________________________________________
void EdbRunAccess::ClearCuts()
{
  for(int i=0; i<3; i++) {
    if(eCuts[i]) eCuts[i]->Clear();
  }
}

///_________________________________________________________________________
void EdbRunAccess::Print()
{
  printf(" EdbRunAccess: eAFID=%d  eCLUST=%d\n",eAFID, eCLUST );
  printf("file: %s\n", eRunFileName.Data());
  for(int i=0; i<3; i++) 
    if(eLayers[i]) eLayers[i]->Print();
 }

///_________________________________________________________________________
bool EdbRunAccess::InitRun(const char *runfile, bool do_update)
{
  if(eRun) SafeDelete(eRun);
  if(runfile) eRunFileName=runfile;
  if(eAFID==2) ReadVAfile();
  if( gSystem->AccessPathName(eRunFileName.Data(), kFileExists) ) {
    Log(1,"EdbRunAccess::InitRun","ERROR! open file: %s\n",eRunFileName.Data());
    return false;
  }
  if(do_update) eRun=new EdbRun(eRunFileName.Data(),"UPDATE");
  else          eRun=new EdbRun(eRunFileName.Data());
  if(!eRun) {Log(1,"EdbRunAccess::InitRun","ERROR! can't open file %s\n",eRunFileName.Data()); return false; }
  if(!(eRun->GetTree())) {Log(1,"EdbRunAccess::InitRun","ERROR! Vews tree is missed %s\n",eRunFileName.Data());return false;}
  for(int i=0; i<3; i++) if(eLayers[i]==0) GetMakeLayer(i);
  if(eAFID==11) {
    eRun->GetMarks()->Print();
    if(eRun->GetMarks()) eAffStage2Abs = eRun->GetMarks()->Stage2Abs();
    if(eAffStage2Abs) eAffStage2Abs->Print();
    else Log(1,"EdbRunAccess::InitRun","ERROR! Stage2Abs do not available (eAFID=11)");
  }
  if(do_update) {
    if(FillVP()<1) return false;
    eVP[1]->Write("views_s1");
    eVP[2]->Write("views_s2");
    eRun->GetTree()->GetCurrentFile()->Purge();
  } else {
    if( eVP[1] ) delete  eVP[1];
    if( eVP[2] ) delete  eVP[2];
    eVP[1] = dynamic_cast<EdbPattern*>(gDirectory->Get("views_s1"));
    eVP[2] = dynamic_cast<EdbPattern*>(gDirectory->Get("views_s2"));
  }

  int nn=0;
  if( eVP[1] && eVP[2] ) nn = eVP[1]->N()+eVP[2]->N();
  if(nn != eRun->GetEntries()) 
    if(FillVP()<1) 
    {
      Log(1,"EdbRunAccess::InitRun","no views accepted");
      return false;
    }
  eVP[1]->Transform(eLayers[1]->GetAffineXY());
  eVP[2]->Transform(eLayers[2]->GetAffineXY());
  eXmin = eVP[1]->Xmin();
  eXmax = eVP[1]->Xmax();
  eYmin = eVP[1]->Ymin();
  eYmax = eVP[1]->Ymax();

  GuessNviewsPerArea();
  if(gEDBDEBUGLEVEL>2) PrintStat();
  for(int i=0; i<3; i++) { eXstep[i]=400;   eYstep[i]=400; }
  
  Log(2,"EdbRunAccess::InitRun","%s AFID=%d",runfile, eAFID);
  return true;
}

///_________________________________________________________________________
void EdbRunAccess::GuessNviewsPerArea()
{
  eNviewsPerArea = 0;
  int nviews=eRun->GetEntries();
  if(!nviews) return;
  
  int nareas=0;
  int idmin = kMaxInt;
  int idmax = kMinInt;
  for(int iview=0; iview<nviews; iview++) {
    int aid = eRun->GetEntry(iview,1,0,0,0,0)->GetAreaID();
    if(aid>idmax) 
    {
      idmax=aid;
      nareas++;
    }
    if(aid<idmin) idmin=aid;
  }
  eNareas=nareas;
  eFirstArea=idmin;
  eLastArea=idmax;
  eNviewsPerArea=nviews/eNareas/2;
  /*
  if( eNviewsPerArea*eNareas*2 != nviews ) 
  {                                           // try if value in run header is correct
    int nva=0;
    if( eRun->GetHeader() && eRun->GetHeader()->GetArea() )
    {
      nva = eRun->GetHeader()->GetArea()->GetN();
      if( nva>0 && nva*eNareas*2==nviews )    eNviewsPerArea=nva;
    }
  }
  */
  if(eNviewsPerArea*eNareas*2 != nviews)
    Log(1,"EdbRunAccess::GuessNviewsPerArea","WARNING: eNviewsPerArea*eNareas*2 != nviews: %d*%d*2 != %d",
      eNviewsPerArea, eNareas, nviews );
}

///_________________________________________________________________________
void EdbRunAccess::ReadVAfile()
{
  int k = eRunFileName.Length();
  TString s(eRunFileName.Data(),k-4); s+="va.root";
  Log( 2,"EdbRunAccess::ReadVAfile","Get view corrections from %s",s.Data() );
  std::unique_ptr<TFile> f(TFile::Open(s.Data(),"READ"));
  if (!f || f->IsZombie())  Log(1,"EdbRunAccess::ReadVAfile","ERROR! can not open file %s",s.Data());
  else                      eViewCorr = (TObjArray*)f->Get("viewcorr");
}

///_________________________________________________________________________
void EdbRunAccess::SetCutLeft(int side, float wmin)
{
   float xmin[5]={ eViewXmin[side]                 , eViewYmin[side], -1, -1, wmin};
   float xmax[5]={ eViewXmin[side] + OverlapX(side), eViewYmax[side],  1,  1, 50};
   ClearCuts(); AddSegmentCut(side,1,xmin,xmax);
}
///_________________________________________________________________________
void EdbRunAccess::SetCutRight(int side, float wmin)
{
   float xmin[5]={ eViewXmax[side] - OverlapX(side), eViewYmin[side], -1, -1, wmin};
   float xmax[5]={ eViewXmax[side]                 , eViewYmax[side],  1,  1, 50};
   ClearCuts(); AddSegmentCut(side,1,xmin,xmax);
}
///_________________________________________________________________________
void EdbRunAccess::SetCutTop(int side, float wmin)
{
   float xmin[5]={ eViewXmin[side] + OverlapX(side), eViewYmax[side] - OverlapY(side), -1, -1, wmin};
   float xmax[5]={ eViewXmax[side] - OverlapX(side), eViewYmax[side]                 ,  1,  1, 50};
   ClearCuts(); AddSegmentCut(side,1,xmin,xmax);
}
///_________________________________________________________________________
void EdbRunAccess::SetCutBottom(int side, float wmin)
{
   float xmin[5]={ eViewXmin[side] + OverlapX(side), eViewYmin[side]                , -1, -1, wmin};
   float xmax[5]={ eViewXmax[side] - OverlapX(side), eViewYmin[side]+ OverlapY(side),  1,  1, 50};
   ClearCuts(); AddSegmentCut(side,1,xmin,xmax);
}

///_________________________________________________________________________
bool EdbRunAccess::InitRunFromRWC(char *rwcname, bool bAddRWD, const char* options)
{
  if(!eRun)
    if( gSystem->AccessPathName(rwcname, kFileExists) ) {
      Log(1,"EdbRunAccess::InitRunFromRWC","ERROR open file: %s\n",rwcname);
      return false;
    }
  eRun=new EdbRun(eRunFileName.Data(),"RECREATE");
  if(!eRun) return false;
  
  AddRWC(eRun, rwcname, bAddRWD, options); 
  
  return true;
}
///_________________________________________________________________________
bool EdbRunAccess::AddRWDToRun(char *rwdname, const char* options)
{
  Log(2,"EdbRunAccess::AddRWDToRun"," %s\n",rwdname);
  if(!eRun) return false;
  if( gSystem->AccessPathName(rwdname, kFileExists) ) {
      Log(1,"EdbRunAccess::AddRWDToRun","ERROR open file: %s\n",rwdname);
      return false;
  }

  int f;//fragment index
  f=eRun->GetHeader()->GetNareas();
  f++;
  if( !AddRWD( eRun, rwdname, f, options) ) return false;  
  eRun->GetHeader()->SetNareas(f);
  return true;
}

///______________________________________________________________________________
EdbLayer *EdbRunAccess::GetMakeLayer(int id)
{
  if(id<0) return 0;  if(id>2) return 0;
  if(!GetLayer(id))  eLayers[id] = new EdbLayer();
  return GetLayer(id);
}

///_________________________________________________________________________
int EdbRunAccess::GetViewsXY( int ud, TArrayI &entr, 
			      float x, float y, float rv)
{
  // return all views in rv-surrounds of the x-y point
  int nv=0;
  if(ud<0||ud>2) return nv;
  EdbPattern *pat=GetVP(ud);
  if(!pat)                             return nv;
  EdbSegP *seg=0;
  float r;
  for( int i=0; i<pat->N(); i++ ) {
    seg = pat->GetSegment(i);
    r = Sqrt( (seg->X()-x)*(seg->X()-x)+(seg->Y()-y)*(seg->Y()-y));
    if( r>rv ) continue;
    entr[nv++] = seg->ID();
  }
  return nv;
}


///_________________________________________________________________________
int EdbRunAccess::GetEntryXY(int ud, float x, float y)
{
  // find view with the center closest to x,y and return it's entry

  int entry = -1;
  EdbPattern *pat=GetVP(ud);
  if(!pat)                             return -1;
  EdbSegP *seg=0;
  float r,rmin=9999999.;
  int imin=-1;
  for( int i=0; i<pat->N(); i++ ) {
    seg = pat->GetSegment(i);
    r = Sqrt( (seg->X()-x)*(seg->X()-x)+(seg->Y()-y)*(seg->Y()-y));
    if( r>rmin ) continue;
    rmin = r;
    entry = seg->ID();
    imin = i;
  }
  //if(imin>0) pat->GetSegment(imin)->Print();
  return entry;
}

///_________________________________________________________________________
int EdbRunAccess::GetVolumeXY(EdbSegP &s, EdbPatternsVolume &vol)
{
  // assuming s as the prediction in plate RS (also Z), fill vol with the segments 
  // of closest up and down views
  // todo: Z of layers must be setted!

  TArrayI entr(2);
  float x,y; 
  x = s.X()+ (eLayers[1]->Z()-s.Z())*s.TX();
  y = s.Y()+ (eLayers[1]->Z()-s.Z())*s.TY();
  entr[0] = GetEntryXY(1,x,y);
  //printf(":GetVolumeXY: 1 %f %f \n",x,y);

  x = s.X()+ (eLayers[2]->Z()-s.Z())*s.TX();
  y = s.Y()+ (eLayers[2]->Z()-s.Z())*s.TY();
  entr[1] = GetEntryXY(2,x,y);
  //printf(":GetVolumeXY: 2 %f %f \n",x,y);

  //printf("EdbRunAccess::GetVolumeXY: entr = %d %d\n",entr[0],entr[1]);
  int nrej=0;
  return GetVolumeData(vol,2,entr,nrej);
}

///_________________________________________________________________________
int EdbRunAccess::GetPatternXYcut(EdbSegP &s0, int side, EdbPattern &pat, float dr, float dt)
{
  Log(2,"EdbRunAccess::GetPatternXYcut","side %d     [%f %f %f %f] +-%5.1f +-%5.3f:", side, s0.X(), s0.Y(), s0.TX(),s0.TY(),dr,dt );

  EdbPattern p0;
  int nseg = GetPatternXY(s0, side, p0, dr);
  Log(2,"EdbRunAccess::GetPatternXYcut","nseg = %d",nseg);
  Log(2,"EdbRunAccess::GetPatternXYcut","side = %d zlayer=%f zseg= %f",side,eLayers[side]->Z(),s0.Z() );
  float x = s0.X()+ (eLayers[side]->Z()-s0.Z())*s0.TX();
  float y = s0.Y()+ (eLayers[side]->Z()-s0.Z())*s0.TY();
  int ic=0;
  for(int i=0; i<nseg; i++)
  {
    EdbSegP *s = p0.GetSegment(i);
    if( Abs(s->X()-x) > dr) continue;
    if( Abs(s->Y()-y) > dr) continue;
    if( Abs(s->TX()-s0.TX()) > dt) continue;
    if( Abs(s->TY()-s0.TY()) > dt) continue;
    pat.AddSegment(*s);
    ic++;
  }
  Log(2,"EdbRunAccess::GetPatternXYcut","selected = %d",ic);
  return ic;
}

///_________________________________________________________________________
int EdbRunAccess::GetPatternXY(EdbSegP &s, int side, EdbPattern &pat, float rmin)
{
  // assuming s as the prediction in plate RS (also Z), fill vol with the segments 
  // of closest up and down views
  // todo: Z of layers must be setted!

  if(side<1||side>2) return 0;
  TArrayI entr(100000);
  float x = s.X()+ (eLayers[side]->Z()-s.Z())*s.TX();
  float y = s.Y()+ (eLayers[side]->Z()-s.Z())*s.TY();
  int nv = GetViewsXY(side,entr,x,y, Max(rmin,(float)300.) );
  pat.SetZ(eLayers[side]->Z());
  int nrej=0;
  return GetPatternData(pat,side,nv,entr,nrej);
}


///_________________________________________________________________________
int EdbRunAccess::GetVolumeArea(EdbPatternsVolume &vol, int area)
{
  // get raw segments as a volume for the given "area"

  int maxA=10000;            //TODO
  TArrayI entr1(maxA);
  float xmin,xmax,ymin,ymax;
  int n1 = GetViewsArea( 1, entr1, area, xmin,xmax,ymin,ymax);
  TArrayI entr2(maxA);
  xmin -= eXstep[1]+eXstep[1]/10.;
  xmax += eXstep[1]+eXstep[1]/10.;
  ymin -= eYstep[1]+eYstep[1]/10.;
  ymax += eYstep[1]+eYstep[1]/10.;
  int n2 = GetViewsArea( 2, entr2, xmin,xmax,ymin,ymax );
  //  int n2 = GetViewsAreaMarg( 2, entr2, area, eXstep+eXstep/10.,eYstep+eYstep/10. );

  printf("n1=%d n2=%d\n",n1,n2);

  TArrayI eall(n1+n2);
  int ic=0;
  for(int i=0; i<n1; i++)    eall[ic++]=entr1[i];
  for(int i=0; i<n2; i++)    eall[ic++]=entr2[i];

  TArrayI ind(ic);
  TArrayI srt(ic);
  TMath::Sort(ic,eall.GetArray(),ind.GetArray(),0);
  for (int i=0; i<ic; i++)    srt[i] = eall[ind[i]];

  int nrej=-1;
  int nseg = GetVolumeData(vol, ic, srt, nrej);

#ifdef _WINDOWS
  Log(2,"GetVolumeArea","Area:%3d (%3.0f %%)  views:%4d/%4d   %6d +%6d segments are accepted; %6d - rejected"
#else
  Log(2,"GetVolumeArea","Area:%3d (%3.0f \%%)  views:%4d/%4d   %6d +%6d segments are accepted; %6d - rejected"
#endif
        ,area, 100.*area/eNareas
        ,n1,n2
	,vol.GetPattern(0)->N()
	,vol.GetPattern(1)->N()
	,nrej );
	 
  if(nseg != vol.GetPattern(0)->N()+vol.GetPattern(1)->N()) 
  Log(2,"GetVolumeArea","WARNING!!! nseg = %d != %d\n", nseg,vol.GetPattern(0)->N()+vol.GetPattern(1)->N());
  return ic;
}

///_________________________________________________________________________
int EdbRunAccess::ViewSide(const EdbView *view) const
{
  int side=0;
  if(     view->GetNframesTop()==0) side=2;         // 2- bottom
  else if(view->GetNframesBot()==0) side=1;         // 1- top
  else return 0;
  if(eInvertSides) side = 3-side;
  return side;
}

///_________________________________________________________________________
int EdbRunAccess::GetVolumeData(EdbPatternsVolume &vol, 
                                int nviews, 
                                TArrayI &srt, 
                                int &nrej)
{
  // get raw segments as a volume for the given array of entries

  EdbSegP    segP;
  EdbView    *view = eRun->GetView();
  int   side;
  int   nseg=0;
  nrej=0;
  int   entry;
  int   nsegV;
  
  for(int iu=0; iu<nviews; iu++) {
    entry = srt[iu];
    if(eCLUST)       {
      view = eRun->GetEntry(entry,1,1,1);
      view->AttachClustersToSegments();
    }
    else view = eRun->GetEntry(entry);

    nsegV = view->Nsegments();
    side = ViewSide(view);

    //printf(" EdbRunAccess::GetVolumeData: %d nsegV = %d side=%d\n",iu, nsegV,side);
    if(side<1||side>2) continue;

    //vol.Print();
    for(int j=0;j<nsegV;j++) {
      if(!AcceptRawSegment(view,j,segP,side,entry)) {
        nrej++;
        continue;
      }
      nseg++;
      vol.GetPattern(side-1)->AddSegment( segP);
    }
  }

  //printf("nseg: %d \t rejected: %d\n", nseg, nrej );

  vol.GetPattern(0)->SetSegmentsZ();
  vol.GetPattern(1)->SetSegmentsZ();
  return nseg;
}

///_________________________________________________________________________
int EdbRunAccess::GetPatternData( EdbPattern &pat, int side,
				  int nviews, 
				  TArrayI &srt, 
				  int &nrej    )
{
  // get raw segments as a pattern for the given array of entries and given side

  EdbSegP    segP;
  EdbView    *view = eRun->GetView();
  int   nseg=0;
  nrej=0;
  int   entry;
  int   nsegV;
  
  for(int iu=0; iu<nviews; iu++) {
    entry = srt[iu];
    if(eCLUST)       {
      view = eRun->GetEntry(entry,1,1,1);
      view->AttachClustersToSegments();
    }
    else view = eRun->GetEntry(entry);

    nsegV = view->Nsegments();
    if( ViewSide(view) != side )   continue;

    for(int j=0;j<nsegV;j++) {
      if(!AcceptRawSegment(view,j,segP,side,entry)) {
	nrej++;
	continue;
      }
      nseg++;
      pat.AddSegment( segP);
    }
  }
  return nseg;
}

///_________________________________________________________________________
int EdbRunAccess::GetPatternView( EdbPattern &pat, int side, int entry, int &nrej )
{
  // get raw segments as a pattern for the given entry in a Views tree
  
  EdbSegP    segP;
  EdbView    *view = eRun->GetView();
  int   nseg=0;
  nrej=0;
  
  if(eCLUST)       {
    view = eRun->GetEntry(entry,1,1,1);
    view->AttachClustersToSegments();
  }
  else 
    view = eRun->GetEntry(entry);
  
  int nsegV = view->Nsegments();
  
  for(int j=0;j<nsegV;j++) {
    if(!AcceptRawSegment(view,j,segP,side,entry)) {
      nrej++;
      continue;
    }
    nseg++;
    pat.AddSegment( segP);
  }
  
  return nseg;
}

///_________________________________________________________________________
int EdbRunAccess::GetPatternDataForPrediction( int id, int side, EdbPattern &pat )
{
  // get raw segments belonging to views having a given id (view->GetHeader()->GetTrack()) and a given side
  // the result is stored into an EdbPattern object
  // this routine is called by EdbRunTracking::FindCandidates and to select microtracks from a raw.root file 
  // acquired for a specific prediction
  // Special case: if id=-1 - accept all views for the given side

  EdbSegP    segP;
  EdbView    *view = eRun->GetView();
  //int nviews = eRun->GetEntries();
  int nviews = eVP[side]->N();
  Log(2,"EdbRunAccess::GetPatternDataForPrediction","%d views are selected for side %d",nviews,side);
  int   nseg=0;
  int   nsegV;
  
  for(int ie=0; ie<nviews; ie++) {
    int entry=eVP[side]->GetSegment(ie)->ID();
    if(eCLUST)       {
      view = eRun->GetEntry(entry,1,1,1);
      view->AttachClustersToSegments();
    }
    else view = eRun->GetEntry(entry);

    if( id>=0 && view->GetHeader()->GetTrack()!=id) continue;

    nsegV = view->Nsegments();
    if( ViewSide(view) != side )   continue;
    //Log(2,"EdbRunAccess::GetPatternDataForPrediction","ie = %d   nsegV = %d",ie,nsegV);

    for(int j=0;j<nsegV;j++) {
      if(!AcceptRawSegment(view,j,segP,side,entry))   continue;
      nseg++;
      segP.SetScanID(pat.ScanID());
      segP.SetSide(pat.Side());
      pat.AddSegment( segP);
    }
  }
  Log(2,"GetPatternDataForPrediction","%d segments for pred %d, side %d from %s",nseg,id,side, eRunFileName.Data());
  return nseg;
}



///_________________________________________________________________________
int EdbRunAccess::GetViewsArea(int ud, TArrayI &entr, int area,
			       float &xmin, float &xmax, float &ymin, float &ymax )
{
  // get all views of "area" for the side "ud"; fill entr with the entries, 
  // calculate area limits

  int nv=0;
  if(ud<0) return nv;  if(ud>2) return nv;

  EdbPattern *pat=GetVP(ud);
  if(!pat)                             return nv;
  xmin=eXmax;  xmax=eXmin;
  ymin=eYmax;  ymax=eYmin;
  EdbSegP *seg=0;
  for( int i=0; i<pat->N(); i++ ) {
    seg = pat->GetSegment(i);
    if( seg->Aid(0) != area ) continue;
    if(xmin>seg->X()) xmin=seg->X();
    if(xmax<seg->X()) xmax=seg->X();
    if(ymin>seg->Y()) ymin=seg->Y();
    if(ymax<seg->Y()) ymax=seg->Y();
    entr[nv++]=seg->ID();
  }  
  return nv;
}

///_________________________________________________________________________
int EdbRunAccess::GetViewsArea(int ud, TArrayI &entr, 
			       float xmin, float xmax, float ymin, float ymax )
{
  // get all views in a given limits

  int nv=0;
  if(ud<0) return nv;  if(ud>2) return nv;

  EdbPattern *pat=GetVP(ud);
  if(!pat)                             return nv;
  EdbSegP *seg=0;
  for( int i=0; i<pat->N(); i++ ) {
    seg = pat->GetSegment(i);
    if(seg->X()<xmin)     continue;
    if(seg->X()>xmax)     continue;
    if(seg->Y()<ymin)     continue;
    if(seg->Y()>ymax)     continue;
    entr[nv++]=seg->ID();
  }
  return nv;
}

///_________________________________________________________________________
int EdbRunAccess::GetViewsAreaMarg(int ud, TArrayI &entr, int area, 
				   float xmarg, float ymarg)
{
  // return area views with margins
  // ud: 0-base, 1-down, 2-up         //TODO - to check this convention!

  int nv=0;
  if(ud<0 || ud>2) return nv;
  EdbPattern *pat=GetVP(ud);
  if(!pat)  return nv;

  EdbPattern ptmp;
  EdbSegP *seg=0;
  int N=pat->N();
  for( int i=0; i<N; i++ ) {
    seg = pat->GetSegment(i);
    if( seg->Aid(0) == area )  ptmp.AddSegment( *seg );
  }
  
  float xmin = ptmp.Xmin();
  float xmax = ptmp.Xmax();
  float ymin = ptmp.Ymin();
  float ymax = ptmp.Ymax();

  for( int i=0; i<N; i++ ) {
    seg = pat->GetSegment(i);
    if( seg->X() < xmin-xmarg ) continue;
    if( seg->X() > xmax+xmarg ) continue;
    if( seg->Y() < ymin-ymarg ) continue;
    if( seg->Y() > ymax+ymarg ) continue;
    entr[nv++] = seg->ID();
  }

  return nv;
}

///_________________________________________________________________________
void EdbRunAccess::PrintStat()
{
  if(!eRun) return;

  printf("\nSome run statistics:\n");
  printf(  "--------------------\n");
  int nviews=eRun->GetEntries();
  EdbRunHeader *header = eRun->GetHeader();
  if(header) eNareas=header->GetNareas();
  printf("entries    : %d\n", nviews );
  printf("areas      : %d   in the range: ( %d : %d )\n", 
	 eNareas, eFirstArea,eLastArea );
  int var=0;
  if(eNareas>0) {
    var =  nviews/eNareas;
    printf("views/area : %d  (%d/side)\n", var, var/2 );
  }
  if(var*eNareas != nviews) printf("\t WARNING: areas are not equal\n");

  for(int ud=0; ud<3; ud++) {
    EdbPattern *pat=GetVP(ud);
    int nempty=0;
    if(pat) {
      printf( "side:%d  xmin,xmax  = %f %f   ymin,ymax = %f %f\n", 
	      ud,
	      pat->Xmin(), 
	      pat->Xmax(),
	      pat->Ymin(),
	      pat->Ymax()
	      );
      nempty = CheckEmptyViews(*pat);
      if(nempty>0) printf("\t WARNING: %d empty views in layer %d\n",nempty,ud);
    }
  }
  printf("\n");
}

///_________________________________________________________________________
void EdbRunAccess::CheckRunLine()
{
  if(!eRun) return;

  cerr << endl << eRunFileName.Data();
  int nviews=eRun->GetEntries();
  cerr<<"\tv: "<< nviews;
  eNareas=0;
  if(eRun->GetHeader()) {
    eNareas =eRun->GetHeader()->GetNareas();
  } else 
    cerr <<"  RUN HEADER is missing!";
  
  cerr<<"  a: "<<eNareas;
  int var=0;
  if(eNareas>0) {
    var =  nviews/eNareas/2;
    cerr<<"   v/a/s: "<<var;
  }
  cerr<<endl;
  
  if(eNareas>0) 
    if(var*2*eNareas != nviews) {
      cerr<<"\tWARNING: areas and views numbers are mismatching!\n";
      for(int ud=0; ud<3; ud++) {
	EdbPattern *pat=GetVP(ud);
	if(pat) {
	  int nv = pat->N();
	  cerr<<"\t\t views("<<ud<<"): "<<nv;
	  if(nv>0) {
	    cerr<<"  Areas in range: ";
	    cerr<<pat->GetSegment(0)->Aid(0)<<" : "<<pat->GetSegment(nv-1)->Aid(0);
	    cerr<<"    Views in range: ";
	    cerr<<pat->GetSegment(0)->Aid(1)<<" : "<<pat->GetSegment(nv-1)->Aid(1);
	  }
	  cerr<<endl;
	}
      }
    }
  
  for(int ud=0; ud<3; ud++) {
    EdbPattern *pat=GetVP(ud);
    int nempty=0;
    if(pat) {
      nempty = CheckEmptyViews(*pat);
      if(nempty>0) cerr<<"\tWARNING: "<<nempty<<" empty views in layer "<<ud<<endl;
      float meanseg = CheckMeanSegsPerView(*pat);
      printf("Mean Segments/view = %10.1f\n",meanseg);
      int thres=1;
      int n0v = Check0Views(*pat, thres);
      if(n0v) printf("%d views with <%d segments!\n",n0v,thres);
    }
  }
}

///_________________________________________________________________________
float EdbRunAccess::CheckMeanSegsPerView(EdbPattern &pat)
{
  Long_t nseg=0, nview=pat.N();
  if(!nview)                       return 0;
  for(int i=0; i<nview; i++)   nseg += pat.GetSegment(i)->Flag();
  return 1.*nseg/nview;
}
///_________________________________________________________________________
int EdbRunAccess::Check0Views(EdbPattern &pat, int thres)
{
  Long_t n0=0, nview=pat.N();
  if(!nview)                       return 0;
  for(int i=0; i<nview; i++)   if( pat.GetSegment(i)->Flag() < thres) n0++;
  return n0;
}

/*  TODO
///_________________________________________________________________________
int EdbRunAccess::CheckViewsOccupancy(EdbPattern &pat)
{
  TH1F h("","", 200,0., 10000.);
  for(int i=0; i<pat.N(); i++)   h.Fill( pat.GetSegment(i)->Flag() );
  
  return nempty;
}
*/

///_________________________________________________________________________
int EdbRunAccess::CheckEmptyViews(EdbPattern &pat)
{
  int nempty=0;
  EdbSegP *s=0;
  for(int i=0; i<pat.N(); i++) {
    s = pat.GetSegment(i);
    if( s->Vid(0)>0 )    continue;
    if( s->Vid(1)>0 )    continue;
    nempty++;
  }
  return nempty;
}

//_________________________________________________________________________
void EdbRunAccess::CheckViewStep()
{
  CheckViewStep(1);
  CheckViewStep(2);
}

///_________________________________________________________________________
void EdbRunAccess::CheckViewStep(int ud)
{
  // assumed that in the tree (and in eVP) views are consecutive in the order of acquisition
  if(ud<0 || ud>2)                        return;
  EdbPattern *pat=GetVP(ud);    if(!pat)  return;

  TH1F hx("viewXstep","viewXstep",3000,200,500);
  TH1F hy("viewYstep","viewYstep",3000,200,500);
  int n=pat->N();               if(n<2)   return;
  for( int i=0; i<n-1; i++ ) {
    EdbSegP *s1 = pat->GetSegment(i);
    EdbSegP *s2 = pat->GetSegment(i+1);
    if( s1->Aid(0) != s1->Aid(0) )  continue;   // skip different areas
    if( s1->Aid(1) == s2->Aid(1) )  continue;   // skip same view
    float dx = Abs(s2->X()-s1->X());
    float dy = Abs(s2->Y()-s1->Y());
    if(dx>dy) 
      if(dx>50&&dx<500&&dy<50)   //check that steps are in reasonable limits
        hx.Fill(dx);
    if(dy>dx) 
      if(dy>50&&dy<500&&dx<50)   //check that steps are in reasonable limits
        hy.Fill(dy);
  }
  //hx.Write(); hy.Write(); f.Close();
  eXstep[ud] = hx.GetMean();
  eYstep[ud] = hy.GetMean();
  Log(2,"EdbRunAccess::CheckViewStep","side %d:   stepX(%d) = %f +- %f    stepY(%d) = %f +- %f",
      ud, (int)(hx.GetEntries()), hx.GetMean(),hx.GetRMS(), (int)(hy.GetEntries()), hy.GetMean(),hy.GetRMS() );
}

///_________________________________________________________________________
int EdbRunAccess::FillVP()
{
  // fill patterns up/down with dummy segments with center of view coordinates 
  // to profit of the pattern tools for the fast views selection
  // the convention is: 
  // segp->ID() - tree entry number
  // segp->X() - coordinates in the external (brick) RS (after transformation)
  // segp->Y() - /
  // segp->Aid(areaID,viewID);

  
  Log(2,"EdbRunAccess::FillVP","Get data from %s with AFID=%d",eRunFileName.Data(),eAFID);
  
  if(!eRun) { Log(1,"EdbRunAccess::FillVP","ERROR: run is not opened\n"); return 0; }
  EdbView        *view = eRun->GetView();
  EdbViewHeader  *head = view->GetHeader();

  int nentr = eRun->GetEntries();

  eRun->GetTree()->Draw(">>lst", eHeaderCut );
  TEventList *lst = (TEventList*)(gDirectory->GetList()->FindObject("lst"));
  if(!lst) {Log(1,"EdbRunAccess::FillVP","ERROR!: events list (lst) did not found!"); return 0;}
  nentr =lst->GetN();
  if(nentr<1) {Log(1,"EdbRunAccess::FillVP","%d views accepted by cut %s", nentr, eHeaderCut.GetTitle()); return 0;}

  if( eVP[1] ) delete  eVP[1];
  if( eVP[2] ) delete  eVP[2];
  eVP[1] = new EdbPattern( 0.,0.,   0., nentr);
  eVP[2] = new EdbPattern( 0.,0., 200., nentr);  //TODO

  EdbSegP segP;

  eFirstArea = 9999999;
  eLastArea  = 0;

  int side;
  int nseg,ncl;
  for(int ie=0; ie<nentr; ie++ ) {
    int iv = lst->GetEntry(ie);
    view = eRun->GetEntry(iv,1,0,0,0,0);
    
    nseg = head->GetNsegments();
    ncl  = head->GetNclusters();
    if(eAFID==0){
      segP.Set( iv,view->GetXview(),view->GetYview(), 0,0,ncl,nseg);
    }
    else if(eAFID==1 || eAFID==2) {
      segP.Set( iv,0,0, 0,0,ncl,nseg);
      segP.Transform( head->GetAffine() );
    }
    else if(eAFID==11) {
      segP.Set( iv,0,0, 0,0,ncl,nseg);
      if(!eAffStage2Abs) Log(1,"","eAffStage2Abs =0");
      segP.Transform( eAffStage2Abs );
    }
    else if(eAFID==20) {  //to remove
      segP.Set( iv,0,0, 0,0,ncl,nseg);
      EdbAffine2D *aff = (EdbAffine2D*)eViewCorr->At(iv);
      if(aff) {
        segP.Transform( aff );
        head->GetAffine()->Print();
        aff->Print();
      } else Log(1,"EdbRunAccess::FillVP", "WARNING! aff is missing for entry %d",iv);
    }
    if( view->GetAreaID() < eFirstArea )  eFirstArea= view->GetAreaID();
    if( view->GetAreaID() > eLastArea )   eLastArea = view->GetAreaID();
    
    segP.SetAid( view->GetAreaID()    , view->GetViewID() );
    segP.SetVid( head->GetNframesTop(), head->GetNframesBot() );

    side = ViewSide(view);
    if(side<1||side>2) continue;
    eVP[side]->AddSegment(segP);
  }
  Log(2,"EdbRunAccess::FillVP"," %d entries in limits X=(%f : %f) Y=(%f %f)\n", nentr,
      eVP[1]->Xmin(), eVP[1]->Xmax(), 
      eVP[1]->Ymin(), eVP[1]->Ymax() );
  return nentr;
}

///______________________________________________________________________________
void EdbRunAccess::ReadImageMatrixCorrection(int side, const char *file)
{
  ReadImageMatrixCorrection(eCorrMap[side], file);
}

///______________________________________________________________________________
void EdbRunAccess::ReadImageMatrixCorrection(EdbCell2 &map, const char *file)
{
  Log(1,"EdbRunAccess::ReadImageMatrixCorrection","from %s ",file);
  FILE *f = fopen(file,"r");
  if(!f) { Log(1,"EdbRunAccess::ReadImageMatrixCorrection","file %s not found!",file); return; }
  char str[256]; 
  fgets(str, 256, f);
  int nxpix, nypix;
  float xpix, ypix;
  if( 4 != sscanf(str,"%d %d  %f %f", &nxpix, &nypix, &xpix, &ypix)) return;
  
  float  mi[2] = { -nxpix*Abs(xpix)/2., -nypix*Abs(ypix)/2. };
  float  ma[2] = {  nxpix*Abs(xpix)/2.,  nypix*Abs(ypix)/2. };
  int    n[2] = { nxpix, nypix };
  map.InitCell(2, n, mi, ma);
  int nc=map.Ncell();
  for( int i=0; i<nc; i++ ) {
    TArrayD *a = new TArrayD(2);
    map.AddObject(i, (TObject*)a );
    map.SetBin(i,0);
  }
  
  int i0,j0;
  float x0,y0,dx0,dy0;
  for(int i=0; i<nxpix; i++) {
    for(int j=0; j<nypix; j++) {
      if(fgets(str, 256, f)==NULL) Log(1,"EdbRunAccess::ReadImageMatrixCorrection","ERROR: input file too short!");
      if( sscanf(str, "%d %d %f %f %f %f", &i0, &j0, &x0, &y0, &dx0, &dy0) == 6 )  
      {
	TArrayD *dxy = (TArrayD *)map.GetObject(x0,y0,0);
	(*dxy)[0]+=dx0;
	(*dxy)[1]+=dy0;
	map.Fill(x0, y0);
      }
    }
  }
  fclose(f);

  nc=map.Ncell();
  for( int i=0; i<nc; i++ ) {
    TArrayD *dxy = (TArrayD*)( map.GetObject(i, 0 ) );    if(!dxy) continue;
    int w = map.Bin(i);                                   if(!w) continue;
    (*dxy)[0] /= w;
    (*dxy)[1] /= w;
  }  
}

///______________________________________________________________________________
void EdbRunAccess::SetImageCorrection(int side, const char *str)
{
  // fedra.link.ImageCorrSideX:                 scaleX scaleY rot
  float scaleX, scaleY, rot;
  if( 3 == sscanf(str,"%f %f %f", &scaleX, &scaleY, &rot) )
  { 
    if(!(side>0&&side<3)) Log(1,"EdbRunAccess::SetImageCorrection:", "side %d out of range!", side);
    else
    {
      eImageCorr[side].ZoomX(scaleX);
      eImageCorr[side].ZoomY(scaleY);
      eImageCorr[side].Rotate(rot);
    }
  }
  Log(2,"EdbRunAccess::SetImageCorrection","%s",str);
}

///______________________________________________________________________________
void EdbRunAccess::SetImageCorrectionAff(int side, const char *str)
{
  // fedra.link.ImageCorrSideX:                 a11 a12 a21 a22 b1 b2
  if(!(side>0&&side<3)) 
    Log(1,"EdbRunAccess::SetImageCorrection:", "side %d out of range!", side);
  else
    eImageCorr[side].Set(str);
  Log(2,"EdbRunAccess::SetImageCorrectionAff","%s",str);
}

///______________________________________________________________________________
float EdbRunAccess::SegmentWeight(const EdbSegment &s)
{
  if     (eWeightAlg==1) return Sqrt(s.GetPuls()*s.GetSigmaY()/2.);
  else if(eWeightAlg==2) return s.GetSigmaX();
  return s.GetPuls();
}

///______________________________________________________________________________
void EdbRunAccess::ApplyImageMatrixCorr(int side, float &x, float &y)
{
  TArrayD *dxy=0;
  if(side>0&&side<3)
    if( dxy = (TArrayD*)eCorrMap[side].GetObject(x,y,0) )
    {
      x -= (*dxy)[0];
      y -= (*dxy)[1];
    }
}

///______________________________________________________________________________
bool EdbRunAccess::AcceptRawSegment(EdbView *view, int id, EdbSegP &segP, int side, int entry)
{
  EdbSegment *seg = view->GetSegment(id);

  if( !PassCuts(side,*seg) )     return false;
  
  if(eTracking>=0)
    if(((Int_t)(seg->GetUniqueID()>>16))!=eTracking) return false;

  float pix, chi2;
  if(eCLUST) {
    if(eUseExternalSurface) {
      SetSegmentAtExternalSurface(seg,side);
    } else {
      pix = GetRawSegmentPix(seg);
      segP.SetVolume( pix );
      chi2 = CalculateSegmentChi2( seg,
				   GetCond(1)->SigmaXgr(),  //TODO: side logic
				   GetCond(1)->SigmaYgr(), 
				   GetCond(1)->SigmaZgr());
      if(chi2>GetCutGR())  return false;
      segP.SetChi2( chi2 );
    }
  }

  if(eDoViewAnalysis) {
    eHViewXY[side].Fill( seg->GetX0(), seg->GetY0());
  }

  EdbLayer  *layer=GetLayer(side);
  const EdbAffine2D *affview = 0;
  if      (eAFID==1) affview = view->GetHeader()->GetAffine();
  else if (eAFID==11) affview = eAffStage2Abs;
  else if(eAFID==2) {
    affview = (EdbAffine2D*)eViewCorr->At(entry);
    if(id==0) {
      view->GetHeader()->GetAffine()->Print();
      affview->Print();
    }
  }
      
  float x,y,z,tx,ty,puls;
  tx   = seg->GetTx()/layer->Shr();
  ty   = seg->GetTy()/layer->Shr();
  x    = seg->GetX0();
  y    = seg->GetY0();

  if(eDoImageCorr) 
  {
    float xt = eImageCorr[side].Xtrans(x,y);
    float yt = eImageCorr[side].Ytrans(x,y);
    x = xt;
    y = yt;  
  }

  if(eDoImageMatrixCorr)  ApplyImageMatrixCorr(side,x,y);

  x  += layer->Zcorr()*tx;
  y  += layer->Zcorr()*ty;
  z    = layer->Z() + layer->Zcorr();
  
  float xx, yy, txr, tyr;
  if( eAFID==0 ) {
    xx = x+view->GetXview();
    yy = y+view->GetYview();
    txr = tx;
    tyr = ty;
  } 
  else if(eAFID==11) {
    x += view->GetXview();
    y += view->GetYview();
    xx  = affview->A11()*x+affview->A12()*y+affview->B1();
    yy  = affview->A21()*x+affview->A22()*y+affview->B2();
    txr = affview->A11()*tx+affview->A12()*ty;                   // rotate also angles
    tyr = affview->A21()*tx+affview->A22()*ty;
  }
  else {
    //seg->Transform( affview );
    xx  = affview->A11()*x+affview->A12()*y+affview->B1();
    yy  = affview->A21()*x+affview->A22()*y+affview->B2();
    txr = affview->A11()*tx+affview->A12()*ty;                   // rotate also angles
    tyr = affview->A21()*tx+affview->A22()*ty;
  }
  
  if(eWeightAlg==2) puls = seg->GetPuls();           // use SigmaX for cuts only
  else              puls = SegmentWeight(*seg);

  EdbAffine2D *aff = layer->GetAffineTXTY();
  float txx = aff->A11()*txr+aff->A12()*tyr+aff->B1();
  float tyy = aff->A21()*txr+aff->A22()*tyr+aff->B2();
  
  segP.Set( seg->GetID(),xx,yy,txx,tyy,puls,0);
  segP.SetZ( z );
  segP.SetDZ( seg->GetDz()*layer->Shr() );
  segP.SetW( puls );
  segP.SetVolume( seg->GetVolume() );
  segP.SetChi2( seg->GetSigmaX() );
  segP.SetProb( seg->GetSigmaY() );       //test: grains density after LASSO
  segP.SetVid(entry,id);
  segP.SetAid(view->GetAreaID(),view->GetViewID(), side);
  
  return true;
}

///______________________________________________________________________________
bool  EdbRunAccess::PassCuts(int id, EdbSegment &seg)
{
  float var[5];
  var[0] = seg.GetX0();
  var[1] = seg.GetY0();
  var[2] = seg.GetTx();
  var[3] = seg.GetTy();
  var[4] = SegmentWeight(seg);

  int nc = NCuts(id);
  for(int i=0; i<nc; i++) {
    if( !(GetCut(id,i)->PassCut(var)) )  return false;
  }
  return true;
}

///______________________________________________________________________________
float  EdbRunAccess::CalculateSegmentChi2( EdbSegment *seg, float sx, float sy, float sz )
{
  //assumed that clusters are attached to segments
  double chi2=0;
  EdbCluster *cl=0;
  TObjArray *clusters = seg->GetElements();
  if(!clusters) return 0;
  int ncl = clusters->GetLast()+1;
  if(ncl<=0)     return 0;

  float xyz1[3], xyz2[3];             // segment line parametrized as 2 points
  float xyz[3];
  bool inside=true;

  xyz1[0] = seg->GetX0() /sx;
  xyz1[1] = seg->GetY0() /sy;
  xyz1[2] = seg->GetZ0() /sz;
  xyz2[0] = (seg->GetX0() + seg->GetDz()*seg->GetTx()) /sx;
  xyz2[1] = (seg->GetY0() + seg->GetDz()*seg->GetTy()) /sy;
  xyz2[2] = (seg->GetZ0() + seg->GetDz())              /sz;

  double d;
  for(int i=0; i<ncl; i++ ) {
    cl = (EdbCluster*)clusters->At(i);
    xyz[0] = cl->GetX()/sx;
    xyz[1] = cl->GetY()/sy;
    xyz[2] = cl->GetZ()/sz;
    d = EdbMath::DistancePointLine3(xyz,xyz1,xyz2, inside);
    chi2 += d*d;
  }

  return TMath::Sqrt(chi2/ncl);
}

///______________________________________________________________________________
EdbSegment  *EdbRunAccess::GetRawSegmentN( int vid, int sid, int rs)
{
  EdbView *view = eRun->GetEntry(vid);
  EdbSegment *s = view->GetSegment(sid);
//  printf("(%d %d)  %d  %f %f\n",vid,sid, s->GetPuls(), view->GetXview(),view->GetYview());
  return s;
}

///______________________________________________________________________________
void EdbRunAccess::ApplyCorrections( const EdbView &view, EdbSegment &s, const int rs )
{
  // rs=0 local view RS
  // rs==1 plate RS (use layer 1&2 information Note that the last emlink
  //       should be done with all corrections off to get the consistent data here
  // rs==1 brick RS - use layer 0 info

  if(rs==0)  return;    // do nothing - remain in local view RS
  if(rs>0)              // use layers information
  {
    if(eAFID) s.Transform( view.GetHeader()->GetAffine() );
    else {
      s.SetX0( s.GetX0()+view.GetXview());
      s.SetY0( s.GetY0()+view.GetYview()); 
    }
    EdbSegP sp(0,s.GetX0(),s.GetY0(),s.GetTx(),s.GetTy());
    sp.SetZ( GetLayer( ViewSide(&view) )->Zxy() ); //forget original Z here...
    GetLayer( ViewSide(&view) )->CorrectSeg(sp);   //get to plate RS
    if(rs==2)  {
      GetLayer(0)->CorrectSeg(sp);                 //get to brick RS
      sp.SetZ( sp.Z() + GetLayer(0)->Zxy() );
    }
    s.SetTx( sp.TX() );
    s.SetTy( sp.TY() );
    s.SetX0( sp.X() );
    s.SetY0( sp.Y() );
    s.SetZ0( sp.Z() );
  }
}

///______________________________________________________________________________
EdbSegment  *EdbRunAccess::GetRawSegment( EdbView &view, int sid, int rs )
{
  EdbSegment *s = view.GetSegment(sid);
  ApplyCorrections( view, *s, rs );
  return s;
}

///______________________________________________________________________________
EdbSegment  *EdbRunAccess::GetRawSegment( int vid, int sid, int rs )
{
  EdbView *view = eRun->GetEntry(vid);
  return GetRawSegment( *view, sid, rs );
}

///______________________________________________________________________________
float  EdbRunAccess::GetRawSegmentPix( EdbSegment *seg )
{
  //assumed that clusters are attached to segments
  float pix=0;
  EdbCluster *cl=0;
  TObjArray *clusters = seg->GetElements();
  if(!clusters) return 0;
  int ncl = clusters->GetLast()+1;
  for(int i=0; i<ncl; i++ ) {
    cl = (EdbCluster*)clusters->At(i);
    pix += cl->GetArea();
  }
  return pix;
}

///______________________________________________________________________________
bool  EdbRunAccess::SetSegmentAtExternalSurface( EdbSegment *seg, int side )
{
  // set the segment coordinates correspondent to the cluster coord with maximal Z
  // assumed that clusters are attached to segments
  TObjArray *clusters = seg->GetElements();
  if(!clusters) return false;
  int ncl = clusters->GetLast()+1;
  EdbCluster *cl=0, *clZ=(EdbCluster*)clusters->At(0);
  for(int i=0; i<ncl; i++ ) {
    cl = (EdbCluster*)clusters->At(i);
    if     (side==1)  { if( clZ->Z() < cl->Z() )  clZ=cl; }   //side=1: top     use point of max Z
    else if(side==2)  { if( clZ->Z() > cl->Z() )  clZ=cl; }   //side=2: bottom  use point of min Z    
  }
  seg->SetX0(clZ->X());
  seg->SetY0(clZ->Y());
  seg->SetZ0(clZ->Z());
  return true;
}

///______________________________________________________________________________
void EdbRunAccess::AddSegmentCut(int xi, const char *cutline)
{
  float var[10];
  int onoff=-1;
  if(sscanf(cutline,"%d %f %f %f %f %f %f %f %f %f %f",&onoff, 
      &var[0],&var[1], &var[2],&var[3], &var[4],&var[5], &var[6],&var[7], &var[8], &var[9] )==11) 
      if(onoff>-1)  {
        AddSegmentCut(1,xi,var);
        AddSegmentCut(2,xi,var);
      }
}
///______________________________________________________________________________
void EdbRunAccess::AddSegmentCut(int layer, int xi, float var[10])
{
  if(!eCuts[layer])  eCuts[layer] = new TObjArray();
  eCuts[layer]->Add( new EdbSegmentCut(xi,var) );
}

///______________________________________________________________________________
void EdbRunAccess::AddSegmentCut(int layer, int xi, float min[5], float max[5])
{
  float var[10]={ min[0],max[0], min[1],max[1], min[2],max[2], min[3],max[3], min[4],max[4] };
  AddSegmentCut(layer,xi,var);
}

///______________________________________________________________________________
bool EdbRunAccess::CopyRawDataXY( float x0, float y0, float dR, const char *file )
{
  Log(2,"CopyRawDataXY", "dR = %f around (%f %f) into file %s",dR, x0,y0, file);
  if(!eRun) return false;
  EdbRun r(*eRun, file);
  EdbView    *view = eRun->GetView();
  for(int side=1; side<=2; side++) {
    TArrayI entr(100000);
    int nv = GetViewsXY(side,entr,x0,y0, dR);
    for(int i=0; i<nv; i++) {
      view = eRun->GetEntry(entr[i]);
      r.AddView(view);
    }
    Log(2,"CopyRawDataXY", "%d views copied for side %d",nv, side);
  }
  r.Close();
  return true;
}

///______________________________________________________________________________
void EdbRunAccess::CheckViewSize()
{
  for(int ud=1; ud<3; ud++ ) {
    eViewXmin[ud] = eHViewXY[ud].XminA(0.1);
    eViewXmax[ud] = eHViewXY[ud].XmaxA(0.1);
    eViewYmin[ud] = eHViewXY[ud].YminA(0.1);
    eViewYmax[ud] = eHViewXY[ud].YmaxA(0.1);
    Log(2,"EdbRunAccess::CheckViewSize"," side%d    X:[%f %f]   Y:[%f %f]", ud,
     eViewXmin[ud],eViewXmax[ud],eViewYmin[ud],eViewYmax[ud] );
  }
}

///______________________________________________________________________________
TH2F *EdbRunAccess::CheckUpDownOffsets()
{
  TH2F *h = new TH2F("UpDownOffsets","dXdY offsets between Up and correspondent Down views",100,-10,10,100,-10,10 );
  EdbPattern *p1 = GetVP(1); if(!p1) { Log(1,"EdbRunAccess::CheckUpDownOffsets","Error: GetVP(1) empty"); return h;}
  EdbPattern *p2 = GetVP(2); if(!p1) { Log(1,"EdbRunAccess::CheckUpDownOffsets","Error: GetVP(2) empty"); return h;}
  int  N1 = p1->N(), N2 = p2->N();
  if(N1 != N2) Log(1,"EdbRunAccess::CheckUpDownOffsets","Warning: different number of views: %d -up,  %d - down", p1->N(), p2->N() );
  int n = TMath::Min( N1, N2 );

//segP.SetAid( view->GetAreaID()    , view->GetViewID() );

  TArrayL   viewid1(N1), viewid2(N2);       //view/area id to sincronize the views
  TArrayI      ind1(N1), ind2(N2);
  for(int i=0; i<N1; i++)     viewid1[i] = p1->GetSegment(i)->Aid(0) * 1000000 + p1->GetSegment(i)->Aid(1);
  for(int i=0; i<N2; i++)     viewid2[i] = p2->GetSegment(i)->Aid(0) * 1000000 + p2->GetSegment(i)->Aid(1);
  TMath::Sort(N1,viewid1.GetArray(),ind1.GetArray(),0);
  TMath::Sort(N2,viewid2.GetArray(),ind2.GetArray(),0);

  for(int i=0; i<n; i++) {
    EdbSegP *s1=p1->GetSegment(ind1[i]);
    EdbSegP *s2=p2->GetSegment(ind2[i]);
    h->Fill( s2->X()-s1->X(), s2->Y()-s1->Y() );
  }
  float sx = h->GetRMS(1);
  float sy = h->GetRMS(2);
  if(sx>0.5 || sy>0.5)
      Log(1,"EdbRunAccess::CheckUpDownOffsets","WARNING!!! XY stage possible malfunction (rms>0.5 micron) rmsX = %7.4f  rmsY = %7.4f (by %d views)",
	  h->GetRMS(1), h->GetRMS(2), n );
  else Log(2,"EdbRunAccess::CheckUpDownOffsets","rmsX = %7.4f  rmsY = %7.4f (by %d views)", sx, sy, n );
  return h;
}


///_________________________________________________________________________
void EdbRunAccess::FillDZMaps()
{
  // fill DZ maps

  int eNsegMin=100;
  
  Log(2,"EdbRunAccess::FillDZmaps","Get data from %s with AFID=%d", eRunFileName.Data(),eAFID);
  if(!eRun) { Log(1,"EdbRunAccess::FillDZMaps","ERROR: run is not opened\n"); return; }
  EdbView        *view = eRun->GetView();
  EdbViewHeader  *head = view->GetHeader();
  int nentries = eRun->GetEntries();

  Double_t *x1   = new Double_t[nentries];
  Double_t *y1   = new Double_t[nentries];
  Double_t *x2   = new Double_t[nentries];
  Double_t *y2   = new Double_t[nentries];
  Double_t *z1   = new Double_t[nentries];
  Double_t *z2   = new Double_t[nentries];
  Double_t *z3   = new Double_t[nentries];
  Double_t *z4   = new Double_t[nentries];
  Int_t cnt_good=0;

  int cnt1=0, cnt2=0;
  for(int ie=0; ie<nentries; ie++ ) {
    view = eRun->GetEntry(ie,1,0,0,0,0);
    if( head->GetNsegments() < eNsegMin ) continue;
    int side = ViewSide(view);
    if( side == 1 ) {
      x1[cnt1] = ViewX(view);
      y1[cnt1] = ViewY(view);
      z1[cnt1] = view->GetZ1();
      z2[cnt1] = view->GetZ2();
      cnt1++;
    }
    else if( side == 2 ) {
      x2[cnt2] = ViewX(view);
      y2[cnt2] = ViewY(view);
      z3[cnt2] = view->GetZ3();
      z4[cnt2] = view->GetZ4();
      cnt2++;
    }
  }

  for(int i=0; i<4; i++) if( eGraphZ[i] )  delete  eGraphZ[i];
  for(int i=0; i<3; i++) if( eGraphDZ[i] ) delete  eGraphDZ[i];

  eGraphZ[0]  = new TGraph2D("graphZ1","graphZ1",cnt1,x1,y1,z1);
  eGraphZ[1]  = new TGraph2D("graphZ2","graphZ2",cnt1,x1,y1,z2);
  eGraphZ[2]  = new TGraph2D("graphZ3","graphZ3",cnt2,x2,y2,z3);
  eGraphZ[3]  = new TGraph2D("graphZ4","graphZ4",cnt2,x2,y2,z4);

  std::unique_ptr<TFile> f(TFile::Open("dz.root","RECREATE"));
  if (!f || f->IsZombie())  Log(1,"EdbRunAccess::FillDZMaps","ERROR! can not open file dz.root");
  else                      for(int i=0; i<4; i++) if( eGraphZ[i] ) eGraphZ[i]->Write();
}
