//-- Author :  Valeri Tioukov   19.05.2002
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbPattern                                                           //
//                                                                      //
// Segments pattern                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TIndexCell.h"
#include "TPolyLine.h"
#include "EdbAffine.h"
#include "EdbPattern.h"
#include "EdbVertex.h"
#include "EdbPhys.h"
#include "EdbLog.h"
#include "vt++/CMatrix.hh"
#include "vt++/VtVector.hh"

ClassImp(EdbSegmentsBox)
ClassImp(EdbTrackP)
ClassImp(EdbPattern)
ClassImp(EdbPatternsVolume)

using namespace MATRIX;
using namespace TMath;

//______________________________________________________________________________
EdbSegmentsBox::EdbSegmentsBox(int nseg)
{
  if(nseg>0) eSegments = new TClonesArray("EdbSegP",nseg);
  else       eSegments = new TClonesArray("EdbSegP");
  Set0();
}
 
//______________________________________________________________________________
EdbSegmentsBox::EdbSegmentsBox(float x0, float y0, float z0, int nseg)
{
  if(nseg>0) eSegments = new TClonesArray("EdbSegP",nseg);
  else       eSegments = new TClonesArray("EdbSegP");
  eX=x0;  eY=y0;  eZ=z0;
  eDZkeep=0;
}
 
//______________________________________________________________________________
EdbSegmentsBox::~EdbSegmentsBox( )
{
  if(eSegments) {
    eSegments->Delete();
    delete eSegments;
    eSegments=0;
  }
}
 
//______________________________________________________________________________
void EdbSegmentsBox::Set0()
{
  eX=0;  eY=0;  eZ=0;
  eDZkeep=0;
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment(int id, float x, float y, float tx, float ty, 
			    float w, int flag)
{
  return new((*eSegments)[N()])  EdbSegP( id,x,y,tx,ty,w,flag );
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment( int i, EdbSegP &s )
{
  return new((*eSegments)[i])  EdbSegP( s );
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment( EdbSegP &s )
{
  return new((*eSegments)[N()])  EdbSegP( s );
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegment( EdbSegP &s1, EdbSegP &s2 )
{
  EdbSegP *s = new((*eSegments)[N()])  EdbSegP( s1 );
  s->MergeTo(s2);
  return s;
}
 
//______________________________________________________________________________
EdbSegP *EdbSegmentsBox::AddSegmentNoDuplicate( EdbSegP &s )
{
  // This function will add the segment s to the current list of
  // segments checking for duplicates. If any segment in the list is found,
  // that is compatible to this, the segment will _not_ be added,
  // instead a warning message will be printed (according to    
  // verbosity level)
  // This function might be slow. 
  
  //cout << endl  << endl;
  //cout <<"EdbSegmentsBox::AddSegmentNoDuplicate( EdbSegP &s )" << endl;
  //s.PrintNice();
  Int_t eSegmentsN = GetN();
  //cout << "Segments already in the eSegments list: " << eSegmentsN << endl;
  //cout << "Check now if segment s is already in the list..." << endl;
  
  EdbSegP* sComparison; 
  
  for (int i=0; i<eSegmentsN; ++i) {
    
    sComparison=(EdbSegP* )eSegments->At(i);
    //sComparison->PrintNice();
    //Bool_t EdbSegP::IsEqual(const TObject *obj) const
    //cout << "Evaluating   s.IsEqual(sComparison) = " << s.IsEqual(sComparison) << endl;
    if (s.IsEqual(sComparison)==kTRUE) {
      //cout << "Tested Segment is witin one sigma equal with one already in the list. DO NOT ADD THIS Segment..." << endl;
      return NULL;
    } 
    
    //bool EdbSegP::IsCompatible( EdbSegP &s, float nsigx, float nsigt ) const
    /*
     cout << "Evaluating   s.IsCompatible(*sComparison,1,1) = " << s.IsCompatible(*sComparison,1,1) << endl;
     if (s.IsCompatible(*sComparison,1,1)==kTRUE) {
      cout << "Tested Segment is witin one sigma compatible with one already in the list. DO NOT ADD THIS Segment..." << endl;
       return NULL;
     }
     */
    
  }

  //cout <<" EdbSegmentsBox::AddSegmentNoDuplicate( EdbSegP &s )   Segment s added." << endl;
  return new((*eSegments)[N()])  EdbSegP( s );
}
 
 
 
 
//______________________________________________________________________________
void EdbSegmentsBox::Reset()
{
  Set0();
  if(eSegments) eSegments->Clear();
}
 
//______________________________________________________________________________
float EdbSegmentsBox::DiffAff(EdbAffine2D *aff)
{
  EdbSegmentsBox p,p1;
  p.AddSegment(0,Xmin(),Ymin(),0.,0.);
  p.AddSegment(0,Xmax(),Ymin(),0.,0.);
  p.AddSegment(0,Xmin(),Ymax(),0.,0.);
  p.AddSegment(0,Xmax(),Ymax(),0.,0.);

  p1.AddSegment(0,Xmin(),Ymin(),0.,0.);
  p1.AddSegment(0,Xmax(),Ymin(),0.,0.);
  p1.AddSegment(0,Xmin(),Ymax(),0.,0.);
  p1.AddSegment(0,Xmax(),Ymax(),0.,0.);

  p1.Transform(aff);
  return p.Diff(p1);
}
 
//______________________________________________________________________________
float EdbSegmentsBox::Diff(EdbSegmentsBox &p)
{
  // return the mean difference beteween pattern elements
  int nseg = TMath::Min( N(), p.N() );

  if(nseg<1) return 0;

  EdbSegP *s1=0, *s2=0;

  float dx=0, dy=0;
  double sdx=0;
  for(int i=0; i<nseg; i++ ) {
    s1 =   GetSegment(i);
    s2 = p.GetSegment(i);
    dx = s2->X() - s1->X();
    dy = s2->Y() - s1->Y();
    sdx += TMath::Sqrt( dx*dx + dy*dy);
  }
  return sdx/nseg;
}

//______________________________________________________________________________
float EdbSegmentsBox::GetSize(Int_t XorY)
{
  // return the size of the segments box in X//Y
  // simple approach using max/min values of segments
  // works only if 2 or more segemts are in it.
  int nseg = N();
  if(nseg<2) return 0;
  EdbSegP *s1=0;

  float minX=0; float maxX=0;
  float minY=0; float maxY=0;
  for(int i=0; i<nseg; i++ ) {
    s1 =   GetSegment(i);
    if (i==0) { 
      minX=s1->X(); maxX=s1->X(); 
      minY=s1->Y(); maxY=s1->Y();
    }
    maxX=TMath::Max(maxX,s1->X());
    minX=TMath::Min(minX,s1->X());
    maxY=TMath::Max(maxY,s1->Y());
    minY=TMath::Min(minY,s1->Y());
  }
  float sizeX=TMath::Abs(maxX-minX);
  float sizeY=TMath::Abs(maxY-minY);
  
  if (XorY==1) { return sizeY; }
  else { return sizeX; }
}

//______________________________________________________________________________
float EdbSegmentsBox::GetSizeXY()
{
  return (GetSize(0)*GetSize(1));
}

//______________________________________________________________________________
float EdbSegmentsBox::GetTrackDensity()
{
  // Returns number of tracks per area (in microns usually)
  int nseg = N();
  if(nseg<2) return 0;
  float size=GetSizeXY();
  return float(nseg)/size;
}

//______________________________________________________________________________
float EdbSegmentsBox::GetTrackDensitymm2()
{
  // Returns number of tracks per area in millimeter squared
  // Original area is in microns x microns
  return GetTrackDensity()*1e6;
}
 
//______________________________________________________________________________
void EdbSegmentsBox::RemovePosMargins(float  x0, float  y0, float  dx, float  dy)
{
  int nseg = N();
  for(int i=nseg-1; i>=0; i-- )    {
    EdbSegP *s = GetSegment(i);
    if( s->eX > x0+dx || s->eX < x0-dx || s->eY > y0+dy || s->eY < y0-dy ) {
      eSegments->RemoveAt(i);
    }
  }
  eSegments->Compress();
  Log(2,"EdbSegmentsBox::RemovePosMargins","%d out of %d segments left", N(), nseg );
}

//______________________________________________________________________________
void EdbSegmentsBox::RemoveAngMargins(float  tx0, float  ty0, float  dtx, float  dty)
{
  int nseg = N();
  for(int i=nseg-1; i>=0; i-- )    {
    EdbSegP *s = GetSegment(i);
    if( s->eTX > tx0+dtx || s->eTX < tx0-dtx || s->eTY > ty0+dty || s->eTY < ty0-dty ) eSegments->RemoveAt(i);
  }
  eSegments->Compress();
  Log(2,"EdbSegmentsBox::RemoveAngMargins","%d out of %d segments left", N(), nseg );
}

//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsPlate(int plate)
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetPlate( plate );
  }
}
 
//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsZ()
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetZ( Z() );
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsDZ(float dz)
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetDZ( dz );
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsW(float w)
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetW( w );
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::SetSegmentsFlag(int flag)
{
  int nseg = N();
  for(int i=0; i<nseg; i++ )    {
    GetSegment(i)->SetFlag( flag );
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::ProjectTo(const float dz)
{
  eZ += dz;  eDZkeep += dz;
  int nseg = N();
  for(int i=0; i<nseg; i++ ) GetSegment(i)->PropagateToDZ(dz);
/*  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);
    p->SetX( p->X() + p->TX()*dz );
    p->SetY( p->Y() + p->TY()*dz );
    p->SetZ( Z() );
  }*/
}
 
//______________________________________________________________________________
int EdbSegmentsBox::CalculateXY( EdbSegmentsBox *pat, EdbAffine2D *aff )
{
  int n=N();
  if( n>pat->N() )  n=pat->N();
  if(n<2) return 0;

  aff->Calculate(this,pat);
  return 1;
}

//______________________________________________________________________________
int EdbSegmentsBox::CalculateAXAY( EdbSegmentsBox *pat, EdbAffine2D *aff )
{
  int n=N();
  if( n>pat->N() )  n=pat->N();
  if(n<2) return 0;

  float *ax1 = new float[n];
  float *ay1 = new float[n];
  float *ax2 = new float[n];
  float *ay2 = new float[n];

  EdbSegP *p1,*p2;
  for(int i=0; i<n; i++ ) {
    p1 = GetSegment(i);
    p2 = pat->GetSegment(i);
    ax1[i] = p1->TX();
    ay1[i] = p1->TY();
    ax2[i] = p2->TX();
    ay2[i] = p2->TY();
  }
  aff->Calculate(n,ax1,ay1,ax2,ay2);

  delete ax1;
  delete ay1;
  delete ax2;
  delete ay2;

  return 1;
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformA( const EdbAffine2D *aff )
{
  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = aff->A11()*p->TX() + aff->A12()*p->TY() + aff->B1();
    ty = aff->A21()*p->TX() + aff->A22()*p->TY() + aff->B2();
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformShr( float shr )
{
  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = p->TX()/shr;
    ty = p->TY()/shr;
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::TransformARot( const EdbAffine2D *aff )
{
  // apply to the angles only rotation members of transformation

  EdbSegP *p;
  float tx,ty;

  int nseg = N();
  for(int i=0; i<nseg; i++ ) {
    p = GetSegment(i);

    tx = aff->A11()*p->TX() + aff->A12()*p->TY();
    ty = aff->A21()*p->TX() + aff->A22()*p->TY();
    p->SetTX(tx);
    p->SetTY(ty);
  }
}

//______________________________________________________________________________
void EdbSegmentsBox::Print(Option_t *opt) const
{
  int nseg=GetN();
  printf("EdbSegmentsBox: %d segments\n", nseg );
  for(int i=0; i<nseg; i++) GetSegment(i)->Print();
} 

//______________________________________________________________________________
//______________________________________________________________________________
EdbTrackP::EdbTrackP(int nseg)
{
  Set0();
  if(nseg>0) eS  = new TSortedList();
  if(nseg>0) { eSF = new TSortedList();    eSF->SetOwner(); }
}

//______________________________________________________________________________
EdbTrackP::EdbTrackP(EdbSegP &seg) : EdbSegP(seg)
{
  //create empty track with segment parameters
  eS=0;
  eSF=0;
  eM=0;
  eDE=0;
  ePDG=-999;
  eVTAS = 0;
  eVTAE = 0;
  eNpl=0;
  eN0=0;
}
 
//______________________________________________________________________________
EdbTrackP::EdbTrackP(EdbSegP *seg, float m) : EdbSegP( *seg )
{
  eS=0;
  eSF=0;
  eM=0;
  eDE=0;
  ePDG=-999;
  eVTAS = 0;
  eVTAE = 0;
  eNpl=0;
  eN0=0;
  AddSegment(seg);
  SetM(m);
}
 
//______________________________________________________________________________
EdbTrackP::~EdbTrackP()
{
  if(eS)    { eS->Clear();  delete eS;  eS=0;  }
  if(eSF)   { eSF->Clear(); delete eSF; eSF=0; }
}

//______________________________________________________________________________
void EdbTrackP::Set0()
{
  ((EdbSegP*)this)->Set0();
  eS=0;
  eSF=0;
  eM=0;
  eDE=0;
  ePDG=-999;
  eVTAS = 0;
  eVTAE = 0;
  eNpl=0;
  eN0=0;
  ePerrUp=0;
  ePerrDown=0;
  return;
}

//______________________________________________________________________________
void EdbTrackP::AddVTA(EdbVTA *vta)
{
  if(vta) {
    if(vta->Zpos()==0)      eVTAE=vta;
    else if(vta->Zpos()==1) eVTAS=vta;
  }
}

//______________________________________________________________________________
void EdbTrackP::ClearVTA()
{
  eVTAE=0;
  eVTAS=0;
}

//______________________________________________________________________________
void EdbTrackP::ClearVTA(EdbVTA *vta)
{
  if(vta->Zpos()==0)      eVTAE=0;
  else if(vta->Zpos()==1) eVTAS=0;
}

//______________________________________________________________________________
void EdbTrackP::Copy(const EdbTrackP &tr)
{
  // do the physical copy of segments
  Clear();
  ((EdbSegP*)(this))->Copy( *((EdbSegP*)(&tr)) );
  SetM(tr.M());
  SetPDG(tr.PDG());
  SetNpl(tr.Npl());  
  SetN0(tr.N0());
  SetDE(tr.DE());
  AddVTA(tr.VTAS());
  AddVTA(tr.VTAE());

  int nseg=tr.N();
  for(int i=0; i<nseg; i++)
    AddSegment(new EdbSegP(*tr.GetSegment(i)));
  nseg=tr.NF();
  for(int i=0; i<nseg; i++)
    AddSegmentF(new EdbSegP(*tr.GetSegmentF(i)));
  if (eS) eS->SetOwner();
  if (eSF) eSF->SetOwner();
}

//______________________________________________________________________________
void EdbTrackP::Transform(const EdbAffine2D &aff)
{
  // apply transformation to all track elements
  ((EdbSegP*)(this))->Transform(&aff);
  for(int i=0; i<N(); i++)  GetSegment(i)->Transform(&aff);
  for(int i=0; i<NF(); i++) GetSegmentF(i)->Transform(&aff);
}

//______________________________________________________________________________
int EdbTrackP::CheckAliasSegments()
{
  int nalias=0;
  for(int i=0; i<N(); i++) 
    if( GetSegment(i)->Track() != ID()) nalias++;
  return nalias;
}

//______________________________________________________________________________
int EdbTrackP::RemoveAliasSegments()
{
  int nalias=0;
  EdbSegP *s=0;
  for(int i=0; i<N(); i++) {
    s = GetSegment(i);
    if( s->Track() != ID()) { 
      nalias++;
      RemoveSegment(s);
    }
  }
  return nalias;
}

//______________________________________________________________________________
int EdbTrackP::CheckMaxGap()
{
  int ngap=0;
  if(N()<2) return 0;
  EdbSegP *s1=0, *s2=0;
  int gap = 0;
  for(int i=0; i<N()-1; i++) {
    s1 = GetSegment(i);
    s2 = GetSegment(i+1);
    gap = TMath::Abs(s2->PID()-s1->PID());
    if( gap > ngap) ngap = gap;
  }
  return ngap;
}
//______________________________________________________________________________
int EdbTrackP::GetSegmentsMCTrack( int& nsegf ) const
{
  // return: ID of the MC track contributed to this reconstructed track 
  //         with the maximal number of segments
  // nsegf - return number of segments originated from the selected MC track

  int nseg = N();
  if (nseg < 2) return -1;
  if (nseg > 300) nseg = 300;
  int count[300];
  int i, f = 0;
  for(i=0; i<nseg; i++)
  {
    count[i] = 0;
    f = GetSegment(i)->MCTrack();
    if (f < 0) continue;
    for (int j=0; j<nseg; j++)
    {
	if ( f == GetSegment(j)->MCTrack()) count[i]++;
    }
  }
  int cmax = 0;
  int mctrackmax = -1;
  for(i=0; i<nseg; i++)
  {
	if (count[i] > cmax)
	{
	    cmax = count[i];
	    mctrackmax = GetSegment(i)->MCTrack();
	}
  }
  nsegf = cmax;
  return mctrackmax;
}

//______________________________________________________________________________
Float_t EdbTrackP::Wmean() const
{
  int n = N();
  double wtot=0;
  for(int i=0; i<n; i++) wtot+=GetSegment(i)->W();
  wtot/=n;
  return (Float_t)wtot;
}

//______________________________________________________________________________
int EdbTrackP::GetSegmentsFlag( int& nsegf ) const
{
  printf(" EdbTrackP::GetSegmentsFlag: obsolete, to be deleted!  now replaced by GetSegmentsMCTrack\n");

  int nseg = N();
  if (nseg < 2) return -1;
  if (nseg > 300) nseg = 300;
  int count[300];
  int i, f = 0;
  for(i=0; i<nseg; i++)
  {
    count[i] = 0;
    f = GetSegment(i)->Flag();
    if (f < 0) continue;
    for (int j=0; j<nseg; j++)
    {
	if ( f == GetSegment(j)->Flag()) count[i]++;
    }
  }
  int cmax = 0;
  int flagmax = -1;
  for(i=0; i<nseg; i++)
  {
	if (count[i] > cmax)
	{
	    cmax = count[i];
	    flagmax = GetSegment(i)->Flag();
	}
  }
  nsegf = cmax;
  return flagmax;
}
//______________________________________________________________________________
int EdbTrackP::GetSegmentsAid( int& nsegf ) const
{
  int nseg = N();
  if (nseg < 2) return -1;
  if (nseg > 300) nseg = 300;
  int count[300];
  int aid[300], ai;
  int i, f = 0;
  for(i=0; i<nseg; i++)
  {
    count[i] = 0;
    aid[i] = -1;
    ai = GetSegment(i)->Aid(1);
    f = ai + (GetSegment(i)->Aid(0))*1000;
    if (f < 0) continue;
    aid[i] = ai;
    for (int j=0; j<nseg; j++)
    {
	if ( f == (GetSegment(j)->Aid(1) + (GetSegment(i)->Aid(0))*1000)) count[i]++;
    }
  }
  int cmax = 0;
  int aidmax = -1;
  for(i=0; i<nseg; i++)
  {
	if (count[i] > cmax)
	{
	    cmax = count[i];
	    aidmax = aid[i];
	}
  }
  nsegf = cmax;
  return aidmax;
}
//______________________________________________________________________________
float EdbTrackP::Wgrains() const
{
  float w=0.;
  int nseg=N();
  for(int i=0; i<nseg; i++)    w+=GetSegment(i)->W();
  return w;
}

//______________________________________________________________________________
void EdbTrackP::AddTrack(const EdbTrackP &tr)
{
  int nseg=tr.N();
  int nsegf=tr.NF();
  for(int i=0; i<nseg; i++)
    AddSegment(tr.GetSegment(i));
  for(int i=0; i<nsegf; i++)
    AddSegmentF(new EdbSegP(*(tr.GetSegmentF(i))));
}

//______________________________________________________________________________
EdbSegP *EdbTrackP::GetSegmentWithClosestZ( float z, float dzMax )
{
  float dzmin=dzMax;
  EdbSegP *sbest=0;
  int nseg=N();
  for(int i=0; i<nseg; i++) {
     EdbSegP *s = GetSegment(i);
     float dz = Abs(s->eZ-z);
     if( dz < dzmin ) {
       dzmin=dz;
       sbest=s;
     }
   }
  return sbest;
}

//______________________________________________________________________________
void EdbTrackP::FitTrack()
{
  // track fit by averaging of segments parameters

  int nseg=N();
  float x=0,y=0,z=0,tx=0,ty=0,w=0;
  EdbSegP *seg=0;
  for(int i=0; i<nseg; i++) {
    seg = GetSegment(i);
    x  += seg->X();
    y  += seg->Y();
    z  += seg->Z();
    tx += seg->TX();
    ty += seg->TY();
    w  += seg->W();
  }
  x  /= nseg;
  y  /= nseg;
  z  /= nseg;
  tx /= nseg;
  ty /= nseg;
  Set(ID(),x,y,tx,ty,w,Flag());
  SetZ(z);
}

//______________________________________________________________________________
int  EdbTrackP::FitTrackKFS( bool zmax, float X0, int design )
{
  // if (zmax==true)  track parameters are calculated at segment with max Z
  // if (zmax==false) track parameters are calculated at segment with min Z
  // TODO: track parameters??
  //
  // Note: momentum for track must be defined with SetP() prior to call this routine!
  // Note: momentum dispersion for track should be defined with SetErrorP()
  //	   prior to call this routine - it is necessary for vertex fit
  //
  // on the output: Chi2: the full chi-square (not divided by NDF); NDF = 4
  //                Prob: is Chi2 probability (area of the tail of Chi2-distribution)
  //                      If we accept events with Prob >= ProbMin then ProbMin is the 
  //                      probability to reject the good event


  int nseg=N();
  if(nseg<=0) return 0;
  if(NF()) {
    eSF->Delete();
  }

  if(SP()<0.000001) SetErrorP(1.);    // TODO: razobratsia s etimi impul'sami!
  //printf("%d segments to fit\n",N());

  //TODO - eliminate constants!!!

  float dPb = 0., p = 0., m = 0.13957, e = 0.13957, de = 0., pa = 0., pn = 0.;
  float e0 = e;
  float eTPb = 1000./1300.;
  float pcut = 0.050;
  double teta0sq;
  double dz, ptx, pty;
  int step;
  int istart, iend;

  VtVector *par[260], *parpred[260], *pars[260], *meas[260];
  VtSqMatrix *pred[260];
  VtSymMatrix *cov[260], *covpred[260], *covpredinv[260], *covs[260], *dmeas[260];
 
  int i=0;
  if(nseg == 1) {
    EdbSegP *s = GetSegment(0);
    EdbSegP segf(*s);

    segf.Set(s->ID(),s->X(),s->Y(),s->TX(),s->TY(),1.,s->Flag());
    segf.SetZ(s->Z());
    segf.SetErrorP ( SP() );
    segf.SetChi2(0.);
    segf.SetProb(1.);
    segf.SetP( P() );
    segf.SetPID( s->PID() );
    segf.SetDZ(s->DZ());

    AddSegmentF(new EdbSegP(segf));
    return 0;
  }

  if(nseg>259)   return -1;

  EdbSegP *seg0=0;
  EdbSegP *seg=0;

  if(  GetSegment(N()-1)->Z() <  GetSegment(0)->Z() )
  {
    if (zmax)
    {
	step=-1;
	istart=nseg-1;
	iend=0;
    }
    else
    {
	step=1;
	istart=0;
	iend=nseg-1;
    }
  }
  else
  {
    if (!zmax)
    {
	step=-1;
	istart=nseg-1;
	iend=0;
    }
    else
    {
	step=1;
	istart=0;
	iend=nseg-1;
    }

  }

  seg0 = GetSegment(istart);

  par[istart] = new VtVector(  (double)(seg0->X()), 
			       (double)(seg0->Y()),  
			       (double)(seg0->TX()), 
			       (double)(seg0->TY()) );
  meas[istart] = new VtVector(*par[istart]);
  pred[iend] = new VtSqMatrix(4);
//  (*pred[istart]).clear();
//  (*pred[istart])(0,0) = 1.;
//  (*pred[istart])(1,1) = 1.;
//  (*pred[istart])(2,2) = 1.;
//  (*pred[istart])(3,3) = 1.;
  cov[istart] = new VtSymMatrix(4);             // covariance matrix for seg0
  for(int k=0; k<4; k++) 
    for(int l=0; l<4; l++) (*cov[istart])(k,l) = (seg0->COV())(k,l);
  dmeas[istart] = new VtSymMatrix(*cov[istart]);             // covariance matrix for seg0

  Double_t chi2=0.; 

  i = istart;
  p = P();
  m = M();
  if (p < pcut) p = pcut; 
  e = TMath::Sqrt((double)p*(double)p + (double)m*(double)m);
  e0 = e;
  while( (i+=step) != iend+step ) {

    seg = GetSegment(i);
 
    VtSymMatrix dms(4);   // multiple scattering matrix
    dms.clear();

    dz = seg->Z()-seg0->Z();
    ptx = (*par[i-step])(2);                        // previous tx
    pty = (*par[i-step])(3);                        // previous ty
    dPb = dz*TMath::Sqrt(1.+ptx*ptx+pty*pty); // thickness of the Pb+emulsion cell in microns
    if ((design != 0) && (p > pcut))
    {
	de = EdbPhysics::DeAveragePb(p, m, TMath::Abs(eTPb*dPb));
	if (de < 0.) de = 0.;
	if (design < 0) de = -de;
	if (dz >= 0.)
	    e = e - de;
	else
	    e = e + de;
	if (e < m) e = m;
	pn = TMath::Sqrt((double)e*(double)e - (double)m*(double)m);
	if (pn <= pcut) pn = pcut;
	pa = 0.5*(p + pn);
	p  = pn;
    }
    else
    {
	pa = p;
    }
    teta0sq = EdbPhysics::ThetaMS2( pa, m, dPb, X0 );
    dms(0,0) = teta0sq*dz*dz/3.;
    dms(1,1) = dms(0,0);
    dms(2,2) = teta0sq;
    dms(3,3) = dms(2,2);
//    dms(2,0) = teta0sq*TMath::Abs(dz)/2.;
    dms(2,0) = teta0sq*dz/2.;
    dms(3,1) = dms(2,0);
    dms(0,2) = dms(2,0);
    dms(1,3) = dms(2,0);

    pred[i-step] = new VtSqMatrix(4);        //propagation matrix for track parameters (x,y,tx,ty)
    pred[i-step]->clear();

    (*pred[i-step])(0,0) = 1.;
    (*pred[i-step])(1,1) = 1.;
    (*pred[i-step])(2,2) = 1.;
    (*pred[i-step])(3,3) = 1.;
    (*pred[i-step])(0,2) = dz;
    (*pred[i-step])(1,3) = dz;

    parpred[i] = new VtVector(4);            // prediction from seg0 to seg
    *parpred[i] = (*pred[i-step])*(*par[i-step]);

    covpred[i] = new VtSymMatrix(4);         // covariation matrix for prediction
    *covpred[i] = (*pred[i-step])*((*cov[i-step])*((*pred[i-step]).T()))+dms;

    dmeas[i] = new VtSymMatrix(4);           // original covariation  matrix for seg
    for(int k=0; k<4; k++) 
      for(int l=0; l<4; l++) (*dmeas[i])(k,l) = (seg->COV())(k,l);

    covpredinv[i] = new VtSymMatrix(4);
    (*covpredinv[i]) = (*covpred[i]).dsinv();
    VtSymMatrix dmeasinv(4);
    dmeasinv  = (*dmeas[i]).dsinv();
    cov[i] = new VtSymMatrix(4);
    (*cov[i]) = (*covpredinv[i]) + dmeasinv;
    (*cov[i]) = (*cov[i]).dsinv();

    meas[i] = new VtVector( (double)(seg->X()), 
			    (double)(seg->Y()),  
			    (double)(seg->TX()), 
			    (double)(seg->TY()) );

    par[i] = new VtVector(4);
    (*par[i]) = (*cov[i])*((*covpredinv[i])*(*parpred[i]) + dmeasinv*(*meas[i]));   // new parameters for seg

//    chi2 += ((*par[i])-(*parpred[i]))*((*covpredinv[i])*((*par[i])-(*parpred[i]))) + 
//	    ((*par[i])-(*meas[i]))*(dmeasinv*((*par[i])-(*meas[i])));

    VtSymMatrix dresid(4);
    dresid = (*dmeas[i]) - (*cov[i]);
    dresid = dresid.dsinv();

    chi2 += ((*par[i])-(*meas[i]))*(dresid*((*par[i])-(*meas[i])));

    seg0 = seg;
  }

  Set(ID(),(float)(*par[iend])(0),(float)(*par[iend])(1),
	   (float)(*par[iend])(2),(float)(*par[iend])(3),1.,Flag());
  SetZ(GetSegment(iend)->Z());
  SetPID(GetSegment(iend)->PID());
  SetCOV( (*cov[iend]).array(), 4 );

  //SetChi2((float)chi2);
  //SetProb( (float)TMath::Prob(chi2,nseg*4));

// Smoothing

  double chi2p=0;

  pars[iend] = new VtVector(*par[iend]);
  covs[iend] = new VtSymMatrix(*cov[iend]);
  VtSymMatrix dresid(4);
  dresid = (*dmeas[iend]) - (*covs[iend]);
  dresid = dresid.dsinv();
  chi2p = ((*pars[iend])-(*meas[iend]))*(dresid*((*pars[iend])-(*meas[iend])));

  EdbSegP segf;
  segf.Set(ID(),(float)(*pars[iend])(0),(float)(*pars[iend])(1),
	   (float)(*pars[iend])(2),(float)(*pars[iend])(3),1.,Flag());
  segf.SetZ(GetSegment(iend)->Z());
  segf.SetCOV( (*covs[iend]).array(), 4 );
  segf.SetErrorP ( SP() );
  segf.SetChi2((float)chi2p);
  segf.SetProb( (float)TMath::Prob(chi2p,4));
  segf.SetW( (float)nseg );
  segf.SetP( P() );
  segf.SetPID( GetSegment(iend)->PID() );
  segf.SetDZ(seg->DZ());
  
  AddSegmentF(new EdbSegP(segf));

  i=iend; 
  double DE=0.;
  EdbPhysics::DeAveragePbFastSet(P(), M());
   
  while( (i-=step) != istart-step ) {
	VtSqMatrix BackTr(4);
	BackTr = (*cov[i])*(((*pred[i]).T())*(*covpredinv[i+step]));
	pars[i] = new VtVector(4);
	covs[i] = new VtSymMatrix(4);
	(*pars[i]) = (*par[i]) + BackTr*((*pars[i+step])-(*parpred[i+step]));
	(*covs[i]) = (*cov[i]) + BackTr*(((*covs[i+step])-(*covpred[i+step]))*BackTr.T());
	dresid = (*dmeas[i]) - (*covs[i]);
	dresid = dresid.dsinv();
	chi2p = ((*pars[i])-(*meas[i]))*(dresid*((*pars[i])-(*meas[i])));
//	chi2 += chi2p;

	segf.Set(ID(),(float)(*pars[i])(0),(float)(*pars[i])(1),
		 (float)(*pars[i])(2),(float)(*pars[i])(3),1.,Flag());
	segf.SetZ(GetSegment(i)->Z());
	segf.SetCOV( (*covs[i]).array(), 4 );
	segf.SetErrorP ( SP() );
	segf.SetChi2((float)chi2p);
	segf.SetProb( (float)TMath::Prob(chi2p,4));
	segf.SetW( (float)nseg );
	segf.SetP( P() );
	segf.SetPID( GetSegment(i)->PID() );
	segf.SetDZ(seg->DZ());
	AddSegmentF(new EdbSegP(segf));
	dz = eTPb*TMath::Abs(GetSegment(i)->Z() - GetSegment(i+step)->Z()); 
        dPb = dz*TMath::Sqrt(1.+(*pars[i])(2)*(*pars[i])(2)+(*pars[i])(3)*(*pars[i])(3)); // thickness of the Pb+emulsion cell in microns
	DE += EdbPhysics::DeAveragePbFast(P(),M(),TMath::Abs(dPb));
  }
  SetChi2((float)chi2);
  SetProb((float)TMath::Prob(chi2,4*(nseg-1)));
  SetW( (float)nseg );
  SetDE( (float)DE );

//  DEBUG
//  printf(" e0 - e = %f, de = %f\n", e0-e, DE);

// Delete matrixes and vectors

  delete par[istart];
  par[istart] = 0;
  delete cov[istart];
  cov[istart] = 0;
  delete meas[istart];
  meas[istart] = 0;
  delete dmeas[istart];
  dmeas[istart] = 0;
  delete pred[istart];
  pred[istart] = 0;
  delete pars[istart];
  pars[istart] = 0;
  delete covs[istart];
  covs[istart] = 0;
  i=istart; 
  while( (i+=step) != iend+step ) {
    delete pred[i];
    pred[i] = 0;
    delete parpred[i];
    parpred[i] = 0;
    delete covpred[i];
    covpred[i] = 0;
    delete covpredinv[i];
    covpredinv[i] = 0;
    delete par[i];
    par[i] = 0;
    delete cov[i];
    cov[i] = 0;
    delete meas[i];
    meas[i] = 0;
    delete dmeas[i];
    dmeas[i] = 0;
    delete pars[i];
    pars[i] = 0;
    delete covs[i];
    covs[i] = 0;
  }
  return 0;
}

///______________________________________________________________________________
int EdbTrackP::EstimatePositionAt( Float_t z, EdbSegP &ss )
{
  // use coordinates of 2 nearest to z points for track extrapolation or interpolation
  // TODO: dz=0: mean?

  float x1,y1,tx,ty,dz;
  if( N()<2 )
  {
    x1=X();
    y1=Y();
    tx=TX();
    ty=TY();
    dz=z-GetSegment(0)->Z();
  }
  else
  {
    EdbSegP *s1=0,*s2=0;
    float   dz1,dz2;   dz1=dz2=kMaxInt;
    for(int i=0; i<N(); i++)
    {
      EdbSegP *s = GetSegment(i);
      float dz = Abs(s->Z()-z);
      if( dz < dz1 ) { dz1 = dz; s1=s; }
    }
    for(int i=0; i<N(); i++)
    {
      EdbSegP *s = GetSegment(i);
      float dz = Abs(s->Z()-z);
      if( dz < dz2 && s != s1 ) { dz2 = dz; s2=s; }
    }
    float dz0 = s2->Z()-s1->Z();
    if(abs(dz0)<0.000000001) {printf("dz0==0\n"); return 0;}
    float dx0 = s2->X()-s1->X();
    float dy0 = s2->Y()-s1->Y();
    tx  = dx0/dz0;
    ty  = dy0/dz0;
    dz  = z-s1->Z();
    x1  = s1->X();
    y1  = s1->Y();
  }
  ss.SetX( x1 + dz*tx );
  ss.SetY( y1 + dz*ty );
  ss.SetZ( z );
  ss.SetDZ( dz );           // keep dz distance
  ss.SetTX( tx );
  ss.SetTY( ty );
  return 1;
}

///______________________________________________________________________________
float   EdbTrackP::MakePredictionTo( Float_t z, EdbSegP &ss )
{
  float dz = Zmax()-Zmin();
  const EdbSegP *tr=0;
  if     ( z <= Zmin() )      tr = TrackZmin();   //TODO: fitted - not fitted
  else if( z >= Zmax() )      tr = TrackZmax();
  else {
    for(int i=0; i<N(); i++)       {
      EdbSegP *s = GetSegment(i);
      if( Abs(s->Z()-z) < dz ) { dz = Abs(s->Z()-z); tr=s; }     // select nearest segment (TODO: correct interpolation)
    }
  }
  if(!tr) // no segments in this track: use track body
    tr = this;
  dz = Abs(tr->Z()-z);
  ss.Copy(*tr);
  ss.PropagateTo(z);
  return dz;
}

///______________________________________________________________________________
int EdbTrackP::MakeSelector( EdbSegP &ss, bool followZ )
{
  if(N()<2) return 0;
  ss.SetCOV( GetSegment(0)->COV() );             // TODO ?
  const EdbSegP *tr;
  if (NF())
  {
    if( followZ ) tr = TrackZmax();
    else  tr = TrackZmin();
  }
  else
  {
    if( followZ ) tr = GetSegmentLast();
    else  tr = GetSegmentFirst();
  }
  ss.SetTX(tr->TX());
  ss.SetTY(tr->TY());
  ss.SetX(tr->X());
  ss.SetY(tr->Y());
  ss.SetZ(tr->Z());
  ss.SetPID(tr->PID());

  if( tr->PID() > GetSegment(0)->PID() )     return 1;
  if( tr->PID() > GetSegment(N()-1)->PID() ) return 1;
  return -1;
}

//______________________________________________________________________________
float EdbTrackP::CHI2()
{
  double dtx=0,dty=0,chi2=0;
  EdbSegP *seg=0;
  int    nseg=N();
  for(int i=0; i<nseg; i++) {
    seg = GetSegment(i);
    dtx = seg->TX()-TX();
    dty = seg->TY()-TY();
    chi2 += TMath::Sqrt( dtx*dtx/seg->STX() + 
			 dty*dty/seg->STY() );
  }
  chi2  /= nseg;
  return chi2;
}

//______________________________________________________________________________
float EdbTrackP::CHI2F()
{
  double dtx=0,dty=0,chi2=0;
  EdbSegP *s=0, *sf=0;
  int    nseg=N();
  for(int i=0; i<nseg; i++) {
    s  = GetSegment(i);
    sf = GetSegmentF(i);
    dtx = s->TX() - sf->TX();
    dty = s->TY() - sf->TY();
    chi2 += TMath::Sqrt( dtx*dtx/s->STX() + 
			 dty*dty/s->STY() );
  }
  chi2  /= nseg;
  return chi2;
}

//______________________________________________________________________________
Double_t EdbTrackP::GetBTEfficiency()
{
  // Returns Basetrack Scanning efficiency, estimated
  // by BT in the Tracks and the number of holes.
  // Formula: efficiency_per_BT = ( N()-2 ) / ( Npl()-2 ) 
  // This makes sense when having a track with 
  // at least 3 plates crossing.
  Double_t nseg=(Double_t)N();
  Double_t npl=(Double_t)Npl();
  
  if (nseg==2 && npl==2) return 1;
  Double_t BTEfficiency=(nseg-2)/(npl-2);
  
  return BTEfficiency;
}


//______________________________________________________________________________
void EdbTrackP::Print() 
{
  int nseg=0, nsegf=0;
  nseg = N();
  nsegf = NF();

  printf("EdbTrackP with %d segments and %d fitted segments\n", nseg, nsegf );
  printf("particle mass = %f\n", M() );
  ((EdbSegP*)this)->Print(); 

  if(nseg) 
    for(int i=0; i<nseg; i++)
      ((EdbSegP*)(eS->At(i)))->Print(); 

}

//______________________________________________________________________________
void EdbTrackP::PrintNice() 
{
  int nseg=0, nsegf=0;
  nseg = N();
  nsegf = NF();

  printf("EdbTrackP with %d segments and %d fitted segments. Estimated BT Scanning efficiency = %.02f \n", nseg, nsegf, GetBTEfficiency() );
  printf("  N  PID     ID          X             Y              Z        TX       TY     W      P      Flag     MC     Track    Chi2    Prob     Mass\n");
  printf("    %3d %8d %13.2f %13.2f %13.2f %7.4f  %7.4f %5.1f  %7.2f %7d %7d %7d  %7.4f  %7.4f  %5.3f\n",
	   PID(), ID(),X(),Y(),Z(),     TX(),   TY(),  W(),  P(), Flag(), MCEvt(),   Track(), Chi2(), Prob(),   M());

  EdbSegP *s=0;
  if(nseg) 
    for(int i=0; i<nseg; i++) {
      s = GetSegment(i);
      printf("%3d %3d %8d %13.2f %13.2f %13.2f %7.4f  %7.4f %5.1f  %7.2f %7d %7d %7d  %7.4f  %7.4f\n",
	     i, s->PID(), s->ID(),s->X(),s->Y(),s->Z(),  s->TX(),s->TY(),  s->W(),  s->P(), s->Flag(), s->MCEvt(), s->Track(),s->Chi2(),s->Prob());
    }
}


//______________________________________________________________________________
EdbVertex  *EdbTrackP::VertexS()
{
    if(eVTAS) return eVTAS->GetVertex();
    return 0;
}

//______________________________________________________________________________
EdbVertex  *EdbTrackP::VertexE()
{
    if(eVTAE) return eVTAE->GetVertex();
    return 0;
}

//______________________________________________________________________________
Float_t EdbTrackP::Zmin() const
{
  Float_t zmin = Z(); 
  if (N() && GetSegmentFirst()->Z() < zmin) zmin = GetSegmentFirst()->Z();
  return zmin;
}

//______________________________________________________________________________
Float_t EdbTrackP::Zmax() const
{
  Float_t zmax = Z(); 
  if (N() && GetSegmentLast()->Z() > zmax) zmax = GetSegmentLast()->Z();
  return zmax;
}
//______________________________________________________________________________
void EdbTrackP::SetPerrUp(Float_t perrUp)
{
  ePerrUp=perrUp;
}
//______________________________________________________________________________
void EdbTrackP::SetPerrDown(Float_t perrDown)
{
  ePerrDown=perrDown;
}
//______________________________________________________________________________
void EdbTrackP::SetPerr(Float_t perrDown, Float_t perrUp)
{
  ePerrUp=perrUp;
  ePerrDown=perrDown;
}

//______________________________________________________________________________
const EdbSegP  *EdbTrackP::TrackStart() const 
{ 
  if(!N())    return (EdbSegP*)this;
  if(Dir()<0) return GetSegmentLast();
  if(Dir()>0) return GetSegmentFirst();
  return (EdbSegP*)this;
}

//______________________________________________________________________________
const EdbSegP  *EdbTrackP::TrackEnd() const 
{ 
  if(!N())    return (EdbSegP*)this;
  if(Dir()>0) return GetSegmentLast();
  if(Dir()<0) return GetSegmentFirst();
  return (EdbSegP*)this;
}

//==============================================================================
EdbPattern::EdbPattern() : EdbSegmentsBox()
{
  eCell     = new TIndexCell();
  Set0();
}

//______________________________________________________________________________
EdbPattern::EdbPattern(float x0, float y0, float z0, int n) : EdbSegmentsBox(x0,y0,z0,n) 
{
  Set0();
  eCell     = new TIndexCell();
}

//______________________________________________________________________________
EdbPattern::EdbPattern(EdbPattern &p):EdbSegmentsBox(p)
{
  eID = p.eID;
  ePID = p.ePID;
  eSide = p.eSide;
  eScanID = p.eScanID;
  eFlag = p.eFlag;
  eNAff = p.eNAff;
  eStepX = p.eStepX;  eStepY = p.eStepY;
  eStepTX = p.eStepTX;  eStepTY = p.eStepTY;
  for(int i=0; i<4; i++) eSigma[i] = p.eSigma[i];
  if(p.eCell) eCell = new TIndexCell(*(p.eCell));
}
 
//______________________________________________________________________________
EdbPattern::~EdbPattern()
{
  if(eCell)  { delete eCell;  eCell=0; }
}

//______________________________________________________________________________
void EdbPattern::Set0()
{
  eID   = 0;
  ePID  = 0;
  eFlag = 0;
  eNAff = 0;
  eStepX=eStepY=eStepTX=eStepTY=0;
  for(int i=0; i<4; i++) eSigma[i]=0;
  eSide=0;
}

//______________________________________________________________________________
EdbSegP *EdbPattern::FindSegment(int id)
{
  for(int i=0; i<N(); i++) if(GetSegment(i)->ID()==id) return GetSegment(i);
  return 0;
}

//______________________________________________________________________________
float EdbPattern::SummaryPath()
{
  // calculate the microscope path of the eventual prediction scan of this pattern
  float sum=0,dx,dy;
  if(N()<2) return sum;
  for(int i=0; i<N()-1; i++ ) {
    dx = GetSegment(i+1)->X()-GetSegment(i)->X();
    dy = GetSegment(i+1)->Y()-GetSegment(i)->Y();
    sum += Sqrt(dx*dx+dy*dy);
  }
  return sum;
}

//______________________________________________________________________________
void EdbPattern::FillCell( float stepx, float stepy, float steptx, float stepty )
{
  // fill cells with fixed size at z=zPat

  TIndexCell *cell = Cell();
  if(cell) cell->Drop();
  SetStep(stepx,stepy,steptx,stepty);

  float x,y,tx,ty,dz;
  Long_t  val[5];  // x,y,ax,ay,i
  EdbSegP *p;
  int npat = N();
  for(int i=0; i<npat; i++ ) {
    p = GetSegment(i);
    dz = Z() - p->Z();
    tx = p->TX();
    ty = p->TY();
    x  = p->X() + tx*dz;
    y  = p->Y() + ty*dz;
    val[0]= (Long_t)( x / stepx  );
    val[1]= (Long_t)( y / stepy  );
    val[2]= (Long_t)( tx/ steptx );
    val[3]= (Long_t)( ty/ stepty );
    val[4]= (Long_t)(i);
    cell->Add(5,val);
  }
  cell->Sort();

}

//______________________________________________________________________________
int EdbPattern::FindCompliments(EdbSegP &s, TObjArray &arr, float nsigx, float nsigt)
{
  // return the array of segments compatible with the
  // prediction segment s with the accuracy of nsig (in sigmas)

  float dz = Z()-s.Z();
  float x = s.X() + s.TX()*dz;
  float y = s.Y() + s.TY()*dz;
  float sx = TMath::Sqrt( s.SX() + s.STX()*dz*dz );
  float sy = TMath::Sqrt( s.SY() + s.STY()*dz*dz );

  float stx = s.STX();
  if (stx <= 0.) stx = 1.;
  float sty = s.STY();
  if (sty <= 0.) sty = 1.;

//   printf("Step: %f %f %f %f\n",StepX(),StepY(),StepTX(),StepTY());
//   printf("Sigm: %f %f %f %f\n",sx,sy,TMath::Sqrt(stx),TMath::Sqrt(sty));


  long vcent[4] = { (long)(x/StepX()),
		    (long)(y/StepY()),
		    (long)(s.TX()/StepTX()),
		    (long)(s.TY()/StepTY())  };
  long vdiff[4] = { (long)(sx*nsigx/StepX()+1),
		    (long)(sy*nsigx/StepY()+1),
		    (long)(TMath::Sqrt(stx)*nsigt/StepTX()+1),
		    (long)(TMath::Sqrt(sty)*nsigt/StepTY()+1) };

  sy*=sy;
  sx*=sx;

  float dx,dy,dtx,dty;

  int nseg=0;
  EdbSegP *seg=0;


  //printf("find compliments in pattern %d\n",ID());
  //s.Print();

  long vmin[4],vmax[4];
  for(int i=0; i<4; i++) {
    vmin[i] = vcent[i]-vdiff[i];
    vmax[i] = vcent[i]+vdiff[i];
    //printf("%f +- %f \t",(float)vcent[i],(float)vdiff[i]);
  }
  //printf("\n");

  ///TODO: move this cycle into iterator

  TIndexCell *c1=0,*c2=0,*c3=0,*c4=0;
  for(vcent[0]=vmin[0]; vcent[0]<=vmax[0]; vcent[0]++) {
    c1 = eCell->Find(vcent[0]);
    if(!c1) continue;
    for(vcent[1]=vmin[1]; vcent[1]<=vmax[1]; vcent[1]++) {
      c2 = c1->Find(vcent[1]);
      if(!c2) continue;
      for(vcent[2]=vmin[2]; vcent[2]<=vmax[2]; vcent[2]++) {
	c3 = c2->Find(vcent[2]);
	if(!c3) continue;
	for(vcent[3]=vmin[3]; vcent[3]<=vmax[3]; vcent[3]++) {
	  c4 = c3->Find(vcent[3]);
	  if(!c4) continue;

	  for(int i=0; i<c4->N(); i++) {
	    seg = GetSegment(c4->At(i)->Value());
	    dtx=s.TX()-seg->TX();
	    if( dtx*dtx > stx*nsigt*nsigt )    continue;
	    dty=s.TY()-seg->TY();
	    if( dty*dty > sty*nsigt*nsigt )    continue;
	    dz = seg->Z()-s.Z();
	    dx=s.X()+s.TX()*dz-seg->X();
	    if( dx*dx > sx*nsigx*nsigx )           continue;
	    dy=s.Y()+s.TY()*dz-seg->Y();
	    if( dy*dy > sy*nsigx*nsigx )           continue;
	    //if(!s.IsCompatible(*seg,nsigx,nsigt)) continue;
	    arr.Add(seg);
	    nseg++;
	  }

	}
      }
    }
  }


//    TIndexCellIterV itr( eCell, 4, vcent, vdiff );
//    const TIndexCell *c=0;
//    while( (c=itr.Next()) ) 
//      for(int i=0; i<c->N(); i++) {
//        seg = GetSegment(c->At(i)->Value());
//        if(!s.IsCompatible(*seg,nsigx,nsigt)) continue;
//        arr.Add(seg);
//        nseg++;
//      }

  return nseg;
}

//______________________________________________________________________________
void EdbPattern::SetSegmentsPID()
{
  int nseg = N();
  for(int i=0; i<nseg; i++) {
    GetSegment(i)->SetPID(ID());  //PID of the segment must be ID of the pattern!
  }
}

//______________________________________________________________________________
void EdbPattern::SetSegmentsScanID(EdbID id)
{
  for(int i=0; i<N(); i++)     GetSegment(i)->SetScanID(id);
}

//______________________________________________________________________________
Int_t EdbPattern::AddPattern(EdbPattern &p)
{
  int n=p.N();
  Expand( N()+n );
  for(int i=0; i<n; i++)
    AddSegment( *(p.GetSegment(i)) );
  return n;
}

//______________________________________________________________________________
EdbPattern *EdbPattern::ExtractSubPatternXY(float xmin, float xmax, float ymin, float ymax)
{
  EdbPattern *pat = new EdbPattern( X(), Y(), Z() );
  EdbSegP *s;
  int nseg = N();

  for(int i=0; i<nseg; i++) {
    s = GetSegment(i);
    if(s->X() < xmin || s->X() > xmax || s->Y() < ymin || s->Y() > ymax) continue;
    pat->AddSegment(*s);
  }

  // Set ID() and PID() to have consistent values of this pattern:
  pat->SetID(this->ID());
  pat->SetPID(this->PID());

  return pat;
}

//______________________________________________________________________________
EdbPattern *EdbPattern::ExtractSubPattern(float min[5], float max[5], int MCEvt)
{
  //
  // return the copy of selected segments
  //
  EdbPattern *pat = new EdbPattern( X(), Y(), Z() );
  EdbSegP *s;
  int nseg = N();

  for(int i=0; i<nseg; i++) {
    s = GetSegment(i);
		
    if (s->MCEvt()!=MCEvt && s->MCEvt()>0 && MCEvt>=0) continue;
    // Do not continue not in case MCEvt was not specified at all.
    // This allows backward compability.
    
    if(s->X()  < min[0])   continue;
    if(s->Y()  < min[1])   continue;
    if(s->TX() < min[2])   continue;
    if(s->TY() < min[3])   continue;
    if(s->W()  < min[4])   continue;

    if(s->X()  > max[0])   continue;
    if(s->Y()  > max[1])   continue;
    if(s->TX() > max[2])   continue;
    if(s->TY() > max[3])   continue;
    if(s->W()  > max[4])   continue;

    pat->AddSegment(*s);
  }
  
  // Set ID() and PID() to have consistent values of this pattern:
  pat->SetID(this->ID());
  pat->SetPID(this->PID());
  
  return pat;
}

//______________________________________________________________________________
void EdbPattern::Reset()
{
  ((EdbSegmentsBox *)this)->Reset();
  Set0();
  if(eCell)     eCell->Drop();
}

//______________________________________________________________________________
float EdbPattern::Xmean()
{
  Double_t s=0;
  int n=N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += GetSegment(i)->eX;
  s /= n;
  return (float)s;
}
/*
//______________________________________________________________________________
float EdbPattern::Xmin() const // reload virtual function
{
  int n=N(); if(n<1) return 0;
  float v0=GetSegment(0)->eX;
  if(n>1) for(int i=1; i<n; i++) { float v=GetSegment(i)->eX; 
    if(v>kMinInt) v0=(v<v0)?v:v0; }                               // hack to exclude the infinity values (to check in viewsideal!)
  return v0;
}
//______________________________________________________________________________
float EdbPattern::Xmax() const // reload virtual function
{
  int n=N(); if(n<1) return 0;
  float v0=GetSegment(0)->eX;
  if(n>1) for(int i=1; i<n; i++) { float v=GetSegment(i)->eX; 
    if(v<kMaxInt)  v0=(v>v0)?v:v0; }
  return v0;
}
//______________________________________________________________________________
float EdbPattern::Ymin() const // reload virtual function
{
  int n=N(); if(n<1) return 0;
  float v0=GetSegment(0)->eY;
  if(n>1) for(int i=1; i<n; i++) { float v=GetSegment(i)->eY; 
    if(v>kMinInt)  v0=(v<v0)?v:v0; }
  return v0;
}
//______________________________________________________________________________
float EdbPattern::Ymax() const // reload virtual function
{
  int n=N(); if(n<1) return 0;
  float v0=GetSegment(0)->eY;
  if(n>1) for(int i=1; i<n; i++) { float v=GetSegment(i)->eY; 
    if(v<kMaxInt) v0=(v>v0)?v:v0; }
  return v0;
}
*/
//______________________________________________________________________________
float EdbPattern::Ymean()
{
  Double_t s=0;
  int n=N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += GetSegment(i)->eY;
  s /= n;
  return (float)s;
}

 
////////////////////////////////////////////////////////////////////////////////
//
//   EdbPatternsVolume
//
////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
EdbPatternsVolume::EdbPatternsVolume()
{
  ePatterns   = new TObjArray();
  eTracksCell = 0;
  ePatternsCell = 0;
  Set0();
}

//______________________________________________________________________________
EdbPatternsVolume::EdbPatternsVolume(EdbPatternsVolume &pvol)
{
  ePatterns   = new TObjArray();
  eTracksCell = 0;
  ePatternsCell = 0;
  Set0();

  pvol.PassProperties(*this);
  int npat,nseg;
  npat = Npatterns();
  for(int j=0; j<npat; j++) {
    nseg = GetPattern(j)->N();
    for(int i=0; i<nseg; i++ ) {
      pvol.GetPattern(j)->AddSegment( *(GetPattern(j)->GetSegment(i)) );
    }
  }
}

//______________________________________________________________________________
EdbPatternsVolume::~EdbPatternsVolume()
{
  if(ePatterns) {
    ePatterns->Delete();
    delete ePatterns;
    ePatterns=0;
  }
  if(eTracksCell)   { delete eTracksCell;   eTracksCell=0; }
  if(ePatternsCell) { delete ePatternsCell; ePatternsCell=0; }
}

//______________________________________________________________________________
void EdbPatternsVolume::Set0()
{
  eX=0; eY=0; eZ=0;
  eDescendingZ=0;
}
 
//______________________________________________________________________________
int EdbPatternsVolume::DropCouples()
{
  int count=0;
  int npat=Npatterns();
  for(int i=0; i<npat; i++ )
    count += GetPattern(i)->Cell()->DropCouples(4);
  if(count) printf("%d couples are dropped in volume cells\n",count);
  return count;
}
 
//______________________________________________________________________________
void EdbPatternsVolume::SetPatternsID()
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ )
    GetPattern(i)->SetID(i);
}
 
//______________________________________________________________________________
void EdbPatternsVolume::Transform( const EdbAffine2D *aff )
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ )  {
    GetPattern(i)->Transform(aff);
    GetPattern(i)->TransformARot(aff);
  }
}
 
//______________________________________________________________________________
void EdbPatternsVolume::Centralize()
{
  // find geometrical center (XY) of all patterns and set it as the center of 
  // coordinates  to simplify transformations
  // To be used before any operations on patterns

  float xc=0;
  float yc=0;
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    xc += GetPattern(i)->Xmax() + GetPattern(i)->Xmin();
    yc += GetPattern(i)->Ymax() + GetPattern(i)->Ymin();
  }
  xc = xc/Npatterns()/2;
  yc = yc/Npatterns()/2;

  Centralize(xc,yc);
}

//______________________________________________________________________________
void EdbPatternsVolume::Centralize( float xc, float yc )
{
  eX = xc;  eY=yc;
  Shift(-xc,-yc);
  float npat = Npatterns();
  for(int i=0; i<npat; i++ ) 
    GetPattern(i)->SetKeep(1,0,0,1,0,0);
}

//______________________________________________________________________________
void EdbPatternsVolume::PrintAff() const
{
  EdbAffine2D a;
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    GetPattern(i)->GetKeep(a);
    printf(" %3d  (%5d) Z=%13.3f ",i,GetPattern(i)->NAff(), GetPattern(i)->Z() ); a.Print();
  }
}

//______________________________________________________________________________
void EdbPatternsVolume::PrintStat( Option_t *opt) const
{
  int npat = Npatterns();
  printf("\nVolume statistics for %d patterns\n",npat);

  float dx,dy;
  EdbPattern *pat=0;
  printf("pat# \t segments \t dX \t\tdY \t meanDist \n");
  for(int i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    if(pat) {
    dx = pat->Xmax() - pat->Xmin();
    dy = pat->Ymax()- pat->Ymin();
    printf(" %d\t %d\t %10.2f \t %10.2f \t %10.4f \n", 
	   i, pat->GetN(),dx,dy, TMath::Sqrt(dx*dy/pat->GetN()) );
    }
  }

  npat=Npatterns();
  for(int i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    if(pat) if(pat->Cell()) pat->Cell()->PrintStat();
  }
}
 
//______________________________________________________________________________
void EdbPatternsVolume::PrintStat(EdbPattern &pat) const
{
} 

//______________________________________________________________________________
void EdbPatternsVolume::Print() const
{
  int npat = Npatterns();
  printf("\nEdbPatternsVolume with %d patterns\n",npat);
  EdbPattern *pat=0;
  for(int i=0; i<npat; i++ ) {
    pat = GetPattern(i);
    printf(" id=%3d pid=%3d   x:y:z =  %12.3f %12.3f %12.3f  n= %8d   ScanID = %s \tside=%d\n", 
	   pat->ID(), pat->PID(), pat->X(),pat->Y(),pat->Z(),pat->N(), pat->ScanID().AsString(), pat->Side() );
  }
  printf("\n");
} 

//______________________________________________________________________________
void EdbPatternsVolume::AddPattern( EdbPattern *pat )
{
  ePatterns->Add(pat);
}

//______________________________________________________________________________
void EdbPatternsVolume::AddPatternAt( EdbPattern *pat,int id )
{
  if(id>=0&&id<1000 && pat)
  {
    if(ePatterns->GetSize()<id+1) ePatterns->Expand(id+1);
    pat->SetID(id);
    ePatterns->AddAt(pat,id);
  } else {
    Log(1,"EdbPatternsVolume::AddPatternAt","Error! id out of range (0,1000): id = %d",id);
    if(!pat) Log(1,"EdbPatternsVolume::AddPatternAt","Error! pat is not defined!");
 }
}

//______________________________________________________________________________
EdbPattern *EdbPatternsVolume::GetPattern( int id ) const
{
  if(Npatterns()>id) return (EdbPattern*)ePatterns->At(id);
  else return 0;
}


//______________________________________________________________________________
void EdbPatternsVolume::PassProperties( EdbPatternsVolume &pvol )
{
  pvol.SetXYZ(eX,eY,eZ);
  EdbAffine2D a;
  EdbPattern *p=0;

  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) {
    p = GetPattern(i);
    p->GetKeep(a);
    EdbPattern *psel = new EdbPattern( p->X(),p->Y(),p->Z() );
    psel->SetKeep( a.A11(), a.A12(), a.A21(), a.A22(), a.B1(),a.B2() );
    pvol.AddPattern(psel);
  }
  pvol.SetPatternsID();
}


//______________________________________________________________________________
void EdbPatternsVolume::Shift( float x, float y )
{
  EdbAffine2D aff;
  aff.ShiftX(x);
  aff.ShiftY(y);

  int npat = Npatterns();
  for(int i=0; i<npat; i++)
    GetPattern(i)->Transform(&aff);
}

//______________________________________________________________________________
void EdbPatternsVolume::DropCell()
{
  int npat = Npatterns();
  for(int i=0; i<npat; i++ ) GetPattern(i)->Cell()->Drop();
}

//______________________________________________________________________________
EdbPattern  *EdbPatternsVolume::GetPatternByPlate(int plate, int side)
{
  EdbPattern *p=0;
  for(int i=0; i<Npatterns(); i++) {
    p = GetPattern(i);
    if(p) if( p->Plate()==plate && p->Side()==side ) return p;
  }
  return 0;
}

///______________________________________________________________________________
EdbPattern  *EdbPatternsVolume::InsertPattern(EdbPattern *pat, Bool_t descendingZ)
{
  // insert new pattern using it's Z as a main criteria for the positioning
  EdbPattern *p  = GetPatternByPlate( pat->Plate(), pat->Side() );
  if(p)  {
    Log(1,"EdbPatternsVolume::InsertPattern",
	"ERROR: Attempt to insert new pattern for already existing plate/side z = %f znew = %f", 
	pat->Plate(), pat->Side(), p->Z(), pat->Z());
    return p;
  }
  AddPattern(pat);
  SortPatternsByZ(descendingZ);
  return pat;
}

//______________________________________________________________________________
void EdbPatternsVolume::SortPatternsByZ(Bool_t descendingZ)
{
  int n = Npatterns();
  int   *ind = new int[n];
  float *z   = new float[n];
  for(int i=0; i<Npatterns(); i++) {
    ind[i] = i; 
    z[i] = GetPattern(i)->Z();
  }
  TMath::Sort(n,z,ind,descendingZ);
  TObjArray *pnew = new TObjArray(ePatterns->GetSize());
  for(int i=0; i<n; i++)     pnew->Add(GetPattern(ind[i]));
  ePatterns->Clear();
  delete ePatterns;
  ePatterns = pnew;
  delete ind;
  delete z;
}



//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternZLowestHighest(Bool_t lowestZ) const
{
  // By default it returns the pattern with the lowest Z position. This is
  // NOT necessarily the first pattern of the patterns volume.
  
  if (Npatterns()==0) {
    Log(1,"EdbPatternsVolume::GetPatternZLowestHighest","ERROR: Attempt to get Pattern, but NO patterns are in EdbPatternsVolume!");
    return 0;
  }
  if (gEDBDEBUGLEVEL>3) cout << "EdbPatternsVolume::GetPatternZLowestHighest  lowestZ = " << lowestZ << endl;
  
  EdbPattern *pat=0;
  pat=GetPattern(0);
  if (Npatterns()==1) return pat;
  float zpat0,zpat1;
  zpat0=GetPattern(0)->Z();
  zpat1=GetPattern(Npatterns()-1)->Z();
  if (gEDBDEBUGLEVEL>2) {
  cout << "EdbPatternsVolume  Pattern At Position 0    has Z = " << zpat0 << endl;
  cout << "EdbPatternsVolume  Pattern At Position last has Z = " << zpat1 << endl;
  }
  if (zpat0>zpat1) pat=GetPattern(Npatterns()-1);
  if (kTRUE==lowestZ) return pat;
  if (kFALSE==lowestZ && zpat0<zpat1) return GetPattern(Npatterns()-1);
  return GetPattern(0);
}


//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternSucceding(EdbPattern* pat) const
{
	return GetPatternNext(pat->Z(),1);
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternPreceding(EdbPattern* pat) const
{
	return GetPatternNext(pat->Z(),-1);
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternNext(float z, int dir) const
{
  // Return next pattern in either downstream (dir=+1) (i.e. Z gets bigger)
  //                  or in either   upstream (dir=-1) (i.e. Z gets smaller)
  // direction.
  // If the last/first pattern is already given and one ask therefor 
  // for a pattern that does not exists, then return a NULL pattern.
  
  EdbPattern *pat=0;
  float dz=dir*99999999.;
  float zpat;
  for(int i=0; i<Npatterns(); i++) {
    zpat=GetPattern(i)->Z();
    if(dir>0) if(zpat>z) if(zpat-z<dz) {dz=zpat-z; pat=GetPattern(i);}
    if(dir<0) if(zpat<z) if(zpat-z>dz) {dz=zpat-z; pat=GetPattern(i);}
  }
  if (gEDBDEBUGLEVEL>2) 
		if (pat == 0) cout <<  "EdbPattern* EdbPatternsVolume::GetPatternNext()   WARNING: pattern is NULL ! Check your input!" << endl;
  return pat;
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternByPID(int pid) const
{
  // Return pattern having PID()  == pid.
  // If there is no such pattern, return a NULL pattern.
  //
  // This is important because some codings (libShower for example).
  // rely on getting Patterns by the PIDs.
  
  EdbPattern *pat=0;
	for(int i=0; i<Npatterns(); i++) {
		pat = GetPattern(i);
    if (GetPattern(i)->PID() == pid ) return pat;
	}
	
  if (gEDBDEBUGLEVEL>2) 
		if (pat == 0) cout <<  "EdbPattern* EdbPatternsVolume::GetPatternNext()   WARNING: pattern is NULL ! Check your input!" << endl;
  return NULL;
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::GetPatternByZ(float z, float tolerance) const
{
  // Return pattern having Z()  == z +- tolerance(microns) due to roundings.
  // If there is no such pattern, return a NULL pattern.
  //
  // This is important because some codings (libShower for example).
  // rely on getting Patterns by the Z.
  // Warning: If tolerance is to big, there is not a unique solution anymore
  
  EdbPattern *pat=0;
	for(int i=0; i<Npatterns(); i++) {
		pat = GetPattern(i);
    if ( TMath::Abs(GetPattern(i)->Z() - z)<tolerance ) return pat;
	}
	
  if (gEDBDEBUGLEVEL>2) 
		if (pat == 0) cout <<  "EdbPattern* EdbPatternsVolume::GetPatternNext()   WARNING: pattern is NULL ! Check your input!" << endl;
  return NULL;
}

//______________________________________________________________________________
EdbPattern* EdbPatternsVolume::NextPattern(float z, int dir) const
{
  // Return next pattern in either downstream (dir=+1)
  //                  or in either   upstream (dir=-1)
  // direction.
  
  EdbPattern *pat=0;
  float dz=dir*99999999.;
  float zpat;
  for(int i=0; i<Npatterns(); i++) {
    zpat=GetPattern(i)->Z();
    if(dir>0) if(zpat>z) if(zpat-z<dz) {dz=zpat-z; pat=GetPattern(i);}
    if(dir<0) if(zpat<z) if(zpat-z>dz) {dz=zpat-z; pat=GetPattern(i);}
  }
  return pat;
}

//______________________________________________________________________________
int EdbPatternsVolume::FindComplimentsVol(EdbSegP &ss, TObjArray &arr, float nsig, float nsigt, int Dpat)
{
  EdbPattern *pat  = 0;
  int npat = Npatterns();
  bool p_inverse_z = (GetPattern(1)->Z() - GetPattern(0)->Z()) > 0. ? false : true ;
  int ii = 0;
  int dpat = 1;
  int pid = 0;
  if (p_inverse_z) dpat = -1;
  for (int i = 0; i<npat; i++)
  {
    ii = i;
    if (p_inverse_z) ii = npat-1-i;
    pat = GetPattern(ii);  
    if (ss.Z() < pat->Z())
    {
	pid = ii - dpat;
	if ( pid >= npat )  ss.SetPID(npat-1);
	if ( pid <  0 )     ss.SetPID(0);
	else	            ss.SetPID(pid);
	break;
    }
    else if (ss.Z() == pat->Z())
    {
	ss.SetPID(ii);
	break;
    }
  }

  int p0 = ss.PID();
  int pend   = p0 + Dpat;
  if (pend >= npat) pend = npat - 1;
  int pstart   = p0 - Dpat;
  if (pstart < 0) pstart = 0;

  //arr.Clear();
  int nseg = 0;
  int n = 0;
  for(int i=pstart; i <=pend; i++ ) {
    pat = GetPattern(i);
    if(!pat)                   continue;
    n = pat->FindCompliments(ss,arr,nsig,nsigt);
    nseg += n;
  }

  return nseg;
}

//______________________________________________________________________________
Float_t EdbPatternsVolume::Xmean()
{
	if (Npatterns()==0) return 0;
	EdbPattern* pat = GetPattern(0);
	Double_t s=0;
  int n=pat->N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += pat->GetSegment(i)->eX;
  s /= n;
  return (float)s;
}

//______________________________________________________________________________
Float_t EdbPatternsVolume::Ymean()
{
	if (Npatterns()==0) return 0;
	EdbPattern* pat = GetPattern(0);
	Double_t s=0;
  int n=pat->N(); if(n<1) return 0;
  for(int i=0; i<n; i++) s += pat->GetSegment(i)->eY;
  s /= n;
  return (float)s;
}
