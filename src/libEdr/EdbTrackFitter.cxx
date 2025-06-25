//-- Author :  Valeri Tioukov  8.10.2005

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TGraphErrors.h"
#include "TVector3.h"

#include "EdbLog.h"
#include "EdbMath.h"
#include "EdbPhys.h"
#include "EdbAffine.h"
#include "EdbLayer.h"
#include "EdbPVRec.h"
#include "EdbScanCond.h"
#include "EdbTrackFitter.h"

#include "vt++/CMatrix.hh"
#include "vt++/VtVector.hh"
#include "vt++/VtDistance.hh"
#include "smatrix/Functions.hh"
#include "smatrix/SVector.hh"

ClassImp(EdbTrackFitter);

using namespace TMath;
using namespace MATRIX;
using namespace VERTEX;

//______________________________________________________________________________
EdbTrackFitter::EdbTrackFitter()
{
  SetDefaultBrick();
}

//______________________________________________________________________________
void EdbTrackFitter::SetDefaultBrick()
{
  eX0            = EdbPhysics::kX0_Pb();
  eTPb           = 1000./1300.;
  ePcut          = 0.050;
  eM             = 0.13957;
  eDE_correction = false;
}

//______________________________________________________________________________
void EdbTrackFitter::Print()
{
  printf("EdbTrackFitter seetings:\n");
  printf("eX0      = %f\n",eX0);
  printf("eM       = %f\n",eM);
  printf("\n");
}

//________________________________________________________________________
float EdbTrackFitter::Chi2Seg( EdbSegP *tr, EdbSegP *s)
{
  // Estimation of chi2 using the covariance matrix
  // Return value:        Prob: is Chi2 probability (area of the tail of Chi2-distribution)
  //                      If we accept couples with Prob >= ProbMin then ProbMin is the 
  //                      probability to reject the good couple
  //
  // The mass and momentum of the tr are used for multiple scattering estimation

  double dz;
  float prob;
  VtVector par( (double)(tr->X()), 
		(double)(tr->Y()),  
		(double)(tr->TX()), 
		(double)(tr->TY()) );
  VtSymMatrix cov(4);             // covariance matrix for seg0 (measurements errors)
  for(int k=0; k<4; k++) 
    for(int l=0; l<4; l++) cov(k,l) = (tr->COV())(k,l);

  Double_t chi2=0.; 
  dz = s->Z()-tr->Z();
  VtSqMatrix pred(4);        //propagation matrix for track parameters (x,y,tx,ty)
  pred.clear();
  pred(0,0) = 1.;
  pred(1,1) = 1.;
  pred(2,2) = 1.;
  pred(3,3) = 1.;
  pred(0,2) = dz;
  pred(1,3) = dz;
  VtVector parpred(4);            // prediction from seg0 to seg
  parpred = pred*par;
  VtSymMatrix covpred(4);         // covariance matrix for prediction
  covpred = pred*(cov*pred.T());

  VtSymMatrix dmeas(4);           // original covariance  matrix for seg2
  for(int k=0; k<4; k++) 
    for(int l=0; l<4; l++) dmeas(k,l) = (s->COV())(k,l);
  
  covpred = covpred.dsinv();
  dmeas   = dmeas.dsinv();
  cov = covpred + dmeas;
  cov = cov.dsinv();
  
  VtVector meas( (double)(s->X()), 
		 (double)(s->Y()),  
		 (double)(s->TX()), 
		 (double)(s->TY()) );

  par = cov*(covpred*parpred + dmeas*meas);   // new parameters for seg
  chi2 = (par-parpred)*(covpred*(par-parpred)) + (par-meas)*(dmeas*(par-meas));
  prob = (float)TMath::Prob(chi2,4);

  tr->Set(tr->ID(),(float)par(0),(float)par(1),(float)par(2),(float)par(3),tr->W(),tr->Flag());
  tr->SetCOV( cov.array(), 4 );
  tr->SetChi2((float)chi2);
  tr->SetProb(prob);
  tr->SetZ(s->Z());
  tr->SetW(tr->W()+s->W());
  return TMath::Sqrt(chi2/4.);
}

//______________________________________________________________________________
float EdbTrackFitter::Chi2ACP( EdbSegP s1, EdbSegP s2, EdbScanCond &cond)
{
  // Exact copy of the function traditionally used for the linking of couples
  //
  // fast estimation of chi2 in the special case when the position 
  // errors of segments are negligible in respect to angular errors: 
  // sigmaXY/dz << sigmaTXY
  // application: up/down linking, alignment (offset search)
  //
  // All calculation are done in the track plane which remove the 
  // dependency of the polar angle (phi)

  TVector3 v1,v2,v;
  v1.SetXYZ( s1.TX(), s1.TY() , -1. );
  v2.SetXYZ( s2.TX(), s2.TY() , -1. );
  v.SetXYZ(  -( s2.X() - s1.X() ),
	     -( s2.Y() - s1.Y() ),
	     -( s2.Z() - s1.Z() ) );

  float phi = v.Phi();
  v.RotateZ(  -phi );
  v1.RotateZ( -phi );
  v2.RotateZ( -phi );

  float dz  = v.Z();
  float tx  = v.X()/dz;
  float ty  = v.Y()/dz;
  float stx   = cond.SigmaTX(tx);
  float sty   = cond.SigmaTY(ty);

  float prob1=1., prob2=1.;
  //if(iprob) {
  prob1 = cond.ProbSeg(tx,ty,s1.W());
  prob2 = cond.ProbSeg(tx,ty,s2.W());
    // }
  float dtx1 = (v1.X()-tx)*(v1.X()-tx)/stx/stx/prob1;
  float dty1 = (v1.Y()-ty)*(v1.Y()-ty)/sty/sty/prob1;
  float dtx2 = (v2.X()-tx)*(v2.X()-tx)/stx/stx/prob2;
  float dty2 = (v2.Y()-ty)*(v2.Y()-ty)/sty/sty/prob2;

  float chi2a=TMath::Sqrt(dtx1+dty1+dtx2+dty2)/2.;
  return chi2a;
}

//______________________________________________________________________________
float EdbTrackFitter::Chi2ASegLL( EdbSegP &s1, EdbSegP &s2, EdbSegP &seg, EdbScanCond &cond1, EdbScanCond &cond2)
{
  // fast estimation of chi2 in the special case when the position
  // errors of segments are negligible in respect to angular errors:
  // sigmaXY/dz << sigmaTXY
  // application: up/down linking,
  //
  // All calculation are done in the track plane which remove the
  // dependency of the polar angle (phi)
  // In this function the Likelihood estimated by Andrey and stored in s1.eChi2 and s2.eChi2 are used

  Float_t dz =  s2.Z() - s1.Z();
  seg.Set( 0,
	   0.5*(s1.X() + s2.X()),
	   0.5*(s1.Y() + s2.Y()),
	   (s2.X() - s1.X())/dz,
	   (s2.Y() - s1.Y())/dz,
	   s1.W() + s2.W(),
	   0);
  seg.SetZ( 0.5*(s2.Z()+s1.Z()) );

  Float_t  phi = -seg.Phi();
  Double_t s   = TMath::Sin(phi);
  Double_t c   = TMath::Cos(phi);
  Float_t  tx  = seg.TX()*c-seg.TY()*s;
  Float_t  ty  = seg.TX()*s+seg.TY()*c;
  Float_t  tx1 = s1.TX()*c-s1.TY()*s;
  Float_t  ty1 = s1.TX()*s+s1.TY()*c;
  Float_t  tx2 = s2.TX()*c-s2.TY()*s;
  Float_t  ty2 = s2.TX()*s+s2.TY()*c;

  Float_t stx   = cond1.SigmaTX(tx);
  Float_t sty   = cond1.SigmaTY(ty);
  Float_t prob1 = cond1.ProbLL(tx,ty,s1.Chi2());
  Float_t prob2 = cond2.ProbLL(tx,ty,s2.Chi2());

  Float_t dtx1  = (tx1-tx)*(tx1-tx)/stx/stx/prob1;
  Float_t dty1  = (ty1-ty)*(ty1-ty)/sty/sty/prob1;
  Float_t dtx2  = (tx2-tx)*(tx2-tx)/stx/stx/prob2;
  Float_t dty2  = (ty2-ty)*(ty2-ty)/sty/sty/prob2;
  Float_t chi2a = Sqrt(dtx1+dty1+dtx2+dty2)/2.;

  seg.SetChi2(chi2a);
  return chi2a;
}

///______________________________________________________________________________
float EdbTrackFitter::Chi2ASeg( EdbSegP &s1, EdbSegP &s2, EdbSegP &seg, EdbScanCond &cond1, EdbScanCond &cond2)
{
  // fast estimation of chi2 in the special case when the position
  // errors of segments are negligible in respect to angular errors:
  // sigmaXY/dz << sigmaTXY
  // application: up/down linking,
  //
  // All calculation are done in the track plane which remove the
  // dependency of the polar angle (phi)

  Float_t dz =  s2.Z() - s1.Z();
  seg.Set( 0,
           0.5*(s1.X() + s2.X()),
           0.5*(s1.Y() + s2.Y()),
           (s2.X() - s1.X())/dz,
           (s2.Y() - s1.Y())/dz,
           s1.W() + s2.W(),
           0);
  seg.SetZ( 0.5*(s2.Z()+s1.Z()) );

  Float_t  phi = -seg.Phi();
  Double_t s   = TMath::Sin(phi);
  Double_t c   = TMath::Cos(phi);
  Float_t  tx  = seg.TX()*c-seg.TY()*s;
  Float_t  ty  = seg.TX()*s+seg.TY()*c;
  Float_t  tx1 = s1.TX()*c-s1.TY()*s;
  Float_t  ty1 = s1.TX()*s+s1.TY()*c;
  Float_t  tx2 = s2.TX()*c-s2.TY()*s;
  Float_t  ty2 = s2.TX()*s+s2.TY()*c;

  Float_t stx   = cond1.SigmaTX(tx);
  Float_t sty   = cond1.SigmaTY(ty);
  Float_t prob1 = cond1.ProbSeg(tx,ty,s1.W());
  Float_t prob2 = cond2.ProbSeg(tx,ty,s2.W());

  Float_t dtx1  = (tx1-tx)*(tx1-tx)/stx/stx/prob1;
  Float_t dty1  = (ty1-ty)*(ty1-ty)/sty/sty/prob1;
  Float_t dtx2  = (tx2-tx)*(tx2-tx)/stx/stx/prob2;
  Float_t dty2  = (ty2-ty)*(ty2-ty)/sty/sty/prob2;
  Float_t chi2a = Sqrt(dtx1+dty1+dtx2+dty2)/2.;

  seg.SetChi2(chi2a);
  return chi2a;
}

///______________________________________________________________________________
float EdbTrackFitter::Chi2PSeg( EdbSegP &s1, EdbSegP &s2, EdbSegP &seg, EdbScanCond &cond1, EdbScanCond &cond2)
{
  // fast estimation of chi2 in the special case when only the position
  // errors used - no any angular informtion:
  // application: alignment

  seg.Set( 0,
           0.5*(s1.X() + s2.X()),
           0.5*(s1.Y() + s2.Y()),
           0.5*(s1.TX() + s2.TX()),
           0.5*(s1.TY() + s2.TY()),
           s1.W() + s2.W(),
           0);
  seg.SetZ( 0.5*(s2.Z()+s1.Z()) );

  Double_t dx = s2.X()-s1.X();
  Double_t dy = s2.Y()-s1.Y();
  Float_t  sx = cond1.SigmaX(s1.Theta());  // to test!
  Float_t  sy = cond1.SigmaY(s1.Theta());  // use the same value since MCS is the main source of errors(theta) in SND conditions
  Double_t chi2 = Sqrt( dx*dx/sx/sx + dy*dy/sy/sy );
  seg.SetChi2(chi2);
  return chi2;
}

//______________________________________________________________________________
double EdbTrackFitter::ProbSegMCS(EdbSegP *s1, EdbSegP *s2){
  if(s1->P()<=0)s1->SetP(ePdef);
  return EdbPVRec::ProbeSeg(s1,s2,eX0,eM);
}
//______________________________________________________________________________
float EdbTrackFitter::Chi2SegMCS(const EdbSegP &s1, const EdbSegP &s2){
  
  // double DZemul=s1.DZem()+s2.DZem();
  ///Use hardcoded value because EdbSegP::DZem() contains sum of chi2 for microtracks.
  double DZemul=88;

  double dzem=0.5*(s1.DZ()+s2.DZ())+DZemul;
  double dz=TMath::Abs(s1.Z()-s2.Z());
  double dzpb=dz-dzem;
  if(dzpb<=0)return 0;
  double dist=EdbSegP::Distance(s1,s2);
  dist*=dzpb/dz;
  double mom=s1.P();
  if(mom<=0)mom=ePdef;
  double theta0=EdbPhysics::ThetaMCS(mom,eM,dist,eX0);
  
  double theta=TMath::Sqrt2()*EdbSegP::Angle(s1,s2);
  double chi=theta/theta0;
  Log(5,"Chi2SegMCS",Form("dist=%6.4f (%6.4f), p=%6.2f, th=%4.3g(%4.3g) => chi2=%5.4g",dist,eX0,mom,theta,theta0,chi*chi));
  return chi*chi;
}
//______________________________________________________________________________
float EdbTrackFitter::Chi2SegM( EdbSegP s1, EdbSegP s2, EdbSegP &s, EdbScanCond &cond1, EdbScanCond &cond2)
{
  // full estimation of chi2 without covariance matrix - the result seems to be identical to Chi2Seg 
  // VT: 19-Sep-2007
  //
  // Input: 2 segments passed by value because them will be modified during calculations
  // Return value - the result of the fit - passed by value

  // 1) calcualte the mean direction vector for the 2-seg group

  float dz = s2.Z() - s1.Z(); 
  float tbx=0, tby=0, wbx=0, wby=0;
  if(Abs(dz) > 0.1 ) {
    tbx = (s2.X() - s1.X())/(s2.Z() - s1.Z());
    tby = (s2.Y() - s1.Y())/(s2.Z() - s1.Z());
    wbx = Sqrt( cond1.SigmaX(tbx)*cond1.SigmaX(tbx) + cond2.SigmaX(tbx)*cond2.SigmaX(tbx) ) / dz;
    wby = Sqrt( cond1.SigmaY(tby)*cond1.SigmaY(tby) + cond2.SigmaY(tby)*cond2.SigmaY(tby) ) / dz;
  }
  float w1x = 1./cond1.SigmaTX(s1.TX());  w1x*=w1x;
  float w1y = 1./cond1.SigmaTY(s1.TY());  w1y*=w1y;
  float w2x = 1./cond2.SigmaTX(s2.TX());  w2x*=w2x;
  float w2y = 1./cond2.SigmaTY(s2.TY());  w2y*=w2y;

  float TX = (s1.TX()*w1x + s2.TX()*w2x + tbx*wbx)/(w1x+w2x+wbx);
  float TY = (s1.TY()*w1y + s2.TY()*w2y + tby*wby)/(w1y+w2y+wby);

  // 2) calcualte the COG of the 2-seg group

  w1x = 1./cond1.SigmaX(TX);    w1x *= w1x;
  w1y = 1./cond1.SigmaY(TY);    w1y *= w1y;
  w2x = 1./cond2.SigmaX(TX);    w2x *= w2x;
  w2y = 1./cond2.SigmaY(TY);    w2y *= w2y;
  float Z = (s1.Z()*(w1x+w1y) + s2.Z()*(w2x+w2y))/(w1x+w1y+w2x+w2y);
  float X = (s1.X()*(w1x+w1y) + s2.X()*(w2x+w2y))/(w1x+w1y+w2x+w2y);
  float Y = (s1.Y()*(w1x+w1y) + s2.Y()*(w2x+w2y))/(w1x+w1y+w2x+w2y);
  //printf("COG: x %f  y %f z %f\n",X,Y,Z);

  s.SetX(X);
  s.SetY(Y);
  s.SetTX(TX);
  s.SetTY(TY);
  s.SetW(s1.W()+s2.W());
  s.SetZ(Z);

  float PHI = ATan2(TY,TX);   // angle of the 2-seg group plane
  float T   = Sqrt(TX*TX+TY*TY);

  EdbAffine2D aff;
  aff.Rotate(PHI);
  //aff.ShiftX(X);
  //aff.ShiftY(Y);
  aff.Invert();
  s1.Transform(&aff);
  s2.Transform(&aff);
  s.Transform(&aff);

  float stx1   = cond1.SigmaTX(T), sty1   = cond1.SigmaTY(0);
  float stx2   = cond2.SigmaTX(T), sty2   = cond2.SigmaTY(0);
  w1x = 1./(stx1*stx1);  w1y = 1./(sty1*sty1);
  w2x = 1./(stx2*stx2);  w2y = 1./(sty2*sty2);

  //printf("w1x %f w2x %f w1y %f w2y %f\n",w1x,w2x,w1y,w2y);
  float chi2t =  ( (s1.TX()-T)*(s1.TX()-T)*w1x +
		    s1.TY()*s1.TY()        *w1y +
		   (s2.TX()-T)*(s2.TX()-T)*w2x +
		    s2.TY()*s2.TY()        *w2y )/4.;
  //printf("angular component of chi2 = %f\n",chi2t);

  float sx1   = cond1.SigmaX(T), sy1   = cond1.SigmaY(0);
  float sx2   = cond2.SigmaX(T), sy2   = cond2.SigmaY(0);

  float dx1 = s1.X()-(s.X()+(s1.Z()-s.Z())*s.TX());  
  float dy1 = s1.Y()-(s.Y()+(s1.Z()-s.Z())*s.TY());  
  float dx2 = s2.X()-(s.X()+(s2.Z()-s.Z())*s.TX());  
  float dy2 = s2.Y()-(s.Y()+(s2.Z()-s.Z())*s.TY());  

  float chi2pos = dx1*dx1/sx1/sx1+dy1*dy1/sy1/sy1+dx2*dx2/sx2/sx2+dy2*dy2/sy2/sy2/4.;
  //printf("position component of chi2 = %f\n",chi2pos);
  s.SetChi2( Sqrt(chi2t+chi2pos) );
  aff.Invert();
  s.Transform(&aff);
  //s.Print();

  return s.Chi2();
}

//______________________________________________________________________________
float EdbTrackFitter::MaxChi2Seg(EdbTrackP &t)
{
  // return the maximal seg-to-seg chi2 along the track
  float chimax=0,chi=0;
  if(t.N()<2) return chimax;
  for(int i=0; i<t.N()-1; i++) {
    EdbSegP s(*t.GetSegment(i));
    chi =  Chi2Seg(&s,t.GetSegment(i+1));
    chimax = chi > chimax ? chi : chimax;
  }
  return chimax;
}

//______________________________________________________________________________
float EdbTrackFitter::MeanChi2Seg(EdbTrackP &t)
{
  float meanchi=0;
  if(t.N()<2) return meanchi;
  for(int i=0; i<t.N()-1; i++) {
    EdbSegP s(*t.GetSegment(i));
    meanchi +=  Chi2Seg(&s,t.GetSegment(i+1));
  }
  return meanchi /= (t.N()-1);
}


//______________________________________________________________________________
float EdbTrackFitter::MaxKink(EdbTrackP &t)
{
  // return the maximal seg-to-seg kink theta angle
  float kink=0,theta=0;
  if(t.N()<2) return kink;
  for(int i=0; i<t.N()-1; i++) {
    theta =  Theta(*t.GetSegment(i),*t.GetSegment(i+1));
    kink = theta > kink ? theta : kink;
  }
  return kink;
}

//______________________________________________________________________________
float EdbTrackFitter::MeanKink(EdbTrackP &t)
{
  float meankink=0;
  if(t.N()<2) return 0;
  for(int i=0; i<t.N()-1; i++) 
    meankink +=  Theta(*t.GetSegment(i),*t.GetSegment(i+1));
  return meankink /= (t.N()-1);
}

//______________________________________________________________________________
float EdbTrackFitter::Theta( EdbSegP &s, EdbSegP &s1 )
{
  return Sqrt( (s.TX()-s1.TX())*(s.TX()-s1.TX()) + (s.TY()-s1.TY())*(s.TY()-s1.TY()) );
}

//______________________________________________________________________________
bool EdbTrackFitter::SplitTrack( EdbTrackP &t, EdbTrackP &t1, int isplit )
{
  // split track t in 2 at the point isplit - will be the first segment of of t1
  if(t.N()<isplit) return false;
  for(int i=t.N()-1; i>=isplit; i--) {
    t1.AddSegment(   t.GetSegment(i) );
    t.RemoveSegment( t.GetSegment(i) );
  }
  t.SetCounters();
  //t1.FitTrackKFS(true, X0, 0);
  t1.SetCounters();
  t1.SetM(t.M());
  t1.SetP(t.P());
  //t1.FitTrackKFS(true, X0, 0);

  return true;
}

//______________________________________________________________________________
int EdbTrackFitter::SplitTrackByKink( EdbTrackP *t, TObjArray &tracks, float maxkink )
{
  // split track t in several tracks accourding to maxkink
  // return total number of tracks after splitting; 
  // after splitting all new tracks added to the array "tracks"

  if(t->N()<1) return 0;
  if(t->N()<2) return 1;
  int nsplit=1;
  for(int i=t->N()-1; i>0; i--) 
    if( Theta(*t->GetSegment(i),*t->GetSegment(i-1)) >= maxkink ) {
      EdbTrackP *t1 = new EdbTrackP();
      SplitTrack(*t,*t1, i);
      tracks.Add(t1);
      nsplit++;
    }
  return nsplit;
}

//______________________________________________________________________________
float EdbTrackFitter::PMS_KF(EdbTrackP &t, float p0, float probbest)
{
  // select track momentum in a way to have the given chi2-probablity calculated by KF

  if(t.N()<2) return 0;
  float prob=0;
  float pu=100., pl=0., p=p0;
  int nstep = 0;
  while( Abs(prob-probbest)>0.001 ) {
    nstep++;
    t.SetP(p);
    t.FitTrackKFS(true);
    prob = t.Prob();
    if(prob<probbest) pu=p;
    else pl=p;
    if(nstep>30) break;
    p = (pu+pl)/2.;
  }
  if(nstep>20) printf("Warning in EdbTrackFitter::PMS_KF: nstep=%d     nseg=%d  p=%f  prob=%f\n",nstep,t.N(),p,prob);
  return t.P();
}

//______________________________________________________________________________
int EdbTrackFitter::FitTrackLine(EdbTrackP &tr)
{
  // track fit by averaging of segments parameters and put them as the track parameters
  float x,y,z,tx,ty,w;
  FitTrackLine(tr,x,y,z,tx,ty,w);
  tr.Set(tr.ID(),x,y,tx,ty,w,tr.Flag());
  tr.SetZ(z);
  
  tr.ClearF();
  int nseg=tr.N();
  for(int i=0; i<nseg; i++) {
     EdbSegP *s = tr.GetSegment(i);
     float dz = s->Z()-z;
     EdbSegP *segf = new EdbSegP(*s);

    segf->Set( s->ID(), x + tx*dz, y+ty*dz, tx, ty, 1.,s->Flag());
    segf->SetErrorP ( tr.SP() );
    //segf.SetChi2(0.);
    //segf.SetProb(1.);
    segf->SetP( tr.P() );
    tr.AddSegmentF(segf);
  }
  
  return tr.N();
}

//______________________________________________________________________________
int EdbTrackFitter::FitTrackLine(const EdbTrackP &tr, float &x, float &y, float &z, float &tx, float &ty, float &w)
{
  // track fit by averaging of segments parameters and return the mean values
  
  int nseg=tr.N();
  x=0; y=0; z=0; tx=0; ty=0; w=0;
  EdbSegP *seg=0;
  for(int i=0; i<nseg; i++) {
    seg = tr.GetSegment(i);
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
  return nseg;
}

//______________________________________________________________________________
int EdbTrackFitter::RMSprojXY( EdbTrackP &tr, float &ex, float &ey )
{
  // fit track by straight line using coordinates only and return RMS's
  int nseg=tr.N();
  float x[100], y[100], z[100], w[100];
  float x0,y0,z0,tx,ty;
  for(int i=0; i<nseg; i++) {
    EdbSegP *s = tr.GetSegment(i);
    x[i]=s->X();
    y[i]=s->Y();
    z[i]=s->Z();
    w[i]=s->W();
  }
  return EdbMath::LFIT3(x,y,z,w,nseg,x0,y0,z0,tx,ty,ex,ey);
}

//______________________________________________________________________________
int EdbTrackFitter::Fit3Pos( EdbTrackP &tr )
{
  // fit track by straight line using coordinates only and update track parameters
  int nseg=tr.N();
  float x[100], y[100], z[100], w[100];
  float x0,y0,z0,tx,ty,ex,ey;
  for(int i=0; i<nseg; i++) {
    EdbSegP *s = tr.GetSegment(i);
    x[i]=s->X();
    y[i]=s->Y();
    z[i]=s->Z();
    w[i]=s->W();
  }
  int rv=EdbMath::LFIT3(x,y,z,w,nseg,x0,y0,z0,tx,ty,ex,ey);
  tr.SetX(x0);
  tr.SetY(y0);
  tr.SetZ(z0);
  tr.SetTX(tx);
  tr.SetTY(ty);
//TODO errors
  return rv;
}
