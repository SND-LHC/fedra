//-- Author :  Valeri Tioukov   13.06.2000

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbVirtual                                                           //
//                                                                      //
// Virtual algorithms                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TBuffer.h>
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TClass.h"

#include "EdbVirtual.h"
#include "EdbAffine.h"

ClassImp(EdbPoint)
ClassImp(EdbPoint2D)
ClassImp(EdbPoint3D)
ClassImp(EdbPointsBox2D)
ClassImp(EdbPointsBox3D)
ClassImp(EdbAngle2D)
ClassImp(EdbTrack2D)

//_____________________________________________________________________________
void EdbTrack2D::Print( Option_t *opt ) const
{
  printf("EdbTrack2D: %f %f %f %f \n",X(), Y(), TX(),TY() );
}

//_____________________________________________________________________________
void EdbTrack2D::Transform( const EdbAffine2D *a )
{
  EdbPoint2D::Transform(a);
  EdbAngle2D::Transform(a);
}

//_____________________________________________________________________________
void EdbTrack2D::Substruct( EdbTrack2D *t )
{
  EdbPoint2D::Substruct( (EdbPoint*)t );
  EdbAngle2D::Substruct( (EdbAngle2D*)t );
}

//_____________________________________________________________________________
void EdbAngle2D::Print( Option_t *opt ) const
{
  printf("EdbAngle2D: %f %f \n",TX(),TY() );
}

//_____________________________________________________________________________
void EdbAngle2D::Transform( const EdbAffine2D *a )
{
  Float_t tx =  a->A11()*TX() + a->A12()*TY();
  Float_t ty =  a->A21()*TX() + a->A22()*TY();
  SetTX(tx);
  SetTY(ty);
}

//_____________________________________________________________________________
void EdbAngle2D::Substruct( const EdbAngle2D *p )
{
  Float_t tx =  TX() - p->TX();
  Float_t ty =  TY() - p->TY();
  SetTX(tx);
  SetTY(ty);
}

//_____________________________________________________________________________
void EdbPoint3D::Print( Option_t *opt ) const
{
  printf("EdbPoint3D: %f %f %f\n",X(),Y(),Z() );
}

//_____________________________________________________________________________
void EdbPoint3D::Transform( const EdbAffine3D *a )
{
  Float_t x =  a->A11()*X() + a->A12()*Y() + a->A13()*Z() + a->B1();
  Float_t y =  a->A21()*X() + a->A22()*Y() + a->A23()*Z() + a->B2();
  Float_t z =  a->A31()*X() + a->A32()*Y() + a->A33()*Z() + a->B3();
  SetX(x);
  SetY(y);
  SetZ(z);
}

//_____________________________________________________________________________
void EdbPoint3D::Substruct( EdbPoint *p )
{
  Float_t x =  X() - p->X();
  Float_t y =  Y() - p->Y();
  Float_t z =  Z() - p->Z();
  SetX(x);
  SetY(y);
  SetY(z);
}

//_____________________________________________________________________________
void EdbPoint2D::Print( Option_t *opt ) const
{
  printf("EdbPoint2D: %f %f \n",X(),Y() );
}

//_____________________________________________________________________________
void EdbPoint2D::Transform( const EdbAffine2D *a )
{
  Float_t x =  a->A11()*X() + a->A12()*Y() + a->B1();
  Float_t y =  a->A21()*X() + a->A22()*Y() + a->B2();
  SetX(x);
  SetY(y);
}

//_____________________________________________________________________________
void EdbPoint2D::Substruct( EdbPoint *p )
{
  Float_t x =  X() - p->X();
  Float_t y =  Y() - p->Y();
  SetX(x);
  SetY(y);
}

//_____________________________________________________________________________
EdbPointsBox2D::EdbPointsBox2D()
{
  eKeep = new EdbAffine2D();
}

//_____________________________________________________________________________
EdbPointsBox2D::EdbPointsBox2D(const EdbPointsBox2D &pb)
{
  eKeep = new EdbAffine2D(*pb.GetKeep());
}

//_____________________________________________________________________________
EdbPointsBox2D::~EdbPointsBox2D()
{
  if(eKeep) delete eKeep;
}

//_____________________________________________________________________________
void EdbPointsBox2D::SetKeep(float a11,float a12,float a21,
			     float a22,float b1,float b2)
    {eKeep->Set(a11,a12,a21,a22,b1,b2);}

//_____________________________________________________________________________
void EdbPointsBox2D::GetKeep( EdbAffine2D &a )
{
  a.Set( eKeep->A11(),eKeep->A12(),eKeep->A21(),eKeep->A22(),
	 eKeep->B1(),eKeep->B2() );
}

//_____________________________________________________________________________
void EdbPointsBox2D::Transform( const EdbAffine2D *a )
{
  int n=N();
  for( int i=0; i<n; i++ )  At(i)->Transform(a);
  eKeep->Transform(a);
}

//_____________________________________________________________________________
void EdbPointsBox2D::Rotate( float angle )
{
  EdbAffine2D aff(1,0,0,1,0,0);
  aff.Rotate(angle);
  Transform(&aff);
}

//_____________________________________________________________________________
void EdbPointsBox2D::Retransform()
{
  eKeep->Invert();
  Transform(eKeep);
  eKeep->Reset();
}

//_____________________________________________________________________________
void EdbPointsBox2D::Substruct( EdbPointsBox2D *b )
{
  int n=N();
  for( int i=0; i<n; i++ )  At(i)->Substruct( b->At(i) );
}

//_____________________________________________________________________________
Float_t EdbPointsBox2D::Xmin() const
{
  int n=N();
  if( n<1 )                 return 0;
  Float_t a = At(0)->X();
  for( int i=0; i<n; i++ )  a = TMath::Min(a,At(i)->X());
  return a;
}

//_____________________________________________________________________________
Float_t EdbPointsBox2D::Xmax() const
{
  int n=N();
  if( n<1 )                 return 0;
  Float_t a = At(0)->X();
  for( int i=0; i<n; i++ )  a = TMath::Max(a,At(i)->X());
  return a;
}

//_____________________________________________________________________________
Float_t EdbPointsBox2D::Ymin() const
{
  int n=N();
  if( n<1 )                 return 0;
  Float_t a = At(0)->Y();
  for( int i=0; i<n; i++ )  a = TMath::Min(a,At(i)->Y());
  return a;
}

//_____________________________________________________________________________
Float_t EdbPointsBox2D::Ymax() const
{
  int n=N();
  if( n<1 )                 return 0;
  Float_t a = At(0)->Y();
  for( int i=0; i<n; i++ )  a = TMath::Max(a,At(i)->Y());
  return a;
}

//_____________________________________________________________________________
TH1F *EdbPointsBox2D::Xhist()
{
  TH1F *hist = new TH1F("xhist","X coordinate of the box",100, Xmin(),Xmax() );
  int n=N();
  for( int i=0; i<n; i++ )  hist->Fill( At(i)->X() );
  return hist;
}

//_____________________________________________________________________________
TH1F *EdbPointsBox2D::Yhist()
{
  TH1F *hist = new TH1F("yhist","Y coordinate of the box",100, Ymin(),Ymax() );
  int n=N();
  for( int i=0; i<n; i++ )  hist->Fill( At(i)->Y() );
  return hist;
}

//_____________________________________________________________________________
TH2F *EdbPointsBox2D::XYhist()
{
  TH2F *hist = new TH2F("xyhist","X coordinate of the box",
			100, Xmin(), Xmax(),
			100, Ymin(), Ymax() );
  int n=N();
  for( int i=0; i<n; i++ )  hist->Fill( At(i)->X(),  At(i)->Y() );
  return hist;
}

//_____________________________________________________________________________
void EdbPointsBox2D::DrawPoints( int   style, int   col, float size, bool same )
{
  TGraph  *mark = new TGraph( N() );

  int n=N();
  for( int i=0; i<n; i++ )  
    mark->SetPoint( i, At(i)->X(),  At(i)->Y() );

  printf("EdbPointsBox2D::Draw: %d points\n", N() );

  mark->SetMarkerStyle(style);
  mark->SetMarkerColor(col);
  mark->SetMarkerSize(size);
  mark->SetTitle("Points box");
  if(same) mark->Draw("P");
  else     mark->Draw("AP");
}

//_____________________________________________________________________________
void EdbPointsBox2D::Print( Option_t *opt ) const
{
  printf("EdbPointsBox2D: %d elements \n", N() );
}

//_____________________________________________________________________________
void EdbPoint2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbPoint2D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbPoint2D::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      EdbPoint::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, EdbPoint2D::IsA());
      //====end of old versions
   } else {
     EdbPoint2D::Class()->WriteBuffer(R__b,this);
   }
}

//_____________________________________________________________________________
void EdbPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbPoint.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbPoint::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      R__b.CheckByteCount(R__s, R__c, EdbPoint::IsA());
      //====end of old versions
   } else {
     EdbPoint::Class()->WriteBuffer(R__b,this);
   }
}

//_____________________________________________________________________________
void EdbPointsBox2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbPointsBox2D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbPointsBox2D::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      EdbPoint3D::Streamer(R__b);
      R__b >> eKeep;
      R__b.CheckByteCount(R__s, R__c, EdbPointsBox2D::IsA());
      //====end of old versions
   } else {
     EdbPointsBox2D::Class()->WriteBuffer(R__b,this);
   }
}


//_____________________________________________________________________________
void EdbPointsBox3D::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbPointsBox3D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbPointsBox3D::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      EdbPointsBox2D::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, EdbPointsBox3D::IsA());
      //====end of old versions
   } else {
     EdbPointsBox3D::Class()->WriteBuffer(R__b,this);
   }
}


//_____________________________________________________________________________
void EdbPoint3D::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbPoint3D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbPoint3D::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      EdbPoint2D::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, EdbPoint3D::IsA());
      //====end of old versions
   } else {
     EdbPoint3D::Class()->WriteBuffer(R__b,this);
   }
}


//_____________________________________________________________________________
void EdbAngle2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbAngle2D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbAngle2D::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      R__b.CheckByteCount(R__s, R__c, EdbAngle2D::IsA());
      //====end of old versions
   } else {
     EdbAngle2D::Class()->WriteBuffer(R__b,this);
   }
}


//_____________________________________________________________________________
void EdbTrack2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbTrack2D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbTrack2D::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      EdbPoint2D::Streamer(R__b);
      EdbAngle2D::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, EdbTrack2D::IsA());
      //====end of old versions
   } else {
     EdbTrack2D::Class()->WriteBuffer(R__b,this);
   }
}

//______________________________________________________________________________
