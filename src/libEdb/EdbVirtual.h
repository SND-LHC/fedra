#ifndef ROOT_EdbVirtual
#define ROOT_EdbVirtual
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbVirtual                                                           //
//                                                                      //
// collection of virtual classes and methods                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TH1F;
class TH2F;

class EdbAffine2D;
class EdbAffine3D;

class EdbPoint;
class EdbPoint2D;
class EdbPoint3D;

class EdbAngle2D;
class EdbTrack2D;

//______________________________________________________________________________
class EdbPoint {

  // protocol class for point 

 public:
  virtual ~EdbPoint() {}

  virtual Float_t    X()          const = 0;
  virtual Float_t    Y()          const = 0;
  virtual Float_t    Z()          const = 0;
  virtual void    SetX( float x )       = 0;
  virtual void    SetY( float y )       = 0;
  virtual void    SetZ( float z )       = 0;

  virtual void Transform( const EdbAffine2D *a ) {}
  virtual void Transform( const EdbAffine3D *a ) {}

  virtual void Substruct( EdbPoint *p ) = 0;

  virtual void Test() const { printf("EdbPoint::Test \n"); }

  ClassDef(EdbPoint,2)  // virtual point
};

//______________________________________________________________________________
class EdbAngle2D {

  // protocol class for 2D angle 

 public:
  virtual ~EdbAngle2D() {}

  virtual Float_t    TX()          const = 0;  // tangens = deltaX/deltaZ
  virtual Float_t    TY()          const = 0;  // 
  virtual void    SetTX( float x )       = 0;
  virtual void    SetTY( float y )       = 0;
 
  virtual void Transform( const EdbAffine2D *a );

  virtual void Substruct( const EdbAngle2D *a );

  virtual void Test() const { printf("EdbAngle2D::Test \n"); }

  virtual void Print( Option_t *opt="") const;

  ClassDef(EdbAngle2D,2)  // virtual angle
};
//______________________________________________________________________________
class EdbPoint2D : public EdbPoint {

 public:
  virtual ~EdbPoint2D() {}

  virtual Float_t    X()          const = 0;
  virtual Float_t    Y()          const = 0;
  virtual Float_t    Z()          const  { return 0; }
  virtual void    SetX( float x )       = 0;
  virtual void    SetY( float y )       = 0;
  virtual void    SetZ( float z )       {}

  virtual void Transform( const EdbAffine2D *a );

  virtual void Substruct( EdbPoint *p );

  virtual void Print( Option_t *opt="") const;

  virtual void Test() const { printf("EdbPoint2D::Test \n"); }
  virtual void TestPoint2D() const { printf("EdbPoint2D::TestPoint2D \n"); }

  ClassDef(EdbPoint2D,2)  // virtual 2D point
};

//______________________________________________________________________________
class EdbPoint3D : public EdbPoint2D {

 public:
  virtual ~EdbPoint3D() {}

  virtual Float_t    Z()          const = 0;
  virtual void    SetZ( float z )       = 0;

  virtual void Transform( const EdbAffine3D *a );

  virtual void Substruct( EdbPoint *p );

  virtual void Print( Option_t *opt="") const;

  virtual void Test() const { printf("EdbPoint3D::Test \n"); }
  virtual void TestPoint3D() const { printf("EdbPoint3D::TestPoint3D \n"); }

  ClassDef(EdbPoint3D,2)  // virtual 3D point
};


//______________________________________________________________________________
class EdbTrack2D : public EdbPoint2D, public EdbAngle2D {

  // virtual class for 2D track

 public:
  virtual ~EdbTrack2D() {}

  virtual void Transform( const EdbAffine2D *a );

  virtual void Substruct( EdbTrack2D *t );

  virtual void Test() const { printf("EdbTrack2D::Test \n"); }

  virtual void Print( Option_t *opt="") const;

  ClassDef(EdbTrack2D,2)  // virtual track
};

//______________________________________________________________________________
class EdbPointsBox2D : public EdbPoint3D {

  // abstrackt box with geometrical elements collection inside
  // 2 members are assumed:
  // *collection - pointer to the collection of objects
  // *affine     - affine transformation with readonly access from outside
  //               to keep all transformation done with the box
  //               function Retransform will recalculate original data

  EdbAffine2D  *eKeep;   // affine transformation to keep all operations 
                         // done with data

 public:
  EdbPointsBox2D();
  EdbPointsBox2D(const EdbPointsBox2D &pb);
  virtual ~EdbPointsBox2D();

  virtual void Print( Option_t *opt="") const;

  virtual Int_t       N()            const = 0;
  virtual EdbPoint   *At( int i )     const     = 0;

  virtual Float_t    X()  const        {return 0;};
  virtual Float_t    Y()  const        {return 0;};
  virtual Float_t    Z()  const        {return 0;};
  virtual void    SetX( float x ) {};
  virtual void    SetY( float y ) {};
  virtual void    SetZ( float z ) {};

  virtual void Transform( const EdbAffine2D *a );
  virtual void Substruct( EdbPointsBox2D *b );
  virtual void Retransform();

  virtual const EdbAffine2D *GetKeep() const {return eKeep;}
  virtual void GetKeep( EdbAffine2D &aff );
  virtual void SetKeep(float a11,float a12,float a21,float a22,float b1,float b2);

  virtual Float_t Xmin() const;
  virtual Float_t Xmax() const;
  virtual Float_t Ymin() const;
  virtual Float_t Ymax() const;

  virtual Float_t DeltaX() { return Xmax()-Xmin(); }
  virtual Float_t DeltaY() { return Ymax()-Ymin(); }

  virtual void  ShiftX( float offset ){}
  virtual void  ShiftY( float offset ){}
  virtual void  ScaleX( float scaleFactor ){}
  virtual void  ScaleY( float scaleFactor ){}
  virtual void  Rotate( float angle );
  virtual void  SmearXY( float sigmax, float sigmay ){}


  virtual TH1F  *Xhist();
  virtual TH1F  *Yhist();
  virtual TH2F  *XYhist();

  virtual void  DrawPoints(int style=23, int   col=4, float size=1., bool same=false);

  ClassDef(EdbPointsBox2D,2)  // collection of  2D points
};

//______________________________________________________________________________
class EdbPointsBox3D : public EdbPointsBox2D {

 public:
  EdbPointsBox3D(): EdbPointsBox2D() {}
  virtual ~EdbPointsBox3D() {}

  ClassDef(EdbPointsBox3D,2)  // collection of  3D points
};

#endif /* ROOT_EdbVirtual */
