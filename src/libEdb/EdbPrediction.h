#ifndef ROOT_EdbPrediction
#define ROOT_EdbPrediction
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbPrediction                                                        //
//                                                                      //
// prediction for scanning                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClonesArray.h"
#include "TNtuple.h"

#ifndef ROOT_EdbVirtual
#include "EdbVirtual.h"
#endif


//______________________________________________________________________________
class EdbPredictionDC : public TObject, public EdbTrack2D {

private:

  Int_t      eIDp;        // prediction id
  Int_t      eEvent;      // event id
  Int_t      eFlag;       // prediction flag

  Float_t    eTy;         // tanX: deltaX/deltaZ
  Float_t    eTz;         // tanY: deltaY/deltaZ
  Float_t    eYp;         // |_coordinates of the prediction in
  Float_t    eZp;         // |  absolute coord system

public:
  EdbPredictionDC() { }
  EdbPredictionDC( int id, int event, int flag, 
		   float ty,  float tz,  float y,  float z ):
    eIDp(id), eEvent(event), eFlag(flag), 
     eTy(ty), eTz(tz), eYp(y), eZp(z) {}
  virtual ~EdbPredictionDC() { }

  virtual Int_t      GetID()    const { return eIDp; }
  virtual Int_t      GetEvent() const { return eEvent; }
  virtual Int_t      GetFlag()  const { return eFlag; }
  virtual Float_t    GetTy()    const { return eTy; }
  virtual Float_t    GetTz()    const { return eTz; }
  virtual Float_t    GetYp()    const { return eYp; }
  virtual Float_t    GetZp()    const { return eZp; }

  //  virtual void       Set( float y,  float z,  float ty,  float tz );

  // mandatory virtual fuctions

  virtual Float_t    X()          const     { return GetYp(); }
  virtual Float_t    Y()          const     { return GetZp(); }
  virtual void    SetX( float x )           { eYp = x; }
  virtual void    SetY( float y )           { eZp = y; }
  
  virtual Float_t    TX()           const   { return GetTy(); }
  virtual Float_t    TY()           const   { return GetTz(); }
  virtual void    SetTX( float tx )         { eTy = tx; }
  virtual void    SetTY( float ty )         { eTz = ty; }
  
  virtual void    Print( Option_t *opt=0) const;
  void            WriteDC( FILE *file ) const;
  
  ClassDef(EdbPredictionDC,2)  // one prediction (track) a-la $c
};

//______________________________________________________________________________
class EdbPredictionsBox : public TObject, public EdbPointsBox2D {

 private:

  TClonesArray *ePredictions;  // collection of predictions

 public: 
  EdbPredictionsBox();
  EdbPredictionsBox(const EdbPredictionsBox &pb);
  EdbPredictionsBox( int n );
  virtual ~EdbPredictionsBox();

  //members access functions
  void                 AddPrediction( int id, int event, int flag, 
				        float ty,  float tz, float y,  float z);

  Int_t                 GetN()         const;
  EdbPredictionDC      *GetPrediction(int i) const;
  TClonesArray         *GetPredictions()     const { return ePredictions; }

  // mandatory virtual functions:
  Int_t         N()       const { return GetN(); }
  EdbPoint     *At(int i) const      { return (EdbPoint*)GetPrediction(i); }

  //  void  Draw(int style=23, int   col=4, float size=1. );

  // other finctions
  void       Print( Option_t *opt="") const;

  void      WriteDC( const char *file="predictions.dc") const;
  void      WriteDC1( const char *file="predictions.dc") const;

  TNtuple  *MakeNtuple( const char *name="predictions") const;

  int       ReadDC( const char *file );
  int       ReadDC0( const char *file );

  void      Generate( int n );

  //void      Test();

  ClassDef(EdbPredictionsBox,2)  // Predictions a-la $c
};

#endif /* EdbPrediction */
