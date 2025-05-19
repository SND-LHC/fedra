#ifndef ROOT_EdbSegP
#define ROOT_EdbSegP
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbSegP                                                              //
//                                                                      //
// segment with the attributes useful for processing                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TRefArray.h"
#include "TMatrixD.h"
#include "EdbVirtual.h"
#include "EdbID.h"

//______________________________________________________________________________
class EdbSegP : public TObject, public EdbTrack2D {
 
 private:
  Int_t      ePID;            // mother pattern ID
  Int_t      eID;             // segment id (unique in plate)
  Int_t      eVid[2];         // [0]-view entry in the input tree, [1]-segment entry in the view
  Int_t      eAid[2];         // [0]-AreaID, [1]-ViewID
  Int_t      eFlag;
  Int_t      eTrack;          // id of the track (-1) if no track, for a track object it represents the index of trk.root file (tracks tree)
 public:
  Float_t    eX , eY, eZ;     // coordinates
  Float_t    eTX, eTY;        // direction tangents
 private:
  Float_t    eSZ;             // square of the Z-error
  Float_t    eChi2;           // chi-square 
  Float_t    eProb;           // probability
  Float_t    eW;              // weight
  Float_t    eVolume;         // segment volume
  Float_t    eDZ;             // the length of segment along z-axis
  Float_t    eDZem;           // the length of segment along z-axis in the emulsion
  Float_t    eP;              // momentum of the particle
  Int_t      eMCTrack;        // MC track number
  Int_t      eMCEvt;          // MC event number
  TRefArray *eEMULDigitArray; //! AM+AC 27/07/07
  EdbID      eScanID;         // brick:plate:major:minor
  //  Int_t      ePlate;          // plate id
 protected: 

  TMatrixD   *eCOV;            // covariance matrix of the parameters (x,y,tx,ty,p)

 public:

  EdbSegP() {
    eEMULDigitArray =0;
    Set0();
  }
  EdbSegP(int id, float x, float y, float tx, float ty, float w=0, int flag=0);
  EdbSegP(const EdbSegP &s) {  
    eEMULDigitArray = 0;
    Set0(); Copy(s); 
  }

  virtual ~EdbSegP() { 
    SafeDelete(eCOV);
    SafeDelete(eEMULDigitArray);
  }

  //void Transform(EdbAffine2D &aff) { ((EdbTrack2D*)this)->Transform(&aff); }
  static void LinkMT(const EdbSegP* s1,const EdbSegP* s2, EdbSegP* s);
  void    PropagateToDZ( float dz );
  void    PropagateTo( float z );
  void    PropagateToCOV( float z );
  void    MergeTo( EdbSegP &s );
  Float_t ProbLink( EdbSegP &s1, EdbSegP &s2 );
  bool    IsCompatible(EdbSegP &s, float nsigx, float nsigt) const;
  bool    IsInside(float xmin, float xmax, float ymin, float ymax) const;
  bool    IsInside(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) const;

  void    addEMULDigit(TObject* a) {
    if(!eEMULDigitArray) eEMULDigitArray = new TRefArray();
    eEMULDigitArray->Add(a);
  }

  TRefArray* EMULDigitArray() const  { return eEMULDigitArray;}

  void    Set0();
  void    Copy(const EdbSegP &s);
  void    Clear() { eCOV->Clear(); }
  void    Set(int id, float x, float y, float tx, float ty, float w, int flag)
             { eID=id; eX=x; eY=y; eTX=tx; eTY=ty; eW=w; eFlag=flag; }
  void    SetErrors0();
  void    SetErrors() {SetErrors( 1.,1.,0.,.0001,.0001,1.);}
  void    SetErrors( float sx2, float sy2, float sz2, float stx2, float sty2, float sp2=1. );
  void    SetErrorsCOV( float sx2, float sy2, float sz2, float stx2, float sty2, float sp2=1. );
  
  void    SetErrorP( float sp2 ) {
    if(!eCOV) eCOV = new TMatrixD(5,5);
    (*eCOV)(4,4) = (double)sp2;
  }
  
  void    SetCOV( TMatrixD &cov) { 
    if(eCOV) eCOV->Copy(cov);
    else eCOV = new TMatrixD(cov);
  }

  bool CheckCOV() const{ //Check Covariance Matrix, is there?
    if (eCOV) return true;
    else return false;
  }

  void    ForceCOV( TMatrixD &cov) { //to correctly copy COV matrices in MakeTracksTree
    if(eCOV) *eCOV = cov;
    else eCOV = new TMatrixD(cov);
  }

  void SetCOV( double *array, int dim=5) { 
    if(!array) return;
    if(!eCOV)  eCOV = new TMatrixD(5,5);
    for(int k=0; k<dim; k++) 
      for(int l=0; l<dim; l++) (*eCOV)(k,l) = array[k*dim + l];
  }

  TMatrixD &COV() const {return *eCOV;}
  void    SetSZ( float sz )             { eSZ=sz; }
  void    SetZ( float z )               { eZ=z; }
  void    SetDZ( float dz )             { eDZ=dz; }
  void    SetDZem( float dz )           { eDZem=dz; }
  void    SetID( int id )               { eID=id; }
  void    SetPID( int pid )             { ePID=pid; }
  void    SetFlag( int flag )           { eFlag=flag; }
  void    SetTrack( int trid )          { eTrack=trid; }
  void    SetW( float w )               { eW=w; }
  void    SetP( float p )               { eP=p; }
  void    SetProb( float prob )         { eProb=prob; }
  void    SetChi2( float chi2 )         { eChi2=chi2; }
  void    SetVolume( float w )          { eVolume=w; }
  void    SetVid(int vid, int sid)      { eVid[0]=vid; eVid[1]=sid; }
  void    SetAid(int a,   int v, int side=0)        { eAid[0]=a; eAid[1]=100000*side+v; }
  void    SetSide(int side=0)        { int v=eAid[1]%100000;     eAid[1]=100000*side+v; }
  void    SetProbability( float p )     { eProb=p; }
  void    SetMC( int mEvt, int mTrack ) { eMCEvt=mEvt; eMCTrack=mTrack; }
  void    SetPlate( int plateid ) { eScanID.ePlate = plateid; }
  void    SetScanID( EdbID id ) { eScanID = id; }

  Int_t   MCEvt()   const {return eMCEvt;}
  Int_t   MCTrack() const {return eMCTrack;}
  Int_t   ID()      const {return eID;}
  Int_t   PID()     const {return ePID;}
  Int_t   Flag()    const {return eFlag;}
  Int_t   Track()   const {return eTrack;}
  Float_t W()       const {return eW;}
  Float_t P()       const {return eP;}
  Float_t Z()       const {return eZ;}
  Float_t DZ()      const {return eDZ;}
  Float_t DZem()    const {return eDZem;}
  Float_t Prob()    const {return eProb;}
  Float_t Chi2()    const {return eChi2;}
  Float_t Volume()  const {return eVolume;}
  Int_t   Plate()   const {return eScanID.ePlate;}
  EdbID   ScanID()  const {return eScanID;}

  Float_t SX()  const  { if(!eCOV) return 0; return (Float_t)(*eCOV)(0,0); }
  Float_t SY()  const  { if(!eCOV) return 0; return (Float_t)(*eCOV)(1,1); }
  Float_t STX() const  { if(!eCOV) return 0; return (Float_t)(*eCOV)(2,2); }
  Float_t STY() const  { if(!eCOV) return 0; return (Float_t)(*eCOV)(3,3); }
  Float_t SP()  const  { if(!eCOV) return 0; return (Float_t)(*eCOV)(4,4); }
  Float_t SZ()  const  { return eSZ; }
  Int_t   Vid(int i) const {return (i < 0 || i > 1) ? -1 : eVid[i];}
  Int_t   Aid(int i) const {if(i==0) return eAid[i]; else if(i==1) return eAid[i]%100000; else return -1;}
  Int_t   Side()     const {return eAid[1]/100000;}

  // mandatory virtual functions:
  Float_t X()           const {return eX;}
  Float_t Y()           const {return eY;}
  Float_t TX()          const {return eTX;}
  Float_t TY()          const {return eTY;}
  void    SetX(  Float_t x )  { eX=x; }
  void    SetY(  Float_t y )  { eY=y; }
  void    SetTX( Float_t tx ) { eTX=tx; }
  void    SetTY( Float_t ty ) { eTY=ty; }
 
  //other functions
  Float_t Phi()      const {return TMath::ATan2(eTY,eTX);}
  Float_t Theta()    const {return TMath::Sqrt(eTY*eTY+eTX*eTX);}
  static Float_t Distance(const EdbSegP &s1,const EdbSegP &s2);
  static Float_t Angle(const EdbSegP &s1,const EdbSegP &s2);

  Float_t DeltaTheta(EdbSegP* seg1)    const {
                    return TMath::Sqrt( TMath::Power(TX()-seg1->TX(),2)+TMath::Power(TY()-seg1->TY(),2) );
  }
  Float_t DeltaR(EdbSegP* seg1)    const {
                    return TMath::Sqrt( TMath::Power(X()-seg1->X(),2)+TMath::Power(Y()-seg1->Y(),2) );
  }
  void    Print( Option_t *opt="") const;
  void    PrintNice() const;

  Bool_t  IsEqual(const TObject *obj) const;
  Bool_t  IsSortable() const { return kTRUE; }
  Int_t   Compare(const TObject *obj) const;

  ClassDef(EdbSegP,18)  // segment
};

#endif /* ROOT_EdbSegP */
