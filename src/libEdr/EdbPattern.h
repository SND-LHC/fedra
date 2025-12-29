#ifndef ROOT_EdbPattern
#define ROOT_EdbPattern
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbPattern                                                           //
//                                                                      //
// pattern of segments                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TSortedList.h"
#include "TArrayL.h"
#include "TClonesArray.h"
#include "EdbSegP.h"

class EdbAffine2D;
class TIndexCell;
class EdbVTA;
class EdbVertex;


//______________________________________________________________________________
class EdbSegmentsBox : public TObject, public EdbPointsBox2D {
 
 private:
 
  Float_t eX,eY,eZ;         // central point

  TClonesArray *eSegments;  // collection of segments (EdbSegP)

  Float_t eDZkeep;          // eDZkeep = eZ - Zoriginal (before any projections)

  //  void (*PF_COUPLE)(EdbSegP *seg1,EdbSegP *seg2); //pointer to the function

  //protected:
 public:
  void       SetZ(float z) { eZ=z; }

 public:
  EdbSegmentsBox(int nseg=0);
  EdbSegmentsBox(EdbSegmentsBox &box)
    {
      eX=box.X();
      eY=box.Y();
      eZ=box.Z();
      int n=box.N();
      eSegments = new TClonesArray("EdbSegP",n);
      for(int i=0; i<n; i++) 
	AddSegment( *(box.GetSegment(i)) );
    }

  EdbSegmentsBox(float x0, float y0, float z0, int nseg=0);
  virtual ~EdbSegmentsBox();
 
  EdbSegP  *AddSegment(int i, EdbSegP &s);
  EdbSegP  *AddSegment(EdbSegP &s);
  EdbSegP  *AddSegment(EdbSegP &s1, EdbSegP &s2); 
  EdbSegP  *AddSegment(int id, float x, float y, float tx, float ty, 
		      float w=0, int flag=0);
  EdbSegP  *AddSegmentNoDuplicate( EdbSegP &s );
  Int_t     GetN() const {return eSegments ? eSegments->GetEntriesFast() : 0;}
  EdbSegP  *GetSegment(int i) const {return (EdbSegP*)eSegments->At(i); }
  EdbSegP  *GetSegmentLast()  const {return GetN() ? (EdbSegP*)eSegments->Last() : 0;}

  TClonesArray *GetSegments()     const { return eSegments; }
  void         *GetSegmentsAddr()       { return &eSegments; }
  void Expand(int n) {eSegments->Expand(n);}

  void      SetSegmentsZ();
  void      SetSegmentsDZ(float dz);
  void      SetSegmentsW(float w);
  void      SetSegmentsFlag(int flag);

  void      SetSegmentsPlate(int plate);

  // mandatory virtual functions:

  void      SetX(float x) { eX=x; }
  void      SetY(float y) { eY=y; }

  Float_t   X()       const {return eX;}
  Float_t   Y()       const {return eY;}
  Float_t   Z()       const {return eZ;}
  Float_t   DZ()      const {return eDZkeep;}
  Int_t     N()       const {return GetN();}
  EdbPoint *At(int i) const {return (EdbPoint*)GetSegment(i);}
 
  // other functions

  void    Print( Option_t *opt="") const;
  void    Set0();
  void    Reset();
  void    ProjectTo(const float dz);
  void    TransformA(    const EdbAffine2D *affA );
  void    TransformARot( const EdbAffine2D *affA );
  void    TransformShr(  const float shr );
  Int_t   CalculateXY(  EdbSegmentsBox *p , EdbAffine2D *aff  );
  Int_t   CalculateAXAY(EdbSegmentsBox *p , EdbAffine2D *affA );
  Float_t DiffAff( EdbAffine2D *aff );
  Float_t Diff( EdbSegmentsBox &p );
  Float_t GetSize(Int_t XorY);
  Float_t GetSizeX() { return GetSize(0); }
  Float_t GetSizeY() { return GetSize(1); }
  Float_t GetSizeXY();
  Float_t GetTrackDensity();
  Float_t GetTrackDensitymm2();
  void RemovePosMargins(float  x0, float  y0, float  dx, float  dy);
  void RemoveAngMargins(float tx0, float ty0, float dtx, float dty);

  ClassDef(EdbSegmentsBox,1)  // collection of segments
};

//______________________________________________________________________________
class EdbTrackP : public EdbSegP {
 
 private:

  TSortedList *eS;     // array of segments
  TSortedList *eSF;    // array of fitted segments
  Int_t        eNpl;   // number of plates passed through
  Int_t        eN0;    // number of holes (if any)
  Float_t      eM;     // invariant mass of the particle
  Float_t      eDE;    // total energy loss of the particle between first and last segments
  Int_t        ePDG;   // particle ID from PDG
  Float_t      ePerrUp;    // error of P() in upper direction, obtained by MCS,or shower-algorithm
  Float_t      ePerrDown;    // error of P() in lower direction, obtained by MCS,or shower-algorithm
  EdbVTA      *eVTAS;  //! vertex track start is attached to
  EdbVTA      *eVTAE;  //! vertex track end is attached to

 public:

  EdbTrackP(int nseg=0);
  EdbTrackP(EdbSegP &seg);
  EdbTrackP(EdbSegP *seg, float m=0.12);
  EdbTrackP(EdbTrackP &track) { Set0(); Copy(track); }

  virtual ~EdbTrackP();

  void       Set0();
  void       SetOwner() { if(eS) eS->SetOwner(true); }
  void       AddVTA(EdbVTA *vta);
  void       ClearVTA();
  void       ClearVTA(EdbVTA *vta);
  EdbVTA    *VTAS() const {return eVTAS;}
  EdbVTA    *VTAE() const {return eVTAE;}
  EdbVertex *VertexS();
  EdbVertex *VertexE();
  EdbVertex *Vertex(int zpos) {return zpos? VertexS(): VertexE();}

  void    SetPDG( int pdg )  { ePDG=pdg; }
  Int_t   PDG()      const {return ePDG;}
  
  Float_t Wmean()      const;

  void    SetM( float m )  { eM=m; }
  Float_t M()      const {return eM;}
  

  //void    SetErrorsCOV(EdbScanCond &cond);
  void    SetCounters() { SetNpl(); SetN0(); }

  void    SetN0( int n0 )  { eN0 = n0; }
  void    SetN0()          { eN0 = eNpl-N(); }
  Int_t   N0()      const  {return eN0;}

  void    SetDE( float de )  { eDE=de; }
  Float_t DE()      const  {return eDE;}

  void    SetNpl( int npl )  { eNpl=npl; }
  void    SetNpl()  
    { if(eS) eNpl = 1+TMath::Abs(GetSegment(0)->PID() - GetSegment(N()-1)->PID()); }
  Int_t   Npl()      const   {return eNpl;}

  //  TList *S()  const { return eS; }
  //  TList *SF() const { return eSF; }

  //int      N() const  { if(eS)  return eS->GetEntries(); else return 0; } //TODO fast
  Int_t    N()  const  { return (eS) ?  eS->GetSize()  : 0; }
  Int_t    NF() const  { return (eSF)?  eSF->GetSize() : 0; }
  
  Double_t  GetBTEfficiency();

  Float_t  Wgrains() const;
  Int_t	   GetSegmentsFlag( int &nseg ) const;
  Int_t	   GetSegmentsAid( int &nseg ) const;
  Int_t	   GetSegmentsMCTrack( int &nseg ) const;

  EdbSegP *GetSegmentWithClosestZ( float z, float dz);
  
  EdbSegP *GetSegmentFirst()  const {return (eS) ? (EdbSegP*)(eS->First()) : 0;}
  EdbSegP *GetSegmentLast()   const {return (eS) ? (EdbSegP*)(eS->Last()) : 0;}

  EdbSegP *GetSegmentFFirst() const {return (eSF) ? (EdbSegP*)(eSF->First()) : 0;}
  EdbSegP *GetSegmentFLast()  const {return (eSF) ? (EdbSegP*)(eSF->Last()) : 0;}

  EdbSegP *GetSegment(int i)  const {return (eS) ? (EdbSegP*)(eS->At(i)) : 0; }
  EdbSegP *GetSegmentF(int i) const {return (eSF) ? (EdbSegP*)(eSF->At(i)) : 0;}

  EdbSegP *TrackZmin(bool usesegpar=false) const { if(usesegpar || (!eSF)) return GetSegmentFirst(); else return GetSegmentFFirst(); }
  EdbSegP *TrackZmax(bool usesegpar=false) const { if(usesegpar || (!eSF)) return GetSegmentLast();  else return GetSegmentFLast(); }

  EdbSegP *TrackExtremity(bool zpos, bool usesegpar=false) 
    const { return zpos? TrackZmin(usesegpar) : TrackZmax(usesegpar); } // 0-end, 1-start (as in vertex class)

  const EdbSegP *TrackStart() const;
  const EdbSegP *TrackEnd()   const;

  Int_t   Dir()    const {return (DZ()<0) ? -1 : 1;}
  Float_t Zmax()   const;
  Float_t Zmin()   const;
  Float_t Zstart() const {return TrackStart()->Z();}
  Float_t Zend()   const {return TrackEnd()->Z();}
  void    AddTrack(const EdbTrackP &tr);

  void  AddSegment(EdbSegP *s)
  { 
    if(!eS) eS = new TSortedList();
    eS->Add(s);
  }
  void  RemoveSegment(EdbSegP *s)
  { 
    if(!eS) return;
    eS->Remove(s);
    s->SetTrack(-1);
    SetCounters();
  }
  void  SubstituteSegment( EdbSegP *sold, EdbSegP *snew )
  { 
    if(!eS) return;
    eS->Remove(sold);
    eS->Add(snew);
    sold->SetTrack(-1);
  }
  void  AddSegmentF(EdbSegP *s)  
  { 
    if(!eSF) 
      {
        eSF = new TSortedList();
	eSF->SetOwner();
      }
    eSF->Add(s);
  }

  int   RemoveAliasSegments();
  int   CheckMaxGap();
  int   CheckAliasSegments();
  int   SetSegmentsP(float p) {for(int i=0; i<N(); i++) GetSegment(i)->SetP(p); return N();}
  int   SetSegmentsTrack(int id) {for(int i=0; i<N(); i++) GetSegment(i)->SetTrack(id); return N();}
  int   SetSegmentsTrack() {return SetSegmentsTrack(ID());}
  int   FitTrackKFS( bool zmax=false, float X0=5810., int design = 0 );
  int   MakeSelector( EdbSegP &ss, bool followZ=true );
  int   EstimatePositionAt( Float_t z, EdbSegP &ss );
  float MakePredictionTo( Float_t z, EdbSegP &ss );
  float CHI2();
  float CHI2F();
  void  FitTrack();
  void  Copy(const EdbTrackP &tr);
  void  Transform(const EdbAffine2D &tr);
  
  void SetPerrUp(Float_t perrUp);
  void SetPerrDown(Float_t perrDown);
  void SetPerr(Float_t perrDown, Float_t perrUp);
  Float_t PerrUp()      const {return ePerrUp;}
  Float_t PerrDown()      const {return ePerrDown;}
  

  void  Clear()  { if(eS)  eS->Clear(); if(eSF) eSF->Clear(); }
  void  ClearF() { if(eSF) eSF->Clear(); }
  void  Print();
  void  PrintNice();

  ClassDef(EdbTrackP,7)  // track consists of segments
};

//______________________________________________________________________________
class EdbPattern : public EdbSegmentsBox {
 
 private:

  Int_t       eID;             // pattern id in the volume
  Int_t       ePID;            // correspond to the piece ID
  TIndexCell *eCell;           //! associated with eSegments
  Float_t     eStepX ,eStepY;  // bin size for the index cell 
  Float_t     eStepTX,eStepTY; // bin size for the index cell 
  Float_t     eSigma[4];       // accuracy in comparison to neibour pattern
  Int_t       eFlag;           // pattern flag
  Int_t       eNAff;           // number of segments selected for affine calculation
  EdbID       eScanID;         //! main scanning ID for this pattern, defined when possible
  Int_t       eSide;           //! emulsion side 0/1/2, defined when possible

 public:

  EdbPattern();
  EdbPattern( EdbPattern &p );
  EdbPattern(float x0, float y0, float z0, int n=0 );
  virtual ~EdbPattern();
 
  void  Set0();
  void  SetScanID(EdbID id) {eScanID=id;}
  void  SetSigma(float sx,float sy,float stx,float sty) 
                {eSigma[0]=sx; eSigma[1]=sy; eSigma[2]=stx; eSigma[3]=sty;}
  void  SetStep(float stepx, float stepy, float steptx, float stepty ) 
                {eStepX=stepx; eStepY=stepy; eStepTX=steptx; eStepTY=stepty;}
  void  Reset();
  float SummaryPath();
  void  FillCell( float stepx, float stepy, float steptx, float stepty );
  int   FindCompliments(EdbSegP &s, TObjArray &arr, float nsig, float nsigt);
  void  SetSegmentsPID();
  void  SetSegmentsScanID(EdbID id);
  EdbPattern *ExtractSubPattern(float min[5], float max[5], int MCevt=-1);
  EdbPattern *ExtractSubPatternXY(float xmin, float xmax, float ymin, float ymax);
  Int_t       AddPattern(EdbPattern &p);

  void  SetID(int id)   {eID=id;}
  void  SetPID(int pid) {ePID=pid;}
  void  SetNAff(int n)  {eNAff=n;}
  void  SetSide(int side)  {eSide=side;}

  Int_t       NAff()   const {return eNAff;}
  float       StepX()  const {return eStepX;}
  float       StepY()  const {return eStepY;}
  float       StepTX() const {return eStepTX;}
  float       StepTY() const {return eStepTY;}
  int         ID()     const {return eID;}
  int         PID()    const {return ePID;}
  TIndexCell *Cell()   const {return eCell;}
  EdbSegP    *FindSegment(int id);
  float      Xmean();
  float      Ymean();
//  float      Xmin() const;
//  float      Xmax() const;
//  float      Ymin() const;
//  float      Ymax() const;

  EdbID      ScanID() const {return eScanID;}
  Int_t      Plate() const {return eScanID.ePlate;}
  Int_t      Brick() const {return eScanID.eBrick;}
  Int_t      Side()  const {return eSide;}

  ClassDef(EdbPattern,2)  // pattern of segments
};

//______________________________________________________________________________
class EdbPatternsVolume : public TObject {
 
 private:
 
  Float_t    eX,eY,eZ;       // central point

 public:

  TObjArray  *ePatterns;      // collection of patterns

  TIndexCell *eTracksCell;    //! "vidt:vids" - connected segments cell
  TIndexCell *ePatternsCell;  //! "pid:id1:chi2:id2" - all found couples

  Bool_t      eDescendingZ;    // if =0 - z increase in the pattrens array; if =1 - decrease 

 public:
  EdbPatternsVolume();
  EdbPatternsVolume(EdbPatternsVolume &pvol);
  virtual ~EdbPatternsVolume();

  void Set0();
  void SetPatternsID();

  void Transform( const EdbAffine2D *aff );
  void Shift(float x, float y);
  void Centralize();
  void Centralize(float xc, float yc);
  void SetXYZ( float x, float y, float z ) { eX=x; eY=y; eZ=z; }

  Float_t X() const {return eX;}
  Float_t Y() const {return eY;}
  Float_t Z() const {return eZ;}
  Int_t   Npatterns() const { if(ePatterns) 
                        return ePatterns->GetEntriesFast();
                        else return 0; }
  Float_t      Xmean();
  Float_t      Ymean();

  void PassProperties(EdbPatternsVolume &pvol);

  void         AddPattern( EdbPattern *pat );
  void         AddPatternAt( EdbPattern *pat, int id );
  EdbPattern  *GetPattern( int id ) const;
  EdbPattern  *GetPatternZLowestHighest(Bool_t lowestZ=kTRUE) const;
	EdbPattern* GetPatternPreceding(EdbPattern* pat) const;
	EdbPattern* GetPatternSucceding(EdbPattern* pat) const;
	EdbPattern* NextPattern(float z, int dir) const;
  EdbPattern* GetPatternNext(float z, int dir) const;
	EdbPattern* GetPatternByPID(int pid) const;
// 	EdbPattern* GetPatternByZ(float z) const;
	EdbPattern* GetPatternByZ(float z, float tolerance=5) const;

  void   DropCell();

  void   PrintStat( Option_t *opt="") const;
  void   PrintStat( EdbPattern &pat ) const;
  void   PrintAff() const;
  int    DropCouples();

  Long_t Vid(int pid, int sid) const { return pid*1000000+sid; }
  Int_t  Pid(Long_t vid)       const { return vid/1000000; }
  Int_t  Sid(Long_t vid)       const { return vid%1000000; }

  EdbSegP *GetSegment(Long_t vid) const 
    {return GetPattern(Pid(vid))->GetSegment( Sid(vid) );}

  int FindComplimentsVol(EdbSegP &s, TObjArray &arr, float nsig, float nsigt, int dpat);

  EdbPattern  *GetPatternByPlate(int plate, int side);
  EdbPattern  *InsertPattern(EdbPattern *pat, Bool_t descendingZ=0);
  void         SortPatternsByZ(Bool_t descendingZ=0);

  void   Print() const;

  ClassDef(EdbPatternsVolume,1)  // patterns nostri
};

#endif /* ROOT_EdbPattern */
