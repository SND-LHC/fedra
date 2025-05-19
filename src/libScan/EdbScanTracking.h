#ifndef ROOT_EdbScanTracking
#define ROOT_EdbScanTracking

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbScanTracking                                                      //
//                                                                      //
// To handle tracking in the scanset (IO)                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TH3F.h>
#include "TEnv.h"
#include "EdbScanProc.h"
#include "EdbTrackFitter.h"

struct PredictionCut
{
  bool  doXY;             // if false - ignore position cut
  float x0,y0, dx,dy ;    // position cut: reference in xy +-d
  bool  doZ;              // if false - ignore this cut
  float z0, zmin, zmax;   // acceptance in z
  bool  doAng;            // if false - ignore this cut
  float tx,ty, dtx,dty;   // angle use x0,y0,z0 as starting point  +-d
};

//_________________________________________________________________________
class EdbTrackAssembler: public TObject {

 private:
  EdbPattern    eSegments;     // all segments of tracks
  TObjArray     eTracks;       // array of tracks (EdbTrackP) (owner of tracks)
  TObjArray     eTrZ;          // "predictions" - tracks extrapolated to the given z (not owner)
  
  Float_t       eZ;            // the z-position
  
  EdbCell2      eTrZMap;       // map of predictions at given eZ
  Float_t       eMapMarg;      // margin for the map creation
  Int_t         eCellN;        // mean cell occupancy
  
  EdbTrackFitter eFitter;
  
  public:
  float         eDTmax;        // angular acceptance for the fast preselection
  float         eDRmax;        // position acceptance for the fast preselection
  float         eDZGapMax;     // maxgap acceptance for the fast preselection
  float         eProbMin;      // min acceptable probability for segments preselection
  int           eDoUseMCS;     //flag to use MultipleScattering addition for chi2 

  int            eCollisionsRate;  //
  EdbScanCond    eCond;
 
  TH1F           *eHistProbBest;    // prob of the best candidate
  TH1F           *eHistProbAll;     // prob of all candidate
  TH1F           *eHistThetaBest;   // theta of the best candidate
  TH1F           *eHistThetaAll;    // theta of all candidate
  TH1F           *eHistNcnd;        // number of candidates after preliminary selection
  TH2F           *eHistXYseg;       // XY overlap of all segments (for showers tag)
  TH3F           *eHistXYPseg;      // XYP overlap of all segments (for showers tag)
  TH3F           *eHistXYPsegNB;    // XYP overlap of all segments (for showers tag) No Beam
  TH2F           *eHistTXTYseg;     // TXTY overlap of all segments
 
 public:
  EdbTrackAssembler();
  virtual ~EdbTrackAssembler();

  void SetMomentum (float p ){eFitter.ePdef=p;}
  void SetRadLength(float x0){eFitter.eX0=x0; eCond.SetRadX0(x0);}
  
  bool        SameSegment( EdbSegP &s1, EdbSegP &s2 );
  void        DoubletsFilterOut(EdbPattern &p);
  void        InitTrZMap( const char *str );
  void        InitTrZMap( int nx, float xmi, float xma, 
                          int ny, float ymi, float yma,  int ncell);
  void        InitTrZMap();
  void        FillTrZMap();
  void        ExtrapolateTracksToZ(float z, int nsegmin=0);
  void        ExtrapolateTracksToZ_mcs3(float z, int nsegmin=0);
  void        AddPattern(EdbPattern &p);
  EdbTrackP   *AddSegment(EdbSegP &s);                  //owner of the segments!!!
  EdbTrackP   *AddSegmentAsTrack(EdbSegP &s);
  float       ProbSeg(EdbSegP &s1, EdbSegP &s2 );
  void        RecalculateSegmentsProb(EdbTrackP &t);
  bool        AcceptDZGap(EdbTrackP &t, float z);
  void        SetSegmentsErrors();
  void        FitTracks();
  void        CombTracks( TObjArray &selected );
  void        FillXYseg(EdbPattern &p);
  
  void CheckPatternAlignment(EdbPattern &p, EdbPlateP &plate, int nsegmin, EdbAffine2D &acorr);
  
  TObjArray   &Tracks() {return eTracks;}
  
 ClassDef(EdbTrackAssembler,1) // generic class for the tracks assembling from segments
};

//_________________________________________________________________________
class EdbScanTracking: public TObject {

 public:
   int          eNsegMin;
   int          eNgapMax;
   EdbScanProc *eSproc;
   bool         eDoRealign;
   PredictionCut ePRC;
  
 public:
   EdbScanTracking();
   virtual ~EdbScanTracking(){}

   void  TrackAli(EdbPVRec &ali, TEnv &env);
   void  TrackSetBT(EdbID id, TEnv &env, Int_t ix = -1, Int_t iy = -1);
   void  SaveHist(EdbID idset, EdbTrackAssembler &etra);
   void  SetPredictionXY( const char *str  );
   void  SetPredictionZ(  const char *str  );
   void  SetPredictionAng( const char *str );
 
   ClassDef(EdbScanTracking,2) // To handle tracking in the scanset
};


#endif /* ROOT_EdbScanTracking */
