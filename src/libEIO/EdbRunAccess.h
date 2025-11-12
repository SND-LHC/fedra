#ifndef ROOT_EdbRunAccess
#define ROOT_EdbRunAccess

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbRunAccess                                                         //
//                                                                      //
// OPERA data Run Access helper class                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TH2F.h"
#include "TGraph2D.h"
#include "TCut.h"
#include "EdbRun.h"
#include "EdbCell2.h"
#include "EdbPattern.h"
#include "EdbLayer.h"
#include "EdbSegmentCut.h"
#include "EdbScanCond.h"

//______________________________________________________________________________
class EdbRunAccess : public TObject {

 public:
  Int_t     eAFID;        // if =1 - use affine transformations of the fiducial marks
  Int_t     eCLUST;       // 1-use clusters, 0 - do not
  Bool_t    eUseExternalSurface;  // if true - set segment position corrisponding to the very external cluster

  Bool_t    eDoViewAnalysis;     // fill or not the histograms for optional view analysis
  EdbH2     eHViewXY[3];         // XY segments distribution in a view local coords
  Bool_t    eInvertSides;        // 0 -do nothing, 1-invert sides
  Int_t     eWeightAlg;          // 0-puls, 1 - density(former eUseDensityAsW), 2-Likelyhood (eSigmaX)
            //eUseDensityAsW;      // in case of LASSO tracking possible to use eSigmaY as eW
            
  Bool_t   eUseDensityAsW;      // in case of LASSO tracking possible to use eSigmaY as eW  
  // re-taken in by FWM 13072016 (otherwithe compilatiion error in libScan EdbRunAccess)
  
  Int_t       eDoImageCorr;      // apply or not pix/mic correction  when read data (default is 0)
  EdbAffine2D eImageCorr[3];     // include pix/micron corr, shift, rotation of microscope image
                                 // side 0(dummy), 1, 2
                                 
  Int_t       eDoImageMatrixCorr;      // apply image correction matrix
  EdbCell2    eCorrMap[3];             // map of corrections for sides 0(dummy), 1,2

  TCut         eHeaderCut;  // header cut to be applied in run initialization
  Int_t        eTracking;   // to test tracking alorithm: -1-ignored(def),0/1 - trackings to accept
  
  TGraph2D    *eGraphZ[4];    //! keep z1/z2/z3/z4 surfaces using eZ1/eZ2/eZ3/eZ4 for each view
  TGraph2D    *eGraphDZ[3];   //! keep the base/layer1/layer2 thickness calculated using eZ1/eZ2/eZ3/eZ4 for each view
   
 private:
  TString  eRunFileName;
  EdbRun   *eRun;        // pointer to the run to be accessed
  
  EdbAffine2D *eAffStage2Abs; // affine transformation extracted from Marks (if AFID=11)

  Int_t eFirstArea;
  Int_t eLastArea;
  Int_t eNareas;
  Int_t eNviewsPerArea;

  EdbLayer    *eLayers[3]; // base(0),up(1),down(2)  layers
  EdbScanCond *eCond[3];   //
  TObjArray   *eCuts[3];   // arrays of cuts to be applied to segments
  EdbPattern  *eVP[3];     // base/up/down side patterns (0,1,2) with 
                           // the views coordinates inside

  Float_t      eCutGR;     // grain cut (chi)

  Float_t eXmin,eXmax,eYmin,eYmax;  //run limits

  Float_t eXstep[3];
  Float_t eYstep[3];
  Float_t eViewXmin[3], eViewXmax[3];
  Float_t eViewYmin[3], eViewYmax[3];
  
  TObjArray *eViewCorr; //! corrections obtained from the views alignment
  Float_t     eScoreCut[2];   // score cut for segments selection (a,b)
 
 public:
  EdbRunAccess();
  EdbRunAccess(EdbRun *run);
  EdbRunAccess(const char *fname);
  virtual ~EdbRunAccess();

  EdbSegment  *GetRawSegmentN( int vid, int sid, int rs=0);
  EdbSegment  *GetRawSegment( int vid, int sid, int rs=0 );
  EdbSegment  *GetRawSegment( EdbView &v, int sid, int rs=0 );
  void ApplyCorrections( const EdbView &view, EdbSegment &s, const int rs );
  void ReadVAfile();
  
  void Set0();
  void SetImageCorrection(int side, const char *str);
  void SetImageCorrectionAff(int side, const char *str);
  void ReadImageMatrixCorrection(int side, const char *file);
  void ReadImageMatrixCorrection(EdbCell2 &map, const char *file);
  void ClearCuts();
  bool InitRun(const char *runfile=0, bool do_update=false);
  bool InitRunFromRWC(char *rwcname, bool bAddRWD=true, const char* options="");
  bool AddRWDToRun(char *rwdname, const char* options="");
  void Print();
  void PrintStat();
  void GuessNviewsPerArea();
  
  void ApplyImageMatrixCorr(int side, float &x, float &y);

  Int_t GetNareas()        { return eNareas;        }
  Int_t GetNviewsPerArea() { return eNviewsPerArea; }
  
  void  CheckRunLine();
  int   Check0Views(EdbPattern &pat, int thres=1);
  int   CheckEmptyViews(EdbPattern &pat);
  float CheckMeanSegsPerView(EdbPattern &pat);
  void  CheckViewStep();
  void  CheckViewStep(int ud);
  TH2F  *CheckUpDownOffsets();
  void  CheckViewSize();
  void SetScoreCut(const char *cut) { float a,b; sscanf(cut,"%f %f",&a,&b); SetScoreCut(a,b); }
  void SetScoreCut(float a, float b) { eScoreCut[0]=a; eScoreCut[1]=b; }
  float GetScoreA() const { return eScoreCut[0]; }
  float GetScoreB() const { return eScoreCut[1]; }

  EdbRun *GetRun() const {return eRun;}
  EdbLayer *GetMakeLayer(int id);
  EdbLayer *GetLayer(int id)
    { if(eLayers[id]) return (EdbLayer *)eLayers[id]; else return 0; }
 
  int FillVP();
  EdbPattern *GetVP(int ud) const { if(ud>-1&&ud<3) return eVP[ud]; else return 0; }

  int FirstArea() const { return eFirstArea;}
  int LastArea()  const { return eLastArea;}

  //return the entries of views in the run tree:

  int GetVolumeArea(EdbPatternsVolume &vol, int area);
  int GetVolumeData(EdbPatternsVolume &vol, int nviews, TArrayI &srt, int &nrej);

  int GetViewsArea(int ud, TArrayI &entr, int area, 
		   float &xmin, float &xmax, float &ymin, float &ymax );
  int GetViewsArea(int ud, TArrayI &entr,  
		   float xmin, float xmax, float ymin, float ymax );
  int GetViewsAreaMarg(int ud, TArrayI &entr, int area, float xmarg, float ymarg);

  int GetViewsXY(int ud, TArrayI &entr, float x, float y, float r=200.);
  int GetEntryXY(int ud, float x, float y);
  int GetVolumeXY(EdbSegP &s, EdbPatternsVolume &vol);

  int GetPatternXYcut(EdbSegP &s, int side, EdbPattern &pat, float dr, float dt);
  int GetPatternXY(EdbSegP &s, int side, EdbPattern &pat, float rmin=200);
  int GetPatternData(EdbPattern &pat, int side, int nviews, TArrayI &srt, int &nrej);
  int GetPatternDataForPrediction( int id, int side, EdbPattern &pat );
  int GetPatternView( EdbPattern &pat, int side, int entry, int &nrej );

  bool  AcceptRawSegment(EdbView *view, int ud, EdbSegP &segP, int side, int entry);

  int   ViewSide(const EdbView *view) const;
  float ViewX(const EdbView *view) const {return view->GetXview();} // todo
  float ViewY(const EdbView *view) const {return view->GetYview();} // todo

  float SegmentWeight(const EdbSegment &s);

  bool  PassCuts(int ud, EdbSegment &seg);
  int  NCuts(int ud) {
    if(!eCuts[ud])  return 0;
    return eCuts[ud]->GetEntriesFast();
  }

  bool  SetSegmentAtExternalSurface( EdbSegment *seg, int side );
  float GetRawSegmentPix( EdbSegment *seg );
  float CalculateSegmentChi2( EdbSegment *seg, float sx, float sy, float sz );

  void           AddSegmentCut(int xi, const char *cutline );
  void           AddSegmentCut(int ud, int xi, float var[10]);
  void           AddSegmentCut(int ud, int xi, float min[5], float max[5]);
  EdbSegmentCut *GetCut(int ud, int i)
    { return (EdbSegmentCut *)(eCuts[ud]->UncheckedAt(i)); }
  float GetCutGR() const {return eCutGR;}

  void SetCond(int ud, EdbScanCond &cond) { 
    if(ud<0) return;
    if(ud>2) return;
    if(eCond[ud]) delete (eCond[ud]);
    (eCond[ud]) = new EdbScanCond(cond);
  }

  EdbScanCond *GetMakeCond(int ud) {return 0;} //TODO?
  EdbScanCond *GetCond(int ud)
    { if(eCond[ud]) return (EdbScanCond *)eCond[ud]; else return 0; }

  float OverlapX(int ud) {return (ud>0&&ud<3)? (eViewXmax[ud]-eViewXmin[ud] - eXstep[ud]): 0.; }
  float OverlapY(int ud) {return (ud>0&&ud<3)? (eViewYmax[ud]-eViewYmin[ud] - eYstep[ud]): 0.; }

  bool CopyRawDataXY( float x0, float y0, float dR, const char *file );

  void SetCutLeft(   int ud, float wmin );
  void SetCutRight(  int ud, float wmin );
  void SetCutTop(    int ud, float wmin );
  void SetCutBottom( int ud, float wmin );
  
  void FillDZMaps();

  ClassDef(EdbRunAccess,2)  // helper class for access to the run data
};

#endif /* ROOT_EdbRunAccess */

