#ifndef ROOT_EdbLinking
#define ROOT_EdbLinking

#include "EdbScanCond.h"
#include "TEnv.h"
#include "EdbAlignmentV.h"
#include "EdbSEQ.h"

//-------------------------------------------------------------------------------------------------
class EdbLinking : public EdbAlignmentV
{
 public:
  
   bool   eDoSaveCouples;     // save couples and report histograms
   bool   eDoCorrectAngles;
   bool   eDoCorrectShrinkage;
   bool   eDoFullLinking;
   bool   eDoDumpDoubletsTree;
   int    eCPRankingAlg;          // couples ranking algorithm (0-default, 1-likelihood used)
   
   float  eDRfull, eDTfull;       // acceptance for the full linking
   float  eCHI2Pmax;              // acceptance to save into couples tree
   float  eChi2Acorr;             // acceptance to save into couples tree
   float  eNsigmaEQshr;           // equalization cut for shrinkage 
   float  eNsigmaEQlnk;           // equalization cut for linking
   float  eShr0;                  // starting value for shrinkage correction search
   float  eDShr;                  // range for the shrinkage correction
   float  eBinOK;                 // mean cell occupancy
   EdbScanCond eCond;             // scanning conditions for couples ranking
  
   int    eNshr[2];               // number of coins found for shrinkage correction
   int    eNcorrMin;              // min number of segments for the correction calculation
 
  TObjArray   eSegCouples;       // segment couples objects to fill couples format tree
  
  EdbLayer eL1,eL2;              // layers with the geometry and corrections
  
  EdbH2 eHdxyShr[2];             // service histograms used for the shrinkage correction
  
  struct RemoveDoublets {
    int   remove;
    float dr;
    float dt;
    int   checkview;
  } eRemoveDoublets;
  
 public:
  EdbLinking();
  virtual ~EdbLinking(){}

  bool     VerifyShrinkageCorr( int side );
  void     GetPar(TEnv &env);
  void     GetDoubletsPar(TEnv &env);
  void     GetPreselectionPar(EdbSEQ &seq, TEnv &env);
  void     SaveCouplesTree(const char *file=0);
  void     FullLinking(TObjArray &p1, TObjArray &p2);
  void     FullLinking(EdbPattern &p1, EdbPattern &p2);
  void     FillCombinationsAtMeanZ(TObjArray &p1, TObjArray &p2);
  void     CorrectShrinkage(TObjArray &p1, TObjArray &p2, float dshr);
  void     CorrectShrinkage(EdbPattern &p1, EdbPattern &p2, float dshr);
  void     CorrectShrinkage(float dshr);
  void     Link(EdbPattern &p1, EdbPattern &p2, EdbLayer &l1, EdbLayer &l2, TEnv &env, Double_t area1=-1, Double_t area2=-1);
  void     CorrectAngles(TObjArray &p1, TObjArray &p2);
  void     WriteShrinkagePlots();
  
  Double_t EstimatePatternArea(EdbPattern &p);
  void     RankCouples( TObjArray &arr1,TObjArray &arr2 );
  
  void     DoubletsFilterOut(TObjArray &p1, TObjArray &p2, bool fillhist=0);
  void     DumpDoubletsTree(EdbAlignmentV &adup, const char *name);
  
  void     CloneCouplesTree( const char *ifile, const char *ofile, EdbAffine2D *aff=0, TCut *cut=0 );
  
  void     SetApplyCorr(bool corr) {eCorr[0].eApplyCorr=corr;eCorr[1].eApplyCorr=corr;}
  void     GetFillPattern(EdbPattern &p, int layer=0);
  
  void ProduceReport();

  ClassDef(EdbLinking,1)  // microtracks linking in one plate
};
#endif /* ROOT_EdbLinking */
