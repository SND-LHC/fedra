#ifndef ROOT_EdbViewMatch
#define ROOT_EdbViewMatch
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbViewMatch - nearby views matching by clusters&grains - the primary purpose is distortions correction  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "TFile.h"
#include "TEnv.h"
#include "TNtuple.h"
#include "EdbCell2.h"
#include "EdbRun.h"


//______________________________________________________________________________
class EdbClMatch : public TObject {
public:
  Float_t  eX,eY,eZ;  // cluster coords in view RS
  Float_t  eXv,eYv;   // view coords
  Int_t    eView;
  Int_t    eFrame;
  Short_t  eIsCenter;   
  
  EdbClMatch(){};
  EdbClMatch(float x,float y,float z, float xv,float yv, int view, int frame, short i, short j)
  { eX=x; eY=y; eZ=z; eXv=xv; eYv=yv; eView=view; eFrame=frame;}
  virtual ~EdbClMatch(){}
  
  ClassDef(EdbClMatch,1)  // service structure for views correction
};

//______________________________________________________________________________
class EdbViewMatch : public TObject {
private:
  int        eNClMin;      // minimal number of clusters inside grain to be used for corrections
  float      eR2CenterMax; // the maximal distance to the matrix center for the reference cluster
  int        eNCenterMin;  // the minimal number of cluster appearance in the central region
  float      eRmax;        // acceptance for clusters matching
  
  TClonesArray  eCl;       // array of EdbClMatch objects
  TClonesArray  eGr;       // array of EdbSegment objects (cluster "clouds" in this case)
  EdbCell2   eGMap;        // map of grains (EdbSegments)
  EdbCell2   eCorrMap;     // map of corrections
  EdbCell2   eCorrMap0;    // input map of corrections
  int        eIC[2];       // index of the central bin used as a position reference
  
  float eXpix,  eYpix;     // pixel size in microns
  int   eNXpix, eNYpix;    // camera matrix size in pixels
  int   eDirX, eDirY;      // microscope steps directions {-1, 0, 1}

  float  eCorrectionMatrixStepX;
  float  eCorrectionMatrixStepY;

public:
  TFile *eOutputFile;
  bool   eDumpGr;                     // if(1) dump grains tree
  float  eAreaMin, eAreaMax;          // clusters Area limits
  float  eVolumeMin, eVolumeMax;      // clusters Volume limits

public:
  EdbViewMatch();
  virtual ~EdbViewMatch(){}

  void InitGMap();
  void InitCorrMap();
  
  void AddCluster( EdbClMatch *c );
  
  void CalculateGrRef();
  void CalculateGrRefMean();
  void CalculateCorr();
  void DrawCorrMap();
  void GenerateCorrectionMatrix(bool do_add);
  void SetPar(TEnv &env);
  void SetPixelSize(float xpix, float ypix) { eXpix=xpix; eYpix=ypix; }
  void CutGrRef();
  void ReadMatrix2Map(const char *file, EdbCell2 &map);
  
  TNtuple *DumpGr(const char *name);

  void MakeDistortionMap(const char *fname, TEnv &env, const char *usefile=0, const char *addfile=0);
  int MakeViewDirectionList( EdbRun &run, int dirx, int diry, TArrayI &entries );
  
  void Print();
  
  ClassDef(EdbViewMatch,1)  // 
};
#endif /* ROOT_EdbViewMatch */
