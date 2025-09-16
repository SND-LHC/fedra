#ifndef ROOT_EdbScanSet
#define ROOT_EdbScanSet

#include "TList.h"
#include "TIndexCell.h"
#include "EdbBrick.h"
#include "EdbID.h"

//---------------------------------------------------------------------------
class EdbScanSet : public TNamed
{
 public:
  EdbBrickP  eB;       // all layers of the brick defined here
  TIndexCell ePID;     // correspondance between index in eB and the plate id
  TObjArray  ePC;      // Plate Couples. Each couple represented as EdbPlateP 
                       // where 2 sides in reality corresponds to 2 plates 
  TList      eIDS;     // list of the identifiers to be processed

  Int_t      eReferencePlate;
  EdbID      eID;      // id of the scanset itself

 public:
  EdbScanSet(EdbID id);
  EdbScanSet(int brickid=0);
  virtual ~EdbScanSet(){}

  void     Copy(EdbScanSet &sc);

  void     MakePIDList();
  void     MakeNominalSet(EdbID id, int from_plate=57, int to_plate=1, 
			  float z0=0, float dz=-1300, float shr=1, float dzbase=210., float dzem=45. );
  void     Print();
  Int_t    AssembleBrickFromPC();
  Int_t    TransformSidesIntoBrickRS();
  Int_t    TransformBrick(EdbAffine2D aff);
  void     TransformBrick(EdbScanSet &ss);
  Int_t    ShiftBrickZ(Float_t z);
  Int_t    MakeParFiles(Int_t piece=0, const char *dir=".");
  Int_t    ReadIDS(const char *file);
  Int_t    WriteIDS(const char *file=0);
  Bool_t   AddID(EdbID *id, Int_t step);
  Bool_t   SetAsReferencePlate(Int_t pid);
  Bool_t   GetAffP2P(Int_t p1, Int_t p2, EdbAffine2D &aff);
  Float_t  GetDZP2P(Int_t p1, Int_t p2);
  EdbID   *FindNextPlateID(Int_t p, Bool_t direction);
  EdbID   *FindPlateID(Int_t p);
  EdbID   *GetID(Int_t i) {return (EdbID *)(eIDS.At(i));}
  void    SetID(EdbID id) {eID=id;}

  const EdbBrickP& Brick()   const {return eB;}
  EdbBrickP&       Brick()         {return eB;}

  bool      ValidSide(int side) const { if(side>-1&&side<3) return 1; return 0;}
  EdbLayer *GetLayer(Int_t p, Int_t side)  { if(GetPlate(p)) if(ValidSide(side)) return GetPlate(p)->GetLayer(side); return 0; }
  float     Zlayer(int plate, int side)  { if(GetLayer(plate,side)) return GetLayer(plate,side)->Z() + GetPlate(plate)->Z(); else return 0; }
  
  EdbPlateP *GetPlate(Int_t p)  {
    if (ePID.Find(p)) 
      return eB.GetPlate(ePID.Find(p)->At(0)->Value()); 
    else return 0; 
  }
  
  int  GetPlateIDNearestToZ(float z, int dir=0);
  void UpdateBrickWithP2P(EdbLayer &la, int plate1, int plate2);
  void RemovePlate(int pid);
  void UpdateGap(float dz, int pid1, int pid2);

  void WriteGeom(const char *fname);
  void ReadGeom( const char *fname);
  void UpdateIDS(int brick, int ma, int mi);

  ClassDef(EdbScanSet,2)  // 
};

#endif /* ROOT_EdbScanSet */
