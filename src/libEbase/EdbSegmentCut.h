#ifndef ROOT_EdbSegmentCut
#define ROOT_EdbSegmentCut
#include "TObject.h"

//______________________________________________________________________________
class EdbSegmentCut : public TObject {

 private:
  Int_t    eXI;         // 0-exclusive; 1-inclusive cut
  Float_t  eMin[5];     // min  x:y:tx:ty:puls
  Float_t  eMax[5];     // max  x:y:tx:ty:puls

 public:
  EdbSegmentCut() {}
  EdbSegmentCut( int xi, float var[10] );
  virtual ~EdbSegmentCut() {}

  void SetXI(int xi)           {eXI=xi;}
  void SetMin(  float min[5] ) { for(int i=0;i<5;i++) eMin[i]=min[i]; }
  void SetMax(  float max[5] ) { for(int i=0;i<5;i++) eMax[i]=max[i]; }
  int   XI() const {return eXI;}
  float Min(int i) const {return eMin[i];}
  float Max(int i) const {return eMax[i];}
  int  PassCut( float x, float y, float tx, float ty, float puls );
  int  PassCut( float var[5] );
  int  PassCutX( float var[5] );
  int  PassCutI( float var[5] );
  void Print();
  const char *CutLine(char *str, int i=0, int j=0) const;

  ClassDef(EdbSegmentCut,1)  // segment cut
};

#endif /* ROOT_EdbSegmentCut */
