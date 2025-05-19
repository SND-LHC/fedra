#ifndef ROOT_EdbFragmentAlignment
#define ROOT_EdbFragmentAlignment

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbFragmentAlignment                                                 //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TFile.h>
#include <TEnv.h>
#include "EdbID.h"
#include "EdbMosaicIO.h"
#include "EdbPattern.h"
#include "EdbPlateAlignment.h"
#include "EdbLayer.h"
#include "EdbViewMap.h"
#include "EdbRunAccess.h"
#include "EdbCouplesTree.h"

struct AlPar
{
  int        NoScaleRot;
  float      SigmaR;
  float      SigmaT;
  float      OffsetMax;
  float      DZ;
  float      DPHI;
  float      Doublets[4];         //dx,dy,dtx,dty
  int        DoFine;
  int        DoSaveCouples;
};

//______________________________________________________________________________
class EdbFragmentAlignment : public TObject {
  private:

    int eN;                 // number of input headers
    TObjArray   eHarr;      // list of headers
    TObjArray   eParr;      // list of patterns
    EdbPattern *eVC0;       //! original centers
    EdbPattern *eVC;        //! centers after alignment
    
    int       eID;        // fragment id
    int       eSide;      // side
    int    eMinPeak;

    EdbAffine2D   eAff;  // recalculeated scaleX scaleY

  public:
    AlPar     eAP;        // alignment parameters
   
  public:
    EdbFragmentAlignment(){ eN=0; eVC0=0; eVC=0; }
    virtual ~EdbFragmentAlignment(){ eParr.Delete(); SafeDelete(eVC); SafeDelete(eVC0);  }
    
    void SetID(int id) {eID=id;}
    void SetSide(int s) {eSide=s;}
    void SetHarr(TObjArray &ha);
    void SetMinPeak(int mp) {eMinPeak=mp;}
    void SetAlPar( const AlPar &ap, EdbPlateAlignment &al);

    int            N() const {return eN;}
    int            Side() const {return eSide;}
    TObjArray      &GetParr() { return eParr; }
    EdbViewHeader *GetHeader(int i)  const { return (EdbViewHeader *)(eHarr.At(i)); }
    EdbPattern    *GetPattern(int i) const { return    (EdbPattern *)(eParr.At(i)); }

    void ApplyAff();
    int  ViewSideAl( EdbPattern &p1, EdbPattern &p2, EdbAffine2D &aff, bool do_shift );
    void AddPatternAt( EdbPattern *p, int i );
    void SetPatternsOwner() { eParr.SetOwner(); }
    float CheckScaleX(float y0);
    float CheckScaleY(float x0);
    void CheckScale( EdbMosaicPath &mp, EdbAffine2D &aff );
    void AlignFragment( EdbPattern &pf );
    void AlignAndShift( EdbMosaicPath &mp );
    void RealignAndShift( EdbMosaicPath &mp );
    void FillVC();
    void FillVDT(EdbCouplesTree &vdt);

    ClassDef(EdbFragmentAlignment,1)  //Mosaic alignment
};

#endif /* ROOT_EdbMosaic */
