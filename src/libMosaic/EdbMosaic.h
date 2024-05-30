#ifndef ROOT_EdbMosaic
#define ROOT_EdbMosaic

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbMosaic                                                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TFile.h>
#include <TEnv.h>
#include "EdbID.h"
#include "EdbMosaicIO.h"
#include "EdbMosaicPath.h"
#include "EdbFragmentAlignment.h"
#include "EdbPattern.h"
#include "EdbLayer.h"
#include "EdbViewMap.h"
#include "EdbRunAccess.h"

//______________________________________________________________________________
class EdbMosaicAl : public TObject {
  private:

    //TObjArray eVH;    // list of all views headers

    EdbID eID;
    EdbRunAccess  eRAW;      // run access for raw data
    EdbCell2      eCF[3];    // grope headers by cells for sides 0,1,2
    EdbMosaicIO   eMIO;      // mosaic output
    
    EdbLayer      *eCorrMap[3]; // corrections for sides 0,1,2
    
    int     eMinPeak;         // minimum coincidences to accept alignment
    AlPar   eAP;              // alignment parameters
   
  public:
    EdbMosaicAl(){eCorrMap[0]=eCorrMap[1]=eCorrMap[2]=0;}
    virtual ~EdbMosaicAl(){}
    
    void FormFragments( float fx, float fy, TObjArray &harr );
    void AlignFragments();
    void AlignFragment( EdbPattern &pf, TObjArray &a);
    int  ViewSideAl0(EdbPattern &p1, EdbPattern &p2);
    void AlignSpot(TObjArray &parr, EdbMosaicPath &mp);    
    void ProcRun( EdbID id, const TEnv &cenv );

    void AlignFragments_new();
    void ReadPatterns( EdbFragmentAlignment &fa );
    void SetAlPar( EdbFragmentAlignment &fa );
    void SetAlPar( const TEnv &env, AlPar &ap );
    
    ClassDef(EdbMosaicAl,1)  //Mosaic alignment
};


#endif /* ROOT_EdbMosaic */

