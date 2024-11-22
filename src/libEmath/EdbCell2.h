#ifndef ROOT_EdbCell2
#define ROOT_EdbCell2

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbCell2                                                             //
//                                                                      //
//  class to group 2-dim objects for the fast access                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "EdbCell1.h"

//__________________________________________________________________________
class EdbH2 : public TObject {

 protected:

  Int_t     eN[2];    // divisions
  Float_t eMin[2];    // min
  Float_t eMax[2];    // max
  Float_t eBin[2];    // bin size

  Int_t   eNcell;     // eNx*eNy
  Int_t     *eNC;     //[eNcell] number of objects/cell

 public:
  EdbH2();
  EdbH2(int nx, float minx, float maxx, int ny, float miny, float maxy) { Set0(); InitH2(nx, minx, maxx, ny, miny, maxy); }
  EdbH2( const EdbH2 &h );
  ~EdbH2();

  void Set0();
  void Copy( const EdbH2 &h );

  int   InitH2( const EdbH2 &h );
  int   InitH2(int n[2], float min[2], float max[2]);
  int   InitH2(int nx, float minx, float maxx, int ny, float miny, float maxy);
  void  CleanCells();
  void  PrintStat();
  void  Delete();

  int DiscardHighCells(int nmax);

  int Ncell() const {return eN[0]*eN[1];}
  int NX() const {return eN[0];}
  int NY() const {return eN[1];}
  int IX(float x)             const {return (int)((x-eMin[0])/eBin[0]); }
  int IY(float y)             const {return (int)((y-eMin[1])/eBin[1]); }
  int IX(int jcell)           const {return jcell%eN[0]; }
  int IY(int jcell)           const {return jcell/eN[0]; }
  //int Jcell(int ix, int iy)   const {if(ix>0&&ix<eN[0]&&iy>0&&iy<eN[1]) return iy*eN[0]+ix; else return -1;}
  int Jcell(int ix, int iy)   const {if(ix>=0&&ix<eN[0]&&iy>=0&&iy<eN[1]) return iy*eN[0]+ix; else return -1;}
  int Jcell(float x, float y) const {return Jcell( IX(x), IY(y)); }
  int Jcell(float v[2])       const {return Jcell(IX(v[0]),IY(v[1])); }
  float X(int ix)             const {return eMin[0]+eBin[0]*(ix+0.5);}
  float Y(int iy)             const {return eMin[1]+eBin[1]*(iy+0.5);}
  float Xj(int j)             const {return X(IX(j));}
  float Yj(int j)             const {return Y(IY(j));}
  float Xmin()                const {return eMin[0];}
  float Xmax()                const {return eMax[0];}
  float Ymin()                const {return eMin[1];}
  float Ymax()                const {return eMax[1];}

  float XminA(float level=0);
  float XmaxA(float level=0);
  float YminA(float level=0);
  float YmaxA(float level=0);
  
  EdbH1 *ProjectionX();
  EdbH1 *ProjectionY();

  float Xbin()                const {return eBin[0];}
  float Ybin()                const {return eBin[1];}

  int   Bin(float x, float y) const {return (Jcell(x,y)>-1)? eNC[Jcell(x,y)] : 0; }
  int   Bin(int ix, int iy)   const {return (Jcell(ix,iy)>-1)? eNC[Jcell(ix,iy)] : 0; }
  int   Bin(int iv[2])        const {return Bin(iv[0],iv[1]);}
  int   Bin(int j)            const {if(j>-1||j<eNcell) return eNC[j]; else return 0; }

  int   MaxBin();

  void  AddBin(int jcell, int n) { if(jcell>=0&&jcell<eNcell) eNC[jcell]+=n; }
  int   Fill(float x, float y) { return Fill(x,y,1); }
  int   Fill(float x, float y, int n);
  int   Fill(int j, int n=1);
  void  SetBin(int ix, int iy, int n) { if(Jcell(ix,iy)>-1) eNC[Jcell(ix,iy)] = n; }
  void  SetBin(int j, int n) { if(j>-1||j<eNcell) eNC[j] = n; }

  Long_t  Integral();
  Long_t  Integral(int iv[2], int ir[2]);
  Float_t Mean() { return 1.*Integral()/eNcell; }
  TH1I   *DrawSpectrum( const char *name="plot1d", const char *title="EdbH2 DrawSpectrum");
  TH2F   *DrawH2(       const char *name="plot2d", const char *title="EdbH2plot2D");

  ClassDef(EdbH2,2)  // fast 2-dim histogram class (used as a basis for EdbCell2)
};

//__________________________________________________________________________
class EdbPeak2 : public EdbH2 {

 public:
  Int_t    eNpeaks;  // number of found peaks
  TArrayF  ePeak;    //
  TArrayF  eMean3;   //
  TArrayF  eMean;    //
  Float_t  eNorm;    // the norm-factor to be applied in case of smoothing
  TArrayF  eXpeak;   //
  TArrayF  eYpeak;   //

 public:
  EdbPeak2() { InitPeaks(10); }
  EdbPeak2( const EdbH2 &h ) : EdbH2( h ) { InitPeaks(10); }
  ~EdbPeak2() {}

  void    Delete();
  void    Init(const EdbH2 &h, int npeaks=10);
  void    Print();
  void    InitPeaks( int npeaks);
  int     FindPeak(  int iv[2]);
  int     FindPeak(  float v[2]);
  int     FindPeak(  float &x, float &y);
  float   FindPeak9( float &x, float &y);
  float   FindGlobalPeak(float &x, float &y, float ratio=0.1);
  float   ProbPeak();
  float   ProbPeak(float &x, float &y);
  float   ProbPeak(int iv[2], int ir[2]);
  int     WipePeak(int iv[2], int ir[2]);
  int     ProbPeaks(int npeak);
  int     ProbPeaks(int npeak, int ir[2]);
  float   EstimatePeakVolume(int ipeak);
  float   EstimatePeakVolumeSafe(int ipeak);
  float   EstimatePeakMeanPosition(int iv[2], int ir[2], float &x, float &y, int mean=0);
  float   Xmean();
  float   Ymean();
  float   Peak(int i=0)  const { if(i>=0&&i<eNpeaks) return ePeak[i];  else return -1; }
  float   Mean3(int i=0) const { if(i>=0&&i<eNpeaks) return eMean3[i]; else return -1; }
  float   Mean(int i=0)  const { if(i>=0&&i<eNpeaks) return eMean[i];  else return -1; }
  float   Smooth(Option_t *option="k5a");
  TH1F   *DrawSpectrumN( const char *name="spectrumN", const char *title="EdbH2 DrawSpectrum");
  TH2F   *DrawH2N(       const char *name="xyN", const char *title="EdbH2plot2D");

  ClassDef(EdbPeak2,2)  // peak analyser for EdbH2
};

//__________________________________________________________________________
class EdbCell2 : public EdbH2 {

 private:

  Int_t   eCellLim;   // max number of entries into one cell (for memory allocation)
  Int_t   eCapacity;  //  eNcell*eCellLim
  TObject   **epO;    //[eCapasity] array of the pointers to objects
  TObject  ***epC;    //! [eNcell] pointers to cells

 public:
  EdbCell2();
  EdbCell2( const EdbCell2 &cell ) { Set0(); Copy(cell); }
  ~EdbCell2();

  //virtual void  Streamer(TBuffer &R__b);
  void  Set0();
  void  Copy( const EdbCell2 &cell);
  void  InitEPC();
  int   CellLim() { return eCellLim; }
  int   InitCell( EdbCell2 &c ) { return InitCell(c.NX(),c.Xmin(),c.Xmax(),c.NY(),c.Ymin(),c.Ymax(), c.CellLim() ); }
  int   InitCell(int nx, float minx, float maxx, int ny, float miny, float maxy, int maxpercell );
  int   InitCell(int maxpercell, int n[2], float min[2], float max[2] );
  void  Delete();
  void  Reset() {CleanCells(); Delete();}

  bool  AddObject( float v[2], TObject *obj ) { return AddObject(v[0],v[1],obj); }
  bool  AddObject( float x, float y, TObject *obj );
  bool  AddObject( int ix, int iy, TObject *obj );
  bool  AddObject( int j, TObject *obj );
  int   AddObjectOverlap(float x, float y, float dx, float dy, TObject *obj);
  
  int SelectObjects(TObjArray &arr);
  int SelectObjects(int   min[2], int max[2], TObjArray &arr);
  int SelectObjects(float min[2], float max[2], TObjArray &arr);
  int SelectObjectsCJ(int j, int ir, TObjArray &arr) {
    int iv[2] = { IX(j), IY(j) };
    int irr[2] = { ir, ir };
    return SelectObjectsC(iv, irr, arr);
  }
  int SelectObjectsC(int iv[2], int ir[2], TObjArray &arr);
  int SelectObjectsC(float v[2], int ir[2], TObjArray &arr) {
    int iv[2] = { IX(v[0]), IY(v[1]) }; 
    return SelectObjectsC(iv, ir, arr);
  }
  int SelectObjectsC(float v[2], float r, TObjArray &arr) {
      int ir[2] = { (int)(r/Xbin()) + 1, (int)(r/Ybin()) + 1 };
      return SelectObjectsC(v, ir, arr);
     };
  TObject *GetObject(float x, float y, int ientr) const { return GetObject( Jcell(x,y), ientr); }
  TObject *GetObject(int ix, int iy, int ientr) const { return GetObject( Jcell(ix,iy), ientr); }
  TObject *GetObject(int j, int ientr) const { 
    if(j>=0&&j<eNcell&&ientr>=0&&ientr<eCellLim) return epC[j][ientr];
    else return 0;
  }
  void  PrintStat();

  ClassDef(EdbCell2,2)  // class to group 2-dim objects
};


#endif /* EdbCell2 */
