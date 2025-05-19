#ifndef ROOT_EdbImage
#define ROOT_EdbImage
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbImage                                                             //
//                                                                      //
// CCD Image in bytemap format                     r                    //
//                                     (2x3) :   c 012    012           //
//     convention about rows&columns:              1      345           //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TArrayC.h"
#include "TArrayD.h"

class TH2F;
class TH1F;

//______________________________________________________________________________
class EdbImage : public TObject {
 
private:

  Int_t       eColumns;          // number of columns of the matrix (1 row length)
  Int_t       eRows;             // number of rows of the matrix

  Int_t       eColors;           // maximal height of pixel (normally 256)

  TArrayC    *eBuffer;           // body of the image buffer
  Int_t       eBytes;            // number of significative bytes in the buffer

public:
  EdbImage();
  EdbImage( const char *file, const char *format="PGM" );
  EdbImage( int columns, int rows, const char *image, int col=256 );
  virtual ~EdbImage();


  virtual void   Set0();
  virtual void   Print( Option_t *opt="" ) const;

  void AdoptImage(int columns, int rows, char *image, int col=256);

  void SetColors(int col) {eColors=col;}
  char   *GetBuffer()            const { return eBuffer->GetArray(); }
  Int_t   FillBufferDouble(TArrayD &image)  const;
  int    Pixel(int c, int r) const 
    { if(r<0) return 0; if(r>eRows)    return 0;
      if(c<0) return 0; if(c>eColumns) return 0;
      return Pixel(eColumns*r + c); }
  int    Pixel(int ipix) const 
    { return (int)(((unsigned char*)(eBuffer->GetArray()))[ipix]); }
  TH2F   *GetHist2(int flip=0) const;
  TH1F   *GetHist1() const;

  Int_t   Width()  const {return eColumns;}
  Int_t   Height() const {return eRows;}
  Int_t   LoadPGM( const char *file );
  Int_t   LoadRAW( const char *file );
  Int_t   LoadBMP( const char *file );

  Int_t   DumpPGM( const char *file );

  ClassDef(EdbImage,2)  // CCD Image in bytemap format
};

#endif /* ROOT_EdbImage */

