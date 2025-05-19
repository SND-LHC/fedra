//-- Author :  Valeri Tioukov   6.04.2000
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbImage                                                             //
//                                                                      //
// Implementation of Image class                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TBuffer.h>
#include <TH1.h>
#include <TH2.h>
#include <Riostream.h>
#include <TClass.h>
//#include <iostream>
//#include <fstream>
using namespace std;

#include "EdbImage.h"
 
ClassImp(EdbImage)

//______________________________________________________________________________
EdbImage::EdbImage()
{
  Set0();
}

//______________________________________________________________________________
EdbImage::EdbImage(const char *file, const char *format)
{
  Set0();
  eBuffer  = new TArrayC();

  Int_t counter = 0;

  if         ( !strcmp(format,"PGM") )  counter = LoadPGM(file);
  else if    ( !strcmp(format,"RAW") )  counter = LoadRAW(file);
  else if    ( !strcmp(format,"BMP") )  counter = LoadBMP(file);

  printf("counter = %d\n",counter);
  if(counter != eBytes) 
    printf("ERROR:  EdbImage::EdbImage: counter != eBytes \n");
}

//______________________________________________________________________________
EdbImage::~EdbImage()
{
  if(eBuffer) delete eBuffer;
  Set0();
}

//______________________________________________________________________________
void EdbImage::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbImage.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      if (R__v > 1) {
	EdbImage::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	return;
      }
      //====process old versions before automatic schema evolution
      TObject::Streamer(R__b);
      R__b >> eColumns;
      R__b >> eRows;
      R__b >> eColors;
      R__b >> eBuffer;
      R__b >> eBytes;
      R__b.CheckByteCount(R__s, R__c, EdbImage::IsA());
      //====end of old versions
   } else {
     EdbImage::Class()->WriteBuffer(R__b,this);
   }
}
   
//______________________________________________________________________________
void EdbImage::Set0()
{
  eRows    = 0;
  eColumns = 0;
  eBuffer  = 0;
  eBytes   = 0;
  eColors  = 0;
}

//______________________________________________________________________________
EdbImage::EdbImage( int columns, int rows, const char *image, int col )
{
  eColors  = col;
  eColumns = columns;
  eRows    = rows;
  eBytes   = rows*columns;

  int save = 32;                        // 32 bytes safety margin
  int s = eBytes + save;
  char *buffer = new char[s];

  eBuffer = new TArrayC();

  int size = (int)strlen(image);     printf("strlen(image) = %d\n",size);

  memcpy((void *)buffer,(const void *)image,eBytes);
  //strncpy(buffer,image,eBytes);                // problems in case of '\0'

  eBuffer->Adopt(s,buffer);
}

//______________________________________________________________________________
void EdbImage::AdoptImage(int columns, int rows, char *image, int col)
{
  eColors  = col;
  eColumns = columns;
  eRows    = rows;
  eBytes   = rows*columns;
  if(!eBuffer) eBuffer = new TArrayC();
  eBuffer->Adopt(eBytes,image);
}

//______________________________________________________________________________
Int_t EdbImage::FillBufferDouble( TArrayD &image ) const
{
  int npix = Width()*Height();
  if(npix<=0) return 0;
  image.Set(npix);
  for(int i=0; i<npix; i++) image[i]=(Double_t)(Pixel(i));
  return npix;
}

//______________________________________________________________________________
void EdbImage::Print( Option_t *opt ) const
{
  printf("Image: %d x %d \t :  %d bytes\n", eColumns, eRows, eBytes);
}

//______________________________________________________________________________
TH1F *EdbImage::GetHist1() const
{
  TH1F *hist = new TH1F("pix","Pixels distribution", 
			eColors, 0, eColors);
  int i=0;
  unsigned char *buf = (unsigned char*)(eBuffer->GetArray());

  for (int icl=0; icl<eColumns; icl++) {
    for (int ir=0; ir<eRows; ir++) {
      i = ir*eColumns + icl;
      hist->Fill( (int)(buf[i]) );
    }
  }

  return hist;
}

//______________________________________________________________________________
TH2F *EdbImage::GetHist2(int flip) const
{
  TH2F *hist = new TH2F("ccd","CCD image", 
			eColumns, 0, eColumns, 
			eRows, 0, eRows);
  int i=0;
  unsigned char *buf = (unsigned char*)(eBuffer->GetArray());

  int pix;

  for (int icl=0; icl<eColumns; icl++) {
    for (int ir=0; ir<eRows; ir++) {
      i = ir*eColumns + icl;
      pix = (int)(buf[i]);
      if(flip) pix = eColors-1-pix;
      hist->Fill( icl, ir, pix );
    }
  }

  return hist;
}


//____________________________________________________________________________
Int_t EdbImage::LoadRAW( const char *file )
{
/*
  Loads a .raw format file...
  Author: Alfredo Cocco


  eColumns = 512;
  eRows    = 512;

  fstream raw(file, ios::in);
  char byte[262144];
  raw.read(byte,262144);
  raw.close();

  eBytes = 262144;
  int     save = 32;     //safety margins in bytes

  //  int index;
  //  int k = 0;
 
  //  for(int i=0; i<size; i++) {
  //    for(int j=0; j<size; j++) {
  //    index = j + i * size;
  //    buffer[index] = (int) byte[k++];
  //  }
  // }

  eBuffer->Adopt(eBytes+save,byte);
 return eBytes;
 */
  return 0;
}

//____________________________________________________________________________
Int_t EdbImage::LoadPGM( const char *file )
{
/*
  Loads a PGM format file...
  Author: Rosario Peluso
*/
  int      colors, pix, i, j, p2, p5=0;
  char     Id[10], Buffer[BUFSIZ];
  FILE    *in = fopen( file, "rb" );
  int     NRows, NCols; 
  int     save = 32;     //safety margins in bytes

  fgets( Id, 10, in );
  if( (p2 = strncmp(Id, "P2", 2)) && (p5 = strncmp(Id, "P5", 2)) )
    {
      printf("LoadPGM: Error Id\n");      return 0;
    }
 
  while( (fgets(Buffer, BUFSIZ, in) != NULL) && (*Buffer == '#') )
    ;
 
  if( sscanf(Buffer, "%d %d", &NRows, &NCols) != 2 )
    { printf("LoadPGM: Error NRows, &NCols\n");  return 0; }

  eRows    = NRows;
  eColumns = NCols;
  eBytes   = eRows*eColumns;
  
  if( fscanf(in, "%d\n", &colors) != 1 )
    { printf("LoadPGM: Error\n");  return 0; }

  eColors = colors;
 
  if( colors != 255 )
      printf("LoadPGM: Warning: colors = %d \n", colors);
 
  char *Pixel = new char[NRows * NCols + save];
 
  int counter = 0;
  if( !p2 )
    {
      for( j = 0; j < NRows; j++ )
	for( i = 0; i < NCols; i++ )
	  {
	    if( fscanf(in, "%d", &pix) != 1 )   printf("LoadPGM: Error\n");
	    Pixel[j * NCols + i] = pix;
	    counter++;
	  }
    }
  else if( !p5 && ((counter = fread((void *) Pixel, 1, NRows * NCols, in)) != eBytes) )
    {
      printf( "LoadPGM: Error loading binary file\n" );
      return( 0 );
    }
    
  fclose( in );

  eBuffer->Adopt(eBytes+save,Pixel);
 
  printf("EdbImage::LoadPGM: file %s  is readed...\n", file);
 
  return counter;
}

//____________________________________________________________________________________
Int_t EdbImage::DumpPGM( const char *file )
{
  //VT: 3/05/2004

  FILE    *out = fopen( file, "w" );
  if(!out) return 0;
  fprintf(out,"P2\n");                           // grey-scale magic identifier
  fprintf(out,"%d %d\n",Width(),Height());       // image size
  fprintf(out,"255\n");                          // max pixel value
  int nbytes = Width()*Height();
  int ibyte=0;
  unsigned char *buf = (unsigned char*)(eBuffer->GetArray());
  while(1) {
    for(int istr=0; istr<16; istr++)   {
      fprintf(out,"%d ",(int)(buf[ibyte]));
      ibyte++;
      if(ibyte>=nbytes) break;
    }
    fprintf(out,"\n");
    if(ibyte>nbytes-1) break;
  }
  fclose(out);
  return ibyte;
}

//____________________________________________________________________________________
Int_t EdbImage::LoadBMP( const char *file )
{
/*
  Loads a BMP format file...
  Author: Igor Kreslo
*/
typedef struct {
   unsigned short int type;                 /* Magic identifier            */
   unsigned int size;                       /* File size in bytes          */
   unsigned short int reserved1, reserved2;
   unsigned int offset;                     /* Offset to image data, bytes */
} HEADER;

typedef struct {
   unsigned int size;               /* Header size in bytes      */
   int width,height;                /* Width and height of image */
   unsigned short int planes;       /* Number of colour planes   */
   unsigned short int bits;         /* Bits per pixel            */
   unsigned int compression;        /* Compression type          */
   unsigned int imagesize;          /* Image size in bytes       */
   int xresolution,yresolution;     /* Pixels per meter          */
   unsigned int ncolours;           /* Number of colours         */
   unsigned int importantcolours;   /* Important colours         */
} INFOHEADER;

HEADER hdr;
INFOHEADER ihdr;
 
//printf("hdr: %d    ihdr: %d\n",sizeof(hdr), sizeof(ihdr));
 char* ptr;
ifstream raw(file,ios::in | ios::binary);       // ios::binary added! 24.11.2003 (Stas)
 ptr=(char*)&hdr;
 raw.read(ptr,14);
 ptr=(char*)&ihdr;
 raw.read(ptr,40);
 if(ihdr.bits!=8) { printf("Error: %d bits/pixel bitmap is not supported!\n",ihdr.bits); return -1;}
 int N=ihdr.width*ihdr.height;
 char* byte= new char[N];
 raw.read(byte,1024);
 raw.read(byte,N);
 raw.close();
 eBytes=N;
 eColumns=ihdr.width;
 eRows=ihdr.height;
 eColors=256;

 eBuffer->Adopt(eBytes,byte);

 printf("EdbImage::LoadBMP: file %s  is readed...Image %d x %d pixels.\n", file,eColumns,eRows);

 return N;
}

