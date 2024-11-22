//-- Author :  Valeri Tioukov   28-11-2008
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbCell2                                                             //
//                                                                      //
// 2-dim collection of objects                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <algorithm>
#include "TBuffer.h"
#include "TDirectory.h"
#include "TMath.h"
#include "EdbCell2.h"
#include "EdbLog.h"

ClassImp(EdbH2)
ClassImp(EdbCell2)
ClassImp(EdbPeak2)

  using namespace TMath;

//____________________________________________________________________________
EdbH2::EdbH2()
{
  Set0();
}

//____________________________________________________________________________
void EdbH2::Set0()
{
  for(int i=0; i<2; i++) { eN[i]=0; eMin[i]=0; eMax[i]=0; eBin[i]=0; }
  eNcell=0;
  eNC=0;
}

//____________________________________________________________________________
EdbH2::EdbH2( const EdbH2 &h )
{
  Set0();
  Copy(h);
}

//____________________________________________________________________________
void EdbH2::Copy( const EdbH2 &h )
{
  Delete();
  for(int i=0; i<2; i++)     { 
    eN[i]=h.eN[i];
    eMin[i]=h.eMin[i];
    eMax[i]=h.eMax[i];
    eBin[i]=h.eBin[i];
  }
  eNcell=h.eNcell;
  if(h.eNC) {
    eNC = new Int_t[eNcell];
    memcpy( eNC, h.eNC, sizeof(Int_t)*eNcell );
  }
}

//____________________________________________________________________________
EdbH2::~EdbH2()
{ 
  Delete();
}

//____________________________________________________________________________________
void EdbH2::Delete()
{
  if(eNC)    { delete [] eNC;    eNC = 0;  }
  eNcell=0;
  eN[0] = eN[1] = 0;
  eMin[0] = eMin[1] = 0;
  eMax[0] = eMax[1] = 0;
}

//____________________________________________________________________________
int EdbH2::InitH2(const EdbH2 &h)
{
  int n[2]={h.NX(),h.NY()};
  float min[2]={h.Xmin(),h.Ymin()};
  float max[2]={h.Xmax(),h.Ymax()};
  return InitH2(n,min,max);
}

//____________________________________________________________________________
int EdbH2::InitH2(int nx, float minx, float maxx, int ny, float miny, float maxy)
{
  int n[2]={nx,ny};
  float min[2]={minx,miny};
  float max[2]={maxx,maxy};
  return InitH2(n,min,max);
}

//____________________________________________________________________________
int EdbH2::InitH2(int n[2], float min[2], float max[2])
{ 
  Delete();
  for(int i=0; i<2; i++) {
    if(n[i]<1) return 0;
    eN[i]=n[i];
    eMin[i] = min[i];
    eMax[i] = max[i];
    eBin[i] = (eMax[i]-eMin[i])/eN[i];
    //    if(eBin[i]<0.00000001) return 0;
  }
  eNcell = eN[0]*eN[1];
  eNC = new Int_t[eNcell];
  CleanCells();
  return eNcell;
}

//____________________________________________________________________________________
void EdbH2::CleanCells()
{
  if(eNC) memset(eNC,'\0',eNcell*sizeof(Int_t));
}

//____________________________________________________________________________
int EdbH2::Fill(float x, float y, int n)
{
  if(x>eMin[0]&&x<eMax[0]&&y>eMin[1]&&y<eMax[1]) {
    int j = Jcell(x,y);
    return Fill(j,n);
  }
  return 0;
}

//____________________________________________________________________________
int EdbH2::Fill(int j, int n)
{ 
  if( j>-1 && j<eNcell)
  {
    eNC[j] += n;
    return n;
  }
  return 0;
}

//____________________________________________________________________________________
int EdbH2::DiscardHighCells(int nmax)
{
  int ic=0;
  for(int i=0; i<eNcell; i++)  
    if(eNC[i]>nmax) { eNC[i]=0; ic++; }
  return ic;
}

//____________________________________________________________________________________
EdbH1 *EdbH2::ProjectionX()
{
  EdbH1 *h = new EdbH1(eN[0],eMin[0],eMax[0]);
  for(int i=0; i<eNcell; i++)   h->Fill( X(IX(i)), eNC[i] );
  return h;
}

//____________________________________________________________________________________
EdbH1 *EdbH2::ProjectionY()
{
  EdbH1 *h = new EdbH1(eN[1],eMin[1],eMax[1]);
  for(int i=0; i<eNcell; i++)    h->Fill( Y(IY(i)), eNC[i] );
  return h;
}

//____________________________________________________________________________________
float EdbH2::XminA(float level)
{
  EdbH1 *h = ProjectionX();  return h->XminA(level); SafeDelete(h);
}
//____________________________________________________________________________________
float EdbH2::XmaxA(float level)
{
  EdbH1 *h = ProjectionX();  return h->XmaxA(level); SafeDelete(h);
}
//____________________________________________________________________________________
float EdbH2::YminA(float level)
{
  EdbH1 *h = ProjectionY();  return h->XminA(level); SafeDelete(h);
}
//____________________________________________________________________________________
float EdbH2::YmaxA(float level)
{
  EdbH1 *h = ProjectionY();  return h->XmaxA(level); SafeDelete(h);
}

//____________________________________________________________________________________
TH1I *EdbH2::DrawSpectrum(const char *name, const char *title)
{
  int imax = MaxBin();
  TH1I *h = new TH1I(name, title, imax,0,imax);
  for(int i=0; i<eNcell; i++)    h->Fill( eNC[i] );
  return h;
}

//____________________________________________________________________________________
TH2F *EdbH2::DrawH2(const char *name, const char *title)
{
  TObject *obj=0;
  if((obj=gDirectory->FindObject(name))) delete obj;
  TH2F *h = new TH2F(name, title, eN[0],eMin[0],eMax[0],eN[1],eMin[1],eMax[1]);
  for(int i=0; i<eN[0]; i++)
    for(int j=0; j<eN[1]; j++)
      h->Fill( X(i), Y(j), Bin(i,j) );
  return h;
}

//____________________________________________________________________________________
void EdbH2::PrintStat()
{
  printf("X:%5d cells in range (%10.2f %10.2f) with bin %10.2f \n", eN[0],eMin[0],eMax[0],eBin[0]);
  printf("Y:%5d cells in range (%10.2f %10.2f) with bin %10.2f \n", eN[1],eMin[1],eMax[1],eBin[1]);
  if(!eNC)  return;
  int nmax=0, ntot=0;
  int nmin=kMaxInt;
  for(int i=0; i<eNcell; i++) {
    if( eNC[i]>nmax ) nmax=eNC[i];
    if( eNC[i]<nmin ) nmin=eNC[i];
    ntot += eNC[i];
  }
  printf("%d objects in %d cells filled as: (%d:%d) with mean= %f\n", 
	 ntot, eNcell, nmin, nmax, 1.*ntot/eNcell);
}

//____________________________________________________________________________________
Long_t EdbH2::Integral()
{
  Long_t ntot=0;
  for(int i=0; i<eNcell; i++) ntot+=eNC[i];
  return ntot;
}

//____________________________________________________________________________________
Long_t EdbH2::Integral(int iv[2], int ir[2])
{
  Long_t ntot=0, nbin=0;
  for(int ix=iv[0]-ir[0]; ix<=iv[0]+ir[0]; ix++) 
    for(int iy=iv[1]-ir[1]; iy<=iv[1]+ir[1]; iy++) 
      {ntot+=Bin(ix,iy); nbin++;}
  return ntot;
}

//____________________________________________________________________________
int EdbH2::MaxBin()
{
  int peak=0;
  for(int i=0; i<eNcell; i++)  if(eNC[i]>peak)  peak =eNC[i];
  return peak;
}

//______________________________________________________________________________
float EdbPeak2::Smooth(Option_t *option)
{
   // Smooth bin contents of this 2-d histogram using kernel algorithms
   // similar to the ones used in the raster graphics community.
   // Bin contents in the active range are replaced by their smooth values.
   //
   // This code is copied from root_5.22 mainly to not depend on the newest root versions
   // Adopted to fedra by VT
   // to speedup processing do all calculations as integer, so the histogram do not normalised at output
   // return value - norm-factor of the kernel

   if(NX()<1||NY()<1)  return 1.;

   Int_t k5a[5][5] =  { { 0, 0, 1, 0, 0 }, 
			{ 0, 2, 2, 2, 0 }, 
			{ 1, 2, 5, 2, 1 }, 
			{ 0, 2, 2, 2, 0 }, 
			{ 0, 0, 1, 0, 0 } };   // norm = 25
   Int_t k5b[5][5] =  { { 0, 1, 2, 1, 0 },
			{ 1, 2, 4, 2, 1 },
			{ 2, 4, 8, 4, 2 },
			{ 1, 2, 4, 2, 1 },
			{ 0, 1, 2, 1, 0 } };
   Int_t k3a[3][3] =  { { 0, 1, 0 },
			{ 1, 2, 1 },
			{ 0, 1, 0 } };      // norm = 6
   Int_t k3b[3][3] =  { { 1, 2, 1 },        //gauss smoothing
			{ 2, 4, 2 },
			{ 1, 2, 1 } };      // norm = 16

   Int_t ksize_x, ksize_y;
   TString opt = option;
   opt.ToLower();
   Int_t *kernel = &k5a[0][0];
   if (opt.Contains("k5a")) {
     kernel = &k5a[0][0];
     ksize_x = ksize_y = 5;
   }
   else if (opt.Contains("k5b")) {
     kernel = &k5b[0][0];
     ksize_x = ksize_y = 5;
   }
   else if (opt.Contains("k3a")) {
      kernel = &k3a[0][0];
      ksize_x = ksize_y = 3;
   } 
   else if (opt.Contains("k3b")) {
      kernel = &k3b[0][0];
      ksize_x = ksize_y = 3;
   } 
   else return 1.;

   // Kernel tail sizes (kernel sizes must be odd for this to work!) 
   Int_t x_push = (ksize_x-1)/2;
   Int_t y_push = (ksize_y-1)/2; 

   Int_t norm=0;
   for(int i=0; i<ksize_x; i++)
     for(int j=0; j<ksize_y; j++) norm += kernel[i*ksize_y +j];

   int nx = NX(), ny= NY();

   EdbH2 buf( *((EdbH2*)this) );

   Long_t content;
   int bin=0, k=0;

   for (int i=0; i<nx; i++){             // input hist
     for (int j=0; j<=ny; j++) { 
       
       content = 0;
       for (int n=0; n<ksize_x; n++) {          // kernel
	 for (int m=0; m<ksize_y; m++) { 
	   Int_t xb = i+(n-x_push);
	   Int_t yb = j+(m-y_push); 
	   bin = buf.Bin(xb,yb);                      // out of dimensions - return 0;
	   if (!bin) continue;
	   k = kernel[n*ksize_y +m];
	   if (!k) continue;
	   content += k*bin;
	 }
       }
       //SetBin( i, j, Nint(1.*content/norm) );
       SetBin( i, j, content );
     }
   }
   eNorm = (float)norm; 
   return eNorm;
}

//____________________________________________________________________________________
TH1F *EdbPeak2::DrawSpectrumN(const char *name, const char *title)
{
  int   nbin = MaxBin()+1;
  float max = 1.*MaxBin()/eNorm;
  TH1F *h = new TH1F(name, title, nbin,0,max);
  for(int i=0; i<eNcell; i++)    h->Fill( eNC[i]/eNorm );
  return h;
}

//____________________________________________________________________________________
TH2F *EdbPeak2::DrawH2N(const char *name, const char *title)
{
  TObject *obj=0;
  if((obj=gDirectory->FindObject(name))) delete obj;
  TH2F *h = new TH2F(name, title, eN[0],eMin[0],eMax[0],eN[1],eMin[1],eMax[1]);
  for(int i=0; i<eN[0]; i++)
    for(int j=0; j<eN[1]; j++)
      h->Fill( X(i), Y(j), Bin(i,j)/eNorm );
  return h;
}

//____________________________________________________________________________
void EdbPeak2::Delete()
{
  ((EdbH2*)this)->Delete();
  eNpeaks=0;
}

//____________________________________________________________________________
void EdbPeak2::Init( const EdbH2 &h, int npeaks )
{
  Delete();
  ((EdbH2*)this)->Copy(h);
  InitPeaks(npeaks);
}

//____________________________________________________________________________
void EdbPeak2::InitPeaks( int npeaks )
{
  eNpeaks=0;
  ePeak.Set(npeaks);   ePeak.Reset(0);
  eMean3.Set(npeaks);  eMean3.Reset(0);
  eMean.Set(npeaks);   eMean.Reset(0);
  eXpeak.Set(npeaks);  eXpeak.Reset(0);
  eYpeak.Set(npeaks);  eYpeak.Reset(0);
  eNorm = 1.;  // default norm factor;
}

//____________________________________________________________________________
void EdbPeak2::Print()
{
  printf("%d peaks:\n", eNpeaks);
  for(int i=0; i<eNpeaks; i++) 
	  printf("%d : %.1f %.1f %.1f %.1f %.1f \n", i, eXpeak[i], eYpeak[i], ePeak[i], eMean3[i], eMean[i]);
}

//____________________________________________________________________________
float EdbPeak2::ProbPeak( float &x, float &y )
{
  int   iv[2];
  FindPeak(iv);
  x = X(iv[0]);
  y = Y(iv[1]);
  int   ir[2]={1,1};              // 3x3 neigbouring
  return ProbPeak(iv,ir);
}

//____________________________________________________________________________
float EdbPeak2::ProbPeak()
{
  int   iv[2];
  FindPeak(iv);
  int   ir[2]={1,1};              // 3x3 neigbouring
  return ProbPeak(iv,ir);
}

//____________________________________________________________________________
int EdbPeak2::ProbPeaks( int npeaks )
{
  int ir[2]={1,1};
  return ProbPeaks(npeaks,ir);
}

//____________________________________________________________________________
int EdbPeak2::ProbPeaks( int npeaks, int ir[2])
{
  int   iv[2], ic=0;
  for(int i=0; i<npeaks; i++) {
    FindPeak(iv);
    ProbPeak(iv,ir);
    WipePeak(iv,ir);
    ic++;
  }
  return ic;
}

//____________________________________________________________________________
float EdbPeak2::ProbPeak(int iv[2], int ir[2])
{
  float prob=0;
  int   npeak      = Bin(iv);
  int   nbinPeak   = (2*ir[0]+1)*(2*ir[1]+1);
  float meanNeib   = 1.*(Integral(iv,ir) - npeak)/(nbinPeak-1);
  float meanNoPeak = 1.*(Integral() - Integral(iv,ir))/(Ncell()-nbinPeak);
  float x,y;
  float volume = EstimatePeakMeanPosition(iv, ir, x, y,meanNoPeak);  
  Log(3,"ProbPeak","found at (%3d %3d): %4d  %6.3f  %6.3f", iv[0],iv[1],npeak, meanNeib,meanNoPeak);
  eXpeak[eNpeaks] = x;
  eYpeak[eNpeaks] = y;
  ePeak[eNpeaks]  = volume/eNorm;
  eMean3[eNpeaks] = meanNeib/eNorm;
  eMean[eNpeaks]  = meanNoPeak/eNorm;
  prob = (npeak - meanNoPeak)/eNorm;
  eNpeaks++;
  return prob;
}

//____________________________________________________________________________
int EdbPeak2::WipePeak(int iv[2], int ir[2])
{
  int mean = (int)Mean();
  int nbin=0;
  for(int ix=iv[0]-ir[0]; ix<=iv[0]+ir[0]; ix++) 
    for(int iy=iv[1]-ir[1]; iy<=iv[1]+ir[1]; iy++) 
      if(Jcell(ix,iy)>-1)
	{eNC[Jcell(ix,iy)] = mean; nbin++;}
  return nbin;
}

//____________________________________________________________________________
float EdbPeak2::FindPeak9( float &x, float &y)
{
  // particular case of 3x3 peak neiboring, return peak volume and estimate mean position
  int iv[2], ir[2]={1,1};
  FindPeak(iv);
  return EstimatePeakMeanPosition(iv, ir, x, y);
}

//____________________________________________________________________________
float EdbPeak2::EstimatePeakMeanPosition(int iv[2], int ir[2], float &x, float &y, int mean )
{
  double x0=0, y0=0, volume=0;
  for(int ix=iv[0]-ir[0]; ix<=iv[0]+ir[0]; ix++) 
    for(int iy=iv[1]-ir[1]; iy<=iv[1]+ir[1]; iy++) 
      {
	int bin = Bin(ix,iy)-mean;
	x0     += bin*X(ix);
	y0     += bin*Y(iy);
	volume += bin;
      }
  if(volume>0) { x =  x0/volume;  y = y0/volume; }
  return volume;
}

//____________________________________________________________________________
int EdbPeak2::FindPeak(int iv[2])
{
  int peak=0;
  for(int ix=0; ix<eN[0]; ix++)
    for(int iy=0; iy<eN[1]; iy++) {
      int j = Jcell(ix,iy);
      if(eNC[j]>peak) {
        peak =eNC[j];
        iv[0] = ix;
        iv[1] = iy;
      }
    }
  Log(4,"EdbPeak2::FindPeak","of %d at %d %d", peak, iv[0],iv[1] );
  return peak;
}

//____________________________________________________________________________
int EdbPeak2::FindPeak(float v[2])
{
  int iv[2];
  int peak=FindPeak(iv);
  v[0] = X(iv[0]);
  v[1] = Y(iv[1]);
  return peak;
}

//____________________________________________________________________________
int EdbPeak2::FindPeak(float &x, float &y)
{
  int iv[2];
  int peak=FindPeak(iv);
  x = X(iv[0]);
  y = Y(iv[1]);
  return peak;
}

//____________________________________________________________________________
float EdbPeak2::EstimatePeakVolumeSafe(int ipeak)
{
  // add to the background one standard deviation for the safety
  if(ipeak>(eNpeaks-1)||ipeak<0)  return 0;
  return ePeak[ipeak] + 8*eMean3[ipeak] - 9*(eMean[ipeak] + Sqrt(eMean[ipeak]));
}

//____________________________________________________________________________
float EdbPeak2::EstimatePeakVolume(int ipeak)
{
  if(ipeak>(eNpeaks-1)||ipeak<0)  return 0;
  return ePeak[ipeak] + 8*eMean3[ipeak] - 9*eMean[ipeak];
}

//____________________________________________________________________________
float EdbPeak2::Xmean()
{
  Double_t mean=0, sum=0;
  for(int ix=0; ix<eN[0]; ix++) {
    float x = X(ix);
    for(int iy=0; iy<eN[1]; iy++) {
      int j = Jcell(ix,iy);
      mean += x*eNC[j];
      sum  += eNC[j];
    }
  }
  mean /= sum;
  return mean;
}

//____________________________________________________________________________
float EdbPeak2::Ymean()
{
  Double_t mean=0, sum=0;
  for(int iy=0; iy<eN[1]; iy++) {
    float y = Y(iy);
    for(int ix=0; ix<eN[0]; ix++) {
      int j = Jcell(ix,iy);
      mean += y*eNC[j];
      sum  += eNC[j];
    }
  }
  mean /= sum;
  return mean;
}

//____________________________________________________________________________
float EdbPeak2::FindGlobalPeak(float &x, float &y, float ratio)
{
  // to find the center of wide multibin peak select ratio*Nbin of highest bins and calculate CoG
  int n = Ncell();
  TArrayI ind(n);
  Sort(n, eNC , ind.GetArray(), 0);  // sort encreasing order
  Int_t nwipe = Min( (Int_t)((1.-ratio)*n), n-1 );
  double mean=0;
  for(int i=0; i<nwipe; i++) {
    mean += eNC[ind[i]];
    eNC[ind[i]] = 0;
  }
  mean /= nwipe;
  x = Xmean();
  y = Ymean();
  float peakvolume = Integral() - (n-nwipe)*mean;
  Log(3,"EdbPeak2::FindGlobalPeak","of volume %8.2f at %10.2f %10.2f using %d bins", 
      peakvolume,x,y, n-nwipe );
  return peakvolume;
}

//____________________________________________________________________________
EdbCell2::EdbCell2()
{
  Set0();
}

//____________________________________________________________________________
void EdbCell2::Set0()
{
  eCellLim=0;
  eCapacity=0;
  epO=0;
  epC=0;
}

//____________________________________________________________________________
EdbCell2::~EdbCell2()
{ 
  Delete();
}

//____________________________________________________________________________________
void EdbCell2::Delete()
{
  if(epC)    { delete [] epC;    epC = 0;  }
  if(epO)    { delete [] epO;    epO = 0;  }
}

//____________________________________________________________________________
int EdbCell2::InitCell(int nx, float minx, float maxx, int ny, float miny, float maxy, int maxpercell )
{
  int n[2] = {nx,ny};
  float min[2] = {minx,miny};
  float max[2] = {maxx,maxy};
  return InitCell(maxpercell, n, min, max);
}

//____________________________________________________________________________
int EdbCell2::InitCell(int maxpercell, int n[2], float min[2], float max[2])
{
  Log(3,"EdbCell2::InitCell","%d %f %f %d %f %f", n[0],min[0],max[0],n[1],min[1],max[1] );
  Delete();
  InitH2(n, min, max);
  Log(3,"EdbCell2::InitCell","reserve space for maxpercell(%d) * eNcell(%d) = %d objects", 
      maxpercell, eNcell, maxpercell*eNcell);
  eCellLim = maxpercell;
  eCapacity = eNcell*eCellLim;
  epO     = new TObject*[eCapacity];          //- pointers to objects
  epC     = new TObject**[eNcell];            //- pointers to cells
  InitEPC();
  return eNcell;
}

//____________________________________________________________________________
void EdbCell2::InitEPC()
{
  TObject **po = epO;
  epC[0]=po;
  for(int i=1; i<eNcell; i++) {po+=eCellLim; epC[i]=po;}
}

//____________________________________________________________________________
void EdbCell2::Copy( const EdbCell2 &cell )
{
  EdbH2::Copy(cell);
  eCellLim  = cell.eCellLim;
  eCapacity = cell.eCapacity;
  for(int i=0; i<eNcell; i++) 
    for(int j=0; i<eCellLim; i++)  AddObject(i, GetObject(i,j));
}

//______________________________________________________________________________
void EdbCell2::Streamer(TBuffer &R__b)
{
   // Stream an object of class EdbCell2.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      EdbH2::Streamer(R__b);
      R__b >> eCellLim;
      R__b >> eCapacity;
      //Log(1,"EdbCell2::Streamer","Read:  eCellLim =%d eCapacity=%d eNcell=%d ", eCellLim, eCapacity, eNcell);
      
      //R__b.ReadArray(epO);
      epO     = new TObject*[eCapacity]; 
      TObject *obj;
      for (Int_t i = 0; i < eCapacity; i++) {
         obj = (TObject*)R__b.ReadObjectAny(TObject::Class());
         if (obj)  epO[i] = obj;
         else epO[i]=0;
      }
      epC     = new TObject**[eNcell];            //- pointers to cells
      InitEPC();                       //R__b.ReadArray(epC);
      R__b.CheckByteCount(R__s, R__c, EdbCell2::IsA());
   } else {
      R__c = R__b.WriteVersion(EdbCell2::IsA(), kTRUE);
      EdbH2::Streamer(R__b);
      R__b << eCellLim;
      R__b << eCapacity;
      for (Int_t i = 0; i < eCapacity; i++) {
         R__b << epO[i];
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//____________________________________________________________________________________
void  EdbCell2::PrintStat() 
{
  printf("EdbCell2 limit = %4d\n",eCellLim);
  ((EdbH2*)this)->PrintStat();
}

//____________________________________________________________________________
bool EdbCell2::AddObject(int ix, int iy, TObject *obj)
{ 
  int j = Jcell(ix,iy);
  return AddObject(j, obj);
}

//____________________________________________________________________________
bool EdbCell2::AddObject(float x, float y, TObject *obj)
{ 
  int j = Jcell(x,y);
  return AddObject(j, obj);
}

//____________________________________________________________________________
int EdbCell2::AddObjectOverlap(float x, float y, float dx, float dy, TObject *obj)
{
  // works correctly if overlap less then the cell size
  int cnt=0;
  std::vector<int> J;
  J.push_back( Jcell(x,y) );
  J.push_back( Jcell(x+dx,y) );
  J.push_back( Jcell(x-dx,y) );
  J.push_back( Jcell(x,y+dy) );
  J.push_back( Jcell(x,y-dy) );
  std::sort(J.begin(), J.end());
  std::vector<int>::iterator ip = std::unique( J.begin(), J.end() );
  J.resize( std::distance(J.begin(), ip) );
  for (ip=J.begin(); ip!=J.end(); ++ip)
  {
    AddObject(*ip,obj); 
    cnt++;
  }
  return cnt;
}

//____________________________________________________________________________
bool EdbCell2::AddObject(int j, TObject *obj)
{ 
  if(!obj)                    return 0;
  if(j>=eNcell)               return 0;
  if(j<0)                     return 0;
  if( eNC[j] >= eCellLim )    return 0;
  epC[j][eNC[j]] = obj;
  eNC[j]++;
  return 1;
}

//____________________________________________________________________________
int EdbCell2::SelectObjectsC(int vcent[2], int vdiff[2], TObjArray &arr)
{
  int min[2] = {vcent[0]-vdiff[0],vcent[1]-vdiff[1]};
  int max[2] = {vcent[0]+vdiff[0],vcent[1]+vdiff[1]};
  return SelectObjects(min,max, arr);
}

//____________________________________________________________________________
int EdbCell2::SelectObjects(float min[2], float max[2], TObjArray &arr)
{
  int mi[2] = { IX(min[0]), IY(min[1]) };
  int ma[2] = { IX(max[0]), IY(max[1]) };
  return SelectObjects(mi,ma, arr);
}

//____________________________________________________________________________
int EdbCell2::SelectObjects(int min[2], int max[2], TObjArray &arr)
{
  int nobj  = 0;
  int jcell = 0;
  for(int ix = Max(0,min[0]); ix<=Min(eN[0]-1,max[0]); ix++)
    for(int iy = Max(0,min[1]); iy<=Min(eN[1]-1,max[1]); iy++) 
      {
	jcell = Jcell(ix,iy);
	if(eNC[jcell]<1)      continue;
	for(int k=0; k<eNC[jcell]; k++)   { arr.Add( epC[jcell][k] );  nobj++; }
      }
  return nobj;
}

//____________________________________________________________________________
int EdbCell2::SelectObjects(TObjArray &arr)
{
  // select all objects
  int nobj  = 0;
  for(int jcell=0; jcell<eNcell; jcell++) {
    if(eNC[jcell]<1)      continue;
    for(int k=0; k<eNC[jcell]; k++)   { arr.Add( epC[jcell][k] );  nobj++; }
  }
  return nobj;
}
