//-- Author :  Valeri Tioukov   15/02/2024
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbAttachPath - attachment path generation algorithm                 //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include "EdbLog.h"
#include "EdbAttachPath.h"

using namespace TMath;
ClassImp(EdbAttachPath);

//-----------------------------------------------------------------------
void EdbAttachPath::SetPoint(const int i, const float x, const float y, const int id)
{
  eX[i]=x;
  eY[i]=y;
  eID[i]=id;
}

//-----------------------------------------------------------------------
void EdbAttachPath::SetStartingPosition(const float x, const float y) 
{
  int i0 = -1;
  float d0=1000000000;
  for(int i=0; i<eN0; i++)
  {
    float dx = eX[i]-x;
    float dy = eY[i]-y;
    float d = Sqrt(dx*dx+dy*dy);
    if(d<d0) {d0=d; i0=i;}
  }  
  eX0=eX[i0]; 
  eY0=eY[i0];
  Log(2,"EdbAttachPath::SetStartingPosition","%f %f",eX0,eY0);
}
 
//-----------------------------------------------------------------------
float EdbAttachPath::GetMin( TArrayF &a ) 
{
  int n = a.GetSize();
  float a0 = a[0];
  for(int i=1; i<n; i++) if(a[i]<a0) a0=a[i];
  return a0;
}
 
//-----------------------------------------------------------------------
float EdbAttachPath::GetMax( TArrayF &a ) 
{
  int n = a.GetSize();
  float a0 = a[0];
  for(int i=1; i<n; i++) if(a[i]>a0) a0=a[i];
  return a0;
}

//-----------------------------------------------------------------------
void EdbAttachPath::SetStartingAtCenter() 
{
  float x = (Xmax()+Xmin())/2.;
  float y = (Ymax()+Ymin())/2.;
  SetStartingPosition(x,y);
}
 
//-----------------------------------------------------------------------
void EdbAttachPath::OrderPointsRadial()
{
  for(int i=0; i<eN0; i++)
  {
    float dx = eX[i]-eX0;
    float dy = eY[i]-eY0;
    eDist[i] = Sqrt(dx*dx+dy*dy);
  }
  TMath::Sort(eN0,eDist.GetArray(),eInd.GetArray(),0);
  eN=eN0;
}

//-----------------------------------------------------------------------
int EdbAttachPath::GetAlignedNeighbours(const int i0, TArrayI &list) const
{
  // for a given view return list of already aligned neighbours
  int n=0;
  for(int i=0; i<eN0; i++)
  {
    if(eOK[i]) {
      float  dx = eX[i] - eX[i0];
      float  dy = eY[i] - eY[i0];
      float  d  = Sqrt(dx*dx+dy*dy);
      if( d<eR && d>0.001 ) 
      { 
	list[n] = i; n++;
	printf("i = %d  d = %f \n",i,d);
      }
    }
  }
  return n;
}


/*

//-----------------------------------------------------------------------
float EdbAttachPath::InitLineX( const TObjArray &harr, const float y0, const float dy0 )
{
  eN0 = harr.GetEntries();
  int cnt=0;
  for(int i=0; i<eN0; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    eHarr.Add(h);
    if( Abs( h->GetYview()-y0 ) < dy0 ) {
      eDist[i] = h->GetXview();
      cnt++;
    } else    eDist[i] = kMaxLong;
  }
  eN  = cnt;
  TMath::Sort(eN0,eDist.GetArray(),eInd.GetArray(),0);
  eY0 = y0;
  eX0 = ((EdbViewHeader *)(eHarr.At(I(0))))->GetXview();
  return ((EdbViewHeader *)(eHarr.At(I(eN-1))))->GetXview() - eX0;
}

//-----------------------------------------------------------------------
float EdbAttachPath::InitLineY( const TObjArray &harr, const float x0, const float dx0 )
{
  eN0 = harr.GetEntries();
  int cnt=0;
  for(int i=0; i<eN0; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    eHarr.Add(h);
    if( Abs( h->GetXview()-x0 ) < dx0 ) {
      eDist[i] = h->GetYview();
      cnt++;
    } else    eDist[i] = kMaxLong;    
  }
  eN  = cnt;
  TMath::Sort(eN0,eDist.GetArray(),eInd.GetArray(),0);
  eX0 = x0;
  eY0 = ((EdbViewHeader *)(eHarr.At(I(0))))->GetYview();
  return ((EdbViewHeader *)(eHarr.At(I(eN-1))))->GetYview() - eY0;
}

//-----------------------------------------------------------------------
EdbViewHeader *EdbAttachPath::FindNearest( const TObjArray &harr, const float x0, const float y0 )
{
  int n = harr.GetEntries();
  float r2min=kMaxLong;
  EdbViewHeader *h0=0;
  for(int i=0; i<n; i++)
  {
    EdbViewHeader *h=(EdbViewHeader *)(harr.At(i));
    float dx = h->GetXview()-x0;
    float dy = h->GetYview()-y0;
    float r2 = dx*dx+dy*dy;
    if(r2<r2min)
    {
      r2min=r2;
      h0=h;
    }
  }
  return h0;
}

*/