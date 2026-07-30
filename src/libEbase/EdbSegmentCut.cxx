#include "EdbSegmentCut.h"

ClassImp(EdbSegmentCut)

///==============================================================================
EdbSegmentCut::EdbSegmentCut(int xi, float var[10])
{
  eXI=xi; 
  for(int i=0;i<5;i++) {
    eMin[i]=var[i*2]; 
    eMax[i]=var[i*2+1];
  }
}

///______________________________________________________________________________
int EdbSegmentCut::PassCut(float x, float y, float tx, float ty, float puls)
{
  float var[5];
  var[0]=x; var[1]=y; var[2]=tx; var[3]=ty; var[4]=puls;
  return PassCut(var);
}

///______________________________________________________________________________
int EdbSegmentCut::PassCut(float var[5])
{
  if     (eXI==0)  return PassCutX(var);
  else if(eXI==1)  return PassCutI(var);
  return 0;
}

///______________________________________________________________________________
int EdbSegmentCut::PassCutX(float var[5])
{
  // exclusive cut: if var is inside cut volume - return 0

  for(int i=0; i<5; i++) {
    if(var[i]<eMin[i])  return 1;
    if(var[i]>eMax[i])  return 1;
  }
  return 0;
}

///______________________________________________________________________________
const char *EdbSegmentCut::CutLine(char *str, int i, int j) const
{
  if(eXI==0) {

    sprintf(str,
	  "  TCut x%1d%2.2d = \"!(eX0>%f&&eX0<%f && eY0>%f&&eY0<%f && eTx>%f&&eTx<%f && eTy>%f&&eTy<%f && ePuls>%f&&ePuls<%f)\";\n"
	  ,i,j,eMin[0],eMax[0], eMin[1],eMax[1], eMin[2],eMax[2], eMin[3],eMax[3], eMin[4],eMax[4]);

  } else if(eXI==1) {

    sprintf(str,
	  "  TCut i%1d%2.2d = \"(eX0>%f&&eX0<%f && eY0>%f&&eY0<%f && eTx>%f&&eTx<%f && eTy>%f&&eTy<%f && ePuls>%f&&ePuls<%f)\";\n"
	  ,i,j,eMin[0],eMax[0],eMin[1],eMax[1],eMin[2],eMax[2],eMin[3],eMax[3],eMin[4],eMax[4]);
  }
  return str;
}

///______________________________________________________________________________
int EdbSegmentCut::PassCutI(float var[5])
{
  // inclusive cut: if var is inside cut volume - return 1

  for(int i=0; i<5; i++) {
    if(var[i]<eMin[i])  return 0;
    if(var[i]>eMax[i])  return 0;
  }
  return 1;
}

///______________________________________________________________________________
void EdbSegmentCut::Print()
{
  printf("min: %f %f %f %f %f\n",eMin[0], eMin[1], eMin[2], eMin[3], eMin[4] );
  printf("max: %f %f %f %f %f\n",eMax[0], eMax[1], eMax[2], eMax[3], eMax[4] );
}
