//-- Author :  Valeri Tioukov & Daria Morozova   27/06/2025
// flip and rotate raw.root file
#include <iostream>
#include<string>
#include <TEnv.h>
#include <TMath.h>
#include <EdbRun.h>
#include <EdbView.h>
#include "EdbLog.h"
using namespace std;
using namespace TMath;

void FlipRaw(const char *infile, const char *outfile, TEnv &env);
void ViewTurn(EdbView *v, const EdbAffine2D *aff, bool invertZ);

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout<< "\nUsage: \n\t  emflipraw -input=INFILE [-output=OUTFILE -v=verbosyty] \n";
  cout<< "\t  Flip and rotate raw data file INFILE \n";
  cout<< "\t  \t OUTFILE - result of rotation (default flipraw.raw.root)\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default(TEnv &cenv)
{
  cenv.SetValue("flipraw.outdir"          , "..");
  cenv.SetValue("flipraw.AFF"             , "1 0 0 1 0 0");
  cenv.SetValue("flipraw.do_invertZ"      , 0 );
  cenv.SetValue("flipraw.PositionCut"       , "0 0. 0. 0. 0."  ); // do_cut x y dx dy
  cenv.SetValue("flipraw.EdbDebugLevel"   ,  1 );
}

//----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  TEnv cenv("flipraw");
  set_default(cenv);
  cenv.WriteFile("flipraw.save.rootrc");
  if (argc < 2)   { print_help_message();  return 0; }
  gEDBDEBUGLEVEL        = cenv.GetValue("flipraw.EdbDebugLevel" ,  1  );
  const char *outdir    = cenv.GetValue("flipraw.outdir"        , "..");
  
  const char *infile=0;
  const char *outfile=0;
  
  for(int i=1; i<argc; i++ ) 
  {
    char *key  = argv[i];
    
    if(!strncmp(key,"-input=",7))
    {
      if(strlen(key)>7)  { infile = key+7; }
    }
    else if(!strncmp(key,"-output=",8))
    {
      if(strlen(key)>8)  { outfile = key+8; }
    }
     else if(!strncmp(key,"-v=",3))
      {
	if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
      }
   
  }
  
  cenv.ReadFile( "flipraw.rootrc" ,kEnvLocal);
  cenv.WriteFile("flipraw.save.rootrc");
  
  if(infile)
  {
    FlipRaw(infile, outfile, cenv);
  }
  else print_help_message();
}

//--------------------------------------------------------------------------------------------------
void FlipRaw(const char *infile, const char *outfile, TEnv &env)
{
  int   do_cut=0;
  float xcut, ycut, dxcut, dycut;
  if( 5 != sscanf( env.GetValue("flipraw.PositionCut"       , "0 0. 0. 0. 0." ), "%d %f %f %f %f", &do_cut, &xcut, &ycut, &dxcut, &dycut ) ) 
  {
    Log(1, "flipraw", "Error read flipraw.PositionCut!");
    exit;
  }
  
  bool invertZ =   env.GetValue("flipraw.do_invertZ"      ,  0   );
  EdbAffine2D affTotal( env.GetValue("flipraw.AFF"      ,  "1 0 0 1 0 0" ) );
  
  EdbAffine2D affRotation;
  affRotation.Set( affTotal.A11(), affTotal.A12(), affTotal.A21(),affTotal.A22(), 0 , 0);
  
  EdbRun rin( infile );

  EdbRun *rout=0;
  if(!outfile) rout = new EdbRun("flipraw.raw.root", "RECREATE");
  else         rout = new EdbRun(outfile, "RECREATE");

  rin.SetView(rout->GetView());
  
  int n = rin.GetEntries();
  for(int i=0; i<n; i++) {
    EdbView *v = rin.GetEntry(i,1,0,0,0,0);
    bool inside_cut= do_cut? (Abs(v->GetXview()- xcut)< dxcut) && (Abs(v->GetYview() - ycut)<dycut) : true;
    if(inside_cut) {
      v = rin.GetEntry(i,1,1,1,1,1);
      v->Print();
      ViewTurn(v, &affTotal, invertZ);
      rout->AddView();
    }
  }
  
  rout->Close();
  rin.Close();
}

//---------------------------------------------------------
void ViewTurn(EdbView *v, const EdbAffine2D *aff, bool invertZ)
{
  float xv = v->GetXview();
  float yv = v->GetYview();
  float x_new = aff->A11()*xv + aff->A12()*yv + aff->B1();
  float y_new = aff->A21()*xv + aff->A22()*yv + aff->B2();
  EdbAffine2D aa( *(v->GetHeader()->GetAffine()) );
  v->GetHeader()->SetAffine(aa.A11(), aa.A12(), aa.A21(), aa.A22(), x_new, y_new );
  int top = v->GetNframesTop();
  int bot = v->GetNframesBot();
  float z0 = top? v->GetZ2():  v->GetZ3();
  if(invertZ) 
  {
       
       int top = v->GetNframesTop();
       int bot = v->GetNframesBot();
    
       float z1 = v->GetZ1(), z2 = v->GetZ2(), z3 = v->GetZ3(), z4 = v->GetZ4();
       printf("z1 = %f, z2 = %f, z3 = %f,  z4 = %f\n",z1,z2,z3,z4);
       
       z4 = top? -1*v->GetZ1()+1000:0;   //TODO
       z3 = top? -1*v->GetZ2()+1000:0;
       z2 = top? 0:-1*v->GetZ3()+1000;
       z1 = top? 0:-1*v->GetZ4()+1000;
       printf("after flipping: z1 = %f, z2 = %f, z3 = %f,  z4 = %f\n",z1,z2,z3,z4);
       v->SetNframes(bot,top );
       v->SetCoordZ(z1,z2,z3,z4);
       float  z0 = top? z3: z2;
  }

  v->SetCoordXY(x_new,y_new);
  printf("z0 = %f\n",z0);
  
  if(v->GetSegments())
    {
      int nseg=v->Nsegments();
      for( int i=0; i<nseg; i++ )
        {
          EdbSegment *seg = v->GetSegment(i);
          float x = seg->GetX0();
          float y =  seg->GetY0();
          float xr =  aff->A11()*x + aff->A12()*y;
          float yr =  aff->A21()*x + aff->A22()*y;
          float eTx = seg->GetTx();
          float eTy = seg->GetTy();
          
          float tx_r =  aff->A11()*eTx + aff->A12()*eTy;       // rotate angles
          float ty_r =  aff->A21()*eTx + aff->A22()*eTy;
          if(invertZ) 
          {
            int invert_side = (seg->GetSide() == 0) ? 1 : 0;
            seg->Set(xr,yr,z0,tx_r,ty_r, seg->GetDz(), invert_side, seg->GetPuls(), seg->GetID());
          }
          else   seg->Set(xr,yr,seg->GetZ0(),tx_r,ty_r, seg->GetDz(), seg->GetSide(), seg->GetPuls(), seg->GetID());      
        }
    }
}