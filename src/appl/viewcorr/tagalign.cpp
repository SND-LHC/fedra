//-- Author :  Valeri Tioukov   03/06/2024
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbMosaic.h"
#include "EdbMosaicIO.h"

using namespace std;
using namespace TMath;

void AlignNewNopar(EdbID id1, EdbID id2, int side, TEnv &cenv, EdbAffine2D *aff, float dz, EdbScanProc &sproc);
int  AlignID( EdbID id, TEnv &env );
bool Align2Fragments( EdbPattern &p1, EdbPattern &p2, EdbLayer &l1);
EdbPattern *GetPattern( EdbID id, int side );

void print_help_message()
{
  cout<< "\nUsage: \n";
  cout<< "\t  tagalign  -A=ida -B=idb [-v=DEBUG] \n";
  cout<< "\t  tagalign  -id=ID        [-v=DEBUG] \n";

  cout<< "\t\t  ida   - id of the first pattern\n";
  cout<< "\t\t  idb   - id of the second pattern \n";
  cout<< "\t\t  ID    - id of the mos.root file (align up/down in this case) \n";
  
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (tagalign.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into tagalign.save.rootrc file\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default_par(TEnv &cenv)
{
  //set default parameters
  cenv.SetValue("fedra.align.OffsetMax"   , 10000. );
  cenv.SetValue("fedra.align.SigmaR"      , 50.    );
  cenv.SetValue("fedra.align.SigmaT"      , 0.008  );
  cenv.SetValue("fedra.align.DoFine"      , 1      );
  cenv.SetValue("fedra.align.DZ"          , 0      );
  cenv.GetValue("fedra.align.DPHI"        , 0.008  );
  cenv.GetValue("fedra.align.SaveCouples" , 1      );
}

int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("tagalignenv");
  gEDBDEBUGLEVEL     = cenv.GetValue("tagalign.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("tagalign.env"            , "tagalign.rootrc");
  const char *outdir = cenv.GetValue("tagalign.outdir"         , "..");
  
  bool      do_single   = false;
  bool      do_A        = false;
  bool      do_B        = false;
  bool      do_set      = false;
  bool      do_flip     = false;
  EdbID     id;
  EdbID     ida, idb;
  int       al_side     = 1; // 1, 2, 0-base

  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];

    if(!strncmp(key,"-set=",5))
    {
       if(strlen(key)>5) if(id.Set(key+5))   do_set=true;
    }
    else if(!strncmp(key,"-id=",4))
    {
      if(strlen(key)>4) if(id.Set(key+4))   do_single=true;
    }
    else if(!strncmp(key,"-A=",3))
    {
      if(strlen(key)>3) if(ida.Set(key+3))   do_A=true;
    }
    else if(!strncmp(key,"-B=",3))
    {
      if(strlen(key)>3) if(idb.Set(key+3))   do_B=true;
    }
    else if(!strncmp(key,"-side=",6))
    {
      if(strlen(key)>6)	al_side = atoi(key+6);
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
    else if(!strncmp(key,"-flip",5))
    {
      do_flip=true;
    }
  }

  if(!(do_single||do_set)&&(!(do_A&&do_B)))   { print_help_message(); return 0; }
  if(  do_single&&do_set )   { print_help_message(); return 0; }

  set_default_par(cenv);
  cenv.SetValue("tagalign.env"            , env);
  cenv.ReadFile( cenv.GetValue("tagalign.env"   , "tagalign.rootrc") ,kEnvLocal);
  cenv.SetValue("tagalign.outdir"         , outdir);

  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("tagalign.outdir","..");
  cenv.WriteFile("tagalign.save.rootrc");

  printf("\n----------------------------------------------------------------------------\n");
  printf("tagalign  %s\n"      ,id.AsString()	   );
  printf(  "----------------------------------------------------------------------------\n\n");


  if(do_single)
  {
    AlignID(id, cenv);
  }
  else if(do_A&&do_B)
  {
    EdbID id0=ida; id0.ePlate=0;
    EdbScanSet *ss = sproc.ReadScanSet(id0);
    if(ss) {
      EdbAffine2D aff;
      float dz = 0;
      if(ss->GetAffP2P(ida.ePlate, idb.ePlate, aff))	  dz = ss->GetDZP2P(ida.ePlate, idb.ePlate);
      AlignNewNopar(ida,idb, al_side, cenv,&aff, dz, sproc);
    }
  }
  cenv.WriteFile("tagalign.save.rootrc");
  return 1;
}

//-----------------------------------------------------------------------------------
int AlignID( EdbID id, TEnv &cenv )
{
  EdbMosaicIO mio;
  TString file;
  file.Form("p%3.3d/%d.%d.%d.%d.tag.root",
	      id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
  TFile f( file.Data() );
  EdbPattern *p1 = (EdbPattern*)(f.Get( Form("gpat%d_1_0",id.ePlate) ) );
  EdbPattern *p2 = (EdbPattern*)(f.Get( Form("gpat%d_2_0",id.ePlate) ) );
  p1->SetScanID(id); p1->SetSide(1);
  p2->SetScanID(id); p2->SetSide(2);
  
  EdbLayer la;
  Align2Fragments( *p1, *p2, la);
  
}

//-----------------------------------------------------------------------------------
EdbPattern *GetPattern( EdbID id, int side )
{
  EdbMosaicIO mio;
  TString file;
  file.Form("p%3.3d/%d.%d.%d.%d.tag.root",
	      id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
  TFile f( file.Data() );
  EdbPattern *p = (EdbPattern*)(f.Get( Form("gpat%d_%d_0",id.ePlate, side) ) );
  p->SetScanID(id); p->SetSide(side);
  p->SetSegmentsW(1.);
  return p;
}

//-----------------------------------------------------------------------
bool Align2Fragments( EdbPattern &p1, EdbPattern &p2, EdbLayer &l1)
{
  bool success=false;
  int eMinPeak=25;
  bool do_transform = true; 
  
  EdbPlateAlignment av;
  av.eNoScale     = 0;         // calculate shift and rotation
  av.eNoScaleRot  = 1;         // calculate shift only
  av.eOffsetMax   = 10000.;
  av.eDZ          = 0.;
  av.eDPHI        = 0.02;
  av.eDoFine      = 1;
  av.eSaveCouples = 1;
  av.SetSigma( 50, 0.1 );
  av.eDoublets[0] = av.eDoublets[1]=0.01;
  av.eDoublets[2] = av.eDoublets[3]=0.0001;

  av.eSaveCouples=1;
  av.InitOutputFile( Form( "p%.3d/%d_%d.alt.root", p1.Plate(), p1.ID(), p2.ID() ) ); 
  av.Align( p1, p2, 10); //-190
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();

  /*
  if(av.eNcoins > eMinPeak )
  {
    if(do_transform) p1.Transform( affXY );
    l1.GetAffineXY()->Transform( affXY );
    success=true;
  }
  */
  
  av.CloseOutputFile();
  return success;
}

//-------------------------------------------------------------------
void AlignNewNopar(EdbID id1, EdbID id2, int side, TEnv &cenv, EdbAffine2D *aff, float dz, EdbScanProc &sproc)
{
  // Align 2 patterns. All necessary information should be in the envfile
  // Convension about Z(setted while process): the z of id2 is 0, the z of id1 is (-deltaZ) where
  // deltaZ readed from aff.par file in a way that pattern of id1 projected 
  // to deltaZ correspond to pattern of id2

  EdbPlateAlignment av;
  av.eOffsetMax =   cenv.GetValue("fedra.align.OffsetMax"   , 1000. );
  av.SetSigma(      cenv.GetValue("fedra.align.SigmaR"      , 50.  ), 
	            cenv.GetValue("fedra.align.SigmaT"      , 0.008) );
  av.eDoFine      = cenv.GetValue("fedra.align.DoFine"      , 1);
  av.eDZ          = cenv.GetValue("fedra.align.DZ"          , 0);
  av.eDPHI        = cenv.GetValue("fedra.align.DPHI"        , 0.008 );
  av.eSaveCouples = cenv.GetValue("fedra.align.SaveCouples" , 1);

  EdbPattern *p1 = GetPattern(id1,side);
  EdbPattern *p2 = GetPattern(id2,side);
  if(aff) { aff->Print(); p1->Transform(aff);}
 
  TString alfile;  
  sproc.MakeAffName(alfile,id1,id2,"al.root");
  av.InitOutputFile( alfile );
  printf("dz = %f \n",dz);
  av.Align( *p1, *p2 , dz);
  av.CloseOutputFile();
  sproc.UpdateAFFPar( id1, id2, av.eCorrL[0], aff );
}


