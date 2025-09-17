//-- Author :  Valeri Tioukov   09/02/2023

#include <string.h>
#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbRunAccess.h"
#include "EdbLinking.h"
#include "EdbScanProc.h"
#include "EdbPlateAlignment.h"
#include "EdbMosaic.h"
#include "EdbMosaicIO.h"

using namespace std;
using namespace TMath;

void MosaicAlID(EdbID id, const TEnv &cenv);
void MosaicAlSet(EdbScanSet &ss, const TEnv &cenv);

void print_help_message()
{
  cout<< "\nUsage: \n";
  cout<< "\t  viewsideal  -id=ID  [-v=DEBUG] \n";
  
  cout<< "\t\t  ID      - id of the raw.root file formed as BRICK.PLATE.MAJOR.MINOR \n";
  
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (viewsideal.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into viewsideal.save.rootrc file\n";
  cout<<endl;
}

void set_default(TEnv &cenv)
{
  // default parameters for views side alignment (vsa)
  
  cenv.SetValue("fedra.vsa.DoImageCorr"      , 0 );
  cenv.SetValue("fedra.vsa.ImageCorrSide1"   , "1 0 0 1 0 0");
  cenv.SetValue("fedra.vsa.ImageCorrSide2"   , "1 0 0 1 0 0");
  cenv.SetValue("fedra.vsa.DoImageMatrixCorr"      , 0 );
  cenv.SetValue("fedra.vsa.ImageMatrixCorrSide1"   , "");
  cenv.SetValue("fedra.vsa.ImageMatrixCorrSide2"   , "");
  cenv.SetValue("fedra.vsa.hxy"   , "2000 0. 200000. 2000 0. 200000.");
  cenv.SetValue("fedra.vsa.offsetMax"      , 10. );
  cenv.SetValue("fedra.vsa.DZ"             ,  0. );
  cenv.SetValue("fedra.vsa.DPHI"           ,  0. );
  cenv.SetValue("fedra.vsa.Yfrag"          ,  5000  );
  cenv.SetValue("fedra.vsa.AFID"           ,  1  );
  cenv.SetValue("fedra.vsa.Xfrag"          ,  10000  );
  cenv.SetValue("fedra.vsa.Yfrag"          ,  10000  );
  cenv.SetValue("fedra.vsa.MinPeak"        ,  25   );
  cenv.SetValue("fedra.vsa.R0"             ,  1200.);
  cenv.SetValue("fedra.vsa.HeaderCut"      , "1" );
  cenv.SetValue("fedra.vsa.sigmaR"         , 0.5 );
  cenv.SetValue("fedra.vsa.sigmaT"         , 0.005 );
  cenv.SetValue("fedra.vsa.ICUT"      , "-1     -500. 500.   -500.   500.    -1.   1.      -1.   1.       8.  50.");
  cenv.SetValue("fedra.vsa.XCUT"      , "-1     -300. 300.   -200.   200.    -1.   1.      -1.   1.       0.  0.");
  cenv.SetValue("fedra.vsa.DoFine"      ,  1    );
  cenv.SetValue("fedra.vsa.SaveCouples" ,  1    );
  
  cenv.SetValue("viewsideal.outdir"          , ".."  );
  cenv.SetValue("viewsideal.env"             , "viewsideal.rootrc");
  cenv.SetValue("viewsideal.EdbDebugLevel"   ,  1    );
}

int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("viewsidealenv");
  set_default(cenv);
  gEDBDEBUGLEVEL     = cenv.GetValue("viewsideal.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("viewsideal.env"            , "viewsideal.rootrc");
  const char *outdir = cenv.GetValue("viewsideal.outdir"         , "..");
  
  bool      do_single   = false;
  bool      do_set      = false;
  EdbID     id;
  
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
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
  }
  
  if(!(do_single||do_set))   { print_help_message(); return 0; }
  if(  do_single&&do_set )   { print_help_message(); return 0; }
  
  cenv.SetValue("viewsideal.env"            , env);
  cenv.ReadFile( cenv.GetValue("viewsideal.env"   , "viewsideal.rootrc") ,kEnvLocal);
  cenv.SetValue("viewsideal.outdir"         , outdir);
  
  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("viewsideal.outdir","..");
  cenv.WriteFile("viewsideal.save.rootrc");
  
  printf("\n----------------------------------------------------------------------------\n");
  printf("viewsideal  %s\n"      ,id.AsString()	   );
  printf(  "----------------------------------------------------------------------------\n\n");
  printf("fedra.vsa.ImageCorrSide1: %s\n", cenv.GetValue("fedra.vsa.ImageCorrSide1","") );
  printf("fedra.vsa.ImageCorrSide2: %s\n", cenv.GetValue("fedra.vsa.ImageCorrSide2","") );
  
  if(do_single) 
  {
    MosaicAlID(id, cenv);
  }
  if(do_set) 
  {
    printf("\n----------------------------------------------------------------------------\n");
    printf("viewsideal set with id= %s\n", id.AsString() );
    printf("----------------------------------------------------------------------------\n\n");
    
    EdbScanSet *ss = sproc.ReadScanSet(id);
    if(ss)
    {
      MosaicAlSet(*ss,cenv);
    }
  }  
  cenv.WriteFile("viewsideal.save.rootrc");
  return 1;
}


//______________________________________________________________________________
void MosaicAlID(EdbID id, const TEnv &cenv)
{
  EdbMosaicAl mal;
  mal.ProcRun( id, cenv );
}

//______________________________________________________________________________
void MosaicAlSet(EdbScanSet &sc, const TEnv &cenv)
{
  if(sc.eIDS.GetSize()<1) return;
  for(int i=0; i<sc.eIDS.GetSize(); i++) {
    EdbID *id = (EdbID *)(sc.eIDS.At(i));        if(!id)    continue;
    //EdbPlateP *plate = sc.GetPlate(id->ePlate);  if(!plate) continue;
    EdbMosaicAl mal;
    mal.ProcRun( *id, cenv );
  }
}
