//-- Author :  Valeri Tioukov   11/06/2008

#include <string.h>
#include <iostream>
#include <TEnv.h>
#include <TRint.h>
#include "EdbLog.h"
#include "EdbScanProc.h"

using namespace std;

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout<< "\nUsage: \n\t  emlink -id=ID [-p=NPRE -f=NFULL -o=DATA_DIRECTORY -a -v=DEBUG] \n";
  cout<< "\t  emlink -set=ID [-p=NPRE -f=NFULL -o=DATA_DIRECTORY -a -v=DEBUG] \n\n";
  cout<< "\t\t  ID    - id of the data piece or dataset formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout<< "\t\t  NPRE  - number of the prealignments (default is 0)\n";
  cout<< "\t\t  NFULL - number of the fullalignments (default is 0)\n";
  cout<< "\t\t  -a    - apply the angular correction when do prelinking (default is not)\n";
  cout<< "\t\t  -new  - new linking (renewed 06-04-11)\n";
  cout<< "\t\t  -view - check view overlaps\n";
  cout<< "\t\t  DEBUG - verbosity level: 0-print nothing, 1-errors only, 2-normal, 3-print all messages\n";
  cout<< "\t\t  -ix   - x cell bin\n";
  cout<< "\t\t  -iy   - y cell bin\n";
  cout<< "\nExample: \n";
  cout<< "\t  emlink -id=4554.10.1.0 -o=/scratch/BRICKS \n";
  cout<< "\n If the data location directory is not explicitly defined\n";
  cout<< " the current directory assumed to be the brick directory \n";
  cout<< "\n If the parameters file (link.rootrc) is not presented - the default \n";
  cout<< " parameters are used. After the execution them will be saved into link.save.rootrc\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default(TEnv &cenv)
{
  // default parameters for the new linking
  cenv.SetValue("fedra.link.AFID"                ,  1   );   // 1 is usually fine for scanned data; for the db-read data use 0!
  cenv.SetValue("fedra.link.PixelCorr"           , "0 1. 1.");
  cenv.SetValue("fedra.link.CheckUpDownOffset"   ,  1   );   // check dXdY offsets between up and correspondent down views
  cenv.SetValue("fedra.link.BinOK"               , 6.   );
  cenv.SetValue("fedra.link.NcorrMin"            , 100  );
  cenv.SetValue("fedra.link.DoCorrectShrinkage"  , true );
  cenv.SetValue("fedra.link.read.InvertSides"    , 0  );
  cenv.SetValue("fedra.link.read.HeaderCut"      , "1"  );
  cenv.SetValue("fedra.link.read.UseDensityAsW"      , false  );
  cenv.SetValue("fedra.link.read.ICUT"           , "-1     -500. 500.   -500.   500.    -1.   1.      -1.   1.       0.  50.");
  cenv.SetValue("fedra.link.RemoveDoublets"      , "1    5. .1   1");  //yes/no   dr  dt  checkview(0,1,2)
  cenv.SetValue("fedra.link.DumpDoubletsTree"    , false );
  cenv.SetValue("fedra.link.shr.NsigmaEQ"        , 7.5  );
  cenv.SetValue("fedra.link.shr.Shr0"            , .9  );
  cenv.SetValue("fedra.link.shr.DShr"            , .3   );
  cenv.SetValue("fedra.link.shr.ThetaLimits"     ,"0.005  1.");
  cenv.SetValue("fedra.link.DoCorrectAngles"     , true );
  cenv.SetValue("fedra.link.ang.Chi2max"         , 1.5  );
  cenv.SetValue("fedra.link.DoFullLinking"       , true );
  cenv.SetValue("fedra.link.full.NsigmaEQ"       , 5.5  );
  cenv.SetValue("fedra.link.full.DR"             , 30.  );
  cenv.SetValue("fedra.link.full.DT"             , 0.1  );
  cenv.SetValue("fedra.link.full.CHI2Pmax"       , 3.   );
  cenv.SetValue("fedra.link.DoSaveCouples"       , true );
  cenv.SetValue("fedra.link.DumpDoubletsTree"    , false);
  cenv.SetValue("fedra.link.Sigma0"         , "1 1 0.013 0.013");
  cenv.SetValue("fedra.link.PulsRamp0"      , "6 9");
  cenv.SetValue("fedra.link.PulsRamp04"     , "6 9");
  cenv.SetValue("fedra.link.Degrad"         ,  5   );
  
  cenv.SetValue("fedra.link.LLfunction"     , "0.256336-0.16489*x+2.11098*x*x" );
  cenv.SetValue("fedra.link.CPRankingAlg"   , 0 );

  cenv.SetValue("emlink.reportfileformat"   , "pdf" );

  cenv.SetValue("emlink.outdir"          , "..");
  cenv.SetValue("emlink.env"             , "link.rootrc");
  cenv.SetValue("emlink.EdbDebugLevel"   , 1);
}

//----------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }

  TEnv cenv("linkenv");
  set_default(cenv);
  gEDBDEBUGLEVEL        = cenv.GetValue("emlink.EdbDebugLevel" , 1);
  const char *env       = cenv.GetValue("emlink.env"           , "link.rootrc");
  const char *outdir    = cenv.GetValue("emlink.outdir"        , "..");

  bool      do_id      = false;
  bool      do_set     = false;
  bool      do_new     = false;
  bool      do_check   = false;
  bool      do_check_view   = false;
  Int_t     brick=0, plate=0, major=0, minor=0;
  Int_t     npre=0,  nfull=0;
  Int_t	    ix = -1, iy = -1; //cell numbers
  Int_t     correct_ang=0;

  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];

    if     (!strncmp(key,"-id=",4)) 
      {
	if(strlen(key)>4)	sscanf(key+4,"%d.%d.%d.%d",&brick,&plate,&major,&minor);
	do_id=true;
      }
    else if(!strncmp(key,"-set=",5))
      {
	if(strlen(key)>5)	sscanf(key+5,"%d.%d.%d.%d",&brick,&plate,&major,&minor);
	do_set=true;
      }
    else if(!strncmp(key,"-o=",3)) 
      {
	if(strlen(key)>3)	outdir=key+3;
      }
    else if(!strncmp(key,"-p=",3))
      {
	if(strlen(key)>3)	npre = atoi(key+3);
      }
    else if(!strncmp(key,"-f=",3))
      {
	if(strlen(key)>3)	nfull = atoi(key+3);
      }
    else if(!strncmp(key,"-a",2))
      {
	correct_ang=1;
      }
    else if(!strncmp(key,"-new",4))
      {
	do_new=true;
      }
    else if(!strncmp(key,"-check",6))
      {
	do_check=true;
      }
    else if(!strncmp(key,"-view",5))
      {
	do_check_view=true;
      }
    else if(!strncmp(key,"-v=",3))
      {
	if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
      }
    else if(!strncmp(key,"-ix=",4))
      {
	if(strlen(key)>4)	ix = atoi(key+4);
      }
    else if(!strncmp(key,"-iy=",4))
      {
	if(strlen(key)>4)	iy = atoi(key+4);
      }
  }

  if(!(do_id||do_set))   { print_help_message(); return 0; }

  cenv.SetValue("emlink.env"            , env);
  cenv.ReadFile( cenv.GetValue("emlink.env"   , "link.rootrc") ,kEnvLocal);

  const char* reportfileformat = cenv.GetValue("emlink.reportfileformat","pdf");
  cenv.SetValue("emlink.outdir"         , outdir);
  cenv.WriteFile("link.save.rootrc");


  EdbScanProc sproc;
  sproc.eProcDirClient=outdir;
  if(do_id) {
    printf("\n----------------------------------------------------------------------------\n");
    printf("link run with id= %d.%d.%d.%d\n", brick,plate, major,minor);
    printf("----------------------------------------------------------------------------\n\n");

    EdbID id(brick,plate,major,minor);
    if(do_new) {
      EdbID id0=id; id0.ePlate=0;
      EdbScanSet *ss = sproc.ReadScanSet(id0);
      if(ss) {
        EdbPlateP *plate = ss->GetPlate(id.ePlate);
        printf("test0, %d, %d",ix,iy);
        if(plate) sproc.LinkRunTest(id, *plate, cenv, ix, iy);
     }
    }
    else if(do_check_view) sproc.CheckViewOverlaps(id, cenv);
    else if(npre+nfull>0) sproc.LinkRunAll(id, npre, nfull, correct_ang);
  }
  if(do_set) {
    printf("\n----------------------------------------------------------------------------\n");
    printf("link set with id= %d.%d.%d.%d\n", brick,plate, major,minor);
    printf("----------------------------------------------------------------------------\n\n");

    EdbID id(brick,plate,major,minor);
    EdbScanSet *ss = sproc.ReadScanSet(id);
    if(ss) {
      if(do_new) {
        sproc.LinkSetNewTest(*ss, cenv, ix, iy);
        sproc.MakeLinkSetSummary(id, reportfileformat);
      }
      else  if(npre+nfull>0)  sproc.LinkSet(*ss, npre, nfull, correct_ang);
    }
    
    if(do_check) 
       {
	 //int argc2=1;
	 //char *argv2[]={"-l"};
	 //TRint app("APP",&argc2, argv2);
	 sproc.MakeLinkSetSummary(id, reportfileformat);
	 //app.Run();
       }
  }
  cenv.WriteFile("link.save.rootrc");

  return 1;
}
