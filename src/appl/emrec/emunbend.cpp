//-- Author :  Valeri Tioukov   30/07/2024
// read tracks tree with given cuts and unbend dataset
#include<stdlib.h>
#include<string>
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbUnbender.h"

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout<< "\nUsage: \n\t  emunbend -set=ID [ -alg3a -alg3g -alg5g -corraff -v=DEBUG] \n";
  cout<< "\t  Unbend dataset and update set.root\n";
  cout<< "\t\t -alg3a      -  triplets simple (default) \n";
  cout<< "\t\t -alg3g      -  triplets generic \n";
  cout<< "\t\t -alg5g      -  quintets generic \n";
  cout<< "\t\t -corraff    -  save corrections to set.root file\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default(TEnv &cenv)
{
  cenv.SetValue("unbend.outdir"          , "..");
  cenv.SetValue("unbend.read_cut"        , "1");
  cenv.SetValue("unbend.save_hist"       ,  0);
  cenv.SetValue("unbend.save_tree"       ,  0);
//  cenv.SetValue("unbend.NCPmin"          , "50");
  cenv.SetValue("unbend.EdbDebugLevel"   ,   1 );
}

//----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  TEnv cenv("unbendenv");
  set_default(cenv);
  gEDBDEBUGLEVEL        = cenv.GetValue("unbend.EdbDebugLevel" ,  1  );
  const char *outdir    = cenv.GetValue("unbend.outdir"        , "..");
  
  bool do_set     = false;
  bool do_3a      = false;
  bool do_3g      = false;
  bool do_5g      = false;
  bool do_corraff = false;
  EdbID idset;

  for(int i=1; i<argc; i++ ) 
  {
    char *key  = argv[i];
    if(!strncmp(key,"-set=",5))
    {
      if(strlen(key)>5)	if(idset.Set(key+5))   do_set=true;
    }
    if(!strncmp(key,"-corraff",8))  do_corraff = true;
    if(!strncmp(key,"-alg3a",6))    do_3a      = true;
    if(!strncmp(key,"-alg3g",6))    do_3g      = true;
    if(!strncmp(key,"-alg5g",6))    do_5g      = true;
  }

  cenv.ReadFile( "unbend.rootrc" ,kEnvLocal);
  cenv.WriteFile("unbend.save.rootrc");

  if(do_set)
  {
    EdbScanProc sproc;
    sproc.eProcDirClient=outdir;
    EdbPVRec ali;
    int ntrMin = cenv.GetValue("unbend.ntrMin"        , 100);
    const char *cut = cenv.GetValue("unbend.read_cut"        , "1");
    int ntr = sproc.ReadTracksTree( idset,ali, cut );
    Log(1,"unbend","read %d tracks of set %s with cut %s",ntr,idset.AsString(),cut);
    if(ntr<ntrMin) {
      Log(1,"unbend","Warning: tracks not enough (%d < %d) - exit",ntr,ntrMin);
    } else {
      EdbScanSet *ss = sproc.ReadScanSet(idset);
      EdbUnbender ub;
      ub.do_save_hist = cenv.GetValue("unbend.save_hist"        , 0);
      ub.do_save_tree = cenv.GetValue("unbend.save_tree"        , 0);
      
      if     (do_3g) ub.Unbend3g(ali,*ss, cenv);
      else if(do_5g) ub.Unbend5g(ali,*ss, cenv);
      else           ub.Unbend3a(ali,*ss, cenv);
      
      if(do_corraff) sproc.WriteScanSet(idset,*ss);
    }
  }
}
