//-- Author :  Valeri Tioukov   09/02/2023

#include <string.h>
#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include <TChain.h>
#include "EdbLog.h"
#include "EdbRunAccess.h"
#include "EdbLinking.h"
#include "EdbScanProc.h"
#include "EdbPlateAlignment.h"
#include "EdbMosaic.h"
#include "EdbMosaicIO.h"

using namespace std;
using namespace TMath;
void LinkPlate( EdbID id, int from, int nfrag, TEnv &env );
void MergePlate( EdbID id, int from, int nfrag );
void AlignWithBeam( const float beam[4], 
			EdbPattern &p1, EdbPattern &p2 );

void print_help_message()
{
  cout<< "\nUsage: \n";
  cout<< "\t  moslink  -id=ID   [-from=frag0 -nfrag=N  -merge -v=DEBUG] \n";
  cout<< "\t  moslink  -set=ID  [-from=frag0 -nfrag=N  -merge -v=DEBUG] \n";

  cout<< "\t\t  ID    - id of the data piece or data set formed as BRICK.PLATE.MAJOR.MINOR \n";
  cout<< "\t\t  frag0 - the first fragment (default: 0) \n";
  cout<< "\t\t  N     - number of fragments to be processed (default: upto 1000000, stop at first empty) \n";
  cout<< "\t\t  merge - merge all fragments into one cp file \n";
  
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (moslink.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into moslink.save.rootrc file\n";
  cout<<endl;
}

//----------------------------------------------------------------------------------------
void set_default_link(TEnv &cenv)
{
  // default parameters for the new linking
  cenv.SetValue("fedra.link.AFID"                ,  1   );   // 1 is usually fine for scanned data; for the db-read data use 0!
  cenv.SetValue("fedra.link.DoImageCorr"           , 0  );
  cenv.SetValue("fedra.link.ImageCorrSide1"           , "1. 1. 0.");
  cenv.SetValue("fedra.link.ImageCorrSide2"           , "1. 1. 0.");
  
  cenv.SetValue("fedra.link.DoImageMatrixCorr"              , 0  );
  cenv.SetValue("fedra.link.ImageMatrixCorrSide1"           , "");
  cenv.SetValue("fedra.link.ImageMatrixCorrSide2"           , "");
  
  cenv.SetValue("fedra.link.CheckUpDownOffset"   ,  1   );   // check dXdY offsets between up and correspondent down views
  cenv.SetValue("fedra.link.BinOK"               , 6.   );
  cenv.SetValue("fedra.link.NcorrMin"            , 100  );
  cenv.SetValue("fedra.link.DoCorrectShrinkage"  , false );
  cenv.SetValue("fedra.link.read.InvertSides"    , 0  );
  cenv.SetValue("fedra.link.read.HeaderCut"      , "1"  );
  cenv.SetValue("fedra.link.read.UseDensityAsW"      , false  );
  cenv.SetValue("fedra.link.read.ICUT"           , "-1     -500. 500.   -500.   500.    -1.   1.      -1.   1.       0.  50.");
  cenv.SetValue("fedra.link.RemoveDoublets"      , "1    2. .01   1");  //yes/no   dr  dt  checkview(0,1,2)
  cenv.SetValue("fedra.link.DumpDoubletsTree"    , true );
  cenv.SetValue("fedra.link.shr.NsigmaEQ"        , 7.5  );
  cenv.SetValue("fedra.link.shr.Shr0"            , .990  );
  cenv.SetValue("fedra.link.shr.DShr"            , .3   );
  cenv.SetValue("fedra.link.shr.ThetaLimits"     ,"0.0  1.");
  cenv.SetValue("fedra.link.DoCorrectAngles"     , true );
  cenv.SetValue("fedra.link.ang.Chi2max"         , 1.5  );
  cenv.SetValue("fedra.link.DoFullLinking"       , true );
  cenv.SetValue("fedra.link.full.NsigmaEQ"       , 5.5  );
  cenv.SetValue("fedra.link.full.DR"             , 20.  );
  cenv.SetValue("fedra.link.full.DT"             , 0.1  );
  cenv.SetValue("fedra.link.full.CHI2Pmax"       , 3.   );
  cenv.SetValue("fedra.link.DoSaveCouples"       , true );
  cenv.SetValue("fedra.link.Sigma0"         , "1 1 0.007 0.007");
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

int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("moslinkenv");
  gEDBDEBUGLEVEL     = cenv.GetValue("moslink.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("moslink.env"            , "moslink.rootrc");
  const char *outdir = cenv.GetValue("moslink.outdir"         , "..");
  
  bool      do_single   = false;
  bool      do_set      = false;
  bool      do_merge    = false;
  EdbID     id;
  int       from_fragment=0;
  int       n_fragments=0;

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
    else if(!strncmp(key,"-from=",6))
    {
      if(strlen(key)>6) from_fragment = atoi(key+6);
    }
    else if(!strncmp(key,"-nfrag=",7))
    {
      if(strlen(key)>7) n_fragments = atoi(key+7);
    }
    else if(!strncmp(key,"-merge",6))
    {
      do_merge=true;
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
  }

  if(!(do_single||do_set))   { print_help_message(); return 0; }
  if(  do_single&&do_set )   { print_help_message(); return 0; }

  set_default_link(cenv);
  cenv.SetValue("moslink.env"            , env);
  cenv.ReadFile( cenv.GetValue("moslink.env"   , "moslink.rootrc") ,kEnvLocal);
  cenv.SetValue("moslink.outdir"         , outdir);

  EdbScanProc sproc;
  sproc.eProcDirClient = cenv.GetValue("moslink.outdir","..");
  cenv.WriteFile("moslink.save.rootrc");

  printf("\n----------------------------------------------------------------------------\n");
  printf("moslink  %s\n"      ,id.AsString()	   );
  printf(  "----------------------------------------------------------------------------\n\n");


  if(do_single) 
  {
    if(do_merge) MergePlate( id, from_fragment, n_fragments );
    else         LinkPlate(  id, from_fragment, n_fragments, cenv);
  }
  else if(do_set)
  {
    EdbScanSet *ss = sproc.ReadScanSet(id);
    if(ss) {
      int n = ss->eIDS.GetSize();
      for(int i=0; i<n; i++) {
	EdbID *id_pl   = ss->GetID(i);
	if(id_pl) 
	{
	  if(do_merge) MergePlate( *id_pl, from_fragment, n_fragments );
	  else         LinkPlate(  *id_pl, from_fragment, n_fragments, cenv);
	}
      }
    }
  }

  cenv.WriteFile("moslink.save.rootrc");
  return 1;
}

void MergePlate( EdbID id, int mfrom, int nfrag )
{
  if(nfrag==0)
  {
    TString file;
    file.Form("p%3.3d/%d.%d.%d.%d.mos.root",
	      id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor);
    EdbMosaicIO mio;
    mio.Init( file.Data() );
    EdbLayer *mapside = mio.GetCorrMap( id.ePlate, 1 );
    nfrag=mapside->Map().Ncell();
  }

  Log(2,"MergePlate", "%s %d:%d",id.AsString(), mfrom, mfrom+nfrag);  
  TChain allcp("couples");

  for( int mid=mfrom; mid<mfrom+nfrag; mid++ ) {
    allcp.Add( Form("p%3.3d/%d.%d.%d.%d.%d.cp.root", 
		    id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor, mid));
  }
  allcp.Merge( Form("p%3.3d/%d.%d.%d.%d.cp.root", 
		    id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor));
}


void LinkPlate( EdbID id, int from, int nfrag, TEnv &cenv )
{
  EdbScanProc sproc;
  sproc.eProcDirClient="..";
  
  EdbID id0=id; id0.ePlate=0;
  EdbScanSet *ss = sproc.ReadScanSet(id0);
  EdbPlateP *plate = ss->GetPlate(id.ePlate);
  
  EdbMosaicIO mio;
  mio.Init( Form("p%3.3d/%d.%d.%d.%d.mos.root",
		 id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor) );
  
  EdbLayer *mapside = mio.GetCorrMap( id.ePlate, 1 );  
  int nc=mapside->Map().Ncell();
  if(nfrag==0) nfrag=nc;
  
  for( int i=from; i<nfrag; i++ )
  {
    EdbPattern *p1 = mio.GetFragment( id.ePlate, 1, i, true); // (plate,side,id)
    EdbPattern *p2 = mio.GetFragment( id.ePlate, 2, i, true);
    Log(1,"moslink::LinkPlate","fragment %d (of %d) with %d & %d mt", i,nfrag, p1->N(), p2->N() );
    if(p1&&p2) 
    {
      p1->SetScanID(id);
      p2->SetScanID(id);
 
      EdbLinking link;
      link.InitOutputFile( Form("p%3.3d/%d.%d.%d.%d.%d.cp.root", 
			   id.ePlate, id.eBrick, id.ePlate, id.eMajor, id.eMinor, i) );
      link.Link( *p1, *p2, *(plate->GetLayer(2)), *(plate->GetLayer(1)), cenv );
      link.CloseOutputFile();
      delete p1;
      delete p2;
    }
  }
}

//-----------------------------------------------------------------------
void AlignWithBeam( const float beam[4], 
			EdbPattern &p1, EdbPattern &p2 )
{
  //int deb=gEDBDEBUGLEVEL;
  //gEDBDEBUGLEVEL=3;
  EdbPlateAlignment av;
  av.eNoScaleRot  = 0;         // calculate shift only
  av.eOffsetMax   = 10.;
  av.eDZ          = 0.;
  av.eDPHI        = 0.;
  av.eDoFine      = 1;
  av.eSaveCouples = 1;
  av.SetSigma( 1, 0.007 );
  av.eDoublets[0] = av.eDoublets[1]=1;
  av.eDoublets[2] = av.eDoublets[3]=0.005;

  av.eSaveCouples=1;
  av.InitOutputFile( Form( "p%.3d/%d_%d.al.avb.root", p1.ScanID().ePlate, p1.ID(), p2.ID() ) ); 
  av.Align( p1, p2, 0); //-190
  EdbAffine2D *affXY = av.eCorrL[0].GetAffineXY();
  EdbAffine2D *affTXTY = av.eCorrL[0].GetAffineTXTY();
//  if(av.eNcoins > eMinPeak ) 
//  {
//    if(do_shift) p1.Transform( affXY );
//    aff.Transform( affXY );
//  }
  av.CloseOutputFile();
//  return av.eNcoins; 
    //gEDBDEBUGLEVEL=deb;
}

