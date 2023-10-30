//-- Author :  Valeri Tioukov   23/10/2023
#include <assert.h>
#include <string.h>
#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbCouplesTree.h"

using namespace std;

EdbScanProc sproc;
void DoTraceBackCP( EdbID id, const char *treename, TEnv &cenv );
void FillRawPat( const EdbRunAccess &run, 
		 const EdbPattern &pp, 
		       EdbPattern &ppr, 
		       EdbPattern &ppv, 
		 const TArrayI &viewid,
		 const TArrayI &ind
	       );


void GroupViews( EdbPattern &p, EdbCell1 &c );
void GroupViews( EdbID &cpid, const EdbPattern &pp1, const EdbPattern &pp2, EdbCell1 &cv1, EdbCell1 &cv2 );
void WriteRaw( EdbRunAccess &run, EdbRun &runout, EdbCell1 &cv, int rs );

void print_help_message()
{
  cout<< "\n For a given file *.cp.root extract raw segments for all couples in doublets tree \n";
  cout<<   " and create a parallel tree with them\n";
  cout<<   "\n";
  cout<< "\nUsage: \n\t  emdbl_tb -id=ID [-v=DEBUG] \n";
  cout<< "\nExample: \n";
  cout<< "\t  emdbl_tb -id=4554.10.1.1 -v=2\n";
  cout<< "\n If the data location directory if not explicitly defined\n";
  cout<< " the current directory will be assumed to be the brick directory \n";
  cout<< "\n If the parameters file (dbl_tb.rootrc) is not presented - the default \n";
  cout<< " parameters will be used. After the execution them are saved into dbl_tb.save.rootrc file\n";
  cout<<endl;
}

void set_default(TEnv &cenv)
{
  // default parameters
  cenv.SetValue("emdbl_tb.outdir"          , ".."  );
  cenv.SetValue("emdbl_tb.env"             , "traceback.rootrc");
  cenv.SetValue("emdbl_tb.EdbDebugLevel"   ,  1    );
}


int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  
  TEnv cenv("emdbl_tb");
  set_default(cenv);
  gEDBDEBUGLEVEL     = cenv.GetValue("emdbl_tb.EdbDebugLevel" , 1);
  const char *env    = cenv.GetValue("emdbl_tb.env"            , "dbl_tb.rootrc");
  const char *outdir = cenv.GetValue("emdbl_tb.outdir"         , "..");

  bool do_id=false;
  EdbID     id;

  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];
  
    if     (!strncmp(key,"-id=",4))
    {
      if(strlen(key)>4) if(id.Set(key+4))   do_id=true;
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3) gEDBDEBUGLEVEL = atoi(key+3);
    }
  }
  
  if(!(do_id))   { print_help_message(); return 0; }

  cenv.SetValue("emdbl_tb.env"            , env);
  cenv.ReadFile( cenv.GetValue("emdbl_tb.env"   , "traceback.rootrc") ,kEnvLocal);
  cenv.SetValue("emdbl_tb.outdir"         , outdir);
  sproc.eProcDirClient = cenv.GetValue("emdbl_tb.outdir","..");
  cenv.WriteFile("dbl_tb.save.rootrc");
  

  printf("\n----------------------------------------------------------------------------\n");
  printf("Traceback doublets from %s\n"   ,id.AsString());
  printf("----------------------------------------------------------------------------\n\n");

  DoTraceBackCP( id, "doublets1", cenv );
  DoTraceBackCP( id, "doublets2", cenv );

  return 0;
}

//----------------------------------------------------------------
void DoTraceBackCP( EdbID id, const char *treename, TEnv &cenv )
{
  // find raw segments corresponding to the couples tree and create new file .dbl.raw.root 
  // another with parallel couples tree in raw coordinates 
  // (to be joined as a fiends trees for analysis)
    
  TString cpfile;
  sproc.MakeFileName(cpfile, id,"cp.root");
  EdbCouplesTree cptree_in;
  cptree_in.InitCouplesTree(treename,cpfile.Data());
  EdbPattern pp,pp1,pp2;
  cptree_in.GetCPData( &pp,&pp1,&pp2 );
  
  int nseg = pp1.N();
  TArrayI   viewid1(nseg);
  TArrayI   ind1(nseg);
  TArrayI   viewid2(nseg);
  TArrayI   ind2(nseg);

  for(int isg=0; isg<nseg; isg++)
  {
    viewid1[isg] = pp1.GetSegment(isg)->Vid(0);
    viewid2[isg] = pp2.GetSegment(isg)->Vid(0);
  }
  TMath::Sort(nseg, viewid1.GetArray(), ind1.GetArray(),0);
  TMath::Sort(nseg, viewid2.GetArray(), ind2.GetArray(),0);
 
  EdbRunAccess run;
  sproc.InitRunAccess(run, id);
  EdbPattern pp1r(0,0,0,nseg);
  EdbPattern pp2r(0,0,0,nseg);
  EdbPattern pp1v(0,0,0,nseg);
  EdbPattern pp2v(0,0,0,nseg);
 
  FillRawPat(run, pp1, pp1r, pp1v, viewid1, ind1);
  FillRawPat(run, pp2, pp2r, pp2v, viewid2, ind2);

  EdbCouplesTree cptree_out;
  sproc.MakeFileName(cpfile, id,"dbl.cp.root");
  cptree_out.InitCouplesTree(treename,cpfile.Data(),"UPDATE");
  
  for(int i=0; i<nseg; i++)
  {
    EdbSegP *s1 = pp1r.GetSegment(i);
    EdbSegP *s2 = pp2r.GetSegment(i);
    EdbSegP *sv1 = pp1v.GetSegment(i);
    EdbSegP *sv2 = pp2v.GetSegment(i);
    sv1->SetTX(sv2->X());
    sv1->SetTY(sv2->Y());
   cptree_out.Fill( s1, s2, sv1);
  }
  cptree_out.WriteTree();
  cptree_out.Close();
  cptree_in.Close();
}

void FillRawPat( const EdbRunAccess &run, 
		 const EdbPattern &pp, 
		 EdbPattern &ppr, 
		 EdbPattern &ppv, 
		 const TArrayI &viewid,
		 const TArrayI &ind
	       )
{
  int nseg = pp.N();
  int iv=-1;
  EdbView *vin  =0;
  EdbSegP s;
  EdbSegP sv;
  for(int isg=0; isg<nseg; isg++)
  {
    int k =  ind[isg];
    if(viewid[k]!=iv) {
      iv=viewid[k];
      vin  = run.GetRun()->GetEntry(iv);
      printf("isg=%d    iv= %d\n",isg,iv);
    }
    EdbSegP    *sp   = pp.GetSegment(k);
    EdbSegment *sraw = vin->GetSegment(sp->Vid(1));
    s.Set(sraw->GetID(), sraw->GetX0(),sraw->GetY0(), sraw->GetTx(), sraw->GetTy(), sraw->GetPuls(),0 );
    s.SetZ(sraw->GetZ0());
    s.SetVid(sp->Vid(0),sp->Vid(1));
    s.SetAid(sp->Aid(0),sp->Aid(1),sp->Aid(2));
    ppr.AddSegment(k,s);
    
    sv.Set(vin->GetViewID(), vin->GetXview(), vin->GetYview(), 0, 0, 0, 0);
    ppv.AddSegment(k,sv);
  } 
 

}

//-------------------------------------------------------------------------------------------------
void WriteRaw(  EdbRunAccess &run, EdbRun &runout, EdbCell1 &cv, int rs )
{
  int nv = cv.Ncell();
  for( int iv=0; iv<nv; iv++ ) {
    int nsv=cv.Bin(iv);
    if(!nsv)  continue;

    printf("iv = %d nsv = %d\n",iv,nsv);
    
    EdbView *vin  = run.GetRun()->GetEntry(iv);
    EdbView *vout = runout.GetView();
    vout->GetHeader()->Copy( vin->GetHeader() );
    
    for( int is=0; is<nsv; is++ ) {
      EdbSegP *sp = (EdbSegP *)(cv.GetObject(iv,is));
      
      EdbSegment *sraw = vin->GetSegment(sp->Vid(1));
      if(rs>0) run.ApplyCorrections( *vin, *sraw, rs );

      if(gEDBDEBUGLEVEL>3) {
        printf("\n");
        sp->PrintNice();
        sraw->Print();
      }

      vout->AddSegment(sraw);
      
      EdbTrack t( sp->X(), sp->Y(), sp->Z(), sp->TX(), sp->TY(),sp->DZ(), sp->ID() ); 
      vout->AddTrack(&t);  // save sp for crosscheck. Note that the last emlink should be done with all corrections off for consistency

    }
    runout.AddView(vout);
  }
}
