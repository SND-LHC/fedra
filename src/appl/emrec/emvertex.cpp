#include <iostream>
#include "TRint.h"
#include "TStyle.h"
#include "TArrayL64.h"
#include "TMath.h"
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbProcPars.h"
#include "EdbVertex.h"
#include "EdbDisplay.h"
#include "EdbCombGen.h"
#include "EdbVertexComb.h"

using namespace std;
using namespace TMath;

//-----------------------------------------------------------------------------
EdbScanCond  gCond;
EdbID        idset;
EdbPVRec     gAli;
EdbScanProc  gSproc;
EdbVertexRec gEVR;

void VertexRec(EdbID id, TEnv &cenv);
void ReadVertex(EdbID id,TEnv &env);
void MakeScanCondBT(EdbScanCond &cond, TEnv &env);
void SetTracksErrors(TObjArray &tracks, EdbScanCond &cond, float p, float m);
void do_vertex(TEnv &env);
void AddCompatibleTracks(EdbPVRec &v_trk, EdbPVRec &v_vtx);
bool IsCompatible(EdbVertex &v, EdbTrackP &t);
void Display( const char *dsname,  EdbVertexRec *evr, TEnv &env );

//----------------------------------------------------------------------------------------
void print_help_message()
{
  cout<< "\n Vertex reconstruction in the volume. Input *.trk.root, output *.vtx.root\n";
  
  cout<< "\nUsage: \n\t  emvertex -set=ID [-v=DEBUG] \n";
  cout<< "\n\t  emvertex -set=ID [-r -display -v=DEBUG]  \n";
  cout<< "\t\t  r       - read found vertices from *.vtx.root\n";
  cout<< "\t\t  display - start interactive event display\n";
  cout<< "\t\t  DEBUG   - verbosity level: 0-print nothing, 1-errors only, 2-normal, 3-print all messages\n";
  
  cout<< "\n If the parameters file (vertex.rootrc) is not presented - the default \n";
  cout<< " parameters are used. After the execution them will be saved into vertex.save.rootrc\n";
  cout<<endl;
}

//---------------------------------------------------------------------
void set_default(TEnv &env)
{
  // default parameters

  env.SetValue("emvertex.vtx.DZmax"         , 3000.);
  env.SetValue("emvertex.vtx.ProbMinV"      , 0.001);
  env.SetValue("emvertex.vtx.ImpMax"        , 10.);
  env.SetValue("emvertex.vtx.UseMom"        , false);
  env.SetValue("emvertex.vtx.UseSegPar"     , false);
  env.SetValue("emvertex.vtx.QualityMode"   , 0);  // (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
  env.SetValue("emvertex.vtx.cutvtx"        , "(flag==0||flag==3)&&n>4");
  env.SetValue("emvertex.vtx.cuttr"         , "nseg>4&&npl<50");

  env.SetValue("emvertex.addtr.doit"        ,  0 );
  env.SetValue("emvertex.addtr.cuttr"       , "1");

  env.SetValue("emvertex.edd.ajustseg"      ,  0);

  env.SetValue("emvertex.trfit.doit"     ,  1 );
  env.SetValue("emvertex.trfit.P"        , 10 );
  env.SetValue("emvertex.trfit.M"        ,  0.139);
  env.SetValue("emvertex.bt.Sigma0", "0.2 0.2 0.002 0.002" );
  env.SetValue("emvertex.bt.Degrad", 5. );
}

//---------------------------------------------------------------------
void AjustSegmentsDisplay( TObjArray &tarr )
{
  int n=tarr.GetEntries();
  for(int i=0; i<n; i++)
  {
    EdbTrackP *t = (EdbTrackP*)(tarr.At(i));
    int nseg = t->N();
    for(int j=0; j<nseg; j++)
    {
      EdbSegP *s=t->GetSegment(j);
      s->SetDZ(300);
      s->SetW(10);
    }    
  } 
}

//---------------------------------------------------------------------
void Display( const char *dsname,  EdbVertexRec *evr, TEnv &env )
{
  TObjArray *varr = new TObjArray();
  TObjArray *tarr = new TObjArray();
  
  EdbVertex *v=0;
  EdbTrackP *t=0;
  
  int nv = evr->Nvtx();
  printf("nv=%d\n",nv);
  if(nv<1) return;
  
  for(int i=0; i<nv; i++) {
    v = (EdbVertex *)(evr->eVTX->At(i));
    varr->Add(v);
    v->PrintGeom();
//    v->SaveGeom();
    for(int j=0; j<v->N(); j++) {
      EdbTrackP *t = v->GetTrack(j);
      tarr->Add( t );
    }
  }
  
  EdbPVRec *pvr = evr->ePVR;
  if(pvr) {
    int ntr = pvr->Ntracks();
    for(int i=0; i<ntr; i++) 
    {
      EdbTrackP *t = pvr->GetTrack(i);
      if(t->Flag()==999999) tarr->Add(t);
    }
  }
  
  if( env.GetValue("emvertex.edd.ajustseg"     ,  0 ) ) AjustSegmentsDisplay( *tarr );
    
  gStyle->SetPalette(1);
  
  EdbDisplay *ds = EdbDisplay::EdbDisplayExist(dsname);
  if(!ds)  ds=new EdbDisplay(dsname,-10000.,10000.,-10000.,10000.,-10000., 10000.);
  ds->SetVerRec(evr);
  ds->SetArrTr( tarr );
  printf("%d tracks to display\n", tarr->GetEntries() );
  ds->SetArrV( varr );
  printf("%d vertex to display\n", varr->GetEntries() );
  //ds->SetArrSegG( tsegG );
  //printf("%d primary tracks to display\n", tsegG->GetEntries() );
  ds->SetDrawTracks(env.GetValue("emvertex.edd.DrawTracks"     ,  14));
  ds->SetDrawVertex(env.GetValue("emvertex.edd.DrawVertex"     ,  1));
  //   //ds->SetView(90,180,90);
  
  ds->GuessRange(2000,2000,30000);
  ds->SetStyle(1);
  ds->Draw();
  
  //  float s[3] = {0,0,0 };
  //float e[3] = {Vmc[0],Vmc[1],Vmc[2]+600};
  //ds->DrawRef(Vmc,e);
}


//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if (argc < 2)   { print_help_message();  return 0; }
  TEnv cenv("vertexenv");
  set_default(cenv);
  gEDBDEBUGLEVEL        = cenv.GetValue("emvertex.EdbDebugLevel" ,  1  );
  const char *outdir    = cenv.GetValue("emvertex.outdir"        , "..");
  gSproc.eProcDirClient=outdir;
  cenv.ReadFile( "vertex.rootrc" ,kEnvLocal);
 
  bool        do_set     = false;
  bool        do_display = false;
  bool        do_read    = false;
  
  for(int i=1; i<argc; i++ ) {
    char *key  = argv[i];
    if(!strncmp(key,"-set=",5))
    {
      if(strlen(key)>5)	  if(idset.Set(key+5))   do_set=true;
    }
    else if(!strncmp(key,"-r",2))
    {
      do_read=true;
    }
    else if(!strncmp(key,"-v=",3))
    {
      if(strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
    }
    else if(!strncmp(key,"-display",8))
    {
      do_display=true;
    }
  } 
  cenv.WriteFile("vertex.save.rootrc");
 
  if(do_set) 
  {
    if(do_read)
    {
      ReadVertex(idset,cenv);
    } 
    else 
    {
      Log(1,"vertex","set %s",idset.AsString());
      VertexRec(idset,cenv);
    }
  }
  
  cenv.WriteFile("vertex.save.rootrc");
  
  if(do_display)
  {
    int argc2=1;
    char *argv2[]={"-l"};
    TRint app("APP",&argc2, argv2);
    Display("display",&gEVR, cenv);
    app.Run();
  }
  
  return 0;
}

void ReadVertex(EdbID id, TEnv &env)
{
  MakeScanCondBT(gCond, env);
  gAli.SetScanCond( new EdbScanCond(gCond) );
  gEVR.eEdbTracks = gAli.eTracks;
  gEVR.eVTX       = gAli.eVTX;
  gEVR.SetPVRec(&gAli);

  gEVR.eDZmax      = env.GetValue("emvertex.vtx.DZmax"         , 3000.);
  gEVR.eProbMin    = env.GetValue("emvertex.vtx.ProbMinV"      , 0.001);
  gEVR.eImpMax     = env.GetValue("emvertex.vtx.ImpMax"        , 10.);
  gEVR.eUseMom     = env.GetValue("emvertex.vtx.UseMom"        , false);
  gEVR.eUseSegPar  = env.GetValue("emvertex.vtx.UseSegPar"     , false);
  gEVR.eQualityMode= env.GetValue("emvertex.vtx.QualityMode"   , 0);  // (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
  TCut cutvtx      = env.GetValue("emvertex.vtx.cutvtx"        , "(flag==0||flag==3)&&n>4");

  TObjArray* v_out = new TObjArray();
  
  EdbDataProc *dproc = new EdbDataProc();
  TString name;
  gSproc.MakeFileName(name,id,"vtx.root",false);
  int nvtx = dproc->ReadVertexTree(gEVR, name.Data(), cutvtx);
  if(nvtx) {
    int do_addtracks = env.GetValue("emvertex.addtr.doit"         , 0);
    if(do_addtracks)
    {
      TCut cuttr       = env.GetValue("emvertex.addtr.cuttr"        , "1");
      EdbPVRec *vtr = new EdbPVRec();
      vtr->SetScanCond( new EdbScanCond(gCond) );
      gSproc.ReadTracksTree( idset,*vtr, cuttr);
      AddCompatibleTracks( *vtr, gAli , v_out);  // assign to the vertices of gAli additional tracks from vtr if any
      EdbDataProc::MakeVertexTree(&v_out,"flag0.vtx.root");
    }
  }
}

void VertexRec(EdbID id, TEnv &env)
{
  /*
  float x=105;
  float y=163;
  float dx,dy;
  dx=dy=5000;
  float x0=x*1000;
  float y0=y*1000;
  TCut cutvol("cutvol",Form("abs(t.eX-%f)<%f&&abs(t.eY-%f)<%f",x0,dx+500,y0,dy+500));
  */
  TCut cuttr = env.GetValue("emvertex.vtx.cuttr" , "nseg>4&&npl<50");
//  TCut cut=cutvol&&cuttr;
  TCut cut=cuttr;
 
  MakeScanCondBT(gCond,env);
  gAli.SetScanCond( new EdbScanCond(gCond) );
  gSproc.ReadTracksTree( idset,gAli, cut );
  do_vertex(env);
}

void do_vertex(TEnv &env)
{
  //gAli.PrintSummary();
  bool do_trfit   = env.GetValue("emvertex.trfit.doit"     ,  1 );
  float pfit      = env.GetValue("emvertex.trfit.P"        , 10 );
  float mfit      = env.GetValue("emvertex.trfit.M"        ,  0.139);
  if(do_trfit) {
    SetTracksErrors( *(gAli.eTracks), gCond, pfit,mfit );
    //gAli.FitTracks(pfit,mfit );
  }

  gEVR.eEdbTracks = gAli.eTracks;
  gEVR.eVTX       = gAli.eVTX;
  gEVR.SetPVRec(&gAli);

  gEVR.eDZmax      = env.GetValue("emvertex.vtx.DZmax"         , 3000.);
  gEVR.eProbMin    = env.GetValue("emvertex.vtx.ProbMinV"      , 0.001);
  gEVR.eImpMax     = env.GetValue("emvertex.vtx.ImpMax"        , 10.);
  gEVR.eUseMom     = env.GetValue("emvertex.vtx.UseMom"        , false);
  gEVR.eUseSegPar  = env.GetValue("emvertex.vtx.UseSegPar"     , false);
  gEVR.eQualityMode= env.GetValue("emvertex.vtx.QualityMode"   , 0);  // (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)

  printf("%d tracks for vertexing\n",  gEVR.eEdbTracks->GetEntries() );
  
  int nvtx = gEVR.FindVertex();
  printf("%d 2-track vertexes was found\n",nvtx);

  if(nvtx == 0) return;
  int nadd =  gEVR.ProbVertexN();
  TString name;
  gSproc.MakeFileName(name,idset,"vtx.root",false);
  EdbDataProc::MakeVertexTree(*(gEVR.eVTX),name.Data());
}

void MakeScanCondBT(EdbScanCond &cond, TEnv &env)
{
  cond.SetSigma0( env.GetValue("emvertex.bt.Sigma0", "0.2 0.2 0.002 0.002" ) );
  cond.SetDegrad( env.GetValue("emvertex.bt.Degrad", 5. ) );
  cond.SetBins(3, 3, 3, 3);
  cond.SetPulsRamp0(  12., 18. );
  cond.SetPulsRamp04( 12., 18. );
  cond.SetChi2Max( 6.5 );
  cond.SetChi2PMax( 6.5 );
  cond.SetChi2Mode( 3 );
  cond.SetRadX0( 5810. );
  cond.SetName("SND_basetrack");
}

void AddCompatibleTracks(EdbPVRec &v_trk, EdbPVRec &v_vtx, TObjArray &v_out)
{
  int ntr  = v_trk.Ntracks();
  int nvtx = v_vtx.Nvtx();
  Log(1,"AddCompatibleTracks", "%d tracks, %d vertex", ntr,nvtx );
  for(int iv=0; iv<nvtx; iv++)
  {
    bool flag1 = false;
    EdbVertex *v = v_vtx.GetVertex(iv);
    for(int it=0; it<ntr; it++) 
    {
      EdbTrackP *t = v_trk.GetTrack(it);
      if( IsCompatible(*v,*t) ) {
      flag1 = true;
      t->SetFlag(999999);
      v_vtx.AddTrack(t);
      }
      }
    if (!flag1) {
      v_out.Add(v);
    }
  }
}

bool IsCompatible(EdbVertex &v, EdbTrackP &t)
{
  EdbSegP ss;
  t.EstimatePositionAt(v.VZ(),ss);
  float dx=ss.X()-v.VX();
  float dy=ss.Y()-v.VY();
  float r2 = Sqrt(dx*dx+dy*dy);
  float dz = Abs(ss.DZ());
  if(r2<5&&dz<4000) { printf("r2=%.4f dz=%.2f\n",r2,ss.DZ()); return true;}
  return false;
}

//-----------------------------------------------------------------------------
void SetTracksErrors(TObjArray &tracks, EdbScanCond &cond, float p, float m)
{
  int n = tracks.GetEntries();
  Log(2,"SetTracksErrors","refit %d tracks with a new errors and p=%f m=%f",n,p,m);
  for(int i=0; i<n; i++) {
     EdbTrackP *t = (EdbTrackP*)tracks.At(i);
     int nseg = t->N();
     t->SetSegmentsP(p);
     t->SetM(m);
     for(int j=0; j<nseg; j++) {
       EdbSegP   *s = t->GetSegment(j);
       s->SetErrors0();
       cond.FillErrorsCov( s->TX(),s->TY(), s->COV() );
     }
     t->FitTrackKFS();
  }
}
