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
void AddCompatibleTracks(TEnv &env, EdbPVRec &v_trk, EdbPVRec &v_vtx, TObjArray &v_out, TObjArray &v_out2, TNtuple* outTree);
bool IsCompatible(TEnv &env, EdbVertex &v, EdbTrackP &t, float *r2, float *dz);
void Display( const char *dsname,  EdbVertexRec *evr, TEnv &env );
float distance(const EdbTrackP& t, const EdbVertex& v);

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
  env.SetValue("emvertex.trfit.r2max", 5. );
  env.SetValue("emvertex.trfit.dzmax", 4000. );
  env.SetValue("emvertex.trfit.impmax", 5. );

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
  bool        do_vtxrefit = false;
  
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
    else if(!strncmp(key, "-refit",6))
    {
      do_vtxrefit=true;
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

  TObjArray v_out;
  TObjArray v_out2;
  
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
      TNtuple *outTree = new TNtuple("tracks","Tree of matched tracks","chosen:n:vid:tid:nseg:npl:tx:ty:firstp:lastp:r2:dz");
      AddCompatibleTracks(env, *vtr, gAli , v_out, v_out2, outTree);  // assign to the vertices of gAli additional tracks from vtr if any
      EdbDataProc::MakeVertexTree(v_out,"flag0.vtx.root");
      EdbDataProc::MakeVertexTree(v_out2,"flag1.vtx.root");
      TFile *outFile = new TFile("found_tracks.root","RECREATE");
      outTree->Write();
      outFile->Write();
      outFile->Close();
      delete outTree;
      delete outFile;
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

void AddCompatibleTracks(TEnv &env, EdbPVRec &v_trk, EdbPVRec &v_vtx, TObjArray &v_out, TObjArray &v_out2, TNtuple* outTree)
{
  int ntr  = v_trk.Ntracks();
  int nvtx = v_vtx.Nvtx();
  Log(1,"AddCompatibleTracks", "%d tracks, %d vertex", ntr,nvtx );
  for(int iv=0; iv<nvtx; iv++)
  {
    bool flag1 = false;
    std::vector<int> trackids;
    EdbVertex *v = v_vtx.GetVertex(iv);
    Log(1,"AddCompatibleTracks","Looking for parent tracks of vtx: %i\n",v->ID());
    for(int i=0; i<v->N(); i++){
      EdbTrackP *t = (EdbTrackP*)v->GetTrack(i);
      int trid = t->ID();
      trackids.push_back(trid);
    }
    EdbTrackP *t_chosen = 0;
    float r2, dz, imp, imp2, r2max, dzmax, impmax, impmax2; r2max=dzmax=impmax=impmax2=1e9f;
    int founds=0;
    for(int it=0; it<ntr; it++) 
    {
      EdbTrackP *t = v_trk.GetTrack(it);
      int trid = t->ID();
      if (std::find(trackids.begin(), trackids.end(), trid)!=trackids.end()) continue; //Maybe here can be changed to EdbVertex::TrackInVertex(EdbTrackP *t)
      if( IsCompatible(env, *v,*t, &r2, &dz) ) {
        flag1 = true;
        t->SetFlag(999999);
        if( r2 < r2max ) { r2max=r2; dzmax=dz; impmax=imp; impmax2=imp2; t_chosen=t; }
        // v_vtx.AddTrack(t);
        founds++;
        outTree->Fill(0, 1, v->ID(), t->ID() ,t->N(), t->Npl(), t->TX(), t->TY(), t->GetSegmentFirst()->Plate(), t->GetSegmentLast()->Plate(), r2, dz, imp, imp2);
      }
    }
    if (flag1){
      v_vtx.AddTrack(t_chosen);
      outTree->Fill(1, founds, v->ID(), t_chosen->ID(), t_chosen->N(), t_chosen->Npl(), t_chosen->TX(), t_chosen->TY(), t_chosen->GetSegmentFirst()->Plate(), t_chosen->GetSegmentLast()->Plate(), r2max, dzmax, impmax, impmax2);
      Log(1,"AddCompatibleTracks","Closest track found at r2=%.4f dz=%.2f\n",r2max,dzmax);
      v_out2.Add(v);
    }
    else {
      v_out.Add(v);
    }
  }
}

bool IsCompatible(TEnv &env, EdbVertex &v, EdbTrackP &t, float *r2, float *dz, float *imp, float *imp2)
{
  EdbSegP ss;
  EdbSegP *tst = t.GetSegmentFirst();
  float tz = tst->Z();
  if (tz > v.VZ()) return false;   //only tracks starting upstream of the vtx
  t.EstimatePositionAt(v.VZ(),ss);
  float dx=ss.X()-v.VX();
  float dy=ss.Y()-v.VY();
  *r2 = Sqrt(dx*dx+dy*dy);
  *dz = ss.DZ();
  *imp = distance(t, v);
  *imp2 = distance(ss, v);
  float r2max      = env.GetValue("emvertex.trfit.r2max"        , 5. );
  float dzmax      = env.GetValue("emvertex.trfit.dzmax"        , 4000. );
  float impmax      = env.GetValue("emvertex.trfit.impmax"        , 5. );
  if(*r2<r2max&&*dz<Abs(dzmax)) { printf("r2=%.4f dz=%.2f\n",*r2,ss.DZ()); return true;}
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

float distance(const EdbTrackP& t, const EdbVertex& v) {
  float dx = v.VX() - t.X();
  float dy = v.VY() - t.Y();
  float dz = v.VZ() - t.Z();
  float tx = t.TX();
  float ty = t.TY();
  float nom = std::pow((dx*ty - dy*tx),2) + std::pow((dy - dz*ty),2) + std::pow((dz*tx - dx),2);
  float denom = std::pow((tx),2) + std::pow((ty),2) + 1.;
return sqrt(nom/denom);
}
