//-- Author :  Valeri Tioukov   27/01/2011
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbCouplesTree                                                       //
//                                                                      //
// service class for the couples tree access                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "EdbLog.h"
#include "EdbMath.h"
#include "EdbBrick.h"
#include "EdbPattern.h"
#include "EdbCouplesTree.h"

ClassImp(EdbCouplesTree)

//---------------------------------------------------------------------
EdbCouplesTree::EdbCouplesTree()
{
  eS1=0;
  eS2=0;
  eS=0;
  eCP=0;
  eXv=eYv=0;
  ePid1=ePid2=0;
  eApplyCorrections=0;
  eCut="1";
  eEraseMask=0;
  eAcceptMask=0;
  eTree=0;
}

//---------------------------------------------------------------------
EdbCouplesTree::~EdbCouplesTree()
{
// most of the members pointers can be the property of other objects
// If the file should be closed - use Close()
  //SafeDelete(eTree);
  //SafeDelete(eS);
  //SafeDelete(eS1);
  //SafeDelete(eS2);
 //SafeDelete(eEraseMask);
}

//---------------------------------------------------------------------
void EdbCouplesTree::Print()
{
  Log(1,"EdbCouplesTree::Print","s1 s2 s cp:");
  if(eS1) eS1->PrintNice();
  if(eS2) eS2->PrintNice();
  if(eS)  eS->PrintNice();
  if(eCP) eCP->Print();
}

///______________________________________________________________________________
void EdbCouplesTree::Close()
{
  if (eTree) {
    TFile *f = (TFile*)(eTree->GetDirectory()->GetFile());
    if (f) {
      //f->Print();
      if(f->IsWritable()) eTree->AutoSave();
      f->Close();
      //SafeDelete(eTree);
      //SafeDelete(f);
    }
     //SafeDelete(eTree);
     //eTree=NULL;
  }
}

//---------------------------------------------------------------------
bool EdbCouplesTree::InitCouplesTree(const char *name, const char *fname, Option_t *mode)
{
  Log(3,"EdbCouplesTree::InitCouplesTree","%s at file %s as %s",name,fname,mode);
  if(fname) {                                      // if fname==0 - assumed that file is already opened 
    //TFile *f = new TFile(fname, mode);
    //printf("\n ROOT dir before:\n"); gDirectory->ls();  printf("\n");
    TFile *f = TFile::Open(fname, mode);
    //printf("\n ROOT dir after:\n"); gDirectory->ls();  printf("\n");
    if (!f->IsOpen()) return 0;
    f->cd();
  }
  
  if( (strcmp(mode,"READ")==0)) {
    eTree = (TTree*)gDirectory->Get(name);
    if(!eTree) {Log(1,"EdbCouplesTree::InitCouplesTree","ERROR: can not open tree %s",name); return 0;}
    
    eTree->SetBranchAddress("pid1",  &ePid1);
    eTree->SetBranchAddress("pid2",  &ePid2);
    eTree->SetBranchAddress("xv"  ,  &eXv);
    eTree->SetBranchAddress("yv"  ,  &eYv);
    eTree->SetBranchAddress("cp"  ,  &eCP);
    eTree->SetBranchAddress("s1." ,  &eS1);
    eTree->SetBranchAddress("s2." ,  &eS2);
    eTree->SetBranchAddress("s."  ,  &eS);
    return 1;
    if( (strcmp(mode,"UPDATE")==0) ) {
      eTree = (TTree*)gDirectory->Get(name);
      if(eTree) {
	eTree->SetBranchAddress("pid1",  &ePid1);
	eTree->SetBranchAddress("pid2",  &ePid2);
	eTree->SetBranchAddress("xv"  ,  &eXv);
	eTree->SetBranchAddress("yv"  ,  &eYv);
	eTree->SetBranchAddress("cp"  ,  &eCP);
	eTree->SetBranchAddress("s1." ,  &eS1);
	eTree->SetBranchAddress("s2." ,  &eS2);
	eTree->SetBranchAddress("s."  ,  &eS);
	return 1;
      }
    }
  }
  
  Log(2,"EdbCouplesTree::InitCouplesTree","create new tree %s in %s",name, fname);
  eTree = new TTree(name, "couples tree (s1:s2:s)");
  eTree->SetMaxTreeSize(15000000000LL);   //set 15 Gb file size limit)
  
  eTree->Branch("pid1",  &ePid1,"pid1/I");
  eTree->Branch("pid2",  &ePid2,"pid2/I");
  eTree->Branch("xv"  ,  &eXv,"xv/F");
  eTree->Branch("yv"  ,  &eYv,"yv/F");
  eTree->Branch("cp"  ,  &eCP,32000,99);
  eTree->Branch("s1." ,  &eS1,32000,99);
  eTree->Branch("s2." ,  &eS2,32000,99);
  eTree->Branch("s."  ,  &eS,32000,99);
  
  return 1;
}

  /*
//---------------------------------------------------------------------
bool EdbCouplesTree::InitTreeNew(const char *name)
{
  // assumed that the correct file is already opened
  
  eTree = new TTree(name,"couples tree (s1:s2:s)");
  
  eTree->Branch("cp",   "EdbSegCouple",&eCP,32000,99);
  eTree->Branch("s1.",  "EdbSegP",&eS1,32000,99);
  eTree->Branch("s2.",  "EdbSegP",&eS2,32000,99);
  eTree->Branch("s." ,  "EdbSegP",&eS ,32000,99);
  eTree->Branch("dz1",  &eDZ1,"dz1/F");
  eTree->Branch("dz2",  &eDZ2,"dz2/F");
  eTree->Branch("flag", &eFlag,"flag/I");
  eTree->Branch("chi2", &eChi2,"chi2/F");

  eTree->SetAlias("x1","s1.eX+dz1*s1.eTX");
  eTree->SetAlias("x2","s2.eX+dz2*s2.eTX");
  eTree->SetAlias("y1","s1.eY+dz1*s1.eTY");
  eTree->SetAlias("y2","s2.eY+dz2*s2.eTY");
  
eTree->SetAlias( "tx1c" , Form("s1.eTX/(%f)+(%f)",ePC1.eShr,ePC1.eDTX) );
  eTree->SetAlias( "ty1c" , Form("s1.eTY/(%f)+(%f)",ePC1.eShr,ePC1.eDTY) );
  eTree->SetAlias( "tx2c" , Form("s2.eTX/(%f)+(%f)",ePC2.eShr,ePC2.eDTX) );
  eTree->SetAlias( "ty2c" , Form("s2.eTY/(%f)+(%f)",ePC2.eShr,ePC2.eDTY) );

  eTree->SetAlias("x1c" ,  Form("(%f) + s1.eX + tx1c*(%f)",ePC1.eDX,ePC1.eDZ) );
  eTree->SetAlias("y1c" ,  Form("(%f) + s1.eY + ty1c*(%f)",ePC1.eDY,ePC1.eDZ) );
  eTree->SetAlias("x2c" ,  Form("(%f) + s2.eX + tx2c*(%f)",ePC2.eDX,ePC2.eDZ) );
  eTree->SetAlias("y2c" ,  Form("(%f) + s2.eY + ty2c*(%f)",ePC2.eDY,ePC2.eDZ) );

  eTree->SetAlias("dx" , "x2c-x1c" );
  eTree->SetAlias("dy" , "y2c-y1c" );
  eTree->SetAlias("dtx" , "tx2c-tx1c" );
  eTree->SetAlias("dty" , "ty2c-ty1c" );

  eTree->SetAlias("dr" , "sqrt(dx*dx+dy*dy)" );
  eTree->SetAlias("dt" , "sqrt(dtx*dtx+dty*dty)" );
  return 1;
}
*/

//---------------------------------------------------------------------
bool EdbCouplesTree::WriteTree()
{
  if(!eTree)    return 0;
  eTree->Write();
  SafeDelete(eTree);
  return 1;
}

//---------------------------------------------------------------------
Int_t EdbCouplesTree::Fill( EdbSegP *s1, EdbSegP *s2, EdbSegP *s, EdbSegCouple *cp, float xv, float yv, int pid1, int pid2)
{
  eS1=s1;    eS2=s2;  eS=s; eCP=cp; 
  eXv=xv; eYv=yv; ePid1=pid1; ePid2=pid2;
   //eDZ1 = dz1; eDZ2 = dz2; eFlag=flag;
  return Fill();
}

//---------------------------------------------------------------------
Int_t EdbCouplesTree::Fill()
{
  return eTree->Fill();
}

///______________________________________________________________________________
int   EdbCouplesTree::GetEntry(int i)
{
  int size = eTree->GetEntry(i);
  if(eApplyCorrections)  ApplyCorrections();
  return size;
}

///______________________________________________________________________________
void   EdbCouplesTree::ApplyCorrections()
{
}

///______________________________________________________________________________
int EdbCouplesTree::GetCPData( TObjArray &cparr )
{
 //read segments into array of EdbSegCouple objects
  TEventList *lst = InitCutList();
  int nlst =lst->GetN();
  int nseg = 0, nerase=0;
  int entr=0;
  for(int i=0; i<nlst; i++ ) {
    entr = lst->GetEntry(i);

    if(eEraseMask) if(eEraseMask->At(entr)) {
      Log(3,"EdbCouplesTree::GetCPData","erase %d\n",entr);
      nerase++;
      continue;
    }
  
    GetEntry(entr);
    EdbSegCouple *cp = new EdbSegCouple(*eCP);
    cp->eS = new EdbSegP(*eS);
    cp->eS1 = new EdbSegP(*eS1);
    cp->eS2 = new EdbSegP(*eS2);
    cparr.Add( cp );

    nseg++;
  }
   
  SafeDelete(lst);

  Log(2,"EdbCouplesTree::GetCPData","select %d of %d segments by cut %s ; %d erased by mask",nlst, (int)(eTree->GetEntries()), eCut.GetTitle(), nerase);

  return nseg;
}

///______________________________________________________________________________
void EdbCouplesTree::RankCouples(TObjArray &cparr)
{
  // rank couples if in cparr the same segment's pointer appare in several combinations
  EdbSegCouple::SetSortFlag(0);    // sort by CHI2P
  cparr.UnSort();
  cparr.Sort();

  int ncp = cparr.GetEntries();

  for(int i=0; i<ncp; i++) {
    EdbSegCouple *sc = (EdbSegCouple*)(cparr.UncheckedAt(i));
    sc->eS1->SetFlag(0);
    sc->eS2->SetFlag(0);
  }

  for(int i=0; i<ncp; i++) {
    EdbSegCouple *sc = (EdbSegCouple*)(cparr.UncheckedAt(i));
    sc->eS1->SetFlag( sc->eS1->Flag()+1 );
    sc->eS2->SetFlag( sc->eS2->Flag()+1 );
    sc->SetN1(sc->eS1->Flag());
    sc->SetN2(sc->eS2->Flag());
  }

  for(int i=0; i<ncp; i++) {
    EdbSegCouple *sc = (EdbSegCouple*)(cparr.UncheckedAt(i));
    sc->SetN1tot(sc->eS1->Flag());
    sc->SetN2tot(sc->eS2->Flag());
  }

  Log(2,"EdbCouplesTree::RankCouples","%d couples ok", ncp );
}

///______________________________________________________________________________
TEventList *EdbCouplesTree::InitCutList()
{
  if(!eTree)  return  0;
  int nentr = (int)(eTree->GetEntries());  if(nentr<1) return 0;
  eTree->Draw(">>lst", eCut );
  TEventList *lst = (TEventList*)(gDirectory->GetList()->FindObject("lst"));
  if(!lst) {Log(1,"EdbCouplesTree::GetCPData","ERROR!: events list (lst) did not found! In couples tree %d entries",nentr); return 0;}
  return lst;
}

///______________________________________________________________________________
int EdbCouplesTree::PutCPData( TObjArray &cparr)
{
  int n = cparr.GetEntries();
  for(int i=0; i<n; i++)
  {
    EdbSegCouple *cp = (EdbSegCouple *)cparr.UncheckedAt(i);
    if(cp) {
      eCP = cp;
      eS  = cp->eS;
      eS1 = cp->eS1;
      eS2 = cp->eS2;
      Fill();
    }
  }
  return n;
}

///______________________________________________________________________________
int EdbCouplesTree::GetCPData( EdbPattern *pat, EdbPattern *p1, EdbPattern *p2, TIndex2 *trseg )
{
  TEventList *lst = InitCutList();
  int nlst =lst->GetN();

  int nseg = 0, nerase=0;
  pat->SetID(0);

  int entr=0;
  for(int i=0; i<nlst; i++ ) {
    entr = lst->GetEntry(i);

//     if(trseg) {           //exclude segments participating in tracks
//       if( (trseg->Find(ePlate*1000+ePiece,entr) >= 0) )  continue;
//     }

    if(eEraseMask) if(eEraseMask->At(entr)) {
      Log(3,"EdbCouplesTree::GetCPData","erase %d\n",entr);
      nerase++;
      continue;
    }

    GetEntry(entr);

    if(pat) {
      eS->SetZ( eS->Z() + pat->Z() );   /// TO CHECK !!!
      //ToDo: s->SetVid(ePlate*1000+ePiece,entr);
      eS->SetChi2(eCP->CHI2P());
      pat->AddSegment( *eS  );
      nseg++;
    }
    if(p1)  { 
      eS1->SetZ( eS1->Z() + pat->Z() );
      p1->AddSegment(  *eS1 ); 
      nseg++; 
    }
    if(p2)  { 
      eS2->SetZ( eS2->Z() + pat->Z() );
      p2->AddSegment(  *eS2 ); 
      nseg++; 
    }
  }

  SafeDelete(lst);

  Log(2,"EdbCouplesTree::GetCPData","select %d of %d segments by cut %s ; %d are erased by mask",nlst, (int)(eTree->GetEntries()), eCut.GetTitle(), nerase);

  return nseg;
}

///______________________________________________________________________________
int EdbCouplesTree::GetCPDataAcceptedMask( EdbPattern *pat )
{
  if(!eAcceptMask) return 0;
  int n = eAcceptMask->N();
  pat->SetID(0);
  int nseg=0;
  for(int i=0; i<n; i++ ) {
    if(eAcceptMask->At(i)>0) {
      GetEntry(i);
      if(pat) {
        eS->SetZ( eS->Z() + pat->Z() );   /// TO CHECK !!!
        //ToDo: s->SetVid(ePlate*1000+ePiece,entr);
        eS->SetChi2(eCP->CHI2P());
        pat->AddSegment( *eS  );
        nseg++;
      }
    }
  }

  Log(2,"EdbCouplesTree::GetCPDataAcceptedMask","select %d of %d segments",nseg, (int)(eTree->GetEntries()));

  return nseg;
}
