//-- Author :  Valeri Tioukov   21/03/2006
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbScanSet                                                           //
//                                                                      //
// Scanning data assembler                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "EdbLog.h"
#include "EdbScanSet.h"
#include <iostream>

using namespace std;
using namespace TMath;

ClassImp(EdbID)
ClassImp(EdbScanSet)

//----------------------------------------------------------------
EdbScanSet::EdbScanSet(EdbID id)
{
  eID = id;
  eB.SetID(id.eBrick);
  eReferencePlate=0;
}

//----------------------------------------------------------------
EdbScanSet::EdbScanSet(int id)
{
  eB.SetID(id);
  eID.eBrick=id;
  eReferencePlate=0;
}

//----------------------------------------------------------------
void EdbScanSet::UpdateIDS(int brick, int ma, int mi)
{
  // update all identifiers in dataset with (brik,ma,mi) keeping plateids unchanged
  
  eID.Set(brick, eID.ePlate, ma, mi);
  TIter next( &(eIDS) );
  TObject *obj=0;
  while((obj = next())) {
    EdbID *id=(EdbID*)obj;
    id->Set(brick, id->ePlate, ma, mi);
  }
}

//----------------------------------------------------------------
void EdbScanSet::Copy(EdbScanSet &sc)
{
  // copy from ss to this
  eB.Copy(sc.eB);
  eID = sc.eID;

  int npc = sc.ePC.GetEntriesFast();
  for(int i=0; i<npc; i++) {
    EdbPlateP *p = new EdbPlateP();
    p->Copy( *((EdbPlateP *)sc.ePC.At(i)) );
    ePC.Add(p);
  }
 
  TIter next( &(sc.eIDS) );
  TObject *obj=0;
  while((obj = next()))
    eIDS.Add( new EdbID( *((EdbID*)obj) ) );

  MakePIDList();
  eReferencePlate = sc.eReferencePlate;
}

//----------------------------------------------------------------
void EdbScanSet::MakeNominalSet( EdbID id, Int_t from_plate, Int_t to_plate, 
				 Float_t z0, Float_t dz, float shr, float dzbase, float dzem)
{
  eID=id;
  Log(2,"EdbScanSet::MakeNominalSet","%s",eID.AsString());
  eReferencePlate = from_plate;
  eB.SetID(id.eBrick);
  Int_t step= (from_plate>to_plate)?-1:1;
  Int_t p=from_plate;
  Float_t z=z0;
  do {
    EdbPlateP *plate = new EdbPlateP();
    plate->SetID(p);
    plate->SetZ(z);
    plate->SetPlateLayout(dzbase,dzem,dzem);

//     plate->SetZlayer(z,z-dzbase/2-dzem,z+dzbase/2+dzem);
//     plate->GetLayer(0)->SetZlayer(0,-dzbase/2,dzbase/2);
//     plate->GetLayer(1)->SetZlayer(dzbase/2,dzbase/2,dzbase/2+dzem);
//     plate->GetLayer(2)->SetZlayer(-dzbase/2,-dzbase/2-dzem,-dzbase/2);

    plate->GetLayer(1)->SetShrinkage(shr);
    plate->GetLayer(2)->SetShrinkage(shr);
    eB.AddPlate(plate);
    eIDS.Add(new EdbID(id.eBrick,p,id.eMajor,id.eMinor));
    p += step;
    z += dz;
  } while( p!=(to_plate+step) );
  MakePIDList();
  if(gEDBDEBUGLEVEL>2) Print();
}

//----------------------------------------------------------------
Int_t EdbScanSet::MakeParFiles(Int_t piece, const char *dir)
{
  // prepare par files for the dataset for "recset-like" processing

  TString name;
  Log(2,"MakeParFiles","for brick %d with %d plates:\n", eB.ID(), eB.Npl());
  EdbPlateP *p=0;
  for(Int_t i=0; i<eB.Npl(); i++) {
    p = eB.GetPlate(i);
    name.Form("%s/%2.2d_%3.3d.par",dir,p->ID(),piece);
    Log(3,"MakeParFiles","%s\n",name.Data());
    FILE *f = fopen(name.Data(),"w");
    EdbAffine2D *a=p->GetAffineXY();
    fprintf(f,"INCLUDE default.par\n");
    fprintf(f,"ZLAYER 0  %15.3f 0. 0.\n", p->Z());
    fprintf(f,"AFFXY  0  %15.6f %9.6f %9.6f %9.6f %15.6f %15.6f\n",
	    a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2());
    fclose(f);
  }
  return eB.Npl();
}

//----------------------------------------------------------------
Int_t EdbScanSet::AssembleBrickFromPC()
{
  // input: ePC - couples of layers represented as "plates"
  // output: brick with a first plate as a reference one

  ePID.Delete();
  eB.Clear();
  EdbPlateP *plate=0;

  Float_t dz0=214,dz1=45,dz2=45;  //TODO!
  EdbAffine2D aff;
  Int_t npc = ePC.GetEntriesFast();      if(npc<1) return 0;
  EdbPlateP *pc=0;                    //the couple of plates
  Float_t z;
  printf("\nAssembleBrickFromPC\n");
  for(Int_t i=-1; i<npc; i++) {
    plate = new EdbPlateP();

    if(i==-1) {
      pc = (EdbPlateP *)(ePC.At(0));
      plate->SetID(pc->GetLayer(1)->ID());
      z=0;
    }
    else {
      pc = (EdbPlateP *)(ePC.At(i));
      plate->SetID(pc->GetLayer(2)->ID());
      z = eB.GetPlate(eB.Npl()-1)->Z() - pc->GetLayer(1)->Z();
      aff.Reset();
      aff.Transform( pc->GetLayer(1)->GetAffineXY() );
      aff.Invert();
      aff.Transform( eB.GetPlate(eB.Npl()-1)->GetAffineXY() );
      plate->GetAffineXY()->Transform(&aff);
    }
    plate->SetZlayer(z, z - dz0/2 + dz1, z+dz0/2+dz2);                
    plate->GetLayer(0)->SetZlayer(0,-dz0/2,dz0/2);       // internal plate coord
    plate->GetLayer(2)->SetZlayer(-dz0/2,-dz0/2-dz2,-dz0/2);
    plate->GetLayer(1)->SetZlayer( dz0/2, dz0/2, dz0/2+dz1);

    eB.AddPlate(plate);    
  }
  MakePIDList();

  eReferencePlate = eB.GetPlate(0)->ID();

  return eB.Npl();
}

//----------------------------------------------------------------
void EdbScanSet::MakePIDList()
{
  ePID.Delete();
  Long_t v[2];
  for(Int_t i=0; i<eB.Npl(); i++) {
    v[0]= eB.GetPlate(i)->ID();
    v[1]= i;
    ePID.Add(2,v);
  }
  ePID.Sort();
}

//----------------------------------------------------------------
Bool_t EdbScanSet::SetAsReferencePlate(Int_t pid)
{
  EdbPlateP *p = GetPlate(pid);
  if(!p) return kFALSE;
  eReferencePlate=pid;
  Float_t  z = p->Z();
  EdbAffine2D aff( *(p->GetAffineXY()) );
  aff.Invert();
  TransformBrick(aff);
  ShiftBrickZ(-z);
  return kTRUE;
}

//----------------------------------------------------------------
Int_t EdbScanSet::TransformSidesIntoBrickRS()
{
  // to put layers 1&2 of each plate into brick reference system - to have absolute transformations for microtracks
  Int_t npl=eB.Npl();
  if(npl>0)
    for(Int_t i=0; i<npl; i++) 
       eB.GetPlate(i)->TransformSidesIntoPlateRS();
  return npl;
}

//----------------------------------------------------------------
Int_t EdbScanSet::TransformBrick(EdbAffine2D aff)
{
  Int_t npl=eB.Npl();
  if(npl>0)
    for(Int_t i=0; i<npl; i++)  eB.GetPlate(i)->GetAffineXY()->Transform(aff);
  return npl;
}

//----------------------------------------------------------------
void EdbScanSet::TransformBrick(EdbScanSet &ss)
{
  int n = eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id   = GetID(i);                  if(!id) continue;
    EdbPlateP *p  = GetPlate(id->ePlate);    if(!p)  continue;
    EdbPlateP *p1 = ss.GetPlate(id->ePlate); if(!p1) continue;
    p->GetAffineXY()->Transform( p1->GetAffineXY() );
    p->SetZlayer(p1->Z(),p1->Zmin(),p1->Zmax());
    // ToDo: add also AFFTXTY?
  }
}

//----------------------------------------------------------------
Int_t EdbScanSet::ShiftBrickZ(Float_t z)
{
  Int_t npl=eB.Npl();
  if(npl>0)  for(Int_t i=0; i<npl; i++) eB.GetPlate(i)->ShiftZ(z);
  return npl;
}

//----------------------------------------------------------------
Bool_t EdbScanSet::GetAffP2P(Int_t p1, Int_t p2, EdbAffine2D &aff)
{
  // input:  p1 - id of the "from-plate"
  //         p2 - id of the "to-plate"
  // output: aff - the transformation - applied to "from" segment became "to" segment
  Int_t ip1,ip2;
  if(ePID.Find(p1)) ip1 = ePID.Find(p1)->At(0)->Value();  else return false;
  if(ePID.Find(p2)) ip2 = ePID.Find(p2)->At(0)->Value();  else return false;
  aff.Reset();
  aff.Transform( eB.GetPlate(ip2)->GetAffineXY() );
  aff.Invert();
  aff.Transform( eB.GetPlate(ip1)->GetAffineXY() );
  return true;
}

//----------------------------------------------------------------
Float_t EdbScanSet::GetDZP2P(Int_t p1, Int_t p2)
{
  // input:  p1 - id of the "from-plate"
  //         p2 - id of the "to-plate"
  // return: dz - segment of from-plate propagated to dz become in to-plate 
  Int_t ip1,ip2;
  if(ePID.Find(p1)) ip1 = ePID.Find(p1)->At(0)->Value();  else return false;
  if(ePID.Find(p2)) ip2 = ePID.Find(p2)->At(0)->Value();  else return false;
  return  eB.GetPlate(ip2)->Z() - eB.GetPlate(ip1)->Z();
}

//----------------------------------------------------------------
void  EdbScanSet::UpdateBrickWithP2P(EdbLayer &la, int plate1, int plate2)
{
  // update the brick geometry with aff & dz defined in layer (plate1->plate2)
  // in this version assumed that plate1 and plate2 are consequtive in the brick
  // Note eReferencePlate should be correctly defined

  int pfrom = TMath::Abs(plate1-eReferencePlate)<TMath::Abs(plate2-eReferencePlate)? plate2 : plate1;
  int step  = ((pfrom-eReferencePlate)>0)? 1: -1;  //the corrections to be propagated in the opposite direction from ref
  
  Log(3,"UpdateBrickWithP2P","from,step: %d %d   ref: %d  Z(%d)=%f",pfrom,step,eReferencePlate,pfrom,GetPlate(pfrom)->Z());

  float dz           = (GetDZP2P(plate1,plate2) - la.Zcorr());

  EdbAffine2D *affxy   = la.GetAffineXY();
  if(pfrom==plate2) { 
    dz *= -1.;
    affxy->Invert();
  }
  //Log(1,"UpdateBrickWithP2P","%d->%d  dz= %f",plate1,plate2, dz);
  int ipl=pfrom;
  for(int i=0; i<eB.Npl(); i++) {
    EdbPlateP *p = GetPlate(ipl);    if(!p) continue;
    p->SetZlayer( p->Z() + dz, p->Zmin(), p->Zmax() );
    p->GetAffineXY()->Transform(affxy);
    ipl+=step;
  }
  
  EdbPlateP *plate = GetPlate(plate1);
  plate->SetShrinkage( GetPlate(plate1)->Shr() * la.Shr() );
  //plate->GetAffineTXTY()->Print();
  plate->GetAffineTXTY()->Transform( la.GetAffineTXTY() );
  //plate->GetAffineTXTY()->Print();
  
  plate->ApplyCorrectionsLocal(la.Map());
  
  //Log( 1, "UpdateBrickWithP2P","after: Z(%d)=%f", pfrom, GetPlate(pfrom)->Z() );

}

//----------------------------------------------------------------
int EdbScanSet::GetPlateIDNearestToZ(float z, int dir)
{
  int   plate=-1;
  float dz0=kMaxInt;
  for(Int_t i=0; i<eIDS.GetSize(); i++)
  {
    EdbID *id = (EdbID*)eIDS.At(i);
    EdbPlateP *p = GetPlate(id->ePlate);
    float dz = (p->Z()-z);
    if(dir==0||dz==0) if(         Abs(dz)  < dz0) { dz0 = dz; plate = id->ePlate; }
    else if(dir<0)    if( dz>0 && Abs(dz)  < dz0) { dz0 = dz; plate = id->ePlate; }
    else if(dir>0)    if( dz<0 && Abs(dz)  < dz0) { dz0 = dz; plate = id->ePlate; }
  }
  return plate;
}

//----------------------------------------------------------------
void EdbScanSet::Print()
{
  EdbPlateP *p=0;
  printf("EdbScanSet: %s\n", eID.AsString() );
//   int npc = ePC.GetEntriesFast();
//   printf("%d couples\n",npc);
//   for(Int_t i=0; i<npc; i++) {
//     p = (EdbPlateP *)(ePC.At(i));
//     p->GetLayer(1)->Print();
//     p->GetLayer(2)->Print();
//  }
  printf("Brick %d with %d plates:\n", eB.ID(), eB.Npl());
  for(Int_t i=0; i<eB.Npl(); i++) { 
    p = eB.GetPlate(i);
    if(!p) continue;
    EdbAffine2D *a=p->GetAffineXY();
    if(!a) continue;
    printf("%3d  %12.2f       %9.6f %9.6f %9.6f %9.6f %15.6f %15.6f\n",
	   p->ID(), p->Z(), a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2());
  }
  printf("for this brick %d identifiers are defined:\n", eIDS.GetSize());
  for(Int_t i=0; i<eIDS.GetSize(); i++)  
  {
    EdbID *id = (EdbID*)eIDS.At(i);
    TString str = id->AsString();
    EdbPlateP *p = GetPlate(id->ePlate);
    if(p) str += Form("  shr1/base/shr2 = %.3f / %.1f / %.3f", p->GetLayer(1)->Shr(),p->GetLayer(0)->DZ(), p->GetLayer(2)->Shr());
    printf("%s\n",str.Data() );
  }
}

//----------------------------------------------------------------
Int_t EdbScanSet::ReadIDS(const char *file)
{
  char line[256];
  FILE *f = fopen (file, "rt");
  Int_t b,p,mi,ma;
  Int_t cnt=0;
  while(fgets(line, 256, f) != NULL)
    {
      if( sscanf(line, "%d, %d, %d, %d", &b, &p, &mi, &ma) != 4 )  break;
      eIDS.Add(new EdbID(b,p,mi,ma));
      cnt++;
    }
  fclose(f);
  return cnt;
}

//----------------------------------------------------------------
Int_t EdbScanSet::WriteIDS(const char *file)
{
  FILE *f = 0;
  if(file) f = fopen(file, "w");
  EdbID *id;
  for( Int_t i=0; i<eIDS.GetSize(); i++ ) {
    id = (EdbID *)eIDS.At(i);
    if(file) fprintf(f,"%d.%d.%d.%d\n", id->eBrick, id->ePlate, id->eMajor, id->eMinor );
    else      printf("%d.%d.%d.%d\n", id->eBrick, id->ePlate, id->eMajor, id->eMinor );
  }

  if(file) fclose(f);
  return eIDS.GetSize();
}

//----------------------------------------------------------------
void EdbScanSet::ReadGeom(const char *file)
{
  //todo

  char line[256];
  FILE *f = fopen (file, "r"); 
  if(!f) {Log(1,"EdbScanSet::ReadGeom","Can not open file %s", file); return; }

  int b,p,mi,ma;
  int n=0;

  if(fgets(line, 256, f) != NULL) if( sscanf(line, "%d", &n) != 1 ) return;
  printf("%d\n",n);
  if(!n) return;

  ePID.Delete();
  eB.Clear();

  EdbID *id1=0, *id2=0;
  float dz=0;
  for(int i=0; i<n; i++) {
    if(i==0) {
      if( fgets(line, 256, f) == NULL )  return;
      if( sscanf(line, "%d.%d.%d.%d", &b, &p, &mi, &ma) != 4 )  break;
      id1 = new EdbID(b,p,mi,ma);
      if(n==1) eIDS.Add(id1);
    }
    if( n>1 && i<n-1) {
      if( fgets(line, 256, f) == NULL )  return;
      if( sscanf(line, "%f", &dz) != 1 ) break;
      if( fgets(line, 256, f) == NULL )  return;
      if( sscanf(line, "%d.%d.%d.%d", &b, &p, &mi, &ma) != 4 )  break;
      id2 = new EdbID(b,p,mi,ma);
      eIDS.Add(id2);
    }
  }
  fclose(f);
}

//----------------------------------------------------------------
void EdbScanSet::WriteGeom(const char *file)
{
  FILE *f = 0;
  if(file) f = fopen(file, "w");
  if(!f) {Log(1,"EdbScanSet::WriteGeom","Can not open file %s", file); return; }
  EdbID *id1=0, *id2=0;
  int n = eIDS.GetSize();
  fprintf(f,"%d\n", n);
  for( Int_t i=0; i<eIDS.GetSize(); i++ ) {
    id1 = (EdbID *)eIDS.At(i);
    if(i==0) fprintf(f,"%d.%d.%d.%d\n", id1->eBrick, id1->ePlate, id1->eMajor, id1->eMinor );
    if( n>1 && i<n-1) id2 = (EdbID *)eIDS.At(i+1);    else continue;
    fprintf(f,"%f\n", GetDZP2P(id1->ePlate, id2->ePlate));
    fprintf(f,"%d.%d.%d.%d\n", id2->eBrick, id2->ePlate, id2->eMajor, id2->eMinor );
  }

  if(file) fclose(f);
}

//----------------------------------------------------------------
EdbID *EdbScanSet::FindNextPlateID(Int_t plate, Bool_t positive_direction)
{
  // find next valid identifier in respect to plate in positive or negative direction (+-1)
  Log(2,"FindNextPlateID","plate = %d  direction = %d", plate, positive_direction);
  EdbID *id = 0;
  Int_t plmin=kMaxInt, plmax=kMinInt;
  for (Int_t i = 0; i < eIDS.GetSize(); i++) {
    id = (EdbID *)eIDS.At(i);
    if (id->ePlate>plmax) plmax=id->ePlate;
    if (id->ePlate<plmin) plmin=id->ePlate;
  }
  int step = positive_direction? 1 : -1;
  int first = plate +step;
  int last = positive_direction? plmax:plmin;   last += step;
  printf("first=%d, last = %d step = %d\n",first,last, step);
  if( positive_direction &&  first>=last)  return 0;
  if( !positive_direction && first<=last)  return 0;
  for(int i=first; i!=last; i+=step) {
    id = FindPlateID(i);
    if(id) return id;
  }
  return 0;
}

//----------------------------------------------------------------
EdbID *EdbScanSet::FindPlateID(Int_t plate)
{
  EdbID *id;
  for (Int_t i = 0; i < eIDS.GetSize(); i++) {
    id = (EdbID *)eIDS.At(i);
    if (id->ePlate == plate) return id;
  }
  return 0;
}

//----------------------------------------------------------------------------
Bool_t EdbScanSet::AddID(EdbID *id, Int_t step)
{
  if (FindPlateID(id->GetPlate())) return kFALSE;
  for (EdbID *it = (EdbID*)eIDS.First(); it ; it = (EdbID*)eIDS.After(it)) {
    if (step < 0 && id->GetPlate() > it->GetPlate()) {
      eIDS.AddBefore(it, id);
      return kTRUE;
    }
    if (step > 0 && id->GetPlate() < it->GetPlate()) {
      eIDS.AddBefore(it, id);
      return kTRUE;
    }
  }
  eIDS.Add(id);
  return kTRUE;
}

//----------------------------------------------------------------------------
void EdbScanSet::RemovePlate(int pid)
{
  // Remove both from eB and from eIDS all data corresponding to plate id
  // note that gaps and affine transformations may become 
  // inconsistent and had to be updated after this operation
  eB.RemovePlate(pid);
  EdbID   *id = FindPlateID(pid);
  if(id) {
    eIDS.Remove(id);
    MakePIDList();
  }
}

//----------------------------------------------------------------------------
void EdbScanSet::UpdateGap(float dz, int pid1, int pid2)
{
  // update gap between plates pid1 and pid2. All updates will be 
  // done in a way to have the reference plate unchanged
  EdbLayer la; la.SetZcorr(dz);
  UpdateBrickWithP2P(la,pid1,pid2);
}
