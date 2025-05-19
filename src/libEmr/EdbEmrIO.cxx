//-- Author :  Valeri Tioukov, Michele Pozzato 3/06/2008

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  Interface from Edb to EmuRec objects                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "EdbVertex.h"
#include "EdbEmrIO.h"
#include "EdbLog.h"

ClassImp(EdbEmrIO);
ClassImp(EdbEmrFileAccess)

//________________________________________________________________________________________
int EdbEmrFileAccess::FirstPlate()
{
  int size = eSides.GetSize(), i=0;
  EdbPattern *p=0;
  for(i=0; i<size; i++) {
    p = (EdbPattern*)eSides.At(i);
    if(p) break;
  }
  return PlateID(i);
}

//________________________________________________________________________________________
int EdbEmrFileAccess::LastPlate()
{
  int size = eSides.GetSize(), i=0;
  EdbPattern *p=0;
  for(i=size-1; i>=0; i--) {
    p = (EdbPattern*)eSides.At(i);
    if(p) break;
  }
  return PlateID(i);
}

//________________________________________________________________________________________
void EdbEmrFileAccess::PrintStat()
{
  int size = eSides.GetSize();
  EdbPattern *p=0;
  for(int i=0; i<size; i++) {
    p = (EdbPattern*)eSides.At(i);
    if(p) printf("%d views in plate %d side %d\n", p->N(), PlateID(i), SideID(i) );
  }
}

//________________________________________________________________________________________
EdbSegP *EdbEmrFileAccess::AddEntry(int plate, int side, int zone, int view, int entry, float x, float y, float z)
{
  int iside = SideID(plate,side);
  if(iside<0||iside>1000)    {
    Log(1,"EdbEmrFileAccess::AddEntry","Error: iside = %d is out of range! (%d %d)",iside,0,1000);
    return 0;
  }
  if(iside+1>eSides.GetSize()) eSides.Expand(iside+1);
  if(!eSides.At(iside)) eSides.AddAt(new EdbPattern, iside);
  EdbPattern *p = (EdbPattern*)eSides.At(iside);
  EdbSegP s(entry, x, y, 0, 0,  0,0);
  s.SetAid(zone, view,side);
  s.SetZ(z);
  return p->AddSegment(s);
}

//________________________________________________________________________________________
int EdbEmrFileAccess::GetPlateEntries(int plate, int side, TArrayI &entries)
{
  int iside = SideID(plate,side);
  EdbPattern *p = (EdbPattern*)eSides.At(iside);  
  if(!p) return 0;
  int n = p->N();
  entries.Set(n);
  EdbSegP *s;
  for(int i=0; i<n; i++) {
    s = p->GetSegment(i);
    entries[i] = s->ID();
  }
  return n;
}

//________________________________________________________________________________________
void EdbEmrIO::Print()
{
  printf("EdbEmrIO ok\n");
}

//________________________________________________________________________________________
void EdbEmrIO::TrackEdb2Emr(EdbTrackP &tedb)
{
  // Input:  EdbTrackP tedb
  // Output: EmrTrack  temr
  tedb.Print();
}

//________________________________________________________________________________________
void EdbEmrIO::TrackEmr2Edb(EdbTrackP &tedb)
{
  // Input:  EmrTrack  temr
  // Output: EdbTrackP tedb
  tedb.Print();
}

//________________________________________________________________________________________
void EdbEmrIO::VertexEdb2Emr(EdbVertex &vedb)
{
  // Input:  EdbVertex vedb
  // Output: EmrVertex vemr
  vedb.Print();
}

//________________________________________________________________________________________
void EdbEmrIO::VertexEmr2Edb(EdbVertex &vedb)
{
  // Input:  EmrVertex vemr
  // Output: EdbVertex vedb
  vedb.Print();
}
