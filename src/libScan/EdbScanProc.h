#ifndef ROOT_EdbScanProc
#define ROOT_EdbScanProc

#include "EdbRunAccess.h"
#include "EdbDataSet.h"
#include "EdbScanClient.h"
#include "EdbScanSet.h"

class EdbScanProc : public TNamed
{
 public:

  TString eProcDirClient;    // directory path for root data

public: 
  EdbScanProc();
  virtual ~EdbScanProc(){}

  bool    CheckProcDir(int id[4], bool create=true);
  void    MakeFileName(TString &s, int id[4], const char *suffix);
  void    MakeAffName(TString &s, int id1[4], int id2[4], const char *suffix="aff.par");
  bool    GetMap(int brick, TString &map);
  bool    AddParLine(const char *file, const char *line);
  bool    MakeInPar(int id[4], const char *option);
  int     CopyFile(int id1[4], int id2[4], const char *suffix, bool overwrite);
  int     CopyPar(int id1[4], int id2[4], bool overwrite=true) {return CopyFile(id1,id2,"par",overwrite);}
  int     CopyPred(int id1[4],int id2[4], bool overwrite=true) {return CopyFile(id1,id2,"pred.root",overwrite);}
  int     CopyAFFPar(int id1c[4],int id2c[4], int id1p[4], int id2p[4], bool overwrite=true);

  int     ReadPatTXT(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     WritePatTXT(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     ReadPatRoot(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     WritePatRoot(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     ReadPred(EdbPattern &pred, int id[4], int flag=-1) {return ReadPatRoot(pred,id,"pred.root",flag);}
  int     WritePred(EdbPattern &pred, int id[4], int flag=-1) {return WritePatRoot(pred,id,"pred.root",flag);}
  int     ReadFound(EdbPattern &pred, int id[4], int flag=-1) {return ReadPatRoot(pred,id,"found.root",flag);}
  int     WriteFound(EdbPattern &pred, int id[4], int flag=-1) {return WritePatRoot(pred,id,"found.root",flag);}

  EdbRun *InitRun(int id[4]);
  bool    FlashRawDir(EdbScanClient &scan, int id[4]);
  int     LoadPlate(EdbScanClient &scan, int id[4], int attempts=1);
  int     ScanAreas(EdbScanClient &scan, int id[4], int flag=-1, const char *opt="NOCLCLFRAMESUM");

  bool    InitPiece(EdbDataPiece &piece, int id[4]);
  int     ReadPiece(EdbDataPiece &piece, EdbPattern &pat);
  int     ReadPatCP(EdbPattern &pat, int id[4]);
  bool    ApplyAffZ(EdbPattern &pat,int id1[4],int id2[4]);
  bool    GetAffZ(EdbAffine2D &aff, float &z,int id1[4],int id2[4]);
  bool    SetAFFDZ(int id1[4], int id2[4], float dz);
  bool    SetAFF0(int id1[4], int id2[4]);

  int     ConvertAreas(EdbScanClient &scan, int id[4], int flag=-1, const char *opt="NOCLCLFRAMESUM");
  int     LinkRun(int id[4], int noUpdate=1);
  int     LinkRunAll(int id[4], int npre=3, int nfull=1);
  int     Align(int id1[4], int id2[4], const char *option="");
  int     AlignAll(int id1[4], int id2[4], int npre=1, int nfull=3, const char *opt="-z");
  bool    CorrectAffWithPred(int id1[4], int id2[4], const char *opt="-z", int patmin=6);
  bool    ProjectFound(int id1[4],int id2[4]);
  int     FindPredictions(EdbPattern &pred, int id[4], EdbPattern &found, int maxholes=3);
  int     FindPredictions(int id[4], int flag=-1, int maxholes=3);
  void    OptimizeScanPath(EdbPattern &pin, EdbPattern &pout,int brick);
  int     RemoveDublets(EdbPattern &pin, EdbPattern &pout,int brick);

  bool    InitRunAccess(EdbRunAccess &ra, int id[4]);
  //int     FindPredictionsRaw(EdbPattern &pred, EdbPattern &found, EdbRunAccess &ra);

  bool    AddAFFtoScanSet(EdbScanSet &sc, int id1[4], int id2[4]);
  bool    AddAFFtoScanSet(EdbScanSet &sc, int b1, int p1, int s1, int e1,int b2, int p2, int s2, int e2);

  void    PrepareVolumesPred(int id[4], EdbPattern &points, int before=5, int after=5, 
			     int pmin=1, int pmax=57, EdbScanSet *sc=0);

  int     TestAl(int id1[4], int id2[4]);

  void    LogPrint(int brick, const char *rout, const char *fmt, ...);
  void    Print();
  ClassDef(EdbScanProc,1)  // scanned data processing
};
#endif /* ROOT_EdbScanProc */