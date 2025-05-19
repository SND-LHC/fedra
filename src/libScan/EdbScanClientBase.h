#ifndef ROOT_EdbScanClientBase
#define ROOT_EdbScanClientBase

#include "TSocket.h"
#include "TString.h"
#include "EdbRun.h"
#include "EdbPattern.h"

#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? : (X) : (Y))

#define CMD_LEN 256
class EdbScanClientBase
{
public:
  enum ScanType{ stPred, stVolume};
	EdbScanClientBase(){};
	virtual ~EdbScanClientBase(){};
	
	//common part
  virtual int    InitializeSocket()=0;

  virtual int    ScanPreloadAreaS( int id1, int id2, int id3, int id4,
							float x1, float x2, float y1, float y2,	const char *fname, 
							float x1n, float x2n, float y1n, float y2n)=0;
  virtual int    ScanAreas(ScanType st, int id[4], EdbPattern &areas, EdbRun *run, const char* options="")=0;
//  virtual int    ScanAreaS( int id1, int id2, int id3, int id4,	
//														float x1, float y1, float x2, float y2,	const char *fname);
//  virtual int   ScanAreasAsync(int id[4], EdbPattern &areas, EdbRun &run, const char* options="");
  virtual int     AddRWC_(EdbRun* run, char* rwcname, int bAddRWD=true, const char* options="")=0;
	// "short" version of brickID to be passed into sysal command line via BernScanDriver
	virtual Short_t  ShortBrick(Int_t brick)=0;
	virtual int   ConvertAreas(int id[4], EdbPattern &areas, EdbRun &run, const char* options="")=0;

	
	//protocol-specific part
	virtual int    UnloadPlate()=0;
	virtual int    LoadPlate(int BRICK, int PLATE, const char *mapext, int nAttempts=1)=0;
  virtual void   SetParameter(const char* Object, const char* Parameter, const char* Value)=0;
  virtual void   SetClusterThresholds(int TOP, int BOT)=0;
  virtual void   SetOdysseyThresholds(int itop, int ibottom, int size, int TOP, int BOT)=0;
  virtual int    SetFragmentSize(int X, int Y)=0;
  virtual void   AsyncScanAreaS( int id1, int id2, int id3, int id4,	
						 float x1, float x2, float y1, float y2,	const char *fname)=0;
  virtual void   AsyncScanPreloadAreaS( int id1, int id2, int id3, int id4,	
								float x1, float x2, float y1, float y2,	const char *fname, 
								float x1n, float x2n, float y1n, float y2n)=0;
  virtual int    AsyncWaitForScanResult()=0;
	virtual void   Print()=0;
  ClassDef(EdbScanClientBase,1)  // remote scanning from Bern
};

//-------------------------------------------------------------------
class EdbScanClientCommon: public EdbScanClientBase
{
protected:
  TSocket*  eSock;          // socket for connection to scanning mashine
  char      eCMD[CMD_LEN];      // command line
  char      eMess[CMD_LEN];     // message line
  char      eMess1[CMD_LEN];    // message line

  float       eNXview, eNYview;   // the fragment size in view
  float       eXstep, eYstep;     // between views in the fragment - to calculate the fragment size

  int       eMAXSCANATTEMPTS;
  int       eMAXFAILS;          // the max number for predictions failed to scan before break
  int       ePORT;              // port number where ScanDriver expect connection (default 1777)
  TString   eServer;            // address of the scanning mashine

  TString eRawDirServer;     // directory path for raw data files visible from the Scan Server (i.e. "o:/MIC5/tmp_SYSAL")
  TString eRawDirClient;     // directory path for raw data files visible from processing comp (i.e. "./raw/MIC5/tmp_SYSAL")
  TString eProcTgtServer;     // path for raw data files visible from processing comp (i.e. "o:/MIC5/b123456/p012/123456.12.1.1000.raw.root")
  TString eProcPthServer;     // path for raw data files visible from processing comp (i.e. "o:/MIC5/")

  TString eServerCreatedRunName;

public: 
  EdbScanClientCommon();
  virtual ~EdbScanClientCommon();
	void CloseSocket();
	//common part
  virtual int    InitializeSocket();

  virtual int    ScanPreloadAreaS( int id1, int id2, int id3, int id4,
							float x1, float x2, float y1, float y2,	const char *fname, 
							float x1n, float x2n, float y1n, float y2n);
//  virtual int    ScanAreaS( int id1, int id2, int id3, int id4,	
//														float x1, float y1, float x2, float y2,	const char *fname);
//  virtual int   ScanAreasAsync(int id[4], EdbPattern &areas, EdbRun &run, const char* options="");
  virtual int     AddRWC_(EdbRun* run, char* rwcname, int bAddRWD=true, const char* options="");
	//static int     AddRWC__(EdbRun* run, char* rwcname, int bAddRWD=true, const char* options="");
	// "short" version of brickID to be passed into sysal command line via BernScanDriver
	virtual Short_t  ShortBrick(Int_t brick){ return brick%10000; };
	virtual int   ConvertAreas(int id[4], EdbPattern &areas, EdbRun &run, const char* options="");

	
	//protocol-specific part
  virtual int   ScanAreas(ScanType st, int id[4], EdbPattern &areas, EdbRun *run, const char* options="")=0;
	virtual int    UnloadPlate()=0;
	virtual int    LoadPlate(int BRICK, int PLATE, const char *mapext, int nAttempts=1)=0;
  virtual void   SetParameter(const char* Object, const char* Parameter, const char* Value)=0;
  virtual void   SetClusterThresholds(int TOP, int BOT)=0;
  virtual void   SetOdysseyThresholds(int itop, int ibottom, int size, int TOP, int BOT)=0;
  virtual int    SetFragmentSize(int X, int Y)=0;
  virtual void   AsyncScanAreaS( int id1, int id2, int id3, int id4,	
						 float x1, float x2, float y1, float y2,	const char *fname)=0;
  virtual void   AsyncScanPreloadAreaS( int id1, int id2, int id3, int id4,	
								float x1, float x2, float y1, float y2,	const char *fname, 
								float x1n, float x2n, float y1n, float y2n)=0;
  virtual int    AsyncWaitForScanResult()=0;
	virtual void   Print()=0;
  
  virtual bool ServerCreatesTarget()=0;

	const char * GetCmd(){return eCMD; };
	void SetCmd(const char* cmd_){
		int len = strlen(cmd_);
		if(len > CMD_LEN-1)  len = CMD_LEN-1;
		memcpy(eCMD, cmd_, len+1);//strcpy_s(eCMD, 256, cmd_); windows only :(
		eCMD[len]=0;
	};

	const char * GetMess(){return eMess; };
	void SetMess(const char* mess_){
		int len = strlen(mess_);
		if(len > CMD_LEN-1) len = CMD_LEN-1;
		memcpy(eMess, mess_, len+1);
		eMess[len]=0;
	};

	const char * GetMess1(){return eMess1; };
	void SetMess1(const char* mess_){
		int len = strlen(mess_);
		if(len > CMD_LEN-1) len = CMD_LEN-1;
		memcpy(eMess1, mess_, len+1);
		eMess1[len]=0;
	};

	float GetNXView(){return eNXview;};
	void SetNXview(float nxview_){eNXview = nxview_;};

	float GetNYView(){return eNYview;};
	void SetNYview(float nyview_){eNYview = nyview_;};

	float GetXstep(){return eXstep;};
	void SetXstep(float xstep_){eXstep = xstep_;};

	float GetYstep(){return eYstep;};
	void SetYstep(float ystep_){eYstep = ystep_;};

	int GetMaxAttempts(){return eMAXSCANATTEMPTS;};
	void SetMaxAttempts(int attempts_){eMAXSCANATTEMPTS = attempts_;};

	int GetMaxFails(){return eMAXFAILS;};
	void SetMaxFails(int fails_){eMAXFAILS = fails_;};

	int GetPort(){return ePORT;};
	void SetPort(int port_){ePORT = port_;};

	const char * GetServer(){return eServer.Data(); };
	void SetServer(const char* server_){
		eServer = server_;
	};

	const char * GetRawDirServer(){return eRawDirServer.Data(); };
	void SetRawDirServer(const char* rawdirserver_){
		eRawDirServer = rawdirserver_;
	};

	const char * GetRawDirClient(){return eRawDirClient.Data(); };
	void SetRawDirClient(const char* rawdirclient_){
		eRawDirClient = rawdirclient_;
	};

	const char * GetProcTgtServer(){return eProcTgtServer.Data(); };
	void SetProcTgtServer(const char* prctgtserver_){
		eProcTgtServer = prctgtserver_;
	};

	const char * GetProcPthServer(){return eProcPthServer.Data(); };
	void SetProcPthServer(const char* prcpthserver_){
		eProcPthServer = prcpthserver_;
	};

	const char * GetServerCreatedRunName(){return eServerCreatedRunName.Data();};

  ClassDef(EdbScanClientBase,1)  // remote scanning from Bern
};

//-------------------------------------------------------------------

#endif /* ROOT_EdbScanClientBase */
