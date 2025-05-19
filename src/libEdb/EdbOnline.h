#ifndef ROOT_EdbOnline
#define ROOT_EdbOnline
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbOnline                                                            //
//                                                                      //
// head Edb class for online scanning                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"

#ifndef ROOT_EdbRun
#include "EdbRun.h"
#endif


//______________________________________________________________________________
class EdbOnline : TROOT 
{

 private:

  EdbRun        *eRun;                   // pointer to the run object

  TString       eOnlineDirectory;        // directory for online parameters
  TString       eDataDirectory;          // directory for scanning data

  Int_t         eRunID;                  // current run id

 public:
  EdbOnline( const char *name="OnlineRoot", 
	     const char *title="Root object for online scanning" );
  virtual ~EdbOnline();

  EdbRun  *GetRun() const { return eRun; }

  void    Print( Option_t *opt="") const;

  void     SetDataDirectory(   const char *dir) { eDataDirectory=dir;  }
  void     SetOnlineDirectory( const char *dir) { eOnlineDirectory=dir;  }

  int     ReadParameters();
  int     ReadParametersDB();
  int     ReadParameters( const char *file );
  void    WriteParametersDB();

  int     GetRunID() { return eRunID; }

  void     Init();

  //  ClassDef(EdbOnline,1)  // Online application class
};

#endif /* ROOT_EdbOnline */
