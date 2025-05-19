#ifndef ROOT_EdbRunHeader
#define ROOT_EdbRunHeader

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbRunHeader                                                         //
//                                                                      //
// header data assocoated with 1 run                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObjString.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TStopwatch.h"

#ifndef ROOT_EdbStage
#include "EdbStage.h"
#endif

#ifndef ROOT_EdbFiducial
#include "EdbFiducial.h"
#endif

//______________________________________________________________________________
class EdbPlate : public TNamed {

private:

  Int_t      eID;           // plate ID

  Float_t    eUp;           // Up side   | nominal emulsion thickness before shrinkage (on the beam)
  Float_t    eBase;         // base thickness
  Float_t    eDown;         // Down side |  (normally should be the same)
 
  Float_t    eShrinkageU;   // srinkage of the Up side (measured in the test run)
  Float_t    eShrinkageD;   // srinkage of the Down side (measured in the test run)

public:
  EdbPlate( int id=0, float u=0, float b=0, float d=0, float su=0, float sd=0 ): 
    eID(id), eUp(u), eBase(b), eDown(d), eShrinkageU(su), eShrinkageD(sd) {}
  virtual ~EdbPlate() {}

  void Print() const;
  void Set(int id, float u, float b, float d, float su, float sd);

  Int_t     GetID()         const { return eID;         }
  Float_t   GetBase()       const { return eBase;       }
  Float_t   GetUp()         const { return eUp;         }
  Float_t   GetDown()       const { return eDown;       }
  Float_t   GetShrinkageU() const { return eShrinkageU; }
  Float_t   GetShrinkageD() const { return eShrinkageD; }

  Float_t   GetEmulsion()  const;

  void      SetEmulsion( float thickness );
  void      SetID( int id)           { eID = id;          }
  void      SetBase( float base)     { eBase = base;      }
  void      SetUp( float up )        { eUp = up;          }
  void      SetDown( float down )    { eDown = down;      }
  void      SetSrinkageU( float su)  { eShrinkageU = su;  }
  void      SetSrinkageD( float sd)  { eShrinkageD = sd;  }

  int  ReadParameters( const char *str );

  ClassDef(EdbPlate,2)  // emulsion plate specification
};

//______________________________________________________________________________
class EdbCamera : public TNamed {

private:
  Float_t    eWidth;     // dimension of the image along rows
  Float_t    eHeight;    // dimension of the image along columns
  Int_t      eRows;      // number of rows in the matrix
  Int_t      eColumns;   // number of columns in the matrix

public:
  EdbCamera( float w=0, float h=0, int r=0, int c=0 ): 
    eWidth(w), eHeight(h),  eRows(r), eColumns(c) {}
  virtual ~EdbCamera() { }

  void Print() const;
  void SetCamera(float w=0, float h=0, int r=0, int c=0);

  Float_t   GetWidth()   const { return eWidth;   }
  Float_t   GetHeight()  const { return eHeight;  }
  Int_t     GetRows()    const { return eRows;    }
  Int_t     GetColumns() const { return eColumns; }

  ClassDef(EdbCamera,2)  // CCD camera specification (obsolete)
};

//______________________________________________________________________________
class EdbRunHeader : public TNamed {

private:

  Int_t        eRunID;            //  Run identifier

  Int_t        eFlag[10];         //  customize run flags:
                                  //    eFlag[0] = 1  - UTS data
                                  //    eFlag[0] = 2  - SySal data
                                  //    eFlag[1] = 1  - Stage coordinates
                                  //    eFlag[1] = 2  - Absolute (fiducial) coordinates
                                  //    eFlag[2] = 1  - real (stage) coordinates for clusters
                                  //    eFlag[2] = 2  - pixels coordinates for clusters

  TDatime      eStartTime;        //  Date and Time of run starting
  TDatime      eFinishTime;       //  Date and Time of run finishing
  Long_t       eSystemStartTime;  //  System time of run start in msec
  Double_t     eCPU;              //  CPU time of run production

  TString      eComment;          //  Comment to define this scanning session

  EdbArea     *eArea;             //  definition of scanning area (views path)

  EdbPlate    *ePlate;            //  Emulsion plate specification

  EdbStage    *eStage;            //  stage specification

  Float_t      eXmin, eXmax;      // Run scanning limits
  Float_t      eYmin, eYmax;      //
  Int_t        eNareas;           // total number of areas (fragments) in the run

public:
  EdbRunHeader( int n=0 );
  virtual ~EdbRunHeader() { }

  void         SetRunID( int id )          { eRunID = id; }
  void         SetFlag(int i, int value)   { if(i>=0&&i<10) eFlag[i]=value; }
  void         SetLimits(float xmin, float xmax, float ymin, float ymax) 
    { eXmin=xmin; eXmax=xmax; eYmin=ymin; eYmax=ymax; }
  void         SetNareas(int n)            {eNareas=n;}
  void         SetComment(const char *com) { eComment=com; }
  void         SetCPU(Double_t cpu)        { eCPU = cpu; }

  Int_t        GetFlag(int i)       const { if(i>=0&&i<10) return eFlag[i]; return -999; }
  EdbArea     *GetArea()            const { return eArea;      }
  EdbStage    *GetStage()           const { return eStage;     }
  EdbPlate    *GetPlate()           const { return ePlate;     }
  Int_t        GetRunID()           const { return eRunID;     }
  TString      GetComment()         const { return eComment;   }
  TDatime     *GetStartTime()             { return &eStartTime; }
  TDatime     *GetFinishTime()            { return &eFinishTime; }
  Long_t       GetSystemStartTime() const { return eSystemStartTime; }
  Double_t     GetCPU()             const { return eCPU; }
  Int_t        GetNareas()          const { return eNareas; }

  Float_t      GetXmin()            const { return eXmin; }
  Float_t      GetXmax()            const { return eXmax; }
  Float_t      GetYmin()            const { return eYmin; }
  Float_t      GetYmax()            const { return eYmax; }

  // plate access functions usefull for analysis

  Float_t      EmulsionUp()      const { return ePlate->GetUp(); }
  Float_t      EmulsionBase()    const { return ePlate->GetBase(); }
  Float_t      EmulsionDown()    const { return ePlate->GetDown(); }
  Float_t      ShrinkageUp()     const { return ePlate->GetShrinkageU(); }
  Float_t      ShrinkageDown()   const { return ePlate->GetShrinkageD(); }

  // stage access functions usefull for analysis

  Int_t        NpixelsX()     const { return eStage->GetCCD()->GetNx(); }
  Int_t        NpixelsY()     const { return eStage->GetCCD()->GetNy(); }

  Float_t      PixelX()       const { return eStage->PixelX(); }
  Float_t      PixelY()       const { return eStage->PixelX(); }

  Float_t      CameraWidth()  const { return PixelX()*NpixelsX(); }
  Float_t      CameraHeight() const { return PixelY()*NpixelsY(); }

  EdbAffine2D  *CCD2Stage()   const { return eStage->CCD2Stage(); }


  // stage setting functions

  void  SetCCD( int nx=0, 
		int ny=0, 
		float px=0, 
		float py=0,
		const char *name="",
		const char *title="")
    { 
      eStage->GetCCD()->Set( nx, ny, px, py ); 
      eStage->GetCCD()->SetName( name );
      eStage->GetCCD()->SetTitle( title );
    }

  void SetObjective( float mag, 
	             float a, float b, float c, float d, float p, float q,
		     const char *name="",
		     const char *title="")
    { eStage->GetObjective()->Set( mag, a,b,c,d,p,q ); 
      eStage->GetObjective()->SetName( name );
      eStage->GetObjective()->SetTitle( title );
    }

  void SetPlate( int id, 
		 float up, float base, float down, float shru, float shrd,
		 const char *name="",
		 const char *title="")
    { ePlate->Set( id, up, base, down, shru, shrd );
      ePlate->SetName(name);
      ePlate->SetTitle(title);
    }

  void SetArea( int N, 
		float stepx, float stepy,
		int   ft,    int fb , int path )
    { eArea->Set(N, stepx, stepy, ft, fb, path); }

  // other routines

  void         Print();

  ClassDef(EdbRunHeader,2)  // Run specific information
};

#endif /* ROOT_EdbRunHeader */
