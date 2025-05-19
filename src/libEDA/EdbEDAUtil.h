#ifndef __EDAUtil_H__
#define __EDAUtil_H__

#include "RVersion.h"
#if ROOT_VERSION_CODE < ROOT_VERSION(5,22,0)
#define __BEFORE_ROOT_5_22__
#endif

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>

#include<TROOT.h>
#include<TSystem.h>
#include<TRint.h>
#include<TStyle.h>
#include<TObject.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TText.h>
#include<TLine.h>
#include<TH2I.h>
#include<TColor.h>
#include<TString.h>
#include<TApplication.h>
#include<TGLViewer.h>
#include<TGeoManager.h>
#include<TGedEditor.h>
#include<TEveManager.h>
#include<TEveEventManager.h>
#include<TEveElement.h>
#include<TEveLine.h>
#include<TEveText.h>
#include<TEveCompound.h>
#include<TEveBrowser.h>
#include<TEveSelection.h> // To enforce auto-loading of libEve.
#include<TEveBoxSet.h>
#include<TEveTrans.h>
#include<TEveWindow.h>
#include<TEveViewer.h>
#include<TGTab.h>
#ifndef __BEFORE_ROOT_5_22__
#include<TEveWindowManager.h>
#endif

#include<TRootEmbeddedCanvas.h>
#include<TGLFontManager.h>
#include<TGNumberEntry.h>
#include<TGFileDialog.h>
#include<TGSlider.h>
#include<TGLabel.h>
#include<TGListBox.h>
#include<TGTextEditor.h>
#include<TGMsgBox.h>
#include<TGInputDialog.h>
#include<TGColorSelect.h>
#include<TGComboBox.h>

//#ifndef __CINT__
#include<EdbVertex.h>
#include<EdbPattern.h>
#include<EdbDataSet.h>
#include<EdbPVRec.h>
#include<EdbMomentumEstimator.h>
#include<EdbScanSet.h>
#include<EdbScanProc.h>
#include<EdbBrick.h>

//#endif


#define CHECK printf("check %15s %5d\n",__FILE__,__LINE__)


namespace EdbEDAUtil{

	enum EDACOLOR{
		kCOLOR_BY_PLATE,
		kCOLOR_BY_PH,
		kCOLOR_BY_ID,
		kCOLOR_BY_MCID,
		kCOLOR_BY_PARTICLE, // to be done
		kBLACKWHITE
	};
	
	enum EDAEXTENDMODE{
		kExtendAuto,
		kExtendUpDown,
		kExtendDown,
		kExtendUp,
		kExtendNo};
	
	int IsVertex(TObject *o);
	int IsVertex(TEveElement *e);
	int IsTrack(TObject *o);
	int IsTrack(TEveElement *e);
	int IsSegment(TObject *o);
	int IsSegment(TEveElement *e);
	int IsIncludeCouples(EdbPattern *pat);
	TEveCompound * GetTrackElement ( EdbTrackP * t);
	TEvePointSet * GetVertexElement ( EdbVertex * v);
	double CalcIP(EdbSegP *s, double x, double y, double z);
	double CalcIP(EdbSegP *s, EdbVertex *v);
	double CalcDmin(EdbSegP *seg1, EdbSegP *seg2, double *dminz = NULL);
	
	double CalcMinimumKinkAngle( TVector3 vertex, TVector3 daupos,    TVector3 daumom);
	double CalcMinimumKinkAngle( EdbVertex *v1ry, EdbSegP *tdaughter, int z_is_middle_of_base=1);
	double CalcMinimumKinkAngle( EdbSegP *t1ry,   EdbSegP *tdaughter, int z_is_middle_of_base=1);
	double CalcPtmin( TVector3 vertex, TVector3 daupos,    TVector3 daumom);
	double CalcPtmin( EdbVertex *v1ry, EdbSegP *tdaughter, int z_is_middle_of_base=1);
	double CalcPtmin( EdbSegP *t1ry,   EdbSegP *tdaughter, int z_is_middle_of_base=1);
	double CalcKinkAngle(EdbSegP *tparent, EdbSegP *tdaughter);
	double CalcPt(EdbSegP *tparent, EdbSegP *tdaughter);
	double CalcEMCSelectron(EdbTrackP *t);
	
	EdbVertex * CalcVertex(TObjArray *segments); // calc vertex from the segments array (EdbSegP*)
	EdbMomentumEstimator *CalcP(EdbTrackP *t, double& p, double& pmin, double& pmax, bool print=kTRUE);
//	void   CalcP2(EdbTrackP *t, double& p, double& pmin, double& pmax);
//	void   CalcP3(EdbTrackP *t, double& p, double& pmin, double& pmax);
	void   CalcPPartial(EdbTrackP *t, EdbSegP *s1st, EdbSegP *slast, double& p, double& pmin, double& pmax, bool print=kTRUE);
	double CalcDistance(EdbSegP *s1,EdbSegP *s2, double z);
	
	void CalcDTTransLongi(EdbSegP *s1, EdbSegP *s2, double *dtTransverse, double *dtLongitudinal);
	void CalcDTTransLongi(double tx1, double ty1, double tx2, double ty2, double *dtTransverse, double *dtLongitudinal);
	void CalcDXTransLongi(EdbSegP *s1, EdbSegP *s2, double *dxt, double *dxl);
	
	double DTRMS(EdbTrackP *t);
	double DTRMSelectron(EdbTrackP *t);
	double DTRMS1Kink(EdbTrackP *t, int *NKinkAngleUsed=NULL);
	double DTRMSTL(EdbTrackP *t, double *rmsspace, double *rmstransverse, double *rmslongitudinal, int *ndata=NULL);
	double DTRMSTL1Kink(EdbTrackP *t, double *rmsspace, double *rmstransverse, double *rmslongitudinal, int *NKinkAngleUsed=NULL);
	double DTRMSTLGiven1Kink(EdbTrackP *t, int iKink, double *rmsspace, double *rmstransverse, double *rmslongitudinal, int *NKinkAngleUsed=NULL);
	EdbTrackP * CleanTrack(EdbTrackP *t);
	
	void   ErrorMessage(const char *title, const char *message);
	void   ErrorMessage(const char *message);
	bool   AskYesNo(const char *message);
	int    InputNumberInteger(const char *message, int idefault=0);
	double InputNumberReal(const char *message, double default_num=0.0, TGNumberFormat::EStyle es=TGNumberFormat::kNESReal);
	int    InputID(const char *message, EdbID &id);
	
	char * FindProcDirClient();
	int FindBrickIDFromPath();
	
	EdbPVRec * ReadFeedbackPVR(const char *filename = NULL);
	EdbPVRec * ReadMxxPVR(const char *filename = NULL);
	void WritePVRMxx(EdbPVRec *pvr, const char *filename = NULL);
	void WriteTracksMxx(TObjArray *pvr, const char *filename = NULL);
	
	void MakePVRFromTracksArray(TObjArray *tracks_or_segments, EdbPVRec& pvr);
	void FillTracksFromPatterns(EdbPVRec *pvr);
};





class EdbIDDialog : public TGTransientFrame {
	private:
	TGNumberEntryField *eNumberEntry1;
	TGNumberEntryField *eNumberEntry2;
	TGNumberEntryField *eNumberEntry3;
	TGNumberEntryField *eNumberEntry4;
	EdbID& eID;
	int *eRetVal;
	public:
	
	EdbIDDialog(const char *message, EdbID& id, int *retval = NULL) : 
		TGTransientFrame(gClient->GetRoot(), gEve?gEve->GetMainWindow():0,200,100,kMainFrame | kVerticalFrame), eID(id), eRetVal(retval){
		SetLayoutBroken(kTRUE);
		SetWindowName("EdbIDDialog");
		SetCleanup(kDeepCleanup);
		
		CenterOnParent();
		
		unsigned int wsize = 250;
		if(strlen(message)*6>wsize) wsize=strlen(message)*6;
		
		int posy=10;
		
		TGLabel *fLabel = new TGLabel(this, message);
		fLabel->MoveResize(20,posy,wsize-40,20);
		AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
		posy+=30;
		int posx=10;
		int dx = 50;
		eNumberEntry1 = new TGNumberEntryField(this, -1, eID.eBrick, TGNumberFormat::kNESInteger);
		               AddFrame(eNumberEntry1, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
		               eNumberEntry1->MoveResize(posx,posy,dx,20);
		eNumberEntry2 = new TGNumberEntryField(this, -1, eID.ePlate, TGNumberFormat::kNESInteger);
		               AddFrame(eNumberEntry2, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
		               eNumberEntry2->MoveResize(posx+=dx+5,posy,dx,20);
		eNumberEntry3 = new TGNumberEntryField(this, -1, eID.eMajor, TGNumberFormat::kNESInteger);
		               AddFrame(eNumberEntry3, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
		               eNumberEntry3->MoveResize(posx+=dx+5,posy,dx,20);
		eNumberEntry4 = new TGNumberEntryField(this, -1, eID.eMinor, TGNumberFormat::kNESInteger);
		               AddFrame(eNumberEntry4, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
		               eNumberEntry4->MoveResize(posx+=dx+5,posy,dx,20);
		
		posy+=30;
		TGTextButton *fb = new TGTextButton(this,"OK");
		AddFrame(fb, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
		fb->MoveResize((wsize-50)/2,posy,50,20);

		SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
		MapSubwindows();
		MapWindow();
		Resize(wsize,100);
		eNumberEntry4->RequestFocus();
		gClient->WaitFor(this);
	}
	virtual ~EdbIDDialog(){
		Cleanup();
	}
	
	virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t){
		// Process message dialog box event.
		switch (GET_MSG(msg)) {
			case kC_COMMAND:
				switch (GET_SUBMSG(msg)) {
					case kCM_BUTTON:
						//printf("%d\n", msg);
						//if (fRetCode) *fRetCode = (Int_t) parm1;
						eID.eBrick = eNumberEntry1->GetIntNumber();
						eID.ePlate = eNumberEntry2->GetIntNumber();
						eID.eMajor = eNumberEntry3->GetIntNumber();
						eID.eMinor = eNumberEntry4->GetIntNumber();
						if(eRetVal) *eRetVal = kTRUE;
						delete this;
						break;
					default:
						break;
				}
				break;
			default:
				break;
		}
		return kTRUE;
	}
	
	virtual void CloseWindow() { 
		if(eRetVal) *eRetVal = kFALSE;
		delete this;
	}
	
	ClassDef(EdbIDDialog, 0) // Input GUI for integer number.
};


#endif //__EDAUtil_H__
