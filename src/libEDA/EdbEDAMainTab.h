#ifndef __EDA_MAINTAB_H__
#define __EDA_MAINTAB_H__

#include"EdbEDAUtil.h"
#include"EdbEDATrackSet.h"
using namespace EdbEDAUtil;


class EdbEDAMainTab{
	private:
	TGCheckButton *eGCBTrackID;
	TGCheckButton *eGCBTrackAngle;
	TGRadioButton *eGRBSelSeg;
	TGRadioButton *eGRBSelTrack;

	public:
	
	enum{ kSelSeg, kSelTrack};

	TGNumberEntry *eGNENseg;
	TGNumberEntry *eGNEMaxPlate;
	TGNumberEntryField *eGNEWCut;
	TGNumberEntryField *eGNEPHDTRMS;
	
	TGCheckButton *eGClearPrevious;

	TGCheckButton *eGSideOut;
	TGNumberEntry *eGSideOutPlate;
	TGCheckButton *eGCBNeighbor;
	TGCheckButton *eGCBIpCut;
	TGCheckButton *eGCBAngleCut;

	TGNumberEntry *eGNETolDT;
	TGNumberEntry *eGNETolDX;
	TGNumberEntryField *eGNEAX;
	TGNumberEntryField *eGNEAY;
	TGNumberEntryField *eGNENeighborDzUp;
	TGNumberEntryField *eGNENeighborDzDown;
	
	// track search
	TGNumberEntryField *eGNESearchTrack;
	TGNumberEntryField *eGNESearchSegmentIPL;
	TGNumberEntryField *eGNESearchSegmentISeg;

	// search couples
	TGNumberEntry *eGNESearchBTIPL;

	TGCheckButton *eGCBExtention;
	TGCheckButton *eGCBExtention_res;
	TGCheckButton *fTextButton_scanback;

	TGHSlider *eGSliderExtentionZ;

	TGNumberEntryField *eGNEExtentionZ;

	TGRadioButton *eGRBColorPH;
	TGRadioButton *eGRBColorPL;
	TGRadioButton *eGRBColorBK;
	TObjArray *eDrawChecks;

	TGCheckButton *eGCBAnimation;
	TGCheckButton *eGCBDrawArea;

	TGMainFrame *eIDWindow;

	EdbEDATrackSet *eTrackSet; // TrackSet for the search.
	
	
	
	EdbEDAMainTab(void){
		eIDWindow=NULL;
		MakeGUIMain();
	}
	virtual ~EdbEDAMainTab(){}
	
	void MakeGUIMain();
	void MakeGUICuts(TGGroupFrame *fGroupFrame1);
	void MakeGUIIDSearch();

	void Redo(void);
	void ExecTrackSearch();
	void ExecSegmentSearch();
	void ExecCuts(int redraw = kTRUE);
	void ApplyParams();
	
	void SearchCouplesIP(int ipl, double x, double y, double z, double ipcut = 50, TObjArray *cps = NULL);
	void SearchCouplesIP(int ipl, EdbVertex *v, double ipcut = 50, TObjArray *cps = NULL);
	void SearchCouplesIP(TObjArray *segsarr = NULL);
	void SearchCouplesDmin(TObjArray *segsarr, int ipl);
	void SearchCouplesDmin(TObjArray *segsarr = NULL);
	
	void ClearTracks(void);
	void ClearVertices();

	void SetPickSegment();
	void SetPickTrack();
	void DrawOnlySelected();
	void RemoveSelected();
	void SelectAll();
	
	void SetExtention();
	void SetDrawText();
	void SetDrawTrackSets(Bool_t on);
	void SetDrawCheckBox();
	void SetExtentionZSlider();
	void SetExtentionZNumber();
	void DrawExtentionResolution();
	int  GetNsegCut(){ return eGNENseg->GetIntNumber();}
	double GetPHCut()  { return eGNEWCut->GetNumber();}
	double GetPHDTRMS(){ return eGNEPHDTRMS->GetNumber();}
	double GetTX(){ return eGNEAX->GetNumber();}
	double GetTY(){ return eGNEAY->GetNumber();}
	double GetDT(){ return eGNETolDT->GetNumber();}
	
	void CalcMinimumDistance(void);
	// calculate Dmin for the all combinations of the selected tracks.
	// original code from Tomoko Ariga
	void CalcVertex(int clear_previous_local=kTRUE);
	void CalcVertexPlus();
	
	void CalcPSelected();
	void CalcIP();
	void SetDrawAreas();
	void PickCenter();
	void SetAnimation();
	void SetColorMode();
	void ReadListFile(const char *filename=NULL, bool clear_previous=kFALSE);
	const char * WriteListFile(const char *filename=NULL);
	void ReadFilteredText();
	void WriteFilteredText();
	void SetProjection();
	
	void WriteFeedback();
	void ReadFeedback();
	void WriteMxxFile();
	
	void SetOutputFileModeButton();
	
	void MicroTrackSearch();
	
	ClassDef(EdbEDAMainTab, 0) // Main GUI
};

#endif // __EDA_MAINTAB_H__
