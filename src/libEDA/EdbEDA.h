
#ifndef __EDA_H__
#define __EDA_H__

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbEDA                                                               //
//                                                                      //
// Mother of EDA (Event Display and Analyser)                           //
//                                                                      //
// it consists of 5 EdbTrackSet (TS,SB,SF,BT,MN) and 1 EdbVertexSet.    //
//                                                                      //
// It inherits following.                                               //
//  EdbEDACamera     : View, Pictures and Animation                     //
//  EdbEDAIO         : I/O of feedback file, filtered text, cs_info.txt, Basetracks. //
//  EdbEDABrickData  : Manage Brick info and store EdbDataProc. PID <-> IPL conversion. //
//  EdbEDASelection  : Manage Selected segments and tracks.             //
//////////////////////////////////////////////////////////////////////////


#include "EdbEDAUtil.h"
#include "EdbEDATrackSet.h"
#include "EdbEDASets.h"
#include "EdbEDAMainTab.h"
#include "EdbEDATabs.h"
#include "EdbEDAPlotTab.h"
#include "EdbEDADecaySearch.h"
#include "EdbEDAShowerTab.h"

//-------------------------------------------------------------------

using namespace EdbEDAUtil;

class EdbEDAObject : public TObject {
	private:
	int eDraw;
	void *eRef;
	public:
	EdbEDAObject():eDraw(kTRUE),eRef(NULL){}
	~EdbEDAObject(){}
	virtual void Draw(Option_t *option=""){ printf("override this");}
	virtual void SetColorMode(int mode){ printf("override this.\n");}
	void SetReference(void *ref) { eRef=ref;}
	void *GetReference(){ return eRef;}
	bool GetDraw(){ return eDraw;}
	void SetDraw(bool b=kTRUE){ eDraw=b;}
	ClassDef(EdbEDAObject,0) // test
};

class EdbEDAText : public EdbEDAObject{
	private:
	float eX, eY, eZ;
	int    eTextColor;
	int    eTextSize;
	int    eTextFont;
	TString eText;
	void *eRef;
	
	public:
	EdbEDAText(float x, float y, float z, TString str, int color=kWhite) : 
		eX(x), eY(y), eZ(z), eText(str){SetAttribute(14,color);}
	~EdbEDAText(){}

	void SetAttribute(int size=14, int color=kWhite, int fontfile=4){
		eTextSize  = size;
		eTextColor = color;
		eTextFont  = fontfile;
	}
	void SetColorMode(int mode){
		if(mode==kBLACKWHITE) SetAttribute(eTextSize, kBlack, eTextFont);
		else SetAttribute(eTextSize, kWhite, eTextFont);
	}
	int GetTextColor() { return eTextColor;}
	int GetTextSize () { return eTextSize;}
	int GetTextFont () { return eTextFont;}
	void SetXYZ(float x, float y, float z){ eX=x; eY=y; eZ=z;}
	void SetText(TString str){ eText=str;}
	void AddText(TString str){ eText+=str;}
	int N(){return eText.Length();}
	void Draw(Option_t *option="");
	
	ClassDef(EdbEDAText, 0) // a comment class
	
};

class EdbEDALine : public EdbEDAObject{
	private:
	float X1,Y1,Z1, X2,Y2,Z2;
	int col, wid;
	TString text;
	
	public:
	EdbEDALine( float x1, float y1, float z1, float x2, float y2, float z2, int color, int width, const char *comment = "") : X1(x1), Y1(y1), Z1(z1), X2(x2), Y2(y2), Z2(z2), col(color), wid(width), text(comment){}
	EdbEDALine( int color, int width, const char *comment = "") : col(color), wid(width), text(comment){}
	~EdbEDALine(){}
	
	void SetPoint1(float x1, float y1, float z1){ X1=x1; Y1=y1; Z1=z1;}
	void SetPoint2(float x2, float y2, float z2){ X2=x2; Y2=y2; Z2=z2;}
	
	void Draw(Option_t *option="");
	ClassDef(EdbEDALine, 0) // Line as EdbEDAObject
};

//-------------------------------------------------------------------
#include<TGLOverlay.h>
class EdbEDAOverlay : public TGLOverlayElement{
	
	private:
	bool eDraw;
	int eXa,eYa, eLength;
	int eXb,eYb;
	int eBarDirection;
	public:
	EdbEDAOverlay();
	~EdbEDAOverlay(){}
	void Render(TGLRnrCtx& ctx);
	void RenderBar(TGLRnrCtx& ctx);
	void RenderAxes(TGLRnrCtx& ctx);
	void SetDraw(bool b) {eDraw=b;}
	void SetXY(int pixel_x=50, int pixel_y=30);
	void SetLength(int pixel_length=120);
	ClassDef(EdbEDAOverlay, 0) // a class for draw overlays, i.e. scale bars.
};


class EdbEDACamera {
	
	private:
	TTimer * eTimer;
	
	int eAnimationMode;
	float eTheta;
	float eAnimationAmplitude;
	EdbEDAOverlay *eOverlay;
	
	public:
	TGLViewer * eViewer;
	
	enum {X,Y,Z,B};
	
	EdbEDACamera():eTimer(NULL),eAnimationMode(0), eTheta(3.14/2), eAnimationAmplitude(0.002), eViewer(NULL){}
//	~EdbEDACamera(){if(NULL==eTimer) {eTimer->Stop(); delete eTimer;}}
	virtual ~EdbEDACamera(){}
	//TGLViewer * GetGLViewer(){ return eViewer;}
	void SetCamera(int projection = Y);
	void Animation();
	void SetAnimationMode();
	void StartAnimation();
	void StopAnimation();
	void SaveAnimation(const char *filename = NULL, int n=100, int interval=10);
	void SavePictures();
	void Snapshot(const char *filename=NULL);
	TGLViewer * GetGLViewer();
	EdbEDAOverlay *SetOverlay() {eOverlay=new EdbEDAOverlay; eViewer->AddOverlayElement(eOverlay); return eOverlay;}
	EdbEDAOverlay *GetOverlay() {return eOverlay;}
	
	
	ClassDef(EdbEDACamera,0) // part of EdbEDA. Manage camera control and pictures.
};


class EdbEDAIO{
	public:
	int eOutputFileMode;
	enum OutputFileMode{
		kBern,
		kOtherLabs
	};
	
	TCut eCutCP;
	
	EdbEDAIO():eOutputFileMode(kBern), eCutCP("eN1==1&&eN2==1&&eCHI2P<s.eW*0.12-1"){}
	virtual ~EdbEDAIO(){}
	EdbPattern * ReadCouples(int ipl, EdbPattern *pat=NULL);
	EdbPattern * ReadCouplesPID(int pid, EdbPattern *pat=NULL);
	EdbPattern * GetPatternIPL(int pid);
	void WriteFeedbackFile(const char *filename = NULL);
	void WriteFeedbackFile2008(EdbVertex *v = NULL, const char *filename = NULL);
	void ReadFeedbackFile(const char *filename = NULL);
	EdbVertex * ReadFeedbackFile2008(const char *filename = NULL);
	int IdMuon(const char *filename = "../cs_info.txt");
	void OpenTextEditor(const char *filename);
	void WriteFilteredText(const char *filename = NULL);
	EdbPVRec * ReadFilteredText(const char *filename = NULL);
	
	//void PrintTrackFeedback(EdbTrackP *t=NULL, EdbVertex *v1=NULL, EdbVertex *v2=NULL, EdbSegP *s1=NULL, FILE *io=stdout, bool print_segs=kTRUE);
	
	
	ClassDef(EdbEDAIO,0) // part of EdbEDA. Manage I/O.
};

class EdbEDABrickData{
	// Data management

	private:
	double eBeamTX, eBeamTY;
	public:

	EdbID eIDSF;
	EdbID eIDSB;
	EdbID eID;
	
	EdbDataSet  *eDataSet;
	EdbScanProc *eScanProc;
	EdbScanSet  *eScanSet;
	EdbPVRec    *ePVR;
	
	int eDataType;
	
	int eFirstPlate; // first plate of brick, 1
	int eLastPlate;  // last plate of brick, 57
	
	EdbEDABrickData() : eBeamTX(0), eBeamTY(0.058), eDataSet(0), eScanSet(0), ePVR(0), eDataType(-1), eFirstPlate(1), eLastPlate(57){
	}
	virtual ~EdbEDABrickData(){}
	void SetIDSF(EdbID id) { eIDSF=id;}
	void SetIDSB(EdbID id) { eIDSB=id;}
	void SetID(EdbID id) { eID=id;}
	EdbID GetID() { return eID;}
	
	void SetScanSet(EdbScanSet *ss) { eScanSet = ss;}
	EdbScanSet *GetScanSet() { return eScanSet;}
	
	void SetPVR(EdbPVRec *pvr){ ePVR=pvr;}
	EdbPVRec *GetPVR(){if(ePVR==NULL) printf("EdbEDA::GetPVR() is only to obtain global Z coordinate.\n"); return ePVR;}
	
	void SetBrick(int ibrick){
		eID.eBrick   = ibrick;
		eIDSF.eBrick = ibrick;
		eIDSF.eBrick = ibrick;
	}
	int GetBrick() { return eID.eBrick;}
	
	int  GetDataType(){ return eDataType;}
	void SetDataType(int datatype){ eDataType=datatype;}
	
	int GetPlateUp()  {return GetIPL(0);} // upstream of data volume
	int GetPlateDown(){return GetIPL(eDataSet->N()-1);} // downstream of data volume
	
	int GetFirstPlate(){ return eFirstPlate;} // first plate of brick
	int GetLastPlate(){ return eLastPlate;}  // last plate of brick
	
	
	EdbDataSet * GetDataSet(){ return eDataSet;}
	
	int GetNPL(){
		if      (eScanSet) return eScanSet->eB.Npl();
		else if (eDataSet) return eDataSet->N();
		else return 0;
	}
	
	int GetIPL(int PID){
		if(eScanSet) { // using ScanSet for calibration.
			EdbPlateP *pl = eScanSet->eB.GetPlate(PID);
			return pl->ID(); 
		}
		else if( eDataSet ){
			if(PID>=eDataSet->N()) return PID;
			EdbDataPiece *piece = eDataSet->GetPiece(PID);
			if(NULL==piece) return PID;
			return piece->Plate();
		}
		
		return PID;
	}
	
	int GetIPLZ(float z){
		// return ipl just downstream of given z.
		for(int ipl=GetFirstPlate(); ipl<=GetLastPlate(); ipl++){
			float zz= GetZ(ipl)+0.1; // 0.1 is for the case that given z is exactly at the plate.
			if(zz>=z) return ipl;
		}
		return GetLastPlate();
	}
	
	int GetPID(int ipl){
		if(eScanSet) { // return the ordinal number in eB.
			EdbBrickP *b = &eScanSet->eB;
			for(int i=0;i<b->Npl();i++){
				EdbPlateP *pl = b->GetPlate(i);
				if(pl->ID()==ipl) return i;
			}
			return ipl; // if not correspondence found.
		}
		else if( eDataSet ){
			// Get PID from iPlate.
			// if iPlate is out of volume scan, return -1;
			for(int i=0;i<eDataSet->N();i++){
				EdbDataPiece *piece = eDataSet->GetPiece(i);
				if(NULL==piece) continue;
				if(ipl==piece->Plate()){
					return i;
				}
			}
			// printf("error at GetPID(%d)\n", ipl);
		}
		return ipl;
	}
	
	double GetZ(int ipl){
		// Get Z value for give plate number.
		if(eScanSet) {
			EdbPlateP *pl = eScanSet->GetPlate(ipl);
			if(pl!=NULL) return pl->Z();
			return (ipl-57)*1300.0;
		}
		else if ( eDataSet ){
			for(int i=0;i<eDataSet->N();i++){
				if( GetIPL(i) == ipl ) return GetZPID(i);
			}
		}
		else if( ePVR ){
			for(int i=0;i<ePVR->Npatterns();i++){
				EdbPattern *pat = ePVR->GetPattern(i);
				if(pat->Plate()==ipl) return pat->Z();
			}
		}
	 	return (ipl-57)*1300.0;
	}
	
	double GetZPID(int pid){
		// Get Z from PID
		// if pid is out of range, set pid within the range.

		if(eScanSet) {
			EdbBrickP& b = eScanSet->Brick();
			EdbPlateP *pl = b.GetPlate(pid);
			return pl->Z();
		}
		else if (eDataSet){
		
			if(pid<0) pid=0;
			if(pid>eDataSet->N()-1) pid=eDataSet->N()-1;
			EdbDataPiece *piece = eDataSet->GetPiece(pid);
			EdbLayer* l = piece->GetLayer(0);
			return l->Z();
		}
		
		return GetZ(pid);
	}

	void SetBeamAngle(double tx, double ty) { eBeamTX=tx; eBeamTY=ty;}
	double GetBeamTX() { return eBeamTX;}
	double GetBeamTY() { return eBeamTY;}
	
	EdbAffine2D * GetAffine(int ipl);

	ClassDef(EdbEDABrickData, 0) // part of EdbEDA. Brick Data Management for EDA
};


class EdbEDASelection{
	private:
	TEveSelection *eSel;
	int eDisable;
	EdbVertex *eSelectedVertex;  // current(last) selected vertex
	TObjArray *eSelectedVertices; // selected vertices
	TObjArray *eSelectedTracks;   // selected tracks
	TObjArray *eSelectedSegments; // selected segments
	TObjArray *eSelected;


	public:

	EdbEDASelection(){
		eSelectedVertices = new TObjArray;
		eSelectedTracks   = new TObjArray;
		eSelectedSegments = new TObjArray;
		eSelected         = new TObjArray;
		eSelectedVertex = NULL;
		eSel = NULL;
		eDisable = 0;
	}

	virtual ~EdbEDASelection(){
		if(eSel!=NULL){
			eSel->Disconnect("SelectionAdded(TEveElement*)", this, "Added(TEveElement*)");
			eSel->Disconnect("SelectionRemoved(TEveElement*)", this, "Removed(TEveElement*)");
			eSel->Disconnect("SelectionCleared()", this, "Cleared()");
		}
	}

	// init Selection, connect to the TEveManager::eSelection = gSel
	void InitSelection(TEveSelection* sel){
		if(eSel!=NULL){
			eSel->Disconnect("SelectionAdded(TEveElement*)", this, "Added(TEveElement*)");
			eSel->Disconnect("SelectionRemoved(TEveElement*)", this, "Removed(TEveElement*)");
			eSel->Disconnect("SelectionCleared()", this, "Cleared()");
		}

		eSel = sel; // supposed to be sel = gEDA->GetSelection()
		eSel->Connect("SelectionAdded(TEveElement*)",   "EdbEDASelection", this, "Added(TEveElement*)");
		eSel->Connect("SelectionRemoved(TEveElement*)", "EdbEDASelection", this, "Removed(TEveElement*)");
		eSel->Connect("SelectionCleared()",             "EdbEDASelection", this, "Cleared()");
	}
	
	void ClearSelected(){
		eSelected->Clear();
		eSelectedTracks->Clear();
	}
	
	EdbTrackP * GetTrack( TEveElement *e);

	TObjArray * SetSelected(TObjArray *selected = NULL);
	TObjArray * GetSelected(void) {return eSelected;}
	EdbSegP * GetSelected(int i) {return (EdbSegP *)eSelected->At(i);}

	void SetSelectedTracks(TObjArray *arr) {
		if(arr) {
			*eSelectedTracks = *arr; 
			*eSelected       = *arr;
		}
	}

	
	TEveElement * GetSelectedElement () { TEveElement::List_i it= gEve->GetSelection()->EndChildren(); return *(--it);}
	TObjArray * GetSelectedTracks(void)  { return eSelectedTracks; }
	
	int NSelectedTracks() { return eSelectedTracks->GetEntriesFast();}
	int NSelected()       { return eSelected->      GetEntriesFast();}
	EdbTrackP * GetSelectedTrack(int i = -1){
		if(NSelectedTracks()==0) return NULL;
		if(i==-1) i=NSelectedTracks()-1; 
		return (EdbTrackP *) eSelectedTracks->At(i); 
	}
	EdbTrackP * GetSelectedTrack(EdbSegP *s){
		for( int i=0; i<NSelectedTracks(); i++){
			EdbTrackP *t = GetSelectedTrack(i);
			for(int j=0; j<t->N(); j++) if(s==t->GetSegment(j)) return t;
		}
		return NULL;
	}
	
	EdbSegP   * GetSelectedSegment(int i=-1){if(i==-1) i=NSelected()-1;         return (EdbSegP *)   eSelected->At(i);}

	void SelectAll();
	void SetSelection ( TObjArray * selected = NULL);

	EdbVertex *GetSelectedVertex(void) { return eSelectedVertex;}
	void SetSelectedVertex(EdbVertex *v) { eSelectedVertex = v; }
	void ClearSelectedVertex() {eSelectedVertex = NULL;}

	void Added(TEveElement* el);
	void Removed(TEveElement* el);
	void Cleared();
	
	int Disable(int force = kFALSE) { return force ? eDisable = 1 : eDisable++;}
	int Enable (int force = kFALSE) { return force ? eDisable = 0 : eDisable--;}

	ClassDef(EdbEDASelection,0)  // part of EdbEDA. Manage Selection.
};

/*	// Set PID as PlateNumber.
	EdbDataSet *dset = dproc.GetDataSet();
	
	printf("hoge %d %d\n", __LINE__, dset);
	for(int j=0;j<ali.Ntracks();j++){
		EdbTrackP *t = ali.GetTrack(j);
		
		for(int k=0;k<t->N();k++){
			EdbSegP *s = t->GetSegment(k);
			int pid = s->PID();
			
			EdbDataPiece *piece = dset->GetPiece(pid);
			int ipl = piece->Plate();
			s->SetPID(ipl);
		}
		t->SetPID(t->GetSegmentFirst()->PID());
	}
	printf("hoge %d\n", __LINE__);
*/

//-------------------------------------------------------------------
class EdbEDA :
	public EdbEDACamera,
	public EdbEDAIO,
	public EdbEDABrickData,
	public EdbEDASelection,
	public TNamed {
	
	private:
	
	TObjArray *eDrawObjects;
	TObjArray *eTrackSets;
	
	EdbEDATrackSet     *eSetTracks;
	EdbEDAVertexSet    *eVertexSet;
	EdbEDAExtentionSet *eExt;
	TEveEventManager   *eEvent;

	EdbEDAMainTab      *eMainTab;
	EdbEDATrackSetTab  *eTrackSetTab;
	EdbEDAPlotTab      *ePlotTab;
	EdbEDAVertexTab    *eVertexTab;
	EdbEDAOperationTab    *eOperationTab;
	EdbEDASelection *eSel;
	EdbEDADecaySearchTab *eDSTab;
	
	double eScaleZ;
	bool eJapanese;
	
	public:
	
	EdbEDA(const char *filename = NULL, int datatype = 100, TCut rcut = "1", bool auto_run=kTRUE) : 
				TNamed("EDA","Event Display and Analyser"){
		// Constructor for EDA
		//  filename : filename for "LinkDef" or "*.set.root(ScanSet)" or "*.root(linked_track.root format)"
		//  datatype : only for "LinkDef". Data type for EdbDataProc::InitVolume(). put -1 for no InitVolume
		//             if datatype= 100 or 1000. linked tracks will be registred to "TS".
		//  rcut     : Track selection when read.
		//  auto_run : Run automatically.
		OpenFile(filename, datatype, rcut);
		if(auto_run) Run();
	}

	EdbEDA(EdbID IDCalib, const char *ProcDirClient=NULL) {
		Init();
		if(ProcDirClient!=NULL) eScanProc->eProcDirClient=ProcDirClient;
		eID=IDCalib;
		eScanSet = eScanProc->ReadScanSet(eID);
	}
	
	EdbEDA(EdbScanSet *ss, const char *ProcDirClient=NULL) : TNamed("EDA","Event Display and Analyser"){
		Init();
		if(ProcDirClient!=NULL) eScanProc->eProcDirClient=ProcDirClient;
		eScanSet = ss;
		eID = *((EdbID *)ss->eIDS.At(0));
		eID.ePlate = 0;
	}
	
	EdbEDA(EdbDataProc *dproc) : TNamed("EDA","Event Display and Analyser") {
		Init();
		eDataSet  = dproc->GetDataSet();
		GetTrackSet("TS")->AddTracksPVR(dproc->PVR());
		GetTrackSet("TS")->SetNsegCut(3);
		GetTrackSet("TS")->DoSelection();
	}

	void OpenFile(const char *filename = NULL, int datatype = 100, TCut rcut = "1");	

	EdbEDA(EdbPVRec *pvr, bool autorun=1){
		Init();
		ePVR = pvr;
		
		if(pvr->Ntracks()==0){
			printf("*******************************************\n");
			printf(" Worning: No tracks in EdbPVRec object! \n");
			printf("*******************************************\n");
			printf("EDA is made to show \"Tracks\". if you want to show segments in PVR, make tracks from the segments. or you may try EdbEDAUtil::FillTracksFromPatterns(EdbPVRec *)\n");
		}
		
		// if plate number is not set. set plate number as PID.
		int flag_plset=0;
		if(ePVR->eTracks==NULL) ePVR->eTracks = new TObjArray;
		for(int i=0;i<ePVR->Ntracks();i++){
			EdbTrackP *t = ePVR->GetTrack(i);
			EdbSegP *s = t->GetSegmentLast();
			if(s->Plate()!=0) flag_plset++; // if plate number is set.
		}
		
		if(flag_plset==0){ // in case plate number is not set, set PID as plate number
			for(int i=0;i<ePVR->Npatterns();i++){
				EdbPattern *p = ePVR->GetPattern(i);
				if(p==NULL) continue;
				p->SetSegmentsPlate(p->ID());
			}
		}
		
		GetTrackSet("TS")->AddTracksPVR(pvr);
		GetVertexSet()->AddVertices(pvr->eVTX);
		if(autorun) Run();
	}
	
	EdbDataProc * PVRec2DataProc(EdbPVRec *pvr);
	
	void Init(void);
	void Run();	
	EdbScanProc *ScanProc() { return eScanProc;}
	
	EdbEDATrackSet *AddTrackSet(const char *name) { EdbEDATrackSet *set = new EdbEDATrackSet(name); eTrackSets->Add(set); return set;}
	
	void AddDrawObject(EdbEDAObject *o) { eDrawObjects->Add(o);}
	void RemoveDrawObject(TObject *o){ eDrawObjects->Remove(o); eDrawObjects->Sort();}
	void ClearDrawObjects(){ eDrawObjects->Clear();}
	int  NDrawObjects(){ return eDrawObjects->GetEntriesFast();}
	EdbEDAObject *GetDrawObject(int i){ return (EdbEDAObject*) eDrawObjects->At(i);}

	TObjArray * GetCurrentTracks(void){
		TEveEventManager *ev = gEve->GetCurrentEvent();
		TObjArray *tracks = new TObjArray;
		for(TEveElement::List_i it = ev->BeginChildren(); it!=ev->EndChildren(); it++){
			TEveElement *e = *it;
			if(IsTrack(e)) tracks->Add((EdbTrackP *)e->GetUserData());
		}
		return tracks;
	}
	
	EdbTrackP * GetTrack( EdbSegP *s){
		for(int i=0;i<NTrackSets();i++){
			EdbTrackP *t = GetTrackSet(i)->GetTrack(s);
			if(t) return t;
		}
		return NULL;
	}
	
	TEveElement * GetEveElement(TObject *t){
		// return TEveElement of the track, vertex. if element is not in the Event, return NULL.
		TEveEventManager *ev = gEve->GetCurrentEvent();
		for(TEveElement::List_i it = ev->BeginChildren(); it!=ev->EndChildren(); it++){
			TEveElement *e = *it;
			if( t == e->GetUserData() ) return e;
		}
		return NULL;
	}
	
	int NTrackSets() { return eTrackSets->GetEntriesFast();}
	EdbEDATrackSet * GetTrackSet(int i){ return (EdbEDATrackSet *) eTrackSets->At(i);}
	EdbEDATrackSet * GetTrackSet(const char *name){ return (EdbEDATrackSet *) eTrackSets->FindObject(name);}
	EdbEDATrackSet * GetTrackSet(EdbTrackP *t){ 
		for(int i=0;i<NTrackSets();i++) 
		if ( GetTrackSet(i)->IsMember(t) ) 
		return (EdbEDATrackSet *) eTrackSets->At(i);
		return NULL;
	}
	EdbEDATrackSet * GetTrackSet(EdbID& ID){ 
		for(int i=0;i<NTrackSets();i++) {
			EdbID& ID2 = GetTrackSet(i)->GetID();
			if(ID.eMajor==ID2.eMajor&&ID.eMinor==ID2.eMinor) return GetTrackSet(i);
		}
		return NULL;
	}
	
	EdbEDATrackSet * CreateTrackSet(const char *name) {
		EdbEDATrackSet *set=new EdbEDATrackSet(name);
		eTrackSets->Add(set);
		return set;
	}
	TObjArray * GetTrackSets(){ return eTrackSets;}
	
	void RemoveTrackBase(EdbTrackP *t) {for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->RemoveTrackBase(t);}
	void RemoveTrack    (EdbTrackP *t) {for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->RemoveTrack(t);}
	
	void RemoveSelected(){ 
		// remove selected tracks or vertices.
		// to remove vertex, select vertex but not select tracks.
		int nselected = NSelected();
		for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->RemoveTracks(GetSelectedTracks());
		if(nselected==0) eVertexSet->RemoveVertex(GetSelectedVertex());
	}

	void SetOnlySelected(){
		for(int i=0;i<NTrackSets();i++) {
			GetTrackSet(i)->ClearTracks();
			GetTrackSet(i)->SetTracks(GetSelectedTracks());
		}
		eVertexSet->ClearVertices();
		eVertexSet->SetVertex(GetSelectedVertex());
	}
	
	EdbEDAVertexSet * GetVertexSet() { return eVertexSet;}
	void AddVertex(EdbVertex *v) { eVertexSet->AddVertex(v);}
	TObjArray *GetVertices() {return eVertexSet->GetVertices();}
	void ClearVertices() { eVertexSet->ClearVertices();}
	void DrawVertices(int redraw = 0) { eVertexSet->Draw(redraw);}
	
	void SetBrick(int ibrick) { 
		// Set brick ID. it sets also for all TrackSets.
		EdbEDABrickData::SetBrick(ibrick);
		for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->GetID().eBrick=ibrick;
		if(gEve) gEve->GetMainWindow()->SetWindowName(
			Form("OPERA Event Display and Analyser. Brick %d", EdbEDABrickData::GetBrick()));		
	}
	
	void SetDrawTrackID(bool draw)    { for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->SetDrawTrackIDG(draw);}
	void SetDrawTrackAngle(bool draw) { for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->SetDrawTrackAngleG(draw);}
	void SetDrawAreaSet(bool draw)    { for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->GetAreaSet()->SetDrawG(draw);}

	void Draw(int redraw = kFALSE);
	
	void Redraw() { if(gEve)Draw(kTRUE);}
	void UpdateScene() { if(gEve) gEve->GetDefaultGLViewer()->UpdateScene();}
	
	void Reset(){ 
		for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->Clear();
		ClearDrawObjects();
		ClearVertices();
	}
	
	void StorePrevious() {
		for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->StorePrevious();
		GetVertexSet()->StorePrevious();
	}
	
	void RestorePrevious() {
		for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->RestorePrevious();
		GetVertexSet()->RestorePrevious();
	}
	
	EdbEDAExtentionSet * GetExtentionSet(){ return eExt;}
	
	void SetColorMode(int mode){
		// Set color mode,
		//	enum EdbEDAUtil::EDACOLOR{
		//		kCOLOR_BY_PLATE,
		//		kCOLOR_BY_PH,
		//		kCOLOR_BY_ID,
		//		kCOLOR_BY_PARTICLE, // to be done
		//		kBLACKWHITE
		//	};

		printf("color mode %d\n", mode);
		if(mode==kBLACKWHITE) eViewer->SetClearColor(kWhite);
		else eViewer->SetClearColor(kBlack);
		for(int i=0;i<NTrackSets();i++)   GetTrackSet(i)  ->SetColorMode(mode);
		GetVertexSet()->SetColorMode(mode);
		for(int i=0;i<NDrawObjects();i++) GetDrawObject(i)->SetColorMode(mode);
		Redraw();
	}
	void SetLayerLength(double layer_length) {
		// Set LayerLength for all TrackSet.
		for(int i=0;i<NTrackSets();i++){ GetTrackSet(i)->SetLayerLength(layer_length);}
	}
	
	EdbEDAMainTab     *GetMainTab()     {return eMainTab;}
	EdbEDATrackSetTab *GetTrackSetTab() {return eTrackSetTab;}
	EdbEDAPlotTab     *GetPlotTab()     {return ePlotTab;}
	EdbEDAVertexTab   *GetVertexTab()   {return eVertexTab;}
	EdbEDAOperationTab   *GetOperationTab()   {return eOperationTab;}
	EdbEDADecaySearchTab *GetDSTab()    {return eDSTab;}
	EdbTrackP * CheckScanback(EdbTrackP *t) {return GetTrackSet("SB")->FindTrackLast(t->GetSegmentLast());}
	
	double GetScaleZ() { return eScaleZ;}
	void SetScaleZ(double scale) { eScaleZ=scale;}
	bool Japanese(){return eJapanese;}
	void SetJapanese(bool b=kTRUE) {eJapanese = b;}
	
	ClassDef(EdbEDA,0); // Base of Display
};

R__EXTERN EdbEDA* gEDA;


#endif // __EDA_H__
