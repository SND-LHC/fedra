#ifndef __EDA_TRACKSET_H__
#define __EDA_TRACKSET_H__

using namespace EdbEDAUtil;

#include <EdbEDASets.h>
#include<EdbRunTracking.h>

class EdbEDATrackComment:public TObject{
	private:
	TString *eComment;
	EdbTrackP *eTrack;
	
	public:
	EdbEDATrackComment(EdbTrackP *t, const char *cmt):eTrack(t){
		eComment = new TString(cmt);
	}
	EdbTrackP * GetTrack(){ return eTrack;}
	const char * GetComment(){ return eComment->Data();}
	ClassDef(EdbEDATrackComment,0) // for comment
};

/*
class EdbEDATrack: public EdbTrackP{
	private:
	TString *eComment;
	public:
	EdbEDATrack (EdbTrackP *t, char *cmt) : EdbTrackP(*t){
		eComment = new TString(cmt);
	}
	const char * GetComment(){ return eComment->Data();}
	
	ClassDef(EdbEDATrack, 0) // Track class for EDA, inherited from EdbTrackP.
}
*/

class EdbEDATrackSelection {
	int    eNsegCut;
	int    ePlateUp;
	int    ePlateDown;
	int    eNeighborSearch;
	int    eImpactSearch;
	int    eAngularCut;
	int    eSideOut;
	int    eSideOutPlate;
	double eTolDX;
	double eTolDTX;
	double eTolDTY;
	double eTX;
	double eTY;
	
	int eClearPrevious;
	double ePHCut;
	double ePHDTRMS;
	double eNeighborDzUp  ;
	double eNeighborDzDown;
	TObjArray *eSelected;
	
	
	public:
	EdbVertex *eVertex;
	EdbRunTracking eRunTracking; // Run Tracking for microtrack search

	EdbEDATrackSelection(){ Reset(); }
	virtual ~EdbEDATrackSelection(){}
	
	void Reset(){
		eNsegCut        = 3;
		ePlateUp        = 0;
		ePlateDown      = 0;
		eNeighborSearch = 0;
		eImpactSearch   = 0;
		eAngularCut     = 0;
		eSideOut        = 0;
		eSideOutPlate   = 3;
		eTolDX          = 200.;
		eTolDTX         = 2.0;
		eTolDTY         = 2.0;
		eTX   = 0.0;
		eTY   = 0.0;
		eClearPrevious  = 1;
		ePHCut          = 0.0;
		ePHDTRMS        = 0.05;
		eNeighborDzUp   = 5000.;
		eNeighborDzDown = 5000.;
		eVertex         = NULL;
		eSelected       = NULL;
		eRunTracking.ePredictionScan = false;
		SetCondTrackingDefault();
	}
	void SetNsegCut        (int nseg)     { eNsegCut = nseg;}
	void SetUpstreamPlate  (int ipl)      { ePlateUp = ipl;}
	//void SetDownstreamPlate(int ipl)    { ePlateUp = ipl;}
	void SetImpactSearch   (bool b=kTRUE, EdbVertex *v = NULL) { eImpactSearch = b; eVertex = v;}
	void SetNeighborSearch (bool b=kTRUE, TObjArray *selected = NULL, double dzup=-1, double dzdown=-1) { 
		eNeighborSearch = b; eSelected=selected;
		if(dzup>=0)   eNeighborDzUp   = dzup;
		if(dzdown>=0) eNeighborDzDown = dzdown;
	}
	void SetSideOut(bool b=kTRUE) { eSideOut = b;}
	void SetSideOutPlate(int npl) { eSideOutPlate = npl;}
	void SetAngularCut     (bool b=kTRUE) { eAngularCut = b;}
	void SetDX   (double dx)  { eTolDX=dx;}
	void SetDT (double dtx, double dty=-1) { eTolDTX=dtx; eTolDTY=dty; if(eTolDTY<0) eTolDTY=eTolDTX;}
	void SetAngle  (double tx, double ty) { eTX=tx; eTY=ty;}
	void SetClearPrevious  (bool b=kTRUE) { eClearPrevious = b;}
	void SetPHCut      (double phcut =0.0) { ePHCut=phcut;}
	void SetPHDTRMS (double slope) {ePHDTRMS=slope;}
	bool GetImpactSearch () { return eImpactSearch;}
	int Neighborhood (EdbSegP *s, double *dmin=NULL);
	int ImpactSearch (EdbTrackP *t);
	int SideOut(EdbTrackP *t);
	int GetUpstreamPlate() { return ePlateUp;}
	int GetNsegCut() { return eNsegCut;}
	double GetPHCut() { return ePHCut;}
	double GetPHDTRMS() { return ePHDTRMS;}
	int GetSideOutPlate() { return eSideOutPlate;}

	void DoSelection (TObjArray *tracksbase, TObjArray *tracks);
	
	void SetRunTracking(EdbRunTracking &rt) {eRunTracking=rt;}
	EdbRunTracking &GetRunTracking(){ return eRunTracking;}
	
	void SetCondTrackingDefault(){
		eRunTracking.eDeltaRview = 400.;
		eRunTracking.eDeltaTheta = 0.2;
		eRunTracking.eDeltaR = 20.;

		eRunTracking.ePreliminaryPulsMinMT = 7;
		eRunTracking.ePreliminaryChi2MaxMT = 3.0;

		eRunTracking.ePulsMinMT       = 9;
		eRunTracking.ePulsMinDegradMT = 1.5;
		eRunTracking.eChi2MaxMT       = 1.6;

		eRunTracking.ePulsMinBT       = 18;
		eRunTracking.ePulsMinDegradBT = 1.5;
		eRunTracking.eChi2MaxBT       = 1.8;

		eRunTracking.eDegradPos   = 5.;
		eRunTracking.eDegradSlope = 0.003;

		SetCondMTDefault( eRunTracking.eCondMT);
		SetCondBTDefault( eRunTracking.eCondBT);
	}

	void SetCondMTDefault(EdbScanCond &cond){
		cond.SetSigma0( 1., 1., 0.010, 0.010 );     // sigma0 "x, y, tx, ty" at zero angle
		cond.SetDegrad( 5. );                       // sigma(tx) = sigma0*(1+degrad*tx)
		cond.SetBins(0, 0, 0, 0);  //???            // bins in [sigma] for checks
		cond.SetPulsRamp0(  5., 5. );               // in range (Pmin:Pmax) Signal/All is nearly linear
		cond.SetPulsRamp04( 5., 5. );
		cond.SetChi2Max( 6.5 );
		cond.SetChi2PMax( 6.5 );
		cond.SetRadX0( 5810. );
		cond.SetName("OPERA_microtrack");
	}

	void SetCondBTDefault(EdbScanCond &cond){
		cond.SetSigma0( 10., 10., 0.007, 0.007 );   // sigma0 "x, y, tx, ty" at zero angle
		cond.SetDegrad( 2. );                       // sigma(tx) = sigma0*(1+degrad*tx)
		cond.SetBins(0, 0, 0, 0);                   // bins in [sigma] for checks
		cond.SetPulsRamp0(  5., 5. );               // in range (Pmin:Pmax) Signal/All is nearly linear
		cond.SetPulsRamp04( 5., 5. );
		cond.SetChi2Max( 6.5 );
		cond.SetChi2PMax( 6.5 );
		cond.SetRadX0( 5810. );
		cond.SetName("OPERA_basetrack");
	}
	
	ClassDef(EdbEDATrackSelection,0)  // Selection of Tracks
};





class EdbEDATrackSet: public TNamed, public EdbEDATrackSelection{
	
	private:
	
//	int    eVer;
//	int    eRun;
	EdbID eID;
	
	int    eColorMode;
	int    eTrackColor;
	int    eTrackWidth;
	double eTrackLength;

	int    eLayerColor;
	int    eLayerWidth;
	double eLayerLength;
	
	int    eTextColor;
	int    eTextSize;
	int    eTextFont;
	int    eTextPosition;
	
	int    eLayerColorMT;
	
	bool   eDraw;            // control of drawing
	bool   eDrawMT;          // control of drawing microtrack
	bool   eDrawTrackID;     // control of Track ID drawing
	bool   eDrawTrackIDG;    // Global control of Track ID drawing by Main tab.
	bool   eDrawTrackAngle;  // control of Track Angle drawing.
	bool   eDrawTrackAngleG; // Global control of Track angle drawing by Main tab.
	
	int    eExtendMode;
	int    eInterpolation;
	
	EdbPVRec   *ePVR;  // pointer for original EdbPVRec. tracks are copied to eTracksBase through AddTracksPVR().
	EdbDataSet *eDataSet;
	EdbScanSet *eScanSet;
	
	TObjArray *eTracks;     // array of tracks to be drawn
	TObjArray *eTracksBase; // array of tracks, base of track search
	TObjArray *eTracksPrevious;     // previous array of tracks
	
	TObjArray *eComments;
	
	EdbEDAAreaSet      *eAreaSet;

	public:
	
	enum { kEnd, kTop, kMiddle};
	
	EdbEDATrackSet(const char *title, TObjArray *EdbTrackPArray=NULL) 
			: eID(0,0,1,0), ePVR(NULL), eDataSet(NULL), eScanSet(NULL){
		SetName(title);
		eTracks         = new TObjArray;
		eTracksBase     = new TObjArray;
		eTracksPrevious = new TObjArray;
		eComments       = new TObjArray;
		eDraw = kTRUE;
		eColorMode = kCOLOR_BY_PH;
		
		eAreaSet = new EdbEDAAreaSet();
		
		SetTrackAttribute();
		SetTextAttribute();
		SetDrawMT(kTRUE);
		SetDrawTrackID(kTRUE);
		SetDrawTrackAngle(kTRUE);
		SetDrawTrackIDG(kFALSE);
		SetDrawTrackAngleG(kFALSE);
		SetExtendMode(kExtendAuto);
		SetInterpolation(kTRUE);
		
		eTextPosition = kEnd;
		
		if(EdbTrackPArray) {
			SetTracks(EdbTrackPArray);
			SetTracksBase(EdbTrackPArray);
		}
	}
	
	void Print(){
		printf("EdbEDATrackSet : %s\n", GetName());
		printf("TracksBase     : %5d tracks\n", NBase());
		printf("Tracks to Draw : %5d tracks\n", N());
		
	}
	
	void SetTracksBase(TObjArray *tracksbase) {
		// Copy tracks to eTracksBase. tracksbase should be an array of EdbTrackP
		*eTracksBase = *tracksbase;
	}
	void AddTracksBase(TObjArray *tracksbase) {
		// Add tracks to eTracksBase. tracksbase should be an array of EdbTrackP
		if(tracksbase==NULL) return;
		for(int i=0;i<tracksbase->GetEntriesFast();i++) eTracksBase->Add(tracksbase->At(i));
	}
	
	TObjArray * GetTracksBase() { return eTracksBase;}
	
	void AddTrack(EdbTrackP *t)    {
		// Add a track to both eTracksBase and eTracks, if it is not yet registered.
		if(eTracksBase->FindObject(t)==NULL) eTracksBase->Add(t);
		if(eTracks->FindObject(t)==NULL) eTracks->Add(t);
	}
	void AddTracks(TObjArray *arr) {
		// Add tracks to both eTracksBase and eTracks, if it is not yet registerd.
		// arr is a TObjArray of EdbTrackP.
		if(NULL==arr) return; for(int i=0;i<arr->GetEntriesFast();i++) AddTrack( (EdbTrackP *) arr->At(i));
	}
	
	void SetTrack(EdbTrackP *t) {
		// Set a track to be drawn.
		// the track is requiered to be already registered in eTracksBase.
		// if you want to simply draw the track, use AddTrack().
		if(eTracksBase->FindObject(t)&&eTracks->FindObject(t)==NULL) eTracks->Add(t);
	}
	
	void SetTrack(int itrack){
		// Set Tracks to draw with the given track id.
		for(int i=0; i<N(); i++) {
			EdbTrackP *t = GetTrack(i);
			if(t->ID()==itrack) SetTrack(t);
		}
	}
	
	void SetTracks(TObjArray *tracks) {
		// Set tracks to be draw.
		// the tracks are requiered to be already resistered in eTracksBase.
		// if you want to simply draw the tracks, use AddTracks().
		if(NULL==tracks) return; for(int i=0; i<tracks->GetEntriesFast(); i++) SetTrack((EdbTrackP *)tracks->At(i));
	}
	
	TObjArray * GetTracks() { return eTracks;}
	
	void AddTracksPVR(EdbPVRec *pvr) {
		// COPY tracks from EdbPVRec.
		// original PVRec is accessible through GetPVRec().
		ePVR = pvr;
		AddTracksBase(ePVR->eTracks);
		SetTracksAll();
		eAreaSet->SetAreas(ePVR);
	}
	EdbPVRec * GetPVRec() { return ePVR;}
	void SetPVRec(EdbPVRec *pvr) { ePVR=pvr;}
	
	void RemoveTrackBase(EdbTrackP *t){ if(eTracksBase->FindObject(t)) {eTracksBase->Remove(t); eTracksBase->Sort();} RemoveTrack(t);}
	void RemoveTrack(EdbTrackP *t){ if(eTracks->FindObject(t)) {eTracks->Remove(t); eTracks->Sort();}}
	void RemoveTracks(TObjArray *tracks) { for(int i=0; i<tracks->GetEntriesFast(); i++)  RemoveTrack((EdbTrackP *) tracks->At(i));}
	
	void ClearTracks(void) {
		// Clear Tracks to be draw.
		eTracks->Clear();
	}
	void ClearTracksBase(void) {
		// Clear eTracksBase.
		eTracksBase->Clear();
	}
	
	void Clear() { 
		// Clear Tracks and TracksBase.
		ClearTracks(); ClearTracksBase(); eAreaSet->Clear();
	}
	
	void SetTracksAll(void) { *eTracks = *eTracksBase;}
	void SetTracksNseg(int nsegcut = 3) {
		eTracks->Clear();
		for(int i=0; i<NBase(); i++) {
			EdbTrackP *t = GetTrackBase(i);
			if(t->N()<nsegcut) continue;
			SetTrack(t);
		}
	}
	
	void RestorePrevious(void) { 
		TObjArray arr = *eTracks;
		*eTracks = *eTracksPrevious;
		*eTracksPrevious = arr;
	}
	void StorePrevious()   { *eTracksPrevious = *eTracks;}

	
	EdbTrackP * SearchCorrespondTrack(EdbTrackP *t);
	
	int MicroTrackSearch(EdbTrackP *t, EdbSegP *pred, int ipl);
	void MicroTrackSearch(EdbTrackP *t);
	int PrepareScanSetForMT();
	
	EdbTrackP * SearchTrack(int itrk) { for(int i=0;i<eTracksBase->GetEntriesFast(); i++) {if(itrk==((EdbTrackP *) eTracksBase->At(i))->ID()) return (EdbTrackP *) eTracksBase->At(i);} return NULL;}
	
	EdbTrackP * SearchSegment(int ipl, int iseg);
	
	EdbTrackP * FindTrack(EdbSegP *s, double dx = 10, double dt = 0.02){
		EdbSegP *ss = new EdbSegP;
		for(int i=0;i<NBase();i++){
			EdbTrackP *t = GetTrackBase(i);
			if(t==NULL) continue;
			*ss = *(t->GetSegmentFirst());
			ss->PropagateTo(s->Z());
			if( fabs(s->X()-ss->X()) < dx){
			if( fabs(s->Y()-ss->Y()) < dx){
			if( fabs(s->TX()-ss->TX()) < dt) {
			if( fabs(s->TY()-ss->TY()) < dt) {
				return t;
			}}}}
			
		}
		delete ss;
		return NULL;
	}

	EdbTrackP * FindTrack(EdbTrackP *t){
		if(eTracks==NULL) return NULL;
		if(eTracks->FindObject(t)) return t;
		return NULL;
	}
	
	EdbTrackP * FindTrackLast(EdbSegP *s, double dx = 20, double dt = 0.03){
		
		EdbSegP ss;
		for(int i=0;i<eTracksBase->GetEntriesFast();i++){
			EdbTrackP *tt = (EdbTrackP *) eTracksBase->At(i);
			for(int j=0;j<tt->N();j++){
				ss.Copy(*tt->GetSegment(j));
				ss.PropagateTo(s->Z());
				if(fabs(s->TX()-ss.TX())<dt){
				if(fabs(s->TY()-ss.TY())<dt){
					if( fabs(s->X()-ss.X())<dx ){
					if( fabs(s->Y()-ss.Y())<dx ){
						return tt;
					}}
				}}
			}
		}
		return NULL;
	}

	void AddSegments(TObjArray *segments){
		if(segments==NULL) return;
		if(segments->GetEntriesFast()==0) return;
		
		for(int i=0;i<segments->GetEntriesFast();i++){
			EdbSegP *s = (EdbSegP *)segments->At(i);
			EdbTrackP *t = new EdbTrackP(s);
			t->SetCounters();
			t->SetPlate(s->Plate());
			AddTrack(t);
		}
		return;
	}
	
	EdbTrackP *GetTrack(EdbSegP *s){
		if(IsTrack(s)) return (EdbTrackP*) s;
		if(IsSegment(s)){
			for(int i=0;i<eTracksBase->GetEntriesFast();i++){
				EdbTrackP *t = (EdbTrackP *) eTracksBase->At(i);
				for(int j=0;j<t->N();j++) if(s==t->GetSegment(j)) return t;
			}
		}
		return NULL;
	}
	
	int N() { return eTracks->GetEntriesFast();}
	int NBase() { return eTracksBase->GetEntriesFast();}
	EdbTrackP * GetTrack(int i){ return (EdbTrackP *) eTracks->At(i);}
	EdbTrackP * GetTrackBase(int i){ return (EdbTrackP *) eTracksBase->At(i);}
	bool IsMember(EdbTrackP *t) {return (eTracksBase->FindObject(t) != NULL) ? kTRUE : kFALSE;}
	
	int  FindUpstreamPlate() {
		if(NBase()==0) return 0;
		int iplup = 1000;
		for(int i=0;i<NBase(); i++) {
			int ipl= GetTrackBase(i)->GetSegmentFirst()->Plate();
			if(ipl<iplup) iplup=ipl;
		}
		return iplup;
	}
	int  FindDownstreamPlate() {
		int ipldw = -1000;
		for(int i=0;i<NBase(); i++) {
			int ipl= GetTrackBase(i)->GetSegmentLast()->Plate();
			if(ipl>ipldw) ipldw=ipl;
		}
		return ipldw;
	}
	
	void SetID(EdbID id) { eID=id;}
	EdbID& GetID() {return  eID;}
//	void SetRun(int run) { eRun = run;}
//	void SetVer(int ver) { eVer = ver;}
//	int  GetRun() { return eRun;}
//	int  GetVer() { return eVer;}
	
	void SetScanSet (EdbScanSet *ss) { eScanSet = ss;}
	EdbScanSet * GetScanSet(){ return eScanSet;}
	void SetDataSet (EdbDataSet *dset) { eDataSet = dset;}
	EdbDataSet * GetDataSet(){ return eDataSet;}
	
	void AddComment(EdbTrackP *t, const char *cmt){ eComments->Add(new EdbEDATrackComment(t,cmt));}
	void ClearComments() {eComments->Clear();}
	const char * GetComment(EdbTrackP *t) { 
		for(int i=0;i<eComments->GetEntriesFast();i++){
			EdbEDATrackComment *cmt = (EdbEDATrackComment *) eComments->At(i);
			if(cmt->GetTrack()==t) {
				return cmt->GetComment();
			}
		}
		return "";
	}
	
	void ReadFile(const char *filename = NULL, int datatype = 100, TCut rcut = "1");
	void ReadTracksTree (const char *scanset_filename, TCut cut="1");
	void ReadTracksTree (EdbID ID, TCut cut = "1"); // for EdbScanProc style
	EdbDataProc *ReadLinkedTracks (const char *lnkdef="lnk.def", TCut cut = "1"); // for EdbDataProc style
	EdbDataProc *ReadLinkedTracksFile (const char *filename="linked_tracks.root", TCut cut = "1"); // for EdbDataProc style
	EdbDataProc *ReadLinkedTracks (EdbDataProc *dproc); // for EdbDataProc style. expect EdbDataProc after InitVolume.

	void ReadPredictionScan(int BRICK, int SBVERSION, int RUN_PRED, int UseMicrotrack=0); // for EdbDataProc style.
	void ReadPredictionScan(EdbID id, bool force_update_setroot=0); // for EdbScanProc style
	void ReadPredictionScan(EdbScanSet *ss); // for EdbScanProc style
	void ReadTextTracks(const char *filename = NULL);
	void ReadMxxFile(const char *filename);
	void ReadListFile  (const char *filename = NULL, bool clear_previous=kTRUE);
	const char *WriteListFile  (const char *filename = NULL, bool append = kFALSE, bool open_editor = kFALSE);

	void DoVertexing(EdbVertexRec *VR = NULL, int ntrkcut = 3);
	void DoSelection() { EdbEDATrackSelection::DoSelection(eTracksBase,eTracks);}
	
	void SetTracksVertex(EdbVertex *v, bool clear_previous = kTRUE){
		// Set eTracks as the tracks which belongs to the vertex.
		if(clear_previous) ClearTracks();
		for(int j=0 ; j < v->N() ; j++){
			EdbTrackP *t = v->GetTrack(j);
			if(t==NULL) continue;
			SetTrack(t);
		}
	}
		
	void SetTracksVertices(TObjArray *vertices, bool clear_previous = kTRUE){
		// Set eTracks as the tracks which belongs to the vertices.
		printf("Set Tracks for %d Vertices\n", vertices->GetEntriesFast());
		if(clear_previous) ClearTracks();
		for(int i=0;i<vertices->GetEntriesFast();i++){
			EdbVertex *v = (EdbVertex *) vertices->At(i);
			if(v==NULL) continue;
			SetTracksVertex(v, kFALSE);
		}
	}
	
	void SetExtendMode(int mode=kExtendAuto) { eExtendMode=mode;} 
	int GetExtendMode() { return eExtendMode;} 
	void SetInterpolation(int interpolation=kTRUE) { eInterpolation = interpolation;} // kTRUE: fill hole with dashed line.
	
	void SetTrackAttribute(
		int    track_color  =  -1, 
		int    track_width  =   1, 
		double track_length =  -1, 
		int    layer_color  = kYellow, 
		int    layer_width  =   2, 
		double layer_length = 150,
		int    layer_color_MT  = kCyan){
		// set track attribute. -1 means automatic setting.
		eTrackColor  = track_color;
		eTrackWidth  = track_width;
		eTrackLength = track_length;
		eLayerColor  = layer_color;
		eLayerWidth  = layer_width>track_width? layer_width: track_width+1;
		eLayerLength = layer_length;
	}
	
	void GetTrackAttribute( int *track_color , int *track_width, double *track_length, int *layer_color, int *layer_width, double *layer_length, int *layer_color_MT){
		*track_color  = eTrackColor ;
		*track_width  = eTrackWidth ;
		*track_length = eTrackLength;
		*layer_color  = eLayerColor ;
		*layer_width  = eLayerWidth ;
		*layer_length = eLayerLength;
		*layer_color_MT = eLayerColorMT;
	}
	void SetLayerLength(double length) { eLayerLength=length;}
	
	void ResetTrackAttribute(void) { SetTrackAttribute(); }
	void SetColorMode(int mode);
	
	int    GetTrackColor (EdbSegP *s);
	int    GetTrackWidth (EdbSegP *s){ return eTrackWidth;}
	double GetTrackLength(EdbSegP *s, int updown = 1); // up = -1, down = 1
	
	int    GetLayerColor (EdbSegP *s){ return s->Side()==0?eLayerColor:eLayerColorMT;}
	int    GetLayerWidth (EdbSegP *s){ return eLayerWidth;}//eLayerWidth>=0 ? eLayerWidth : (int)s->W()/10+1;}
	double GetLayerLength(EdbSegP *s){ return eLayerLength;};

	
	void SetTextAttribute(int size=14, int color=kWhite, int fontfile=4){
		eTextSize  = size;
		eTextColor = color;
		eTextFont  = fontfile;
	}

	int GetTextColor() { return eTextColor;}
	int GetTextSize () { return eTextSize;}
	int GetTextFont () { return eTextFont;}
	
	void SetTextPosition (int pos) { eTextPosition = pos;}
	int  GetTextPosition () { return eTextPosition;}
	
	void SetDraw(bool b=kTRUE) { eDraw = b;} // Setting Draw or not.
	bool GetDraw() { return eDraw;}
	void SetDrawMT        (bool b=kTRUE) { eDrawMT = b;} // Setting Draw or not.
	bool GetDrawMT() { return eDrawMT;}
	void SetDrawTrackID   (bool draw) { eDrawTrackID = draw;}
	void SetDrawTrackAngle(bool draw) { eDrawTrackAngle = draw;}
	void SetDrawTrackIDG   (bool draw) { eDrawTrackIDG = draw;}
	void SetDrawTrackAngleG(bool draw) { eDrawTrackAngleG = draw;}
	bool GetDrawTrackID   () { return eDrawTrackID;}
	bool GetDrawTrackAngle() { return eDrawTrackAngle;}
	

	void Draw(int redraw=kFALSE);
	void DrawSingleTrack(EdbTrackP *t);
	void DrawSingleSegment(EdbSegP *s, const char *elname, TEveCompound *cmp = NULL, int extendup = 1, int extenddown = 1);
	void DrawSingleComment(EdbEDATrackComment *tcmt);

	// Area related
	void SetAreaSet(EdbEDAAreaSet *set) { eAreaSet=set;}
	EdbEDAAreaSet * GetAreaSet(){ return eAreaSet;}
	void SetDrawAreaSet(bool draw = kTRUE) { eAreaSet->SetDraw(draw);}

	ClassDef(EdbEDATrackSet,0)  // Track Set
};

class EdbEDATrackSetTab {
	private:
	bool eEnable;

	EdbEDATrackSet *eSet;

	TGNumberEntryField *fBrick;
	TGNumberEntryField *fRun;
	TGNumberEntryField *fVer;
	TGTextButton *fPredScan;
	TGTextButton *fReadText;
	TGCheckButton *fColorAuto;
	TGCheckButton *fReadMicro;
	TGCheckButton *fTrackID;
	TGCheckButton *fTrackAngle;
	TGComboBox    *fTrackSet;
	TGColorSelect *fTrackColor;
	TGLineWidthComboBox *fTrackWidth;
	TGColorSelect *fLayerColor;
	TGLineWidthComboBox *fLayerWidth;
	TGNumberEntryField *fLayerLength;
	TGComboBox    *fExtendMode;
	TGComboBox    *fFontSize;
	TGComboBox    *fFontFile;
	TGColorSelect *fFontColor;
	TGMainFrame   *fFontWindow;
	TGTextButton  *fScanID; // Label to indicate scan id.
	
	TGLabel *fDataText;
	
	public:
	EdbEDATrackSetTab(void);
	void Select(int id) {fTrackSet->Select(id, kTRUE);}
	void Selected(const char *txt);
	void SetAttributeGUIs(EdbEDATrackSet *set);
	void Enable(bool b=kTRUE) {eEnable = b;}
	void Update();
	void ToggleAuto(bool b) { fTrackColor->Enable(!b); Update();}
	void ReadPredictionScan();
	void ReadTextTracks();
	void ReadFile();
	void ClearTracks();
	void MakeFontWindow();
	void CloseFontWindow();
	void SetFont();
	void ChangeScanID();
};


#endif // __EDA_TRACKSET_H__
