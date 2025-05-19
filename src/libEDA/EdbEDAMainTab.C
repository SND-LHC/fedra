#include"EdbEDA.h"

extern EdbEDA *gEDA;
using namespace EdbEDAUtil;


ClassImp(EdbEDAMainTab)

void EdbEDAMainTab::ApplyParams(){
	eTrackSet = gEDA->GetTrackSet("TS");
	
	eTrackSet->SetUpstreamPlate  (eGNEMaxPlate->GetIntNumber());
	eTrackSet->SetNsegCut        (eGNENseg->GetIntNumber());
	eTrackSet->SetSideOut        (eGSideOut->IsOn());
	eTrackSet->SetSideOutPlate   (eGSideOutPlate->GetIntNumber());
	eTrackSet->SetNeighborSearch (eGCBNeighbor->IsOn(), gEDA->GetSelected(), 
									eGNENeighborDzUp->GetNumber(), eGNENeighborDzDown->GetNumber());
	eTrackSet->SetImpactSearch   (eGCBIpCut->IsOn(), gEDA->GetSelectedVertex());
	eTrackSet->SetAngularCut     (eGCBAngleCut->IsOn());
	eTrackSet->SetClearPrevious  (eGClearPrevious->IsOn());
	eTrackSet->SetDX  ( eGNETolDX->GetNumber());
	eTrackSet->SetDT  (eGNETolDT->GetNumber(), eGNETolDT->GetNumber());
	eTrackSet->SetAngle ( eGNEAX->GetNumber(),   eGNEAY->GetNumber());
	eTrackSet->SetPHCut (eGNEWCut->GetNumber());
	eTrackSet->SetPHDTRMS(eGNEPHDTRMS->GetNumber());
}

void EdbEDAMainTab::ExecCuts(int redraw){
	printf("Execute cuts\n");
	ApplyParams();
	printf("Apply Parameters\n");
	gEDA->StorePrevious();
	
	if(eTrackSet->GetImpactSearch()){
		if(eGClearPrevious->IsOn()) gEDA->ClearVertices();
		gEDA->AddVertex(gEDA->GetSelectedVertex());
	}
	
	// selection of tracks in the tracks-base.
	eTrackSet->DoSelection();
	
	// selected tracks are anyway shown.
	eTrackSet->SetTracks(gEDA->GetSelectedTracks());
	
	if(redraw) gEDA->Redraw();
}

void EdbEDAMainTab::Redo(void) { gEDA->RestorePrevious(); gEDA->Redraw();}

void EdbEDAMainTab::SearchCouplesIP(int ipl, EdbVertex *v, double ipcut, TObjArray *cps){
	if(NULL==v) return;
	if(NULL==cps) return;
	SearchCouplesIP (ipl, v->X(), v->Y(), v->Z(), ipcut, cps);
}


void EdbEDAMainTab::SearchCouplesIP(TObjArray *segsarr){
	// This will be called from Main tab.
	gEDA->StorePrevious();
	if(NULL==segsarr) segsarr = new TObjArray;
	EdbVertex *selected_vertex=gEDA->GetSelectedVertex();
	if(NULL==selected_vertex) return;
	int ipl      = eGNESearchBTIPL->GetIntNumber();
	double ipcut = eGNETolDX->GetIntNumber();
	if(ipl==0) return;
	
	if(ipl==-1){
		for(int i=1; i<60;i++) SearchCouplesIP ( i, selected_vertex, ipcut, segsarr);
	}
	else SearchCouplesIP (ipl, selected_vertex, ipcut, segsarr);
	
	EdbEDATrackSet *TrackSet = gEDA->GetTrackSet("BT");
	if(eGClearPrevious->IsOn()) TrackSet->ClearTracks();

	TrackSet->AddSegments(segsarr);
	gEDA->Redraw();
}

void EdbEDAMainTab::SearchCouplesIP(int ipl, double x, double y, double z, double ipcut, TObjArray *cps){
	if(NULL==cps) return;
	EdbPattern *pat = gEDA->GetPatternIPL(ipl);
	if(pat==NULL) {
		//ErrorMessage(Form("Cannot read couples. pl=%d",ipl));
		//printf("Cannot read couples. pl=%d. stop.\n",ipl);
		return;
	}
	printf("Search BaseTrack in PL %d with IP cut %.1lf micron wrt (%.1lf, %.1lf, %.1lf)\n",
		ipl, ipcut, x, y, z);
	
	EdbEDATrackSet *TS = eTrackSet;
	printf("%2s %6s %5s %2s %8s %8s\n", "pl","id","IP","ph","tx","ty");
	
	int k=0;
	int phcut = (int) GetPHCut();
	for(int i=0;i<pat->N();i++){
		EdbSegP *s = pat->GetSegment(i);
		if(s->W()<phcut) continue;
		double ip = EdbEDAUtil::CalcIP(s, x, y, z);
		if(ipcut < ip) continue;
		
		// already found as tracks to be drawn?
		int flag=0;
		for(int j=0;j<TS->N();j++){
			EdbTrackP *t = TS->GetTrack(j);
			for(int k=0;k<t->N();k++){
				EdbSegP *ss = t->GetSegment(k);
				if(ss->PID()==s->PID()&&ss->ID()==s->ID()){
					flag++;
				}
			}
		}
		if(flag) {
			printf("%2d %6d %5.1lf %2d %8.4f %8.4f already selected as Track\n", 
						ipl, s->ID(), ip, (int)s->W(), s->TX(), s->TY());
			continue;
		}
		else 	printf("%2d %6d %5.1lf %2d %8.4f %8.4f\n", ipl, s->ID(), ip, (int)s->W(), s->TX(), s->TY());

		cps->Add(s);
		k++;
	}
	
	printf("%d basetracks are selected.\n",k);
}


void EdbEDAMainTab::SearchCouplesDmin(TObjArray *segsarr, int ipl){
	
	EdbPattern *pat = gEDA->GetPatternIPL(ipl);
	if(pat==NULL) {
		//ErrorMessage(Form("Cannot read couples. pl=%d",ipl));
		//printf("Cannot read couples. pl=%d. stop.\n",ipl);
		return;
	}
	printf("Search BaseTrack in PL %d with dmin cut %.1lf micron wrt the selected tracks\n", ipl, eGNETolDX->GetNumber());
	
	eTrackSet->SetNeighborSearch(kFALSE, gEDA->GetSelected());
	
	int k=0;
	printf("%2s %7s %5s %2s %8s %8s %4s\n", "PL","ID","Dmin","PH","TX","TY","Chi2");
	double phcut = GetPHCut();
	for(int i=0;i<pat->N();i++){
		EdbSegP *s = pat->GetSegment(i);
		if(s->W()<phcut) continue;

		double dmin;
		if(eTrackSet->Neighborhood(s,&dmin)){
			// already found as tracks to be drawn?
			int flag=0;
			for(int j=0;j<eTrackSet->N();j++){
				EdbTrackP *t = eTrackSet->GetTrack(j);
				for(int k=0;k<t->N();k++){
					EdbSegP *ss = t->GetSegment(k);
					if(ss->PID()==s->PID()&&ss->ID()==s->ID()){
						flag++;
					}
				}
			}
			
			if(flag) {
				printf("%2d %7d %5.1lf %2d %8.4f %8.4f %4.2f already selected as Track\n", 
					ipl, s->ID(), dmin, (int)s->W(), s->TX(), s->TY(), s->Chi2());				continue;
			}

			if(NULL==segsarr->FindObject(s)) {
				segsarr->Add(s);
				k++;
				printf("%2d %7d %5.1lf %2d %8.4f %8.4f %4.2f\n", 
					ipl, s->ID(), dmin, (int)s->W(), s->TX(), s->TY(), s->Chi2());
			}
		}
	}
	printf("%d basetracks are selected\n",k);
	return;
}


void EdbEDAMainTab::SearchCouplesDmin(TObjArray *segsarr){
	// this is called from MainTab.
	gEDA->StorePrevious();
	ApplyParams();
	if(segsarr==NULL) segsarr = new TObjArray;
	int ipl   = eGNESearchBTIPL->GetIntNumber();
	
	if(ipl==-1) {
		printf("Search BaseTrack in all plate with dmin cut %.1lf micron wrt the selected tracks\n", eGNETolDX->GetNumber());
		for(int i=0;i<60;i++) {
			SearchCouplesDmin(segsarr, i);
		}
	}
	else{
		SearchCouplesDmin(segsarr, ipl);
	}

	EdbEDATrackSet *TrackSet = gEDA->GetTrackSet("BT");
	if(eGClearPrevious->IsOn()) TrackSet->ClearTracks();

	TrackSet->AddSegments(segsarr);
	gEDA->Redraw();
}


void EdbEDAMainTab::MicroTrackSearch(){
	// Micro-track search
	// search +- 1 plate from end of track.
	// fill holes by microtrack.

	if(gEDA->NSelectedTracks()==0){
		ErrorMessage("Please select tracks\n");
		return;
	}
	
	for(int i=0;i<gEDA->NSelectedTracks();i++){
		EdbTrackP *t = gEDA->GetSelectedTrack(i);
		EdbEDATrackSet *set = gEDA->GetTrackSet(t);
		if(set!=gEDA->GetTrackSet("TS")){
			ErrorMessage("Currently only MT search for TrackSet:TS is supported.");
			
			int move=AskYesNo("Do you want to move the selected track to TS?");
			if(move==kFALSE) continue;
			else {
				gEDA->GetTrackSet("TS")->AddTrack(t);
				set->RemoveTrack(t);
				set->RemoveTrackBase(t);
				printf("track %5d move from TrackSet %s to %s.\n", t->ID(), set->GetName(), "TS");
			}
		}
		set->MicroTrackSearch(t);
	}
	gEDA->Redraw();
}


void EdbEDAMainTab::ExecTrackSearch(){ 
	gEDA->StorePrevious();
	int itrk = eGNESearchTrack->GetIntNumber();
	if(eGClearPrevious->IsOn()) eTrackSet->ClearTracks();
	printf("Track Search, itrk = %d\n", itrk);
	EdbTrackP *t = eTrackSet->SearchTrack(itrk);
	if(t) eTrackSet->AddTrack(t);
	gEDA->Redraw();
	
	if(eIDWindow!=NULL) eIDWindow->CloseWindow();
	eIDWindow = NULL;

}
void EdbEDAMainTab::ClearTracks(void) { gEDA->StorePrevious(); eTrackSet->ClearTracks(), gEDA->Redraw();}
void EdbEDAMainTab::ClearVertices()   { gEDA->StorePrevious(); gEDA->ClearVertices(); gEDA->Redraw();}

void EdbEDAMainTab::ExecSegmentSearch(){
	gEDA->StorePrevious();
	int ipl = eGNESearchSegmentIPL->GetIntNumber();
	int iseg = eGNESearchSegmentISeg->GetIntNumber();
	printf("Track Search by Segment, ipl = %d iseg = %d\n", ipl, iseg);
	if(eGClearPrevious->IsOn()) eTrackSet->ClearTracks();
	EdbTrackP *t = eTrackSet->SearchSegment(ipl, iseg);
	if(t) eTrackSet->AddTrack(t);
	else  printf("track not found.\n");
	gEDA->Redraw();
	if(eIDWindow!=NULL) eIDWindow->CloseWindow();
	eIDWindow = NULL;

}

void EdbEDAMainTab::MakeGUIIDSearch(){
	if(eIDWindow!=NULL) {
		eIDWindow->RaiseWindow();
		return;
	}
	TGTransientFrame *fMainFrame = new TGTransientFrame(gClient->GetRoot(),gEve?gEve->GetMainWindow():0,10,300,kMainFrame | kVerticalFrame);
	fMainFrame->SetLayoutBroken(kTRUE);

	eIDWindow = fMainFrame;

	// Search
	int posy=10;

	TGLabel *fLabel = new TGLabel(fMainFrame,"Search:     track  id = ");
	fLabel->MoveResize(10,posy,120,18);
	eGNESearchTrack = new TGNumberEntryField(fMainFrame, -1, 0,TGNumberFormat::kNESInteger);
	fMainFrame->AddFrame(eGNESearchTrack, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
	eGNESearchTrack->MoveResize(130,posy,60,18);
	TGTextButton *fb = new TGTextButton(fMainFrame,"Go");
	fb->MoveResize(200,posy+1,30,18);
	fMainFrame->AddFrame(fb, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
	fb->Connect("Clicked()","EdbEDAMainTab",this,"ExecTrackSearch()");

	posy+=20;
	fLabel = new TGLabel(fMainFrame,"pl = ");
	fLabel->MoveResize(18,posy,30,18);
	eGNESearchSegmentIPL = new TGNumberEntryField(fMainFrame, -1, 0,TGNumberFormat::kNESInteger);
	fMainFrame->AddFrame(eGNESearchSegmentIPL, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
	eGNESearchSegmentIPL->MoveResize(55,posy,30,18);

	fLabel = new TGLabel(fMainFrame,"iseg = ");
	fLabel->MoveResize(90,posy,35,18);
	eGNESearchSegmentISeg = new TGNumberEntryField(fMainFrame, -1, 0,TGNumberFormat::kNESInteger);
	fMainFrame->AddFrame(eGNESearchSegmentISeg, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
	eGNESearchSegmentISeg->MoveResize(130,posy,60,18);

	fb = new TGTextButton(fMainFrame,"Go");
	fb->MoveResize(200,posy+1,30,18);
	fMainFrame->AddFrame(fb, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
	fb->Connect("Clicked()","EdbEDAMainTab",this,"ExecSegmentSearch()");

	fMainFrame->SetMWMHints(kMWMDecorAll,
		kMWMFuncAll,
		kMWMInputModeless);
	fMainFrame->MapSubwindows();

	fMainFrame->Resize(fMainFrame->GetDefaultSize());
	fMainFrame->MapWindow();
	fMainFrame->Resize(250,70);

}


void EdbEDAMainTab::MakeGUICuts(TGGroupFrame *fGroupFrame1){
	int posy=20, posx, dx;
	TGLabel *fLabel;
	TGTextButton *fb;
	fLabel = new TGLabel(fGroupFrame1,"Nseg  >=");
	fLabel->MoveResize(10,posy,50,20);

	eGNENseg = new TGNumberEntry(fGroupFrame1, eTrackSet->GetNsegCut(),11,-1,TGNumberFormat::kNESInteger);
	eGNENseg->MoveResize(70,posy,50,20);

	posy+=20;
	fLabel = new TGLabel(fGroupFrame1," iPlate >= ");
	if(gEDA->Japanese()) fLabel->SetTitle(" iPlate <= ");
	fLabel->MoveResize(10,posy,50,20);
	
	if(eTrackSet->GetUpstreamPlate()==0) eTrackSet->SetUpstreamPlate( eTrackSet->FindUpstreamPlate());
	eGNEMaxPlate = new TGNumberEntry(fGroupFrame1, eTrackSet->GetUpstreamPlate(),11,-1,TGNumberFormat::kNESInteger);
	eGNEMaxPlate->MoveResize(70,posy,50,20);

	posy+=22;
	posx=5;
	fLabel = new TGLabel(fGroupFrame1,"ph>=");
	fLabel->MoveResize(posx,posy,dx=27,18);
	
	posx+=dx;
	eGNEWCut = new TGNumberEntryField(fGroupFrame1, -1, (Double_t) eTrackSet->GetPHCut());
	eGNEWCut->SetFormat(TGNumberFormat::kNESRealOne);
	eGNEWCut->MoveResize(posx,posy,dx=33,18);
	eGNEWCut->SetToolTipText("Cut for PH mean.");
	
	posx+=dx+5;
	fLabel = new TGLabel(fGroupFrame1,"dt/ph=");
	fLabel->MoveResize(posx,posy,dx=40,18);

	posx+=dx;
	eGNEPHDTRMS = new TGNumberEntryField(fGroupFrame1, -1, (Double_t) eTrackSet->GetPHDTRMS());
	eGNEPHDTRMS->SetFormat(TGNumberFormat::kNESRealThree);
	eGNEPHDTRMS->SetLimits(TGNumberFormat::kNELLimitMin, 0);
	eGNEPHDTRMS->MoveResize(posx,posy,dx=40,18);

	
	posx+=dx+8;
	eGClearPrevious = new TGCheckButton(fGroupFrame1,"Clear prev");
	eGClearPrevious->MoveResize(posx,posy,dx=80,20);
	eGClearPrevious->SetState(kButtonDown);
	
	posy+=20;
	eGSideOut = new TGCheckButton(fGroupFrame1,"Penetration and Side-out");
	eGSideOut->MoveResize(10,posy,160,20);
	eGSideOut->SetToolTipText(	"Reject penetrating and side-out tracks,\n"
								"by requesting the scanning area upstream\n"
								"of the track for a given number of plates.");

	eGSideOutPlate = new TGNumberEntry(fGroupFrame1, eTrackSet->GetSideOutPlate(),11,-1,TGNumberFormat::kNESInteger);
	eGSideOutPlate->MoveResize(175,posy,35,18);

	fLabel = new TGLabel(fGroupFrame1,"pl");
	fLabel->MoveResize(210,posy,dx=25,20);

	posy+=23;
	// Check Neighbor search
	eGCBNeighbor = new TGCheckButton(fGroupFrame1,"Neighbor search");
	eGCBNeighbor->MoveResize(10,posy,115,20);
	eGCBNeighbor->SetToolTipText("Search neighbors.\nAngular and position tolerance are valid.");
	

	eGCBIpCut = new TGCheckButton(fGroupFrame1,"IP cut");
	eGCBIpCut->MoveResize(128,posy,100,20);
	eGCBIpCut->SetState(kButtonUp);
	eGCBIpCut->SetToolTipText("IP cut against a selected vertex (Red).\nPosition tolerance is valid.");
	

	posy+=20;
	// Check Angular cut
	eGCBAngleCut = new TGCheckButton(fGroupFrame1,"Angular cut");
	eGCBAngleCut->MoveResize(10,posy,90,20);
	eGCBAngleCut->SetToolTipText("Angular cut against a given angle.\nAngular tolerance is valid.");

	eGNEAX = new TGNumberEntryField(fGroupFrame1);
	eGNEAX->MoveResize(100,posy,60,18);
	eGNEAX->SetFormat(TGNumberFormat::kNESRealThree);

	eGNEAY = new TGNumberEntryField(fGroupFrame1);
	eGNEAY->MoveResize(170,posy,60,18);
	eGNEAY->SetFormat(TGNumberFormat::kNESRealThree);

	posy+=25;

	fLabel = new TGLabel(fGroupFrame1,"tol :");
	fLabel->MoveResize(10,posy,25,20);

	eGNETolDT = new TGNumberEntry(fGroupFrame1, 0.05,11,-1,TGNumberFormat::kNESRealTwo);
	eGNETolDT->MoveResize(35,posy,60,20);
//	eGNETolDT->SetToolTipText("Angular tolerance");

	fLabel = new TGLabel(fGroupFrame1,"rad");
	fLabel->MoveResize(100,posy,25,18);

	eGNETolDX = new TGNumberEntry(fGroupFrame1, 250,11,-1,TGNumberFormat::kNESReal);
	eGNETolDX->MoveResize(130,posy,60,20);
//	eGNETolDX->SetToolTipText("Position tolerance");

	fLabel = new TGLabel(fGroupFrame1,"micron");
	fLabel->MoveResize(190,posy,40,18);

	posy+=25;
	fLabel = new TGLabel(fGroupFrame1,"dz range :  up");
	fLabel->MoveResize(12,posy,75,18);
	eGNENeighborDzUp = new TGNumberEntryField(fGroupFrame1, -1,5000);
	eGNENeighborDzUp->MoveResize(90,posy,50,18);
	eGNENeighborDzUp->SetToolTipText("dZ range for neighbor search");


	fLabel = new TGLabel(fGroupFrame1,"down");
	fLabel->MoveResize(145,posy,30,18);
	eGNENeighborDzDown = new TGNumberEntryField(fGroupFrame1, -1,5000);
	eGNENeighborDzDown->MoveResize(180,posy,50,18);
	eGNENeighborDzDown->SetToolTipText("dZ range for neighbor search");


	// Apply
	fb = new TGTextButton(fGroupFrame1,"Apply");
	fb->MoveResize(128,20,74,40);
	fb->Connect("Clicked()","EdbEDAMainTab",this, "ExecCuts()");

	fb = new TGTextButton(fGroupFrame1,"clr");
	fb->MoveResize(205,20,30,20);
	fb->Connect("Clicked()","EdbEDAMainTab",this,"ClearTracks()");
	
	fb = new TGTextButton(fGroupFrame1,"<");
	fb->MoveResize(205,40,30,20);
	fb->Connect("Clicked()","EdbEDAMainTab",this,"Redo()");
	
	posy+=23;
	fb = new TGTextButton(fGroupFrame1,"ID Search");
	fb->MoveResize(5,posy,80,18);
	fb->Connect("Clicked()","EdbEDAMainTab",this,"MakeGUIIDSearch()");

	fb = new TGTextButton(fGroupFrame1,"MT search");
	fb->MoveResize(130,posy,70,18);
	fb->Connect("Clicked()","EdbEDAMainTab", this, "MicroTrackSearch()");

	
	posy+=23;
	fLabel = new TGLabel(fGroupFrame1,"Search BT : pl");
	fLabel->MoveResize(5,posy,80,18);

	eGNESearchBTIPL = new TGNumberEntry(fGroupFrame1, -1, 0,TGNumberFormat::kNESInteger);
	eGNESearchBTIPL->MoveResize(90,posy,35,18);

	fb = new TGTextButton(fGroupFrame1,"Dmin");
	fb->MoveResize(130,posy,35,18);
	fb->Connect("Clicked()","EdbEDAMainTab", this, "SearchCouplesDmin()");
	fb->SetToolTipText("Search BT by Neibor search. \nif ipl==-1, search in all plate.\nPHcut is valid.");

	fb = new TGTextButton(fGroupFrame1,"IP cut");
	fb->MoveResize(170,posy,40,18);
	fb->Connect("Clicked()","EdbEDAMainTab", this, "SearchCouplesIP()");
	fb->SetToolTipText("Search BT by IP cut. \nif ipl==-1, search in all plate.\nPHcut is valid.");

	fb = new TGTextButton(fGroupFrame1,"clr");
	fb->MoveResize(215,posy,20,18);
	fb->Connect("Clicked()","EdbEDAMainTab", this, "ClearTracks()");

	
}

void EdbEDAMainTab::MakeGUIMain(){

	eTrackSet = gEDA->GetTrackSet("TS");
	
	
	TEveBrowser * browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kLeft);

	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
	frmMain->SetWindowName("XX GUI");
	frmMain->SetCleanup(kDeepCleanup);

	//TGLabel *fLabel;
	TGTextButton *fb;

	TGVerticalFrame* hf = new TGVerticalFrame(frmMain);
	{
	// "fGroupFrame600" group frame
		TGGroupFrame *fGroupFrame1 = new TGGroupFrame(hf,"Cuts");
		fGroupFrame1->SetLayoutBroken(kTRUE);
		MakeGUICuts(fGroupFrame1);
		hf->AddFrame(fGroupFrame1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
		fGroupFrame1->MoveResize(0,0,240,255);


	// Group2
		TGGroupFrame *fGroupFrame2 = new TGGroupFrame(hf,"Selection");
		fGroupFrame2->SetLayoutBroken(kTRUE);
		{
			int posy = 15;
			TGLabel *fLabel = new TGLabel(fGroupFrame2,"method ");
			fLabel->MoveResize(10,posy,50,18);

			eGRBSelSeg = new TGRadioButton(fGroupFrame2,"segment");
			eGRBSelSeg->MoveResize(80,posy,80,18);
			eGRBSelSeg->Connect("Clicked()","EdbEDAMainTab", this, "SetPickSegment()");
			eGRBSelSeg->SetState((EButtonState)1);
			
			eGRBSelTrack = new TGRadioButton(fGroupFrame2,"track");
			eGRBSelTrack->MoveResize(170,posy,65,18);
			eGRBSelTrack->Connect("Clicked()","EdbEDAMainTab", this, "SetPickTrack()");

			posy+=23;
			fb = new TGTextButton(fGroupFrame2,"Only selected");
			fb->MoveResize(10,posy,90,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"DrawOnlySelected()");

			fb = new TGTextButton(fGroupFrame2,"Remove");
			fb->MoveResize(105,posy,60,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"RemoveSelected()");
			fb->SetToolTipText("Remove selected tracks.");

			fb = new TGTextButton(fGroupFrame2,"All");
			fb->MoveResize(170,posy,60,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"SelectAll()");
			fb->SetToolTipText("Select all tracks shown in the display.");

			posy+=23;
			eGCBExtention = new TGCheckButton(fGroupFrame2,"Extention");
			eGCBExtention->MoveResize(15,posy,80,18);
			eGCBExtention->Connect("Clicked()","EdbEDAMainTab", this,"SetExtention()");

			fb = new TGTextButton(fGroupFrame2,"res");
			fb->MoveResize(100,posy,40,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"DrawExtentionResolution()");
			fb->SetToolTipText("Draw extention with measurement error + scattering.");

			fLabel = new TGLabel(fGroupFrame2,"z =");
			fLabel->MoveResize(145,posy,20,18);
			eGNEExtentionZ =  new TGNumberEntryField(fGroupFrame2, -1, gEDA->GetExtentionSet()->GetZ(), TGNumberFormat::kNESInteger);
			eGNEExtentionZ -> MoveResize(170,posy,55,18);
			eGNEExtentionZ -> Connect("TextChanged(char *)","EdbEDAMainTab", this,"SetExtentionZNumber()");
			
			
			posy+=20;
			eGSliderExtentionZ = new TGHSlider(fGroupFrame2,220,kSlider1 | kScaleBoth,-1,kHorizontalFrame);
			eGSliderExtentionZ->SetRange(-80000,0);
			eGSliderExtentionZ->SetPosition((int)gEDA->GetExtentionSet()->GetZ());
			eGSliderExtentionZ->MoveResize(5,posy,225,15);
			eGSliderExtentionZ->Connect("PositionChanged(int)","EdbEDAMainTab", this,"SetExtentionZSlider()");
		
		}
		hf->AddFrame(fGroupFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
		fGroupFrame2->MoveResize(0,200,240,105);
		

	// Group3
		TGGroupFrame *fGroupFrame3 = new TGGroupFrame(hf,"Analysis");
		fGroupFrame3->SetLayoutBroken(kTRUE);
		{
			fb = new TGTextButton(fGroupFrame3,"Min distance");
			fb->MoveResize(10,20,110,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"CalcMinimumDistance()");

			fb = new TGTextButton(fGroupFrame3,"Vertex");
			fb->MoveResize(10,40,68,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"CalcVertex()");
			fb->SetToolTipText("Calc vertex using only the selected segments.");

			fb = new TGTextButton(fGroupFrame3,"+");
			fb->MoveResize(78,40,18,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"CalcVertexPlus()");
			fb->SetToolTipText("Calc vertex without clearing existing vertices.");

			fb = new TGTextButton(fGroupFrame3,"clr");
			fb->MoveResize(98,40,22,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"ClearVertices()");
			
			fb = new TGTextButton(fGroupFrame3,"Momentum");
			fb->MoveResize(125,20,110,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"CalcPSelected()");
			
			fb = new TGTextButton(fGroupFrame3,"Calc IP");
			fb->MoveResize(125,40,110,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"CalcIP()");
		}
		hf->AddFrame(fGroupFrame3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
		fGroupFrame3->MoveResize(0,210,240,72);

	// Group4
		TGGroupFrame *fGroupFrame4 = new TGGroupFrame(hf,"Display");
		fGroupFrame4->SetLayoutBroken(kTRUE);
		{
			int posy=20;
//			TGLabel *fLabel = new TGLabel(fGroupFrame4,"Camera");
//			fLabel->MoveResize(10,posy,50,20);
			int posx = 10;
			int dx = 20;
			fb = new TGTextButton(fGroupFrame4,"X", EdbEDACamera::X); 
			fb->MoveResize(posx,posy,dx,18);
			//fb->SetCommand("gEDA->SetCamera(kEDA_CAMERA_X)");
			fb->Connect("Clicked()","EdbEDAMainTab", this,"SetProjection()");

			posx+=dx+5;
			fb = new TGTextButton(fGroupFrame4,"Y", EdbEDACamera::Y);
			fb->MoveResize(posx,posy,dx,18);
//			fb->SetCommand("gEDA->SetCamera(kEDA_CAMERA_Y)");
			fb->Connect("Clicked()","EdbEDAMainTab", this,"SetProjection()");

			posx+=dx+5;
			fb = new TGTextButton(fGroupFrame4,"Z", EdbEDACamera::Z);
			fb->MoveResize(posx,posy,dx,18);
//			fb->SetCommand("gEDA->SetCamera(kEDA_CAMERA_Z)");
			fb->Connect("Clicked()","EdbEDAMainTab", this,"SetProjection()");

			posx+=dx+5;
			fb = new TGTextButton(fGroupFrame4,"Beam", EdbEDACamera::B);
			fb->MoveResize(posx,posy,dx=50,18);
//			fb->SetCommand("gEDA->SetCamera(kEDA_CAMERA_B)");
			fb->Connect("Clicked()","EdbEDAMainTab", this,"SetProjection()");

			posx+=dx+5;
			fb = new TGTextButton(fGroupFrame4,"Pick center");
			fb->MoveResize(posx,posy,dx=80,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this,"PickCenter()");

			posy+=20;

			TGLabel *fLabel = new TGLabel(fGroupFrame4,"Track Color");
			fLabel->MoveResize(10,posy,70,18);
			//fb = new TGTextButton(fGroupFrame4,"Track Color");
			//fb->SetCommand("gEDA->ShowColorTable()");
			//fb->MoveResize(10,posy,70,18);
			
			eGRBColorPL = new TGRadioButton(fGroupFrame4,"Plate", kCOLOR_BY_PLATE);
			eGRBColorPL->MoveResize(85,posy,50,18);
//			eGRBColorPL->SetCommand("gEDA->SetColorMode(kCOLOR_BY_PLATE)");
			eGRBColorPL->Connect("Clicked()","EdbEDAMainTab", this, "SetColorMode()");
			
			eGRBColorPH = new TGRadioButton(fGroupFrame4,"PH", kCOLOR_BY_PH);
			eGRBColorPH->MoveResize(140,posy,40,18);
			//eGRBColorPH->SetCommand("gEDA->SetColorMode(kCOLOR_BY_PH)");
			eGRBColorPH->Connect("Clicked()","EdbEDAMainTab", this, "SetColorMode()");
			eGRBColorPH->SetState((EButtonState)1);

			eGRBColorBK = new TGRadioButton(fGroupFrame4,"Black", kBLACKWHITE);
			eGRBColorBK->MoveResize(185,posy,50,18);
//			eGRBColorBK->SetCommand("gEDA->SetColorMode(kBLACKWHITE)");
			eGRBColorBK->Connect("Clicked()","EdbEDAMainTab", this, "SetColorMode()");

			posy+=22;
			dx=37;
			
			eDrawChecks = new TObjArray;
			posx=10;
			dx = 45;
			for(int i=0;i<gEDA->NTrackSets();i++){
//				posx+=dx+5;
				if(i==5||i==10||i==15){
					posx=10;
					posy+=20;
				}
				EdbEDATrackSet *set = gEDA->GetTrackSet(i);
				TGCheckButton *cb = new TGCheckButton(fGroupFrame4,set->GetName(),i);
				eDrawChecks->Add(cb);
				cb->MoveResize(posx,posy,dx,20);
				cb->SetState((EButtonState)1);
				cb->Connect("Toggled(Bool_t)","EdbEDAMainTab", this, "SetDrawTrackSets(Bool_t)");
				posx+=dx;
			}

			posy+=22;
			posx = 10;
			eGCBDrawArea = new TGCheckButton(fGroupFrame4,"Scan area");
			eGCBDrawArea->MoveResize(posx,posy,dx=85,20);
			eGCBDrawArea->SetState((EButtonState)1);
			eGCBDrawArea->Connect("Clicked()","EdbEDAMainTab", this,"SetDrawAreas()");

			posx+=dx+10;
			eGCBTrackID = new TGCheckButton(fGroupFrame4,"track ID");
			eGCBTrackID->MoveResize(posx,posy,dx=65,20);
			eGCBTrackID->SetState((EButtonState)0);
			eGCBTrackID->Connect("Clicked()","EdbEDAMainTab", this,"SetDrawText()");

			posx+=dx+5;
			eGCBTrackAngle = new TGCheckButton(fGroupFrame4,"angle");
			eGCBTrackAngle->MoveResize(posx,posy,dx=60,20);
			eGCBTrackAngle->SetState((EButtonState)0);
			eGCBTrackAngle->Connect("Clicked()","EdbEDAMainTab", this,"SetDrawText()");


			posy+=22;
			eGCBAnimation = new TGCheckButton(fGroupFrame4,"Animation");
			eGCBAnimation->MoveResize(10,posy,80,20);
			//fb->SetCommand("gEDA->SetCamera(EDA_CAMERA_B)");
			eGCBAnimation->Connect("Clicked()","EdbEDAMainTab", this,"SetAnimation()");

			fb = new TGTextButton(fGroupFrame4,"op");
			fb->MoveResize(90,posy,20,18);
			fb->Connect("Clicked()", "EdbEDA", gEDA, "SetAnimationMode()");

			fb = new TGTextButton(fGroupFrame4,"Save animation");
			fb->MoveResize(125,posy,100,18);
			fb->Connect("Clicked()", "EdbEDA", gEDA, "SaveAnimation()");

			posy+=21;
			fb = new TGTextButton(fGroupFrame4,"Snapshot");
			fb->MoveResize(10,posy,100,18);
			fb->Connect("Clicked()", "EdbEDA", gEDA, "Snapshot()");

			fb = new TGTextButton(fGroupFrame4,"Save projections");
			fb->MoveResize(125,posy,100,18);
			fb->Connect("Clicked()", "EdbEDA", gEDA, "SavePictures()");

			posy+=21;
			fb = new TGTextButton(fGroupFrame4,"Write track list");
			fb->MoveResize(10,posy,100,18);
			fb->Connect("Clicked()","EdbEDAMainTab",this,"WriteListFile()");

			fb = new TGTextButton(fGroupFrame4,"Read track list");
			fb->MoveResize(125,posy,100,18);
			fb->Connect("Clicked()","EdbEDAMainTab",this,"ReadListFile()");
/*
			posy+=21;
			fb = new TGTextButton(fGroupFrame4,"Write filtered text");
			fb->MoveResize(10,posy,100,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this, "WriteFilteredText()");
			

			fb = new TGTextButton(fGroupFrame4,"Read filtered text");
			fb->MoveResize(125,posy,100,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this, "ReadFilteredText()");
*/
			posy+=21;
			fb = new TGTextButton(fGroupFrame4,"Write feedback");
			fb->MoveResize(10,posy,90,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this, "WriteFeedback()");

			fb = new TGTextButton(fGroupFrame4,"op");
			fb->MoveResize(102,posy,20,18);
			fb->Connect("Clicked()", "EdbEDAMainTab", this, "SetOutputFileModeButton()");

			fb = new TGTextButton(fGroupFrame4,"Read feedback");
			fb->MoveResize(125,posy,100,18);
			fb->Connect("Clicked()","EdbEDAMainTab", this, "ReadFeedback()");
			
			if(gEDA->Japanese()){
				
				posy+=21;
				fb = new TGTextButton(fGroupFrame4,"Write Mxx.all");
				fb->MoveResize(10,posy,90,18);
				fb->Connect("Clicked()","EdbEDAMainTab", this, "WriteMxxFile()");
				
			}
			
			

		}
		hf->AddFrame(fGroupFrame4, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
		fGroupFrame4->MoveResize(0,260,240,240);

	}

	frmMain->AddFrame(hf);
	
	frmMain->MapSubwindows();
	frmMain->Resize();
	frmMain->MapWindow();

	browser->StopEmbedding();
	browser->SetTabTitle("Main", 0);
}



void EdbEDAMainTab::CalcMinimumDistance(void){
	printf("/// Calc Dmin //////////////////////////////////////////////////////////////////////////\n");
	TObjArray * selected = gEDA->GetSelected();
	if(selected->GetEntriesFast()<2) {
		printf("Select more than 2 segment!!\n");
		return;
	}
	int i,j;
	int n = selected->GetEntriesFast();

	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			int id1=0,id2=0;
			EdbSegP *seg1 = (EdbSegP*) selected->At(i);
			if(IsSegment(seg1)) {
				id1=seg1->Track();
				if(id1==-1) id1=seg1->ID();
			}
			else id1=seg1->ID();
			
			EdbSegP *seg2 = (EdbSegP*) selected->At(j);
			
			if(IsSegment(seg2)) {
				id2=seg2->Track();
				if(id2==-1) id2=seg2->ID();
			}
			else id2=seg2->ID();
			
			double dz;
			double r = CalcDmin(seg1,seg2, &dz);
			
			printf(" %5d x %5d dmin = %6.1lf um at %6.1lf um upstream from pl %2d. dT=(%7.4lf, %7.4lf).\n", 
					id1,id2, r, dz, seg1->Plate(), seg1->TX()-seg2->TX(), seg1->TY()-seg2->TY());

			gEve->SetStatusLine(Form("min distance = %6.1lf um at %6.1lf um upstream of pl %2d\n", r, dz, seg1->Plate()));
		}
	}
	printf("////////////////////////////////////////////////////////////////////////////////////////\n");
}

void EdbEDAMainTab::CalcVertexPlus(){ CalcVertex(kFALSE);}

void EdbEDAMainTab::CalcVertex(int clear_previous_local){
	// calculate a vertex with the selected tracks (segments).
	// use the data of (the selected object)->X(), Y(), Z(), TX(), TY().
	// means, if the selected object == segment, use the data of the segment. or it == track, the use the fitted data.
	printf("/// Calc Vertex /////////////////////////////////////////////////\n");
	TObjArray *selected = gEDA->GetSelected();
	
	if(selected->GetEntriesFast()<1){
		ErrorMessage("CalcVertex() error : please select more than 1 segment or track\n");
		return;
	}
	
	if(eGClearPrevious->IsOn()&&clear_previous_local) gEDA->GetVertexSet()->Clear();
	
	
	/*
	// this is for FEDRA vertex using libVt++. but it's not convenient. use only topological.
	
	EdbVertex *v0 = EdbEDAUtil::CalcVertex(selected); // by Tomoo func. pre-vertexing. to know vertex z.
	
	EdbVertex *v;
	if(selected->GetEntriesFast()==1){
		EdbSegP *s  = (EdbSegP *)selected->At(0);
		EdbTrackP *t = gEDA->GetSelectedTrack(s);
		
		EdbVTA *vta = new EdbVTA(t, v0);
		vta->SetFlag(2);
		vta->SetImp(EdbEDAUtil::CalcIP(s,v0));
		v0->AddVTA(vta);
		v=v0;
	}
	else {
		// Calc vertex by FEDRA version.
		v = new EdbVertex;
		// register the tracks to vertex.
		for(int i=0;i<selected->GetEntriesFast();i++){
			EdbSegP *s  = (EdbSegP *)selected->At(i);
			EdbTrackP *t = gEDA->GetSelectedTrack(s);
			
			//if(t==NULL) t = new EdbTrackP(s);
			
			double dz1 = abs(t->GetSegmentFirst()->Z()-v0->Z()); // z distance from 1st segment to vertex
			double dz2 = abs(t->GetSegmentLast()->Z()-v0->Z());  // z distance from last segment to vertex
			
			EdbVTA *vta = new EdbVTA(t, v);
			vta->SetZpos(dz1<dz2?1:0);    // if 1st segment is near to vertex, use 1st segment.
			vta->SetFlag(2);
			vta->SetImp(EdbEDAUtil::CalcIP(s,v));
			v->AddVTA(vta);
		}
		// calculate in EdbVertexRec.
		EdbVertexRec VR;
		VR.eImpMax=1e6;
		EdbPVRec PVR; // needed because of X0.
		PVR.SetScanCond(new EdbScanCond);
		VR.SetPVRec(&PVR);
		VR.eUseSegPar=kTRUE; // definitely need for better IP calculation.
		VR.MakeV(*v);
	}
	*/
	// calc vertex point by Tomoko // VTA is not set.
	// old version. upto 2009 10 29
	EdbVertex *v = EdbEDAUtil::CalcVertex(selected);

	static int vid=0;
	if(vid<gEDA->GetVertexSet()->N()) vid=gEDA->GetVertexSet()->N();
	v->SetID(vid++);
	// register the vertex to selected_vertex.
	for(int i=0;i<selected->GetEntriesFast();i++){
		EdbSegP *s  = (EdbSegP *)selected->At(i);
		EdbTrackP *t = gEDA->GetSelectedTrack(s);
		
		if(t==NULL) t = new EdbTrackP(s);
		
		EdbVTA *vta = new EdbVTA(t, v);
		vta->SetZpos(1);
		vta->SetFlag(2);
		vta->SetImp(EdbEDAUtil::CalcIP(s,v));
		v->AddVTA(vta);
	}
	
	// claculate the dZ from the last plate
	double vertexplz = 1e6;
	int vertexipl = 0;
	for(int i=0;i<gEDA->GetNPL();i++){
		if(gEDA->GetZPID(i)<v->Z()) continue;
		if(gEDA->GetZPID(i)<vertexplz) {
			vertexplz = gEDA->GetZPID(i);
			vertexipl = gEDA->GetIPL(i);
		}
	}
	if(gEDA->GetNPL()==0){
		EdbTrackP *t = gEDA->GetSelectedTrack(0);
		EdbSegP *s = t->GetSegmentFirst();
		vertexplz=s->Z();
		vertexipl=s->Plate();
	}

	gEve->SetStatusLine(Form("Vertex position (x,y,z) = (%8.1lf, %8.1lf, %8.1lf). %.1lf micron upstream of pl %d\n", 
							v->X(), v->Y(), v->Z(), vertexplz-v->Z(), vertexipl));
	
	printf(" --------------------------------------------------------------------------\n");
	printf("Vertex position (x,y,z) = (%8.1lf, %8.1lf, %8.1lf). %.1lf micron upstream of pl %d\n",v->X(), v->Y(), v->Z(), vertexplz-v->Z(), vertexipl);
	
	
	gEDA->WriteFeedbackFile2008(v,"tmp.feedback");
	printf("/// Calc Vertex end//////////////////////////////////////////////\n");
	gEDA->GetVertexSet()->AddVertex(v);
	gEDA->SetSelectedVertex(v);
	gEDA->Redraw();
}

void EdbEDAMainTab::CalcPSelected(){
	printf("/// Calc Momentum, temporary version//////////////////////\n");
	TObjArray *selected = gEDA->GetSelected();
	for(int i=0; i<selected->GetEntriesFast(); i++ ){
		EdbTrackP *t=gEDA->GetTrack((EdbSegP *)selected->At(i));
		if(t==NULL) continue;
		double p, pmin, pmax;

//		CalcP(t,p,pmin,pmax);
//		printf(" %5d %7.4lf %7.4lf p pmin pmax = %6.2f %6.2f %6.2f GeV/c (90%%CL)\n", 
//				t->ID(), t->TX(), t->TY(), p, pmin, pmax);
//		gEve->SetStatusLine(Form(" %5d %7.4lf %7.4lf p pmin pmax = %6.2f %6.2f %6.2f GeV/c (90%%CL)", 
//				t->ID(), t->TX(), t->TY(), p, pmin, pmax));

//		printf("Bern conventional (Magali's)\n");
		CalcP(t,p,pmin,pmax);
//		printf("Bern %5d %7.4lf %7.4lf p pmin pmax = %6.2f %6.2f %6.2f GeV/c (90%%CL)\n", 
//				t->ID(), t->TX(), t->TY(), p, pmin, pmax);
		gEve->SetStatusLine(Form(" %5d %7.4lf %7.4lf p pmin pmax = %6.2f %6.2f %6.2f GeV/c (90%%CL)", 
				t->ID(), t->TX(), t->TY(), p, pmin, pmax));
/*
		printf("Andrea conventional\n");
		CalcP2(t,p,pmin,pmax);
		printf("A_g0 %5d %7.4lf %7.4lf p pmin pmax = %6.2f %6.2f %6.2f GeV/c (90%%CL)\n", 
				t->ID(), t->TX(), t->TY(), p, pmin, pmax);

		printf("Andrea new statistical treatment\n");
		CalcP3(t,p,pmin,pmax);
		printf("A_g2 %5d %7.4lf %7.4lf p pmin pmax = %6.2f %6.2f %6.2f GeV/c (90%%CL)\n", 
				t->ID(), t->TX(), t->TY(), p, pmin, pmax);
*/
	}
	
	printf("////////////////////////////////////////////////////////////////\n");
}


void EdbEDAMainTab::CalcIP(){
	// calculate IP for the selected tracks wrt the given vertex.
	// if the vertex is not given, use the selected vertex.
	EdbVertex *v=gEDA->GetSelectedVertex();
	if(NULL==v) {
		printf("select a vertex!\n");
		return;
	}
	printf("/// Calc IP w.r.t. Vertex %d %8.1lf %8.1lf %8.1lf (red vertex) ///\n",v->ID(), v->X(),v->Y(),v->Z());
	
	TObjArray *selected = gEDA->GetSelected();
	if(selected->GetEntriesFast()==0){
		printf("select segments\n");
		return;
	}
	for(int i=0;i<selected->GetEntriesFast();i++){
		EdbSegP *s = (EdbSegP *) selected->At(i);
		double r = EdbEDAUtil::CalcIP(s,v->X(),v->Y(),v->Z());
		
		double dz = fabs(v->Z()-s->Z());
		double minimum_kink_angle = r / (dz-293.0/2);
		
		printf("%6d PL %2d %7.4lf %7.4lf IP = %6.2lf, dZ = %6.1lf, Minimum kink angle = %6.4lf\n", 
						s->Track(), s->Plate(), s->TX(), s->TY(), r, dz, minimum_kink_angle);
		gEve->SetStatusLine(Form("%4d PL %2d %7.4lf %7.4lf IP = %6.2lf, dZ = %6.1lf, Minimum kink angle = %6.4lf\n", 
						s->Track(), s->Plate(), s->TX(), s->TY(), r, dz, minimum_kink_angle));
	}
}


void EdbEDAMainTab::DrawOnlySelected() {
	gEDA->StorePrevious();
	eTrackSet->ClearTracks();
	eTrackSet->SetTracks(gEDA->GetSelectedTracks());
	gEDA->GetVertexSet()->ClearVertices();
	gEDA->GetVertexSet()->SetVertex(gEDA->GetSelectedVertex());
	gEDA->Redraw();
}
void EdbEDAMainTab::RemoveSelected(){
	gEDA->StorePrevious();
	gEDA->RemoveSelected();
	gEDA->Redraw();
}
void EdbEDAMainTab::SelectAll(){
	gEDA->SelectAll();
	gEDA->Redraw();
}

void EdbEDAMainTab::SetExtention(){
	int extention=(int) eGCBExtention->IsOn();
	if(extention){
		gEDA->GetExtentionSet()->Enable();
		gEDA->GetExtentionSet()->Draw(kTRUE);
	} else {
		gEDA->GetExtentionSet()->Disable();
		gEDA->GetExtentionSet()->UnDraw(kTRUE);
	}
	printf("extention %s\n", eGCBExtention->IsOn()? "On": "Off");
}
void EdbEDAMainTab::SetDrawText() {
	if(eGCBTrackID->IsOn()) gEDA->SetDrawTrackID(kTRUE);
	else gEDA->SetDrawTrackID(kFALSE);
	
	if(eGCBTrackAngle->IsOn()) gEDA->SetDrawTrackAngle(kTRUE);
	else gEDA->SetDrawTrackAngle(kFALSE);
	
	gEDA->Redraw();
}

void EdbEDAMainTab::SetDrawTrackSets(Bool_t on) {
	// Control Draw option for each Track sets.

	TGButton *btn = (TGButton *) gTQSender;
	int id = btn->WidgetId();  // id is ordinal of TrackSets. TS=0, SB=1,,,,

	TGCheckButton *cb = (TGCheckButton *) eDrawChecks->At(id);
	EdbEDATrackSet *set = gEDA->GetTrackSet(id);
	if(cb==NULL||set==NULL) return;

	if(set->NBase()) set->SetDraw(cb->GetState()?kTRUE:kFALSE); 
	
	if(set->NBase()==0) {
		gEve->GetBrowser()->GetTabBottom()->SetTab(0); // set TrackSetTab in front.
	}

	gEDA->GetTrackSetTab()->Select(id);
	gEDA->GetTrackSetTab()->Selected(set->GetName());

/*	
	for(int i=0;i<gEDA->NTrackSets();i++){
		TGCheckButton *cb = (TGCheckButton *) eDrawChecks->At(i);
		EdbEDATrackSet *set = gEDA->GetTrackSet(i);
		if(cb==NULL) continue;
		set->SetDraw(cb->GetState()?kTRUE:kFALSE); 
	}
*/	
	gEDA->Redraw();
}





void EdbEDAMainTab::SetDrawCheckBox(){
	for(int i=0;i<gEDA->NTrackSets();i++){
		TGCheckButton *cb = (TGCheckButton *) eDrawChecks->At(i);
		EdbEDATrackSet *set = gEDA->GetTrackSet(i);
		if(cb==NULL) continue;
		//cb->SetState((EButtonState) set->GetDraw()); 
		
		cb->SetState(set->GetDraw()&&set->NBase()?kButtonDown:kButtonUp);
		
		//if(set->NBase()>0){
		//	cb->SetEnabled(kTRUE);
		//	cb->SetState(set->GetDraw()?kButtonDown:kButtonUp);
		//}
		//else cb->SetEnabled(kFALSE);
	}
}

void EdbEDAMainTab::SetExtentionZSlider(){
	eGNEExtentionZ->SetNumber(eGSliderExtentionZ->GetPosition());
}
void EdbEDAMainTab::SetExtentionZNumber(){
		double z = eGNEExtentionZ ->GetNumber();
		eGSliderExtentionZ->SetPosition((int) z);
		gEDA->GetExtentionSet()->SetZ(z);
		gEDA->GetExtentionSet()->Draw(kTRUE);
	}
	
void EdbEDAMainTab::DrawExtentionResolution(){
	// Draw Extention with resolution
	gEDA->GetExtentionSet()->DrawResolution(kTRUE); 
}

void EdbEDAMainTab::SetProjection(){
	// Set Camera Projection
	TGButton *btn = (TGButton *) gTQSender;
	int id = btn->WidgetId();

	gEDA->SetCamera(id);
}



void EdbEDAMainTab::SetDrawAreas() {
	if(eGCBDrawArea->IsOn()) gEDA->SetDrawAreaSet(kTRUE);
	else gEDA->SetDrawAreaSet(kFALSE);
	//gEDA->Redraw();
	gEDA->Redraw();
}


void EdbEDAMainTab::SetPickSegment(){
	printf("selection mode : Segment\n");
	gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_Element);
	eGRBSelSeg->SetState((EButtonState)1);
	eGRBSelTrack->SetState((EButtonState)0);
}
void EdbEDAMainTab::SetPickTrack() {
	gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_Compound);
	eGRBSelSeg->SetState((EButtonState)0);
	eGRBSelTrack->SetState((EButtonState)1);
}

void EdbEDAMainTab::PickCenter(){
	TGLViewer *v = gEDA->GetGLViewer();
	v->PickCameraCenter();
	gEve->SetStatusLine("Pick a segment as the rotation center");
	printf("Pick a segment as the rotation center.\n");
}

void EdbEDAMainTab::SetAnimation(){ 
	if( eGCBAnimation->IsOn() ) gEDA->StartAnimation();
	else gEDA->StopAnimation();
}

void EdbEDAMainTab::SetColorMode(){
	TGButton *btn = (TGButton *) gTQSender;
	int mode = btn->WidgetId();
	if(mode== kCOLOR_BY_PH){
		eGRBColorPH->SetState(kButtonDown);
		eGRBColorPL->SetState(kButtonUp);
		eGRBColorBK->SetState(kButtonUp);
	}else if(mode== kCOLOR_BY_PLATE){
		eGRBColorPH->SetState(kButtonUp);
		eGRBColorPL->SetState(kButtonDown);
		eGRBColorBK->SetState(kButtonUp);
	}else if(mode==kBLACKWHITE) {
		eGRBColorPH->SetState(kButtonUp);
		eGRBColorPL->SetState(kButtonUp);
		eGRBColorBK->SetState(kButtonDown);
	}
	gEDA->SetColorMode(mode);
}

void EdbEDAMainTab::ReadListFile(const char *filename, bool clear_previous){
	// Read list file.
	// search tracks with plateID and segmentID of BT.
	// if corresponding track is not found in TS, it search in BT (couples file).
	//
	// previously EdbEDATrackSet::ReadListFile is used.
	// however, in order to allow to read BT, this has changed.
	// and now this func play a real role.
	
	gEDA->StorePrevious(); 
	
	// previously EdbEDATrackSet::ReadListFile is used.
	// eTrackSet->ReadListFile(NULL); 
	
	if(filename==NULL){
		TGFileInfo *fi = new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "List file", "*.lst", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDOpen, fi);
		filename=fi->fFilename;
	}

	EdbEDATrackSet *setTS = gEDA->GetTrackSet("TS");
	EdbEDATrackSet *setBT = gEDA->GetTrackSet("BT");
	printf("Read List file : %s\n", filename);

	FILE *fp=NULL;
	if(filename!=NULL) fp = fopen(filename,"rt");
	
	if(fp==NULL){
		printf("List file is not available.\n");
		return;
	}

	if(clear_previous) {
		setTS->ClearTracks();
		setTS->ClearComments();
		setBT->ClearTracks();
		setBT->ClearComments();
	}
	
	char buf[256];
	int ipl, iseg;
	char comment[256];
	// read list, setting TObjArray *selected.
	for(int i=0;NULL!=fgets(buf,sizeof(buf),fp);i++){
		if(feof(fp)||ferror(fp)) break;
		sscanf(buf,"%d %d %s", &ipl, &iseg, comment);
		
		EdbTrackP *t=NULL;
		for(int j=0;j<setTS->NBase();j++){
			EdbTrackP *tt = setTS->GetTrackBase(j);
			if(NULL==tt) continue;
			for(int k=0;k<tt->N();k++){
				EdbSegP *s = tt->GetSegment(k);
				if(s->ID()==iseg){
				if(s->Plate()==ipl){
					setTS->SetTrack(tt); // add tt for drawing
					t=tt;
					break;
				}}
			}
		}
		
		EdbEDATrackSet *set = NULL; // set for comments
		if(t) set = setTS;
		else  set = setBT;
		
		if(t==NULL) {
			// no track found in TS
			// search in BT
			
			EdbPattern *pat = gEDA->GetPatternIPL(ipl);
			for(int j=0; j<pat->N(); j++){
				EdbSegP *s = pat->GetSegment(j);
				if(s->ID()==iseg){
					t = new EdbTrackP(s);
					t->SetCounters();
					t->SetPlate(s->Plate());
					setBT->AddTrack(t);
				}
			}
			
		}
		
		
		if(t==NULL) continue;
		
		unsigned int l;
		for(l=0;l<strlen(buf);l++){
			if(buf[l]=='\"') break;
		}
		if(l<strlen(buf) ){
			strcpy(comment, buf+l+1);
			for(l=0;l<strlen(comment);l++){
				if(comment[l]=='\"') comment[l]='\0';
			}
			if(strlen(comment)!=0) {
				set->AddComment(t,comment);
			}
		}
	
	}

	
	
	
	
	
	
	
	
	
	gEDA->Redraw();
}

const char * EdbEDAMainTab::WriteListFile(const char *filename){
	filename = eTrackSet->WriteListFile(filename);
	
	EdbEDATrackSet *setBT = gEDA->GetTrackSet("BT");
	
	if (setBT->N()){
	     gEDA->GetTrackSet("BT")->WriteListFile(filename, kTRUE);
	}
	
	return filename;
	
}

void EdbEDAMainTab::WriteFilteredText(){
	gEDA->WriteFilteredText();
}
void EdbEDAMainTab::ReadFilteredText(){
	gEDA->ReadFilteredText();
	gEDA->Redraw();
}

void EdbEDAMainTab::SetOutputFileModeButton(){
	gEDA->eOutputFileMode = 
		InputNumberInteger(
			"Feeback output mode : 0 = with comment and true-track-id. 1: without.", 
			gEDA->eOutputFileMode);
}
void EdbEDAMainTab::WriteFeedback() {gEDA->WriteFeedbackFile();}
void EdbEDAMainTab::ReadFeedback()  {gEDA->ReadFeedbackFile();}

void EdbEDAMainTab::WriteMxxFile(){
	TObjArray *tracks = gEDA->GetCurrentTracks();
	WriteTracksMxx(tracks);
}

