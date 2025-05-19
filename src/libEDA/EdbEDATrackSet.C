#include"EdbEDA.h"
#include<TGFontDialog.h>
#include<TEveTextEditor.h>
extern EdbEDA *gEDA;

ClassImp(EdbEDATrackSet)
ClassImp(EdbEDATrackSelection)
ClassImp(EdbEDATrackSetTab)

int EdbEDATrackSelection::ImpactSearch (EdbTrackP *t){
	// select nearer segment. vertex or decay (2ry vertex)
	double dz1 = t->GetSegmentFirst()->Z()-eVertex->Z();
	double dz2 = t->GetSegmentLast()->Z()-eVertex->Z();
	EdbSegP *s = fabs(dz1)<fabs(dz2) ? t->GetSegmentFirst() :  t->GetSegmentLast(); // use nearer segment
	double ip = CalcIP(s, eVertex); 
	if( 0.0 < ip && ip < eTolDX) return 1; // ip<0 if the eVertex == NULL
	else return 0;
}

int EdbEDATrackSelection::Neighborhood (EdbSegP *s, double *dmin){
	int flag=0;
	double dminz;
	if(dmin) *dmin=1e5;
	for(int i=0;i<eSelected->GetEntriesFast();i++){
		EdbSegP *ss = (EdbSegP *)eSelected->At(i);
		if(ss==s) return 1;  
		
		if( fabs(s->TX()-ss->TX()) > eTolDTX) continue; // angular cut
		if( fabs(s->TY()-ss->TY()) > eTolDTY) continue;
		
		double d = CalcDmin( ss, s, &dminz);
		if(dmin) if(*dmin>d) *dmin=d;
//		if(d<10) printf("%lf %lf %lf %lf\n", d, dminz, eNeighborDzDown, eNeighborDzUp);

		if( d > eTolDX) continue;
		
		if( -eNeighborDzDown<dminz&&dminz<eNeighborDzUp) flag++;
		else if(dminz<-eNeighborDzDown){
			if( CalcDistance(ss, s, ss->Z()+eNeighborDzDown)<eTolDX) flag++;
		}
		else if(eNeighborDzUp<dminz){
			if( CalcDistance(ss, s, ss->Z()-eNeighborDzUp)<eTolDX) flag++;
		}
	}
	return flag;
}

int EdbEDATrackSelection::SideOut(EdbTrackP *t){
	// Return 1, if the track goes side-out.
	// - if require 3 scanning area upstream from the 1st segment.
	//EdbPVRec *pvr = gEDA->GetPVR();
	//if(pvr==NULL) return 0;
	
	int sideout=0;
	
	// Upstream Side-out check.
	int npl = 0;
	EdbSegP *s1 = t->GetSegmentFirst();
	EdbSegP s=*s1;
	
	// original. but very slow. EdbPointsBox2D::Xmax().....
	/*
	for(int pid=s1->PID()-1; pid>=0; pid--){
		EdbPattern *p = pvr->GetPattern(pid);
		if(p->N()==0) continue;
		s.PropagateTo(p->Z());
		if(s.X()<p->Xmin()) continue;
		if(s.X()>p->Xmax()) continue;
		if(s.Y()<p->Ymin()) continue;
		if(s.Y()>p->Ymax()) continue;
		
		npl++;
	}
	*/
	/*
	for(int pid=s1->PID()-1; pid>=0; pid--){
		EdbEDAArea *a = gEDA->GetAreaSet()->GetArea(pid);
		s.PropagateTo(a->Z());
		if(s.X()<a->Xmin()) continue;
		if(s.X()>a->Xmax()) continue;
		if(s.Y()<a->Ymin()) continue;
		if(s.Y()>a->Ymax()) continue;
		
		npl++;
	}
	*/
	for(int ipl=s1->Plate()-1; ipl>0; ipl--){
		EdbEDAArea *a = gEDA->GetTrackSet("TS")->GetAreaSet()->GetAreaIPL(ipl);
		if(a==NULL) continue;
		s.PropagateTo(a->Z());
		if(s.X()<a->Xmin()) continue;
		if(s.X()>a->Xmax()) continue;
		if(s.Y()<a->Ymin()) continue;
		if(s.Y()>a->Ymax()) continue;
		
		npl++;
	}
	
	
	if(npl<eSideOutPlate) sideout++;
	/*
	// Downstream Side-out check.
	// - if 0 or only 1 scanning area downstream from the last segment. if the propagated position is in the scanning area of the last plate, eceptionally doesn't cut the track.

	npl = 0;
	s1 = t->GetSegmentLast();
	s=*s1;
	if(s1->PID()==pvr->Npatterns()-1) npl=2; 
	for(int pid=s1->PID()+1; pid<pvr->Npatterns(); pid++){
		EdbPattern *p = pvr->GetPattern(pid);
		if(p->N()==0) continue;
		s.PropagateTo(p->Z());
		if(s.X()<p->Xmin()) continue;
		if(s.X()>p->Xmax()) continue;
		if(s.Y()<p->Ymin()) continue;
		if(s.Y()>p->Ymax()) continue;
		
		npl++;
		if(pid==pvr->Npatterns()-1) npl++; // if the propagated position is in the scanning area of last plate, eceptionally doesn't cut the track.
	}
	if(npl<=1) sideout++;
	*/
	return sideout;
}

void EdbEDATrackSelection::DoSelection(TObjArray *tracksbase, TObjArray *tracks){
	if(tracksbase==NULL||tracks==NULL) return;
	if(eClearPrevious) tracks->Clear();
	
	if(eImpactSearch&&eVertex==NULL){
		eVertex=gEDA->GetSelectedVertex();
		if(eVertex==NULL) {
			eImpactSearch=0;
			printf("please select a vertex.\n");
		}
	}
	
	printf("applying cuts, Nseg %d, PlateUp = %d \n",eNsegCut, ePlateUp);
	if( eSelected==NULL ) eNeighborSearch = 0;
	if( eNeighborSearch ){
		printf(	"Neighbor search On. tolerance : dTX = %.3lf rad, dTY = %.3lf rad, dmin = %.0lf micron\n",eTolDTX,eTolDTY, eTolDX);
	}
	if( eAngularCut ){
		printf(	"Angular cut On. tolerance : dTX = %.3lf rad, dTY = %.3lf rad\n", eTolDTX, eTolDTY);
	}
	if( eImpactSearch ){
		printf( "Impact cut On. tolerance : IP < %.0lf micron\n", eTolDX);
	}
	
	if( eSideOut ){
		printf( "Penetration and Side-Out rejection : require %d scanning area upstream\n", eSideOutPlate);
	}
	int nentr = tracksbase->GetEntriesFast();
	for(int i=0;i<nentr;i++){
		if(i%1000==0) printf("%d / %d. %d tracks\r", i, tracksbase->GetEntriesFast(), tracks->GetEntriesFast());
		EdbTrackP *t = (EdbTrackP *) tracksbase->At(i);
		if(NULL==t) continue;
		// Selections
		
		// Nseg cut
		if( t->N() < eNsegCut) continue;

		// Max plate cut
		int ipl_first = t->GetSegmentFirst()->Plate();
		
		if( ipl_first < ePlateUp) continue;
		
		// Reject Side-out tracks
		if( eSideOut ){
			if( SideOut(t) ) continue;
		}
		
		// Neighbor Search. use eTolDX, eTolDT
		if( eNeighborSearch) {
			if( Neighborhood(t) == 0 ) continue;
		}
		
		// Impact Search. use eTolDX
		if( eImpactSearch ){
			if( ImpactSearch (t) == 0 ) continue;
		}
		
		// Angular Cut
		if( eAngularCut ){
			if(	fabs(t->TX() - eTX) > eTolDTX) continue;
			if(	fabs(t->TY() - eTY) > eTolDTY) continue;
		}
		
		// PH Cut
		double ph = t->Wgrains()/t->N();
		if(ePHCut > ph ) continue;
		
		// PH - Angle deviation cut
		double dtrms = DTRMS(t);
		if(dtrms> ePHDTRMS*(ph-ePHCut)) continue;
		
		if(tracks->FindObject(t)==NULL) tracks->AddLast(t);
	}
	printf("%d tracks out of %d tracks\n", tracks->GetEntriesFast(), tracksbase->GetEntriesFast());
}

double EdbEDATrackSet::GetTrackLength(EdbSegP *s, int updown){  
	// Calculate Track extrapolation length.
	// updon : up = -1, down = 1
	if( eTrackLength>0 ) return eTrackLength;
	
	// ScanSet case
	if(eScanSet!=NULL&&eDataSet==NULL){
		EdbPlateP *pl0 = eScanSet->GetPlate(s->Plate());
		EdbPlateP *pl1 = eScanSet->GetPlate(s->Plate()+updown);
		if(pl1==NULL) {
			// in case some plates are missing.
			for(int i=1;i<eScanSet->eB.Npl();i++){
				int dpl=i*updown;
				pl1 = eScanSet->GetPlate(s->Plate()+dpl);
				if(pl1!=NULL) break;
			}
			if(pl1==NULL) pl1=pl0;
		}
		return fabs( pl0->Z() - pl1->Z() )*0.5; 
	}	
	
	// PVRec case
	// check if the next plate exist (only the case ePVR has patterns)
	//  currently SB,SF doesn't have patterns, to be done.
	if(ePVR){
		if(ePVR->Npatterns()!=0){
			for(int i=0;i<ePVR->Npatterns();i++){
				EdbPattern *pat0 = ePVR->GetPattern(i);
				if( pat0==NULL) continue;
				if( pat0->GetSegment(0)==NULL ) continue;
				if( pat0->GetSegment(0)->Plate() == s->Plate()) {
					int pidnext = i+updown;
					if(pidnext<0||pidnext>ePVR->Npatterns()) return 650;
					EdbPattern *pat1 = ePVR->GetPattern(pidnext);
					if(pat1==NULL) return 650;
					else return fabs( pat1->Z()-pat0->Z() )*0.5;
				}
			}
		}
	}
	
	return fabs( gEDA->GetZ(s->Plate()) - gEDA->GetZ(s->Plate()+updown) )*0.5; 
} 


int EdbEDATrackSet::GetTrackColor(EdbSegP *s){
	// return Track color. if eTrackColor==-1, calculate according to eColorMode
	
	if(eTrackColor!=-1) return eTrackColor;

	if(eColorMode==kCOLOR_BY_PLATE)    { int ipl=s->Plate(); if(ipl==0) ipl=s->PID();return ((ipl+5)%13)*4+52; }
	if(eColorMode==kCOLOR_BY_PH)       { int ph = (int)s->W(); if(ph>32) ph=32; return ((ph-8)%25)*2+50;}
//	if(eColorMode==kCOLOR_BY_PH)       { int ph = (int)s->W(); if(ph>16) ph=16; return (ph-4)*4+51;}
	if(eColorMode==kBLACKWHITE)  { return kBlack; };
	if(eColorMode==kCOLOR_BY_PARTICLE) { return s->Flag();}
	if(eColorMode==kCOLOR_BY_ID)       {
		int col = s->Track()%10+1;
		if (col==5) col = kOrange;
		if (col==10) col = kGreen+3;
		return col;
	}
	if(eColorMode==kCOLOR_BY_MCID)     {
		int col = s->MCTrack()%10+1;
		if (col==5) col = kOrange;
		if (col==10) col = kGreen+3;
		if (col==1) col = kOrange+3;
		if (col<0) col = kWhite; //background
		return col;
	}

	return (s->Plate()%13)*4+52;
}


void EdbEDATrackSet::SetColorMode(int mode){

	if(mode==kBLACKWHITE){
		eTrackColor = kBlack;
		eLayerColor = kBlack;
		eLayerColorMT = kBlack;
		//eLayerWidth = -1;
	}
	if(eColorMode==kBLACKWHITE&&mode!=kBLACKWHITE){
		ResetTrackAttribute();
	}
	eColorMode = mode;
}

void EdbEDATrackSet::DrawSingleSegment(EdbSegP *s, const char *elname, TEveCompound *cmp, int extendup, int extenddown){
	// Draw Single Segment
	
	if(s==NULL) return;
	
	if(cmp==NULL) {
		cmp = new TEveCompound;
		gEve->AddElement(cmp);
	}
	cmp->OpenCompound();
	
		// default extention line
	TEveLine *l = new TEveLine;
	l->SetName(elname);
	l->SetLineWidth(GetTrackWidth(s));
	l->SetUserData(s);
	
	l->SetLineColor(GetTrackColor(s));
	if(extendup){
		//Auto calculation of track_line_length
		double dz = -GetTrackLength(s,-1);
		l->SetNextPoint(s->X()+dz*s->TX(), 
						s->Y()+dz*s->TY(), 
						(s->Z()+dz)*gEDA->GetScaleZ());

	} else {
		l->SetNextPoint(s->X(), s->Y(), s->Z()*gEDA->GetScaleZ());
	}
	if(extenddown){
		
		double dz = GetTrackLength(s,1);

		l->SetNextPoint(s->X()+dz*s->TX(), 
						s->Y()+dz*s->TY(), 
						(s->Z()+dz)*gEDA->GetScaleZ());
	} else {
		l->SetNextPoint(s->X(), s->Y(), s->Z()*gEDA->GetScaleZ());
	}
	cmp->AddElement(l);

	// sensitive layer.
	l = new TEveLine;
	l->SetName(elname);
	l->SetUserData(s);
	l->SetLineWidth(GetLayerWidth(s));

	l->SetLineColor(GetLayerColor(s));
	
	double dz = GetLayerLength(s);
	l->SetNextPoint( s->X()-dz*s->TX(), 
					 s->Y()-dz*s->TY(), 
					(s->Z()-dz)*gEDA->GetScaleZ());
	l->SetNextPoint( s->X()+dz*s->TX(), 
					 s->Y()+dz*s->TY(), 
					(s->Z()+dz)*gEDA->GetScaleZ());

	cmp->AddElement(l);
	
	cmp->CloseCompound();
}

void EdbEDATrackSet::DrawSingleTrack(EdbTrackP *t){
	int i;
	
	if(t==NULL) return;
	TEveCompound *cmp = new TEveCompound;
	gEve->AddElement(cmp);
	
	cmp->SetName( Form("%3s trk %4d pl %2d -> %2d %8.1f %8.1f %8.1f %7.4f %7.4f %d",
				GetName(), t->ID(), t->GetSegmentFirst()->Plate(), t->GetSegmentLast()->Plate(), t->X(),t->Y(),t->Z(),t->TX(),t->TY(), t->Flag()));
	cmp->SetUserData(t);

	int nseg = t->N();
	EdbSegP *last_seg=NULL;
	for(i=0;i<nseg;i++){
		EdbSegP *s = t->GetSegment(i);
		
		if(eDrawMT==kFALSE && s->Side()!=0) continue;

		char elname[256];
		// set name
		if(gEDA->Japanese()) sprintf(elname,
			"%3s %4d seg  itrk %4d pl %2d(%2d) id %6d %3s ph %2d %8.1f %8.1f %8.1f %7.4f %7.4f chi2 %4.2f sflag %3d", 
			GetName(), s->ScanID().eMinor, t->ID(), s->Plate(), 58-s->Plate(), (int)s->ID(), s->Side()? s->Side()==1? "MT1" : "MT2" : "BT", (int)s->W(), 
			s->X(), s->Y(), s->Z(), s->TX(), s->TY(), s->Chi2(), s->Flag());
		else sprintf(elname,
			"%3s %4d seg  itrk %4d pl %2d id %6d %3s ph %2d %8.1f %8.1f %8.1f %7.4f %7.4f chi2 %4.2f sflag %3d", 
			GetName(), s->ScanID().eMinor, t->ID(), s->Plate(), (int)s->ID(), s->Side()? s->Side()==1? "MT1" : "MT2" : "BT", (int)s->W(), 
			s->X(), s->Y(), s->Z(), s->TX(), s->TY(), s->Chi2(), s->Flag());
		if(s->MCEvt()>=0) sprintf(elname,"%s MCEvt %d P %7.4f PdgID %6d",elname,s->MCEvt(),s->P(), s->Flag());

		// set extend or not
		int extendup=1, extenddown=1;
		if(eExtendMode==kExtendAuto){
			if(i==0)		extendup=0;
			if(i==nseg-1)	extenddown=0;
			if(t->N()==1)   extenddown=1;
		} 
		else if (eExtendMode == kExtendUpDown) extendup=extenddown=1;
		else if (eExtendMode == kExtendDown) { extendup=0; extenddown=1;}
		else if (eExtendMode == kExtendUp)   { extendup=1; extenddown=0;}
		else if (eExtendMode == kExtendNo)     extendup=extenddown=0;

		DrawSingleSegment(s, elname, cmp, extendup, extenddown);
		
		// in case that there was holes, draw the interpolation line.
		// from end point of extention(down) of last_seg to that of extention(up) of s.
		if(eInterpolation&&last_seg!=NULL){
			// hole = skip of plate number
			if(abs(last_seg->PID()-s->PID())>1){
				TEveLine *l = new TEveLine;
				l->SetLineStyle(7);
				l->SetName(elname);
				l->SetUserData(last_seg);
				l->SetLineWidth(1);
				l->SetLineColor(kGray+1);

				double dz = -GetTrackLength(s,-1); // to upstream

				l->SetNextPoint(s->X()+dz*s->TX(), 
								s->Y()+dz*s->TY(), 
								(s->Z()+dz)*gEDA->GetScaleZ());

				dz = GetTrackLength(last_seg,1); // to downstream
				
				l->SetNextPoint(last_seg->X()+dz*last_seg->TX(), 
								last_seg->Y()+dz*last_seg->TY(), 
								(last_seg->Z()+dz)*gEDA->GetScaleZ());
				cmp->OpenCompound();
				cmp->AddElement(l);
				cmp->CloseCompound();
			}
		}
		last_seg=s;
	}
	
	if((eDrawTrackID&&eDrawTrackIDG)||(eDrawTrackAngle&&eDrawTrackAngleG)){
		char buf[256];
		strcpy(buf,"");
		if(eDrawTrackID&&eDrawTrackIDG)    sprintf(buf,"%d",t->ID());
		if(eDrawTrackAngle&&eDrawTrackAngleG) sprintf(buf,"%s (%.3lf,%.3lf)", buf, t->TX(), t->TY());
		TEveText* tx = new TEveText;
		tx->SetText(buf);
		tx->SetMainColor(eTextColor);
		tx->SetFontMode(TGLFont::kBitmap);
		tx->SetFontFile(eTextFont);
		tx->SetFontSize(eTextSize);
		
		if(eTextPosition==kTop){
			EdbSegP *s = t->GetSegmentFirst();
			tx->PtrMainTrans()->SetPos(	s->X(), s->Y(), s->Z()*gEDA->GetScaleZ());
		}
		else if(eTextPosition==kMiddle){
			EdbSegP *s1 = t->GetSegmentFirst();
			EdbSegP *sl = t->GetSegmentLast();
			tx->PtrMainTrans()->SetPos( (s1->X()+sl->X())/2, (s1->Y()+sl->Y())/2, (s1->Z()+sl->Z())/2);
		}
		else {
			EdbSegP *s = t->GetSegmentLast();
			tx->PtrMainTrans()->SetPos( s->X()+s->TX()*650, s->Y()+s->TY()*650, (s->Z()+650)*gEDA->GetScaleZ());
		}
		cmp->OpenCompound();
		cmp->AddElement(tx);
		cmp->CloseCompound();
	}
}
/*
void EdbEDATrackSet::DrawSingleComment(EdbEDATrackComment *tcmt){
	
	EdbTrackP *t = tcmt->GetTrack();
	
	if(eTracks->FindObject(t)==NULL) return; // if the corresponding track will not be drawn, skip to draw comment.
	TEveText* tx = new TEveText;
	//printf("%s\n", tcmt->GetComment());
	tx->SetText(tcmt->GetComment());
	// draw text 1mm downstream of the last segment.
	// ! to be selected by eTextPosition......
	EdbSegP *s = t->GetSegmentLast();
	tx->PtrMainTrans()->SetPos(
		s->X()+s->TX()*1000, 
		s->Y()+s->TY()*1000,
		(s->Z()+1000)*gEDA->GetScaleZ());
	
	tx->SetFontMode(TGLFont::kBitmap);
	tx->SetFontSize(eTextSize);
	tx->SetFontFile(eTextFont);
	tx->SetMainColor(eTextColor);
	
	gEve->AddElement(tx);
}
*/
void EdbEDATrackSet::DrawSingleComment(EdbEDATrackComment *tcmt){
	
	EdbTrackP *t = tcmt->GetTrack();
	
	// if ID, Angle are draw. use same TEveText.
	TEveCompound *cmp = EdbEDAUtil::GetTrackElement(t);
	if(cmp==NULL) return; // if the corresponding track will not be drawn, skip to draw comment.

	TEveText *tx = NULL;
	for(TEveElement::List_i it=cmp->BeginChildren(); it!=cmp->EndChildren();it++){
		if(strncmp("TEveText",(*it)->GetElementName(),8)==0) {
			tx = (TEveText *) (*it);
			tx->SetText(Form("%s %s", tx->GetText(), tcmt->GetComment()));
		}
	}
	
	// if ID, Angle are not drawn.
	if(tx==NULL) {
		tx=new TEveText();
		tx->SetText(tcmt->GetComment());
		// draw text 1mm downstream of the last segment.
		// ! to be selected by eTextPosition......
		EdbSegP *s = t->GetSegmentLast();
		tx->PtrMainTrans()->SetPos(
			s->X()+s->TX()*1000, 
			s->Y()+s->TY()*1000,
			(s->Z()+1000)*gEDA->GetScaleZ());
		
		tx->SetFontMode(TGLFont::kBitmap);
		tx->SetFontSize(eTextSize);
		tx->SetFontFile(eTextFont);
		tx->SetMainColor(eTextColor);
		gEve->AddElement(tx);
	}
	
}


void EdbEDATrackSet::Draw(int redraw){
	eAreaSet->Draw();

	if(eDraw==kFALSE) {
		printf("TrackSet %s Draw    0 / %4d tracks. Disabled\n", GetName(), NBase());
		return;
	}
	for(int i=0;i<N();i++){
		EdbTrackP *t = GetTrack(i);
		DrawSingleTrack(t);
	}
	
	for(int i=0;i<eComments->GetEntriesFast();i++){
		EdbEDATrackComment *cmt = (EdbEDATrackComment *) eComments->At(i);
		DrawSingleComment(cmt);
	}
	
	
	if(redraw) gEve->Redraw3D();
	printf("TrackSet %s Draw %4d / %4d tracks\n", GetName(), N(), NBase());
}

EdbTrackP * EdbEDATrackSet::SearchCorrespondTrack(EdbTrackP *t){
	// check if the track corresponds to the scanback track. return scanback id, if not return -1.
	if(t==NULL||eTracks==NULL) return NULL;
	EdbSegP *s = t->GetSegmentLast();
	for(int i=0;i<N();i++){
		EdbSegP *ss = new EdbSegP(*((EdbSegP *)GetTrack(i)));
		ss->PropagateTo(s->Z());
		if(fabs(s->TX()-ss->TX())<0.03){
		if(fabs(s->TX()-ss->TX())<0.03){
			if( fabs(s->X()-ss->X())<100 ){
			if( fabs(s->Y()-ss->Y())<100 ){
				printf("      Track: %d corresponds to scanback %d\n", t->ID(), ss->ID());
				return GetTrack(i);
			}}
		}}
		delete ss;
	}
	return NULL;
}

void EdbEDATrackSet::DoVertexing(EdbVertexRec *VR, int ntrkcut){
	// Do vertexing. if EdbVertexRec is given, use the parameter for reconstruction.
	if(VR==NULL){
		VR = new EdbVertexRec;
		VR->eDZmax      = 5000.;  // maximum z-gap in the track-vertex group
		VR->eProbMin    = 0.001;  // minimum acceptable probability for chi2-distance between tracks
		VR->eImpMax     = 5.;     // maximal acceptable impact parameter [microns] (for preliminary check)
		VR->eUseMom     = false;  // use or not track momentum for vertex calculations
		VR->eUseSegPar  = true;   // use only the nearest measured segments for vertex fit (as Neuchatel)
		VR->eQualityMode= 0;      // vertex quality estimation method (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distanc
		VR->eEdbTracks = GetTracks();
	}
	if(ePVR==NULL) ePVR = new EdbPVRec;
	
	if(ePVR->GetScanCond()==NULL) ePVR->SetScanCond(new EdbScanCond);
	
	VR->eVTX = new TObjArray;
	
	VR->SetPVRec(GetPVRec());
	int nvtx=VR->FindVertex();
	if(nvtx)VR->ProbVertexN();
	
	gEDA->GetVertexSet()->Clear();
	gEDA->GetVertexSet()->AddVertices(VR->eVTX);
	gEDA->GetVertexSet()->SetVerticesNtrk(ntrkcut);
	
	ClearTracks();
	SetTracksVertices(gEDA->GetVertexSet()->GetVertices());

}

void EdbEDATrackSet::ReadFile(const char *filename, int datatype, TCut rcut){
	//  filename : filename for "LinkDef" or "*.set.root(ScanSet)" or "*.root(linked_track.root format)"
	//  datatype : only for "LinkDef". Data type for EdbDataProc::InitVolume(). put -1 for no InitVolume
	//            if datatype= 100 or 1000. linked tracks will be registred to "TS".
	//  rcut     : Track selection when read.
	
	// if filename is not given. open file browser.
	
	if(filename==NULL) {
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "LinkDef", "*.def","Linkeda Tracks", "*.root","ScanSet", "*.set.root", "Mxx file","*.all",0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, fi);
		filename = fi->fFilename;
	}
	
	if(filename==NULL){
		printf("Cancel EdbEDATrackSet::ReadFile\n");
		return;
	}
	
	// check the existence of the file.
	if(gSystem->AccessPathName(filename, kReadPermission)!=0) {
		printf("Error : %s is not found\n", filename);
		return;
	}

	// select an action for the file.
	TString s(filename);
	s.ToLower();
	if(s.EndsWith(".set.root")){
		// Scan set file
		printf("Read Scan Set : %s\n", s.Data());
		TFile *f = TFile::Open(filename);
		eScanSet = (EdbScanSet *)f->Get("set");
		eID = *((EdbID *)eScanSet->eIDS.At(0));
		eID.ePlate = 0;
		
		TString s2;
		gEDA->ScanProc()->MakeFileName(s2, eID, "trk.root", kFALSE);
		int IsExist = gSystem->AccessPathName(s2.Data(), kReadPermission);
		if( IsExist==kFALSE) {
			// linked tracks case
			ReadTracksTree(eID,rcut);
		}
		else {
			// prediction scan
			ReadPredictionScan(eScanSet);
		}
	}
	else if(s.EndsWith(".def")){
		// Link def file
		printf("Read LinkDef : %s\n", s.Data());
		EdbDataProc *dproc = new EdbDataProc(filename);
		eDataSet = dproc->GetDataSet();
		if(datatype>=0) {
			dproc->InitVolume(datatype);
			ReadLinkedTracks(dproc);
		}
	}
	else if(s.EndsWith(".root")){
		// Linked tracks file
		printf("Read Linked_track : %s\n", s.Data());
		ReadLinkedTracksFile(filename,rcut);
	}
	else if(s.EndsWith(".all")){
		// Mxx file (Nagoya style ascii)
		ReadMxxFile(filename);
	}
	
	if(NBase()!=0) SetDraw();
}

void EdbEDATrackSet::ReadTracksTree (const char *scanset_filename, TCut cut){
	TFile *f = TFile::Open(scanset_filename);
	if(f==NULL) return;
	eScanSet = (EdbScanSet *)f->Get("set");
	if(eScanSet==0) {
		fprintf(stderr,"Error at ReadTracksTree : cannot find ScanSet in %s\n",scanset_filename);
		return;
	}
	EdbID id = *((EdbID *) (eScanSet->eIDS.At(0)));
	id.ePlate=0;
	ReadTracksTree(id, cut);
	f->Close();
}

void EdbEDATrackSet::ReadTracksTree (EdbID ID, TCut cut){
	EdbScanProc *sproc = gEDA->ScanProc();
	eScanSet = sproc->ReadScanSet(ID);
	ePVR = new EdbPVRec();
	sproc->ReadTracksTree(ID, *ePVR, cut);
	AddTracksPVR(ePVR);
	eID=ID;
}

EdbDataProc *EdbEDATrackSet::ReadLinkedTracksFile (const char *filename, TCut cut){
	// Read linked_tracks.root format file.
	
	// read file
	EdbDataProc *dproc = new EdbDataProc;
	ePVR = new EdbPVRec;
	dproc->ReadTracksTree(*ePVR, filename, cut);
	
	// if plate number is not set. set plate number as PID.
	int flag_plset=0;
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
	AddTracksPVR(ePVR);
	
	return dproc;
}

EdbDataProc *EdbEDATrackSet::ReadLinkedTracks (const char *lnkdef, TCut cut){
	EdbDataProc *dproc = new EdbDataProc(lnkdef);
	dproc->InitVolume(100, cut);
	ReadLinkedTracks(dproc);
	
	return dproc;
}

EdbDataProc *EdbEDATrackSet::ReadLinkedTracks (EdbDataProc *dproc){
	// read tracks from linked_tracks.root for old structure.
	// plate number will be filled from link-list file.
	
	printf("EdbEDATrackSet::ReadLinkedTracks(EdbDataProc *)\n");
	
	eDataSet  = dproc->GetDataSet();
	EdbPVRec    *pvr   = dproc->PVR();
	
	// set plate number to segment.
	for(int j=0;j<pvr->Ntracks();j++){
		EdbTrackP *t = pvr->GetTrack(j);
		
		for(int k=0;k<t->N();k++){
			EdbSegP *s = t->GetSegment(k);
			int pid = s->PID();
			
			int ipl = eDataSet->GetPiece(pid)->Plate();
			s->SetPlate(ipl);
		}
		t->SetPlate(t->GetSegmentFirst()->Plate());
	}
	
	AddTracksPVR(pvr);
	
	EdbID id =pvr->GetTrack(0)->GetSegment(0)->ScanID();
	
	SetID( EdbID(id.eBrick, 0, id.eMajor, id.eMinor));
	
	return dproc;
}


/*
		EdbDataSet *dset= new EdbDataSet("lnk.def");
	int ipl0 = dset->GetPiece(0)->Plate();
	
	for(int j=0;j<pvr->Ntracks();j++){
		EdbTrackP *t = pvr->GetTrack(j);
		
		for(int k=0;k<t->N();k++){
			EdbSegP *s = t->GetSegment(k);
			int pid = s->PID();
			
			int ipl = dset->GetPiece(pid)->Plate();
			s->SetPID(ipl);
		}
		t->SetPID(t->GetSegmentFirst()->PID());
	}
*/


void EdbEDATrackSet::ReadPredictionScan(int BRICK, int SBVERSION, int RUN_PRED, int UseMicrotrack){
	// this is an inferface function from old style usage.
	// reccoment ReadPredictionScan( EdbID ) for non-Bern people.
	// this force to update scanset file = "bxx.set.root".
	
	// if Brick number is not set, show message.
	// give command line input.
	if(BRICK==0){
		new TGMsgBox(gClient->GetRoot(), 0,
		            "EDA", "Follow your shell window");
		system("pwd");
		printf("Read Scanback ---------------- \n");
		printf("  Enter Brick number : ");
		scanf("%d", &BRICK);
		printf("  Enter Run number : ");
		scanf("%d", &RUN_PRED);
	}
	
	printf("Read Scanback ibrick = %d run = %d ver = %d\n", BRICK, RUN_PRED, SBVERSION);

	eID=EdbID(BRICK, 0, SBVERSION, RUN_PRED);
	
	ReadPredictionScan(eID, kTRUE);
	SetDrawMT(UseMicrotrack);

	SetDraw();
}


void EdbEDATrackSet::ReadPredictionScan(EdbID id, bool force_update_setroot){
	// Read Prediction scan
	// if there is "set.root" already exist, ask whether overwrite or not.
	// if force_update_setroot == kTRUE, overwrite without asking.
	
	if(id.eMinor==0) {
		printf("ReadPredictionScan: strange id.eMinor %d, stop\n", id.eMinor);
		ErrorMessage(Form("ReadPredictionScan: strange id.eMinor %d, stop\n", id.eMinor));
		return;
	}
	
	eID=id;
	TString s;
	gEDA->ScanProc()->MakeFileName(s, eID, "set.root", kFALSE);
	int update_setroot=1;
	if(force_update_setroot) update_setroot=1;
	else if(gSystem->AccessPathName(s.Data(), kReadPermission)==0){
		// if set.root already exist, ask whether update or not.
		update_setroot = AskYesNo(Form("%s is already exist. do you want to update?", s.Data()));
	}
	
	if(update_setroot==1){
		EdbScanSet *ss = new EdbScanSet();
		ss->MakeNominalSet(id,60);
		gEDA->ScanProc()->MakeScannedIDList(id, *ss, 60, 0,"found.root"); 
		ss->eB.SetID(id.eBrick);
		if(ss->eIDS.GetSize()==0){
			printf("ID %d.%d.%d.%d is not scanned\n", id.eBrick, id.ePlate, id.eMajor, id.eMinor);
			return;
		}
		/*
		EdbAffine2D *aff = ss->eB.GetPlate(ss->eB.Npl()-1)->GetAffineXY();
		if(aff->A11()==1.0){
			printf("ReadPredictionScan : no affine parameters found as ScanBack. try to read as ScanForth\n");
			delete ss;
			ss = new EdbScanSet();
			ss->MakeNominalSet(id);
			gEDA->ScanProc()->MakeScannedIDList(id, *ss, 0, 60,"found.root"); 
			ss->eB.SetID(id.eBrick);
			gEDA->ScanProc()->AssembleScanSet(*ss);
		}
		*/
		
		// Z shift
		int ipl = ss->eB.GetPlate(0)->ID();
		float zl = ss->eB.GetPlate(0)->Z();
		float zg = gEDA->GetZ(ipl);
		ss->ShiftBrickZ(-zl);
		ss->ShiftBrickZ(zg);
	
		gEDA->ScanProc()->WriteScanSet(id, *ss);
	
		if(gEDA->GetScanSet()) {
			EdbID idmain = gEDA->GetID();
			printf("Update geometry of %d.%d.%d.%d with %d.%d.%d.%d\n", 
				id.eBrick, id.ePlate, id.eMajor, id.eMinor,
				idmain.eBrick, idmain.ePlate, idmain.eMajor, idmain.eMinor);
			gEDA->ScanProc()->UpdateSetWithAff(id, idmain); 
			ss = gEDA->ScanProc()->ReadScanSet(id);
		}else {
			gEDA->ScanProc()->AssembleScanSet(*ss);
		}
	}
	//ReadPredictionScan(ss);
	
	ReadFile(const_cast<char*>(s.Data()),0,"1");
	
}

void EdbEDATrackSet::ReadPredictionScan(EdbScanSet *ss){
	EdbPVRec *ali = new EdbPVRec;
	gEDA->ScanProc()->ReadFoundTracks(*ss, *ali,-11);

	eScanSet = ss;
	eID = *((EdbID *)eScanSet->eIDS.At(0));
	eID.ePlate = 0;
	
	
	// to be compatible also main set by EdbDataSet
	// this is temporary solution.
	if(gEDA->GetDataSet()){
		for(int i=0;i<ali->Ntracks();i++){
			EdbTrackP *t = ali->GetTrack(i);
			for(int j=0;j<t->N();j++){
				if(t->N()==0) break;
				EdbSegP *s = t->GetSegment(j);

				// set Z from global calibration
				double z = gEDA->GetZ(s->Plate());
				s->SetZ(z);
				// reset affine on this ScanSet.
				EdbPlateP *pl =eScanSet->GetPlate(s->Plate());
				if(pl==NULL) continue;
				EdbAffine2D aff = *(pl->GetAffineXY());
				aff.Invert();
				s->Transform(&aff); // --> plate coordinate

				// apply affine of Main ScanSet (or DataSet).
				EdbAffine2D *affg = gEDA->GetAffine(s->Plate());
				if(NULL==affg) continue; // if no affine params are found in the main calibration, just leave it.
				s->Transform(affg); // --> global coordinate
			}
		}
	}
	
	AddTracksPVR(ali);
	SetDraw();
}




void EdbEDATrackSet::ReadTextTracks(const char *filename){
	// Read Text file data
	// the format should be : itrk ipl ph x y ax ay flag_trans
	// flag_trans control Transformation. Affine=plate->global
	// flag_trans==1 (Affine)
	// flag_trans==2 (Affine)+(Bottom->Middle of base). For manualcheck data
	// flag_trans==3 (Affine)+(Top->Middle of base). For manualcheck data
	int i,j;
	char buf[512];

	FILE *fp;

	if(filename==NULL){
		TGFileInfo *fi = new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDOpen, fi);
		filename=fi->fFilename;
	}

	printf(" the format should be : itrk ipl ph x y ax ay flag_trans runid\n");
	printf(" the Z value by ipl will be used\n");
	printf(" if(flag_affine==1) affine transform will be applied from plate->global\n");
	if(filename==NULL) {
		printf("Read tracks. no filename. end.\n");
		return;
	}
	printf("Read tracks from Text : %s -> TrackSet:%s\n", filename, GetName());
	if(filename!=NULL) fp = fopen(filename,"rt");
	
	if(fp==NULL){
		printf("text track is not available.\n");
		return;
	}
	
	eID=gEDA->GetID();

	ePVR = new EdbPVRec;
	
	EdbScanCond cond;
	for(;fgets(buf,sizeof(buf),fp)!=NULL;){
		if(feof(fp)||ferror(fp)) break;
		int itrk;
		double x,y,z,ax,ay;
		int ipl,ph,flag_affine=1;
		int runid;
		printf("%s",buf);
		int n = sscanf(buf,"%d %d %d %lf %lf %lf %lf %d %d", &itrk, &ipl, &ph, &x, &y, &ax, &ay, &flag_affine,&runid);
		if(n<9) runid=-99;
		eID.eMajor=runid%100;
		eID.eMinor=runid-runid%100;

		if(n==0) n=-1;
		int pid = gEDA->GetPID(ipl);
		
		// if the PID is out of volume scan, skip the segment.
		if(pid==-1) continue;
		z=gEDA->GetZ(ipl);

		EdbSegP *s = new EdbSegP;
		s->Set(itrk,x,y,ax,ay,0,0);
		s->SetZ(z);
		s->SetPID(pid);
		s->SetW(ph);
		s->SetTrack(itrk);
		s->SetPlate(ipl);
		s->SetAid(eID.eMajor, eID.eMinor);

		// fill COV for vertexing
		s->SetErrors();
		cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
		
		// add segment to track
		EdbTrackP *t = NULL;
		
		// search track in sbtrk
		for(j=0;j<N();j++){
			EdbTrackP *tt = GetTrack(j);
			if(s->ID()==tt->ID()) t = tt;
		}
		
		// if not track corresponding.
		if( t==NULL){
			t = new EdbTrackP;
			t->SetID(s->ID());
			t->SetFlag(runid);
			AddTrack(t);
		}
		
		if(flag_affine==1){
			EdbAffine2D *traff = gEDA->GetAffine(ipl);
			if(NULL!=traff){
				s->Transform(traff);
			}
		}
		if(flag_affine==2){
			// in case of manualcheck data at bottom
			EdbAffine2D *traff = gEDA->GetAffine(ipl);
			if(NULL!=traff){
				s->Transform(traff);
				s->SetZ(z-102.5);
				s->PropagateTo(z);
			}
		}
		if(flag_affine==3){
			// in case of manualcheck data at top
			EdbAffine2D *traff = gEDA->GetAffine(ipl);
			if(NULL!=traff){
				s->Transform(traff);
				s->SetZ(z+102.5);
				s->PropagateTo(z);
			}
		}
		
		// add segment into track.
		t->AddSegment(s);
	}
	fclose(fp);
	
	for(i=0;i<N();i++){
		EdbTrackP* t = GetTrack(i);
		EdbSegP* s = t->GetSegmentFirst();
		t->Set(s->ID(),s->X(),s->Y(),s->TX(),s->TY(),s->W(),t->Flag());
		t->SetZ(s->Z());
		t->SetCounters();
	}
	printf("%d tracks available from text file\n", N());
}

void EdbEDATrackSet::ReadMxxFile(const char *filename){
	EdbPVRec *pvr = EdbEDAUtil::ReadMxxPVR(filename);
	AddTracksPVR(pvr);
}


void EdbEDATrackSet::ReadListFile(const char *filename, bool clear_previous){
	// Read list file.
	// Search tracks in this track set.
	// see also EdbEDAMainTab::ReadListFile(), which try to search 
	// in BT files.
	
	if(filename==NULL){
		TGFileInfo *fi = new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "List file", "*.lst", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDOpen, fi);
		filename=fi->fFilename;
	}

	printf("Read List file : %s -> TrackSet : %s\n", filename, GetName());
	FILE *fp=NULL;
	if(filename!=NULL) fp = fopen(filename,"rt");
	
	if(fp==NULL){
		printf("List file is not available.\n");
		return;
	}

	if(clear_previous) {
		ClearTracks();
		ClearComments();
	}
	char buf[256];
	int ipl, iseg;
	char comment[256];
	// read list, setting TObjArray *selected.
	for(int i=0;NULL!=fgets(buf,sizeof(buf),fp);i++){
		if(feof(fp)||ferror(fp)) break;
		sscanf(buf,"%d %d %s", &ipl, &iseg, comment);
		
		EdbTrackP *t=NULL;
		for(int j=0;j<NBase();j++){
			EdbTrackP *tt = GetTrackBase(j);
			if(NULL==tt) continue;
			for(int k=0;k<tt->N();k++){
				EdbSegP *s = tt->GetSegment(k);
				if(s->ID()==iseg){
				if(s->Plate()==ipl){
					SetTrack(tt); // add tt for drawing
					t=tt;
					break;
				}}
			}
		}
		
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
					AddTrack(t);
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
				AddComment(t,comment);
			}
		}
	
	}
}

const char * EdbEDATrackSet::WriteListFile(const char *filename, bool append, bool open_editor){
	// Write List file.
	// ouput plate_number and segment_id. (track ID could change everytime after reconstruction, but segment ID doesn't change.)
	// the result of MT search will be ignored.
	if(filename==NULL){
				TGFileInfo *fi = new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "List file", "*.lst", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		filename=fi->fFilename;
	}

	if(NULL==filename) return NULL;
	printf("Write List file : %s (%6s) <- TrackSet:%s\n", filename, append?"append":"create", GetName());
	FILE *fp;
	if(append) fp = fopen(filename,"at");
	else       fp = fopen(filename,"wt");
	
	if(fp==NULL){
		printf("Couldn't open the file.\n");
		return filename;
	}

	int k=0;
	for(int j=0;j<N();j++){
		EdbTrackP *t=GetTrack(j);
		if(NULL==t) continue;
		
		EdbSegP *s = 0;
		
		if(ePVR){
			// select a segment which from original reconstruction, in order to keep information after microtrack search.
			// if the segment is registered in ePVR, it's from original reconstruction.
			// if not, it's found by MT search or any.
			for(int i=0, flag=0;i<t->N();i++){
				s = t->GetSegment(i);
				for(int ipat=0;ipat<ePVR->Npatterns();ipat++){
					EdbPattern *pat = ePVR->GetPattern(ipat);
					TClonesArray *arr = pat->GetSegments();
					if( arr->FindObject(s) ) {flag=1;break;}
				}
				if(flag) break;
			}
		} else s=t->GetSegmentFirst();
		
		fprintf(fp,"%2d %7d %10s itrk %5d %7.4lf %7.4lf\n", s->Plate(), s->ID(), Form("\"%s\"",GetComment(t)), t->ID(), t->TX(), t->TY());
		k++;
	}

	
	fclose(fp);
	printf("Write Selected : %d tracks\n", k); 

	// Open text editer
	if(open_editor) gEDA->OpenTextEditor(filename);
	return filename;
}


EdbTrackP * EdbEDATrackSet::SearchSegment(int ipl, int iseg) { 
	for(int i=0;i<eTracksBase->GetEntriesFast(); i++) {
		EdbTrackP *t = (EdbTrackP *) eTracksBase->At(i);
		for(int j=0;j<t->N();j++){
			EdbSegP *s = t->GetSegment(j);
			if(s->Plate()!=ipl) continue;
			if(s->ID()!=iseg) continue;
			return t;
		}
	}
	return NULL;
}


void EdbEDATrackSet::MicroTrackSearch(EdbTrackP *t){
	//	void SetCondTrackingDefault(){
	//		eRunTracking.eDeltaRview = 400.;
	//		eRunTracking.eDeltaTheta = 0.2;
	//		eRunTracking.eDeltaR = 20.;
	//
	//		eRunTracking.ePreliminaryPulsMinMT = 7;
	//		eRunTracking.ePreliminaryChi2MaxMT = 3.0;
	//
	//		eRunTracking.ePulsMinMT       = 9;
	//		eRunTracking.ePulsMinDegradMT = 1.5;
	//		eRunTracking.eChi2MaxMT       = 1.6;
	//
	//		eRunTracking.ePulsMinBT       = 18;
	//		eRunTracking.ePulsMinDegradBT = 1.5;
	//		eRunTracking.eChi2MaxBT       = 1.8;
	//
	//		eRunTracking.eDegradPos   = 5.;
	//		eRunTracking.eDegradSlope = 0.003;
	//
	//		SetCondMTDefault( eRunTracking.eCondMT);
	//		SetCondBTDefault( eRunTracking.eCondBT);
	//	}
	//
	//	void SetCondMTDefault(EdbScanCond &cond){
	//		cond.SetSigma0( 1., 1., 0.010, 0.010 );     // sigma0 "x, y, tx, ty" at zero angle
	//		cond.SetDegrad( 5. );                       // sigma(tx) = sigma0*(1+degrad*tx)
	//		cond.SetBins(0, 0, 0, 0);  //???            // bins in [sigma] for checks
	//		cond.SetPulsRamp0(  5., 5. );               // in range (Pmin:Pmax) Signal/All is nearly linear
	//		cond.SetPulsRamp04( 5., 5. );
	//		cond.SetChi2Max( 6.5 );
	//		cond.SetChi2PMax( 6.5 );
	//		cond.SetRadX0( 5810. );
	//		cond.SetName("OPERA_microtrack");
	//	}
	//
	//	void SetCondBTDefault(EdbScanCond &cond){
	//		cond.SetSigma0( 10., 10., 0.007, 0.007 );   // sigma0 "x, y, tx, ty" at zero angle
	//		cond.SetDegrad( 2. );                       // sigma(tx) = sigma0*(1+degrad*tx)
	//		cond.SetBins(0, 0, 0, 0);                   // bins in [sigma] for checks
	//		cond.SetPulsRamp0(  5., 5. );               // in range (Pmin:Pmax) Signal/All is nearly linear
	//		cond.SetPulsRamp04( 5., 5. );
	//		cond.SetChi2Max( 6.5 );
	//		cond.SetChi2PMax( 6.5 );
	//		cond.SetRadX0( 5810. );
	//		cond.SetName("OPERA_basetrack");
	//	}
	//
	
	int ret=PrepareScanSetForMT();
	if(ret==-1) {
		printf("ScanSet error. stop.\n");
		return ;
	}
	
	// Upstream segment search
	for(;;){
		EdbSegP *s = t->GetSegmentFirst();
		int ret = MicroTrackSearch(t,s,s->Plate()-1);
		if(ret==-1) {
			ret = MicroTrackSearch(t,s,s->Plate()-2);
		}
		if (ret==-1) break;
	}
	
	// Fill holes
	for(int i=0;i<t->N()-1;i++){
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);

		for(int j=1;j<=2&&i+j<t->N();j++){
			int ipl=s1->Plate()+j;
			if( ipl!=s2->Plate() ) {
				// search MT from upstream
				int ret=MicroTrackSearch(t, s1, ipl);
				if(ret!=-1) break;
				// search MT from downstream
				ret=MicroTrackSearch(t, s2, ipl);
				if(ret!=-1) break;
			}
			else{
				break;
			}
		}
	}
	
	// Downstream segment search
	for(;;){
		EdbSegP *s = t->GetSegmentLast();
		int ret = MicroTrackSearch(t,s,s->Plate()+1);
		if(ret==-1) break;
	}
}


int EdbEDATrackSet::PrepareScanSetForMT(){

// Scan id check
	EdbID& id = GetID();
	if(id.eMinor==0){
		id.eBrick = gEDA->GetBrick();
		InputID("Please input EdbID for TS.",id);
		SetID(id);
	}
	
	// get ScanSet of Total scan.
	// normally it's not filled even if TS are read via ReadTracksTree(id)
	EdbScanSet *sc = GetScanSet(); 
	
	if(sc==NULL){
		// if sc is NULL, read and assemble ScanSet. and set for TrackSet.
		// this will be executed only 1st time.
		EdbScanProc *sproc = gEDA->ScanProc();
		sc = new EdbScanSet();
		sc->MakeNominalSet(id);
		sproc->MakeScannedIDList(id, *sc, 65, -5,"cp.root"); 
		sc->eB.SetID(id.eBrick);
		int nid = sc->eIDS.GetSize();
		if(nid==0) {
			ErrorMessage(Form("Error : Not cp files! stop. ScanID=%d.%d.%d.%d", id.eBrick, id.ePlate, id.eMajor, id.eMinor));
			return -1;
		}
		
		gEDA->ScanProc()->WriteScanSet(id, *sc);
		if(gEDA->GetScanSet()) {
			gEDA->ScanProc()->UpdateSetWithAff(id, gEDA->GetID()); 
			sc = gEDA->ScanProc()->ReadScanSet(id);
		}else {
			gEDA->ScanProc()->AssembleScanSet(*sc);
		}
		
		if(sc->eB.Npl()!=0) SetScanSet(sc); // set Scanset in TrackSet
	}
	if(sc->eB.Npl()==0){
		// if there is not plate, stop.
		ErrorMessage("No plate in ScanSet :");
		id.Print();
		return -1;
	}
	
	return 0;
}

int EdbEDATrackSet::MicroTrackSearch(EdbTrackP *t, EdbSegP *app, int ipl){
	// Micro-track search.
	// Original from Napoli, modified by Frederic, impremented by Aki.
	// 
	// Use EdbScanProc
	// return -1 = not found, 0 = Basetrack found, 1 or 2 = Microtrack found.
	printf("## Micro-track search ##\n");
	
	
	// Check if there is already a segment on the given ipl.
	int seg_already_exist=0;
	for(int i=0;i<t->N();i++) if( t->GetSegment(i)->Plate()==ipl) seg_already_exist++;
	if(seg_already_exist){
		printf(" A segment is already reconstructed in the track. No search is done.\n\n");
		return -1;
	}
	
	// Check if the give plate is out of brick
	if(ipl<gEDA->GetFirstPlate() || ipl>gEDA->GetLastPlate() ){
		printf("searching plate %d is out of brick. stop.\n", ipl);
		return -1;
	}
	
	int ret = PrepareScanSetForMT();
	if(ret==-1) {
		printf("ScanSet error. stop.\n");
		return -1;
	}
	EdbScanProc *sproc = gEDA->ScanProc();
	
	// get ScanSet of Total scan.
	// normally it's not filled even if TS are read via ReadTracksTree(id)
	EdbScanSet *sc = GetScanSet(); 
	
	// Check if the plate is scanned or not.
	int flag_scanned=0;
	for(int i=0;i<sc->eIDS.GetSize();i++){
		EdbID *id = sc->GetID(i);
		if(ipl==id->ePlate) flag_scanned++;
	}
	if(flag_scanned==0) {
		//ErrorMessage(Form("Plate %d is not scanned.\n", ipl));
		printf("Plate %d is not scanned.\n", ipl);
		return -1;
	}
	// Tracking setting --> Moved in constructor
	//eRunTracking.ePredictionScan = false;
	//SetCondTrackingDefault();

	printf("  Segment pl=%d target_pl=%d\n",app->Plate(), ipl);
	
	// Affine parameters.
	// use main data-set affine to unify affine transformation.
	//EdbPlateP *pl0 = sc->GetPlate(app->Plate());
	//EdbPlateP *pl1 = sc->GetPlate(ipl);

	//EdbAffine2D *aff_vid0 = pl0->GetAffineXY();
	//EdbAffine2D *aff_vid1 = pl1->GetAffineXY();
	
	if(gEDA->GetAffine(app->Plate())==NULL||gEDA->GetAffine(ipl)==NULL) {
		printf("No affine parameters found. stop.\n");
		return -1;
	}

	EdbAffine2D aff_vid0 = *gEDA->GetAffine(app->Plate());
	EdbAffine2D aff_vid1 = *gEDA->GetAffine(ipl);
	
	
	aff_vid1.Print();
	aff_vid0.Print();
	
	aff_vid0.Invert();                   // invert affine
	aff_vid1.Invert();
	
	double z = gEDA->GetZ(ipl);           // z at target plate
	
	EdbSegP app_tr(*app);                 // copy of original track
	app_tr.Transform(&aff_vid0);           // from global to local

	EdbSegP spred(*app);                  // prediction from original track
	spred.SetPlate(ipl);
	spred.PropagateTo(z);                 // propagate to target Z
	
	spred.Transform(&aff_vid1);            // from global to local

//	printf("  Track id = %d\n", t->ID());
//	printf("  Original (global)  -> Plate = %d, X = %2.2f, Y = %2.2f, TX = %2.4f, TY = %2.4f\n", 
//				app->Plate(), app->X(), app->Y(), app->TX(), app->TY());
//	printf("  Original (local)   -> Plate = %d, X = %2.2f, Y = %2.2f, TX = %2.4f, TY = %2.4f\n", 
//				app_tr.Plate(), app_tr.X(), app_tr.Y(), app_tr.TX(), app_tr.TY());
//	printf("  Prediction (local) -> Plate = %d, X = %2.2f, Y = %2.2f, TX = %2.4f, TY = %2.4f\n", 
//				spred.Plate(), spred.X(), spred.Y(), spred.TX(), spred.TY());
	
	EdbID& id = GetID();
	EdbSegP fndbt,fnds1,fnds2,snewpred;
	int idp[4]={id.eBrick, ipl, id.eMajor, id.eMinor};
	sproc->InitRunAccess(eRunTracking,idp);
	
	int status = eRunTracking.FindPrediction(spred, fndbt, fnds1, fnds2, snewpred );
	if (status==-1) {
		printf(" *** MICROTRACK SEARCH ***\n");
		printf(" ***** NOTHING FOUND *****\n");
		if(gEve) gEve->SetStatusLine(Form("Nothing found on pl%02d. prediction from pl%02d", 
									ipl, app->Plate()));
	}
	else {
		printf("  \n");
		printf("  Status = %d, Candidate found\n", status);	
		printf("  Prediction (local) -> Plate = %2d, X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f\n", 
					spred.Plate(), spred.X(), spred.Y(), spred.TX(), spred.TY());
		printf("  Found      (local) -> Plate = %2d, X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f\n", 
					ipl, snewpred.X(), snewpred.Y(), snewpred.TX(), snewpred.TY());
		printf("  Residual                              %8.1f,     %8.1f,      %7.4f,      %7.4f\n", 
					spred.X()-snewpred.X(), spred.Y()-snewpred.Y(), spred.TX()-snewpred.TX(), spred.TY()-snewpred.TY());
		
//		printf("  \n");
		aff_vid1.Invert();            //from local to global 
		snewpred.Transform(&aff_vid1);
		printf("  Found     (global) -> Plate = %2d, X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f\n\n", 
			ipl, snewpred.X(), snewpred.Y(), snewpred.TX(), snewpred.TY()); 
	}

	
	if(status==0){
		EdbSegP *snewpred0 = new EdbSegP(snewpred);

		snewpred0->SetPID(gEDA->GetPID(ipl));
		snewpred0->SetZ(z);
		snewpred0->SetErrors(1.,1.,10.,0.002,0.002);
		snewpred0->SetDZ(-214);

		t->AddSegment(snewpred0);
		t->SetCounters();
		t->FitTrackKFS();
		
		printf("## BASETRACK found ##\n");
		printf(" MT Top : X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f, W = %2.0f, Chi2 = %5.2f,\n",
					fnds1.X(),fnds1.Y(),fnds1.TX(),fnds1.TY(),fnds1.W(),fnds1.Chi2());
		printf(" MT Bot : X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f, W = %2.0f, Chi2 = %5.2f,\n",
					fnds2.X(),fnds2.Y(),fnds2.TX(),fnds2.TY(),fnds2.W(),fnds2.Chi2());
		
		if(gEve) gEve->SetStatusLine(Form("Basetrack found on pl %d. (W,Chi2)= Top(%d,%.2f) Bot(%d,%.2f)", 
									ipl, (int)fnds1.W(), fnds1.Chi2(), (int)fnds2.W(), fnds2.Chi2()));
	}
	
	if(status==1){
		EdbSegP *snewpred1 = new EdbSegP(snewpred);
		
		snewpred1->SetPID(gEDA->GetPID(ipl));  
		snewpred1->SetErrors(1.,1.,10.,0.002,0.002);
		snewpred1->SetZ(z);
		snewpred1->SetDZ(-107);
		snewpred1->SetAid(snewpred1->Aid(0), snewpred1->Aid(1), status);
		t->AddSegment(snewpred1);
		t->SetCounters();
		t->FitTrackKFS();

		printf("## MICROTRACK Top found ##\n");
		printf(" MT Top : X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f, W = %2.0f, Chi2 = %5.2f,\n",fnds1.X(),fnds1.Y(),fnds1.TX(),fnds1.TY(),fnds1.W(),fnds1.Chi2());
		if(gEve) gEve->SetStatusLine(Form("Microtrack Top found on pl %d. (W,Chi2)=(%d,%.2f)", 
									ipl,(int)fnds1.W(), fnds1.Chi2()));
	}

	if(status==2){
		EdbSegP *snewpred2 = new EdbSegP(snewpred);

		snewpred2->SetPID(gEDA->GetPID(ipl));  
		snewpred2->SetErrors(1.,1.,10.,0.002,0.002);
		snewpred2->SetZ(z);
		snewpred2->SetDZ(-107);
		snewpred2->SetAid(snewpred2->Aid(0), snewpred2->Aid(1), status);

		t->AddSegment(snewpred2);
		t->SetCounters();
		t->FitTrackKFS();
		printf("## MICROTRACK Bottom found ## \n");
		printf(" MT Bot : X = %8.1f, Y = %8.1f, TX = %7.4f, TY = %7.4f, W = %2.0f, Chi2 = %5.2f\n",fnds2.X(),fnds2.Y(),fnds2.TX(),fnds2.TY(),fnds2.W(),fnds2.Chi2());
		if(gEve) gEve->SetStatusLine(Form("Microtrack Bottom found on pl %d. (W,Chi2)=(%d,%.2f)", ipl, (int)fnds2.W(), fnds2.Chi2()));
	}
	
	printf("\nMicrotrack search end.\n\n");

	return status;
}
void EdbEDATrackSetTab::ChangeScanID(){
	EdbID& id = eSet->GetID();
	InputID(Form("Input Scan ID for TrackSet:%s",eSet->GetName()), id);
	SetAttributeGUIs(eSet);
}

EdbEDATrackSetTab::EdbEDATrackSetTab (void) : eEnable(kTRUE), fFontWindow(NULL){
	TEveBrowser *browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kBottom);

	TGMainFrame *frame = new TGMainFrame(gClient->GetRoot());
	frame->SetWindowName("XX GUI");
	frame->SetCleanup(kDeepCleanup);

	frame->SetLayoutBroken(kTRUE);

	TGLabel *fLabel;
	TGTextButton *fb;

	int posy=10;
	int posx=10;
	int dx=0;
	fLabel = new TGLabel(frame,"TrackSet :");
	fLabel->MoveResize(posx,posy,dx=60,20);

	posx+=dx+10;
	// list box
	fTrackSet = new TGComboBox(frame);
	frame->AddFrame(fTrackSet);
	for(int i=0;i<gEDA->NTrackSets();i++){
		EdbEDATrackSet *set = gEDA->GetTrackSet(i);
		fTrackSet->AddEntry(set->GetName(),i);
	}
	fTrackSet->Select(0);
	fTrackSet->MoveResize(posx,posy,dx=60,20);
	fTrackSet->Connect("Selected(const char *)", "EdbEDATrackSetTab", this, "Selected(const char *)");

	/*
	posx+=dx;
	TGButton *fb;
	fb = new TGTextButton(fGroupFrame4,"Read TS");
	fb->MoveResize(posx,posy,dx=80,18);
	fb->Connect("Clicked()","EdbEDATrackSetTab", this,"ReadTotalScan()");
	*/
	
	posx+=dx+10;
	fScanID = new TGTextButton(frame,"hoge");
	fScanID->MoveResize(posx,posy,dx=150,20);
	fScanID->Connect("Clicked()", "EdbEDATrackSetTab", this, "ChangeScanID()");
	
/*
	fBrick = new TGNumberEntryField(frame, -1, gEDA->GetBrick());
	fBrick->MoveResize(posx,posy,dx=50,20);

	posx+=dx+5;
	fLabel = new TGLabel(frame,"Ver");
	fLabel->MoveResize(posx,posy,dx=25,20);
	posx+=dx+5;
	fVer     = new TGNumberEntryField(frame, -1, 100);
	fVer->MoveResize(posx,posy,dx=25,20);

	posx+=dx+5;
	fLabel = new TGLabel(frame,"Run");
	fLabel->MoveResize(posx,posy,dx=25,20);
	posx+=dx+5;
	fRun     = new TGNumberEntryField(frame, -1, 1);
	fRun->MoveResize(posx,posy,dx=35,20);
*/
	posx+=dx+5;
	fReadMicro = new TGCheckButton(frame,"Microtrk");
	fReadMicro->MoveResize(posx, posy, dx=70,20);
	fReadMicro->Connect("Toggled(bool)", "EdbEDATrackSetTab", this, "Update()");

	posx+=dx+5;
	fb = fPredScan = new TGTextButton(frame,"Read Pred. Scan");
	fb->MoveResize(posx,posy,dx=100,20);
	fb->Connect("Clicked()","EdbEDATrackSetTab", this,"ReadPredictionScan()");

	posx+=dx+5;
	fb = fReadText = new TGTextButton(frame,"Read File");
	fb->MoveResize(posx,posy,dx=90,20);
	fb->Connect("Clicked()","EdbEDATrackSetTab", this,"ReadFile()");

	posx+=dx+5;
	fb = fReadText = new TGTextButton(frame,"Read Text File");
	fb->MoveResize(posx,posy,dx=90,20);
	fb->Connect("Clicked()","EdbEDATrackSetTab", this,"ReadTextTracks()");

	posx+=dx+10;
	fb = new TGTextButton(frame,"Clear");
	fb->MoveResize(posx,posy,dx=50,20);
	fb->Connect("Clicked()","EdbEDATrackSetTab", this,"ClearTracks()");
	
	
	posy = 35;
	posx=10;

	fDataText = new TGLabel(frame);
	fDataText->MoveResize(posx,posy,dx=70,20);

	posx+=dx+10;
	fColorAuto = new TGCheckButton(frame,"Color auto");
	fColorAuto->MoveResize(posx, posy, dx=80,20);
	fColorAuto->Connect("Toggled(bool)", "EdbEDATrackSetTab", this, "ToggleAuto(bool)");

	posx+=dx;
	fTrackColor = new TGColorSelect(frame);
	fTrackColor->MoveResize(posx, posy, dx=35,20);
	fTrackColor->Connect("ColorSelected(Pixel_t)", "EdbEDATrackSetTab", this, "Update()");

	posx+=dx;
	fTrackWidth = new TGLineWidthComboBox(frame);
	frame->AddFrame(fTrackWidth);
	fTrackWidth->MoveResize(posx, posy, dx=60,20);
	fTrackWidth->Connect("Selected(int)", "EdbEDATrackSetTab", this, "Update()");

	posx+=dx+10;
	fLayerColor = new TGColorSelect(frame);
	fLayerColor->MoveResize(posx, posy, dx=35,20);
	fLayerColor->Connect("ColorSelected(Pixel_t)", "EdbEDATrackSetTab", this, "Update()");

	posx+=dx;
	fLayerWidth = new TGLineWidthComboBox(frame);
	frame->AddFrame(fLayerWidth);
	fLayerWidth->MoveResize(posx, posy, dx=60,20);
	fLayerWidth->Connect("Selected(int)", "EdbEDATrackSetTab", this, "Update()");


	posx+=dx+5;
	fLabel = new TGLabel(frame,"Layer length");
	fLabel->MoveResize(posx,posy,dx=70,20);
	posx+=dx+5;
	fLayerLength = new TGNumberEntryField(frame, -1, 150);
	fLayerLength->MoveResize(posx,posy,dx=30,20);
	fLayerLength->Connect("TextChanged(char *)", "EdbEDATrackSetTab", this, "Update()");

	posx+=dx+5;
	fTrackID = new TGCheckButton(frame,"ID");
	fTrackID->MoveResize(posx, posy, dx=35,20);
	fTrackID->Connect("Toggled(bool)", "EdbEDATrackSetTab", this, "Update()");

	posx+=dx+5;
	fTrackAngle = new TGCheckButton(frame,"Angle");
	fTrackAngle->MoveResize(posx, posy, dx=50,20);
	fTrackAngle->Connect("Toggled(bool)", "EdbEDATrackSetTab", this, "Update()");
/*
	posx+=dx+5;
	TGFontTypeComboBox *fFont = new TGFontTypeComboBox(frame);
	frame->AddFrame(fFont);
	fFont->MoveResize(posx, posy, dx=60,20);
	fFont->Connect("Selected(int)", "EdbEDATrackSetTab", this, "Update()");
*/

	posx+=dx+5;
	fb = new TGTextButton(frame,"Font");
	fb->MoveResize(posx, posy, dx=35,20);
	fb->Connect("Clicked()", "EdbEDATrackSetTab", this, "MakeFontWindow()");


	posx+=dx+10;
	// list box
	fExtendMode = new TGComboBox(frame);
	frame->AddFrame(fExtendMode);
	fExtendMode->AddEntry("ExtendAuto",  kExtendAuto);
	fExtendMode->AddEntry("ExtendUpDown",kExtendUpDown);
	fExtendMode->AddEntry("ExtendDown",  kExtendDown);
	fExtendMode->AddEntry("ExtendUp",    kExtendUp);
	fExtendMode->AddEntry("No Extend",   kExtendNo);
	fExtendMode->MoveResize(posx,posy,dx=115,20);
	fExtendMode->Connect("Selected(const char *)", "EdbEDATrackSetTab", this, "Update()");



	eSet = gEDA->GetTrackSet(0);
	SetAttributeGUIs(eSet);
	
	frame->MapSubwindows();
	frame->Resize();
	frame->MapWindow();

	browser->StopEmbedding();
	browser->SetTabTitle("TrackSets", 2);

}

void EdbEDATrackSetTab::CloseFontWindow(){
	if(fFontWindow!=NULL) fFontWindow->CloseWindow();
	fFontWindow=NULL;
}


void EdbEDATrackSetTab::MakeFontWindow(){
	if(fFontWindow!=NULL){
		fFontWindow->RaiseWindow();
		return;
	}
	
	// main frame
	TGMainFrame *fMainFrame = new TGMainFrame(gClient->GetRoot(),10,300,kMainFrame | kVerticalFrame);
	fFontWindow = fMainFrame;
	
	//  fMainFrame->SetLayoutBroken(kTRUE);
	fMainFrame->SetWindowName("Font");

	TGLabel  *fLabel= new TGLabel(fMainFrame, "Size");
	fMainFrame->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
	fLabel->Resize(70,20);

	// Face Size combo
	fFontSize = new TGComboBox(fMainFrame);
	fMainFrame->AddFrame(fFontSize,new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
	fFontSize->Resize(90, 20);

	Int_t* fsp = &TGLFontManager::GetFontSizeArray()->front();
	Int_t  nums = TGLFontManager::GetFontSizeArray()->size();
	for(Int_t i= 0; i< nums; i++){
		fFontSize->AddEntry(Form("%-2d", fsp[i]), fsp[i]);
	}
	fFontSize->Select( eSet->GetTextSize());
	
	fFontSize->Connect("Selected(Int_t)", "EdbEDATrackSetTab", this, "SetFont()");


	fLabel= new TGLabel(fMainFrame, "Font");
	fMainFrame->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
	fLabel->Resize(70,20);

	// Font File combo

	fFontFile = new TGComboBox(fMainFrame);
	fMainFrame->AddFrame(fFontFile, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
	fFontFile->Resize(90, 20);

	TObjArray* farr = TGLFontManager::GetFontFileArray();
	TIter next(farr);
	TObjString* os;
	Int_t cnt = 0;
	while ((os = (TObjString*) next()) != 0)
	{
		fFontFile->AddEntry(Form("%s", os->GetString().Data()), cnt);
		cnt++;
	}
	fFontFile->Select( eSet->GetTextFont());
	fFontFile->Connect("Selected(Int_t)", "EdbEDATrackSetTab", this, "SetFont()");

	// Font Color combo
	fFontColor = new TGColorSelect(fMainFrame);
	fMainFrame->AddFrame(fFontColor, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
	fFontColor->MoveResize(35,20);
	fFontColor->SetColor( TColor::Number2Pixel(eSet->GetTextColor()));
	fFontColor->Connect("ColorSelected(Pixel_t)", "EdbEDATrackSetTab", this, "SetFont()");


	TGTextButton *fb = new TGTextButton(fMainFrame,"OK");
	fMainFrame->AddFrame(fb, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 5));
	fb->Resize(90,20);
	fb->Connect("Clicked()", "EdbEDATrackSetTab", this, "CloseFontWindow()");

	fMainFrame->SetMWMHints(kMWMDecorAll,
							kMWMFuncAll,
							kMWMInputModeless);
	fMainFrame->MapSubwindows();

	fMainFrame->Resize(fMainFrame->GetDefaultSize());
	fMainFrame->MapWindow();
	fMainFrame->Resize(100,170);
}

void EdbEDATrackSetTab::SetFont(){
	printf("SetFont size=%d, col=%d, font=%d\n", fFontSize->GetSelected(),
		TColor::GetColor(fFontColor->GetColor()),
		fFontFile->GetSelected());
	eSet->SetTextAttribute(
		fFontSize->GetSelected(),
		TColor::GetColor(fFontColor->GetColor()),
		fFontFile->GetSelected());
		
	gEDA->Redraw();

}



void EdbEDATrackSetTab::Selected(const char *txt){
	eSet = gEDA->GetTrackSet(const_cast<char*> (txt));
	printf("%s\n",txt);
	SetAttributeGUIs(eSet);
	if(fFontWindow!=NULL) fFontWindow->CloseWindow();
	fFontWindow=NULL;
}

void EdbEDATrackSetTab::Update(){
	if(eEnable==0) return;

	int track_color = TColor::GetColor(fTrackColor->GetColor());
	if(fTrackColor->IsEnabled()==kFALSE) track_color = -1;
	int track_width = fTrackWidth->GetSelected();
	int layer_color = TColor::GetColor(fLayerColor->GetColor());
	int layer_width = fLayerWidth->GetSelected();
	double track_length=-1;
	double layer_length = fLayerLength->GetNumber();
	eSet->SetTrackAttribute(track_color, track_width, track_length, layer_color, layer_width, layer_length);
	
	eSet->SetDrawMT( fReadMicro->GetState()==kButtonDown? kTRUE : kFALSE);
	eSet->SetDrawTrackID   ( fTrackID->GetState()==kButtonDown? kTRUE: kFALSE );
	eSet->SetDrawTrackAngle( fTrackAngle->GetState()==kButtonDown? kTRUE: kFALSE );
	
	gEDA->SetLayerLength(layer_length);
	eSet->SetExtendMode(fExtendMode->GetSelected());
	gEDA->Redraw();
}

void EdbEDATrackSetTab::SetAttributeGUIs(EdbEDATrackSet *set){
	Enable(kFALSE);
	
	EdbID& id = set->GetID();
	fScanID->SetText(Form("ScanID %d.%d.%d.%d", id.eBrick, id.ePlate, id.eMajor, id.eMinor));
	
	/*
	fBrick->SetNumber(set->GetID()->eBrick);
	
	fVer->SetNumber(set->GetID()->eMajor);
	fRun->SetNumber(set->GetID()->eMinor);
	
	fBrick->SetEnabled( strcmp(set->GetName(),"TS")==0 ? kFALSE : kTRUE);
	fVer  ->SetEnabled( strcmp(set->GetName(),"TS")==0 ? kFALSE : kTRUE);
	fRun  ->SetEnabled( strcmp(set->GetName(),"TS")==0 ? kFALSE : kTRUE);
	*/
	
	fReadMicro->SetState(set->GetDrawMT()? kButtonDown : kButtonUp);
	
	fPredScan ->SetEnabled( strcmp(set->GetName(),"TS")==0 ? kFALSE : kTRUE);
	fReadText ->SetEnabled( strcmp(set->GetName(),"TS")==0 ? kFALSE : kTRUE);
	
	fDataText->SetText(Form(" %d / %d", set->N(), set->NBase()));
	
	int    track_color, track_width, layer_color, layer_width, layer_color_MT;
	double track_length, layer_length;
	set->GetTrackAttribute(&track_color, &track_width, &track_length, &layer_color, &layer_width, &layer_length, &layer_color_MT);
	
	fColorAuto->SetState( track_color == -1 ? kButtonDown : kButtonUp);
	if(track_color == -1) fTrackColor->Disable();
	else {
		fTrackColor->Enable();
		fTrackColor->SetColor(TColor::Number2Pixel(track_color));
	}
	
	fTrackWidth->Select(track_width);
	fLayerColor->SetColor(TColor::Number2Pixel(layer_color));
	fLayerWidth->Select(layer_width);
	
	fTrackID->SetState   (set->GetDrawTrackID()?kButtonDown:kButtonUp);
	fTrackAngle->SetState(set->GetDrawTrackAngle()?kButtonDown:kButtonUp);
	fExtendMode->Select(set->GetExtendMode());
	
	Enable(kTRUE);
}

void EdbEDATrackSetTab::ReadPredictionScan(){
	
	
//	int ibrick = fBrick->GetIntNumber();
//	int ver = fVer->GetIntNumber();
//	int run = fRun->GetIntNumber();

	EdbID& id = eSet->GetID();
	int ret=InputID("Input ID for Prediction scan", id);
	if(ret==0){
		// "close" button was pushed. stop.
		printf("ReadPrediction canceled.");
		return;
	}

	bool read_micro = fReadMicro->GetState()==kButtonDown? kTRUE: kFALSE;
	eSet->SetDrawMT(read_micro);
	
	eSet->ReadPredictionScan(id);
//	eSet->ReadPredictionScan(*id);
//	problem when it read SF data.

	gEDA->Redraw();
	SetAttributeGUIs(eSet);
}

void EdbEDATrackSetTab::ReadTextTracks(){
	eSet->ReadTextTracks(NULL);
	gEDA->Redraw();
	SetAttributeGUIs(eSet);
}

void EdbEDATrackSetTab::ReadFile(){
	eSet->ReadFile();
	gEDA->Redraw();
	SetAttributeGUIs(eSet);
}

void EdbEDATrackSetTab::ClearTracks(){
	eSet->Clear();
	gEDA->Redraw();
	SetAttributeGUIs(eSet);
}
