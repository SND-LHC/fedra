#include"EdbEDA.h"

EdbEDA *gEDA=0;

ClassImp(EdbEDAIO)
ClassImp(EdbEDACamera)
ClassImp(EdbEDABrickData)
ClassImp(EdbEDASelection)
ClassImp(EdbEDA)
ClassImp(EdbEDAText)
ClassImp(EdbEDALine)
ClassImp(EdbEDAOverlay)


#include<TEveViewer.h>

using namespace EdbEDAUtil;




TGLViewer * EdbEDACamera::GetGLViewer(){
	// Get GLViewer which is infront.
	// when you make a new GLViewer, need to define same name for GLViewer and Tab.
	
	
	TGTabElement *tab=gEve->GetBrowser()->GetTabRight()->GetCurrentTab();
	
	TEveViewerList *viewers = gEve->GetViewers();

	TGLViewer *v = NULL;
	for(TEveElement::List_i it = viewers->BeginChildren(); it!=viewers->EndChildren(); it++){
		printf("%s\n", ((TEveViewer *)(*it))->GetName());
		if( strcmp(((TEveViewer *)(*it))->GetName(), tab->GetName())==0 ) {
			// compare the name of viewer and name of tab.
			//printf("GLViewer %s\n", tab->GetName());
			v = ((TEveViewer *)(*it))->GetGLViewer();
			printf("GLViewer in current Tab = %s\n", ((TEveViewer *)(*it))->GetName());
			break;
		}
	}
	if(v==NULL) {
		printf("No viewer found in current Tab. Set Default GLViewer.\n");
		v=gEve->GetDefaultGLViewer();
	}
	eViewer = v;
	
	return v;
}
void EdbEDACamera::SetCamera(int projection){
	switch(projection){
		case X: // X projection
			{
				eViewer->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
				TGLOrthoCamera& cam = (TGLOrthoCamera &)(eViewer->CurrentCamera());
				cam.Reset();
				cam.RotateRad( -3.14159*0.5, 0.0);
				cam.SetEnableRotate(1);
				printf("set camera : X projection\n");
			}
			break;
		case Y: // Y projection
			{	
				eViewer->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
				TGLOrthoCamera& cam = (TGLOrthoCamera &)eViewer->CurrentCamera();
				cam.Reset();
				cam.SetEnableRotate(1);
				printf("set camera : Y projection\n");
			}
			break;
		case Z: // Z
			{
				eViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
				TGLOrthoCamera& cam = (TGLOrthoCamera &)eViewer->CurrentCamera();
				cam.Reset();
				cam.SetEnableRotate(1);
				printf("set camera : Y projection\n");
			}
			break;
		case B: // Beam view
			{
				eViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
				TGLOrthoCamera& cam = (TGLOrthoCamera &)eViewer->CurrentCamera();
				cam.Reset();
				double tx = atan(gEDA->GetBeamTY()/gEDA->GetScaleZ());
				double ty = atan(gEDA->GetBeamTX()/gEDA->GetScaleZ());
				cam.RotateRad(-tx,-ty);
				cam.SetEnableRotate(1);
				printf("set camera : Beam projection. 0.058 rad tilted\n");
			}
			break;
		default : break;
	}
	(eViewer)->DoDraw(); // ???? doesn't work...
}

void EdbEDACamera::SetAnimationMode(){
	
	eAnimationMode = InputNumberInteger("Animation Mode? 0: rotation(default), 1: hand move like, 2: left-right.", eAnimationMode);
	eAnimationAmplitude = InputNumberReal("Amplitude?", eAnimationAmplitude);
	
}
void EdbEDACamera::Animation(){
	
	
	TGLViewer *v=eViewer;
	TGLCamera& cam=v->CurrentCamera();
	if(eAnimationMode==0) cam.RotateRad(0.0,0.02);
	if(eAnimationMode==1) cam.RotateRad( eAnimationAmplitude*sin(eTheta), eAnimationAmplitude*cos(eTheta+=0.04));
	if(eAnimationMode==2) cam.RotateRad(0, eAnimationAmplitude*sin(eTheta+=0.04));
	v->DoDraw();
}

void EdbEDACamera::StartAnimation(){
	if(NULL==eTimer){
		eTimer = new TTimer;
		eTimer->Connect("Timeout()","EdbEDA",this,"Animation()");
		eViewer->Connect("Destroyed()", "EdbEDA", this, "StopAnimation()");
	}
	eTimer->Start(50,0);
}
void EdbEDACamera::StopAnimation() { eTimer->TurnOff();}


void EdbEDACamera::SaveAnimation(const char *filename, int n, int interval){
	printf("gEDA->SaveAnimation(const char *filename, int nview=100, int interval=10)\n");
	if(filename==NULL){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Animated Gif", "*.gif", 0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		filename = fi->fFilename;
	}
	if(filename==NULL) return;

	n=InputNumberInteger("Input N frame per rotation.", n);
	
	TGLViewer *v=eViewer;
	TGLCamera& c=v->CurrentCamera();
	v->SavePicture(Form("%s",filename));
	for(int i=0;i<n;i++){
		printf("save animation %d/%d\r",i,n);
		if(eAnimationMode==0) c.RotateRad(0.0,2*3.14/n);
		if(eAnimationMode==1) c.RotateRad(eAnimationAmplitude*sin(eTheta+=2*3.14/n), 
		                                  eAnimationAmplitude*cos(eTheta));
		if(eAnimationMode==2) c.RotateRad(0,
		                                  eAnimationAmplitude*cos(eTheta+=2*3.14/n));
		v->DoDraw();
		v->SavePicture(Form("%s+%d",filename,interval));
	
	}
	v->SavePicture(Form("%s++",filename));
	printf("save animation done!!.\n");
}

void EdbEDACamera::SavePictures(){
	printf("Save Pictures : x.gif, y.gif, z.gif, b.gif\n");
	FILE *fp = fopen("x.gif","rb");
	if(fp) {
		fclose(fp);
		printf("file already exist. delete and re-try.\n");
		return;
	}
	
	TGLViewer *v = eViewer;
	gEDA->SetCamera(EdbEDACamera::X);
	v->SavePicture("x.gif");
	gEDA->SetCamera(EdbEDACamera::Y);
	v->SavePicture("y.gif");
	gEDA->SetCamera(EdbEDACamera::Z);
	v->SavePicture("z.gif");
	gEDA->SetCamera(EdbEDACamera::B);
	v->SavePicture("b.gif");
}




void EdbEDACamera::Snapshot(const char *filename){
	if(filename==NULL) {
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = {	"GIF", "*.gif", "JPG","*.jpg", "PNG","*.png", 0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		printf("Save current view : %s\n", fi->fFilename);
		filename=fi->fFilename;
	}
	if(filename==NULL) return;
	TGLViewer *v = eViewer;
	TString name=filename;
	v->SavePicture(filename);
}


EdbAffine2D * EdbEDABrickData::GetAffine(int ipl){
	int i;
	
	if(eScanSet){
		EdbPlateP *p = eScanSet->GetPlate(ipl);
		if(p!=NULL) return p->GetAffineXY();
		else return NULL;
	}
	else if (eDataSet){
		EdbDataSet *dset = eDataSet;
		for(i=0;i<dset->N();i++){
			EdbDataPiece *piece = dset->GetPiece(i);
			if(piece->Plate()==ipl){
				EdbLayer *l = piece->GetLayer(0);
				return l->GetAffineXY();
			}
		}
	}
	return NULL;
}

EdbPattern * EdbEDAIO::ReadCouples(int ipl, EdbPattern *pat0){
	printf("ReadCouples DS ipl=%d\n", ipl);
	return ReadCouplesPID(gEDA->GetPID(ipl), pat0);
}



EdbPattern * EdbEDAIO::ReadCouplesPID(int pid, EdbPattern *pat0){
			
	EdbDataSet *dset = gEDA->GetDataSet();
	EdbDataPiece *piece = dset->GetPiece(pid);
	EdbLayer *l = piece->GetLayer(0);
	EdbPattern *pat = new EdbPattern( 0.,0., l->Z());
	pat->SetPID(pid);
	piece->InitCouplesTree("READ");
	piece->GetCPData_new(pat);
	piece->CloseCPData();
	pat->SetSegmentsZ();
	pat->SetID(pid);
	pat->SetSegmentsPID();
	pat->Transform(    l->GetAffineXY()   );
	pat->TransformA(   l->GetAffineTXTY() );
	pat->TransformShr( l->Shr()  );
	
	if(pat0==0) return pat;
	
	EdbScanCond cond;
	
	// copy segments.
	for(int i=0;i<pat->N();i++){
		EdbSegP *s = pat->GetSegment(i);
		if(s->Plate()==0) s->SetPlate(gEDA->GetIPL(pid));
	
		// fill COV for vertexing
		s->SetErrors();
		cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
	
		pat0->AddSegment(*s);
	}
	delete pat;
	return pat0;
}


EdbPattern * EdbEDAIO::GetPatternIPL(int ipl){
	// return Pattern with Couples.
	// if Couples are not loaded yet, load it in PVRec for TS.
	
	EdbEDATrackSet *set = gEDA->GetTrackSet("TS");
	EdbPVRec *pvr = set->GetPVRec();

	EdbPattern *pat=NULL;
	
	for(int i=0;i<pvr->Npatterns();i++){
		EdbPattern *p = pvr->GetPattern(i);
		if(p==NULL) continue;
		if(p->N()>0){
			if(p->GetSegment(0)->Plate()==ipl) {
				pat=p;
				break;
			}
		}
	}
	
	if(pat==NULL) {
		printf("No pattern in TS, pl=%d.\n", ipl);
		return NULL;
	}
	
	// Check if couples data were already loaded.
	if( IsIncludeCouples(pat)==0 ) {
		
		if     (set->GetDataSet()!=NULL) gEDA->ReadCouples(ipl,pat);
		else if(set->GetScanSet()!=NULL) {
			printf("ReadCouples SS ipl=%d\n", ipl);
			EdbID id = *(EdbID *)set->GetScanSet()->eIDS.At(0);
			id.ePlate=ipl;
			id.Print();
			float z = pat->Z();
			gEDA->ScanProc()->ReadPatCP(*pat,id,eCutCP);
			pat->SetZ(z);
			pat->SetSegmentsZ();
			printf("z = %lf\n", pat->Z());
			pat->SetSegmentsPlate(ipl);
		}
	}
	return pat;
}


void EdbEDAIO::WriteFeedbackFile(const char *filename){
	// write feedback file
	// if eOutputFileMode == 0 (kBern), output comments with "//" escape sequence (default).
	//                                  also track id is not 1,2,3,4, but id of track.
	// if eOutputFileMode == 1 (kOtherLabs). doesn't output comments.
	//
	// to convert Bern format, use
	// gawk -f ConvertFeedback.awk tmp.feedback
	//
	// -- ConvertFeedback.awk --
	// BEGIN{ itrk=1;}
	// !/\/\//&&NF==27{print itrk++,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27;}
	// !/\/\//&&NF==9{print $1,$2,$3,$4,$5,$6,$7,$8,$9;}
	// !/\/\//&&NF==10{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10;}
	
	
	if( gEDA->GetVertexSet()->N() == 0 ){
		ErrorMessage("No vertex is found. stop.");
		return;
	}
	
	if(filename==NULL){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Feedback file", "*.feedback", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		if(fi->fFilename!=NULL) filename = fi->fFilename;
	}
	
	FILE *fp;
	if(NULL==filename||strlen(filename)==0) {
		return;
	}
	else {
		fp=fopen(filename,"wt");
		if(NULL==fp){
			printf("Couldn't open file : %s\n", filename);
			return;
		}
	}

	printf("--------- Feedback format---------------\n");
	printf("VID        X        Y        Z 1ry charm tau Ndw Nup OutOfBrick\n");
	printf("trkid v1 v2        x        y        z      tx      ty   ip1   ip2 ");
	printf("    p  pmin  pmax mn pa sb dk ou ps ");
	printf(" n RmaxT RmaxL   rmsT   rmsL pl1 pl2 res\n");
	printf("        ipl        x        y        z      tx      ty type irec ngrains\n");
	
	if(eOutputFileMode==kBern){
		fprintf(fp,"// --------- Feedback format---------------\n");
		fprintf(fp,"// VID X Y Z 1ry charm tau Ndw Nup OutOfBrick\n");
		fprintf(fp,"//   OutOfBrick   : 1= vertex in dead material, 0= others\n");
		fprintf(fp,"// \n");
		fprintf(fp,"// Trkid v1 v2 x y z tx ty ip1 ip2 p pmin pmax mn pa sb dk ou ps n RmaxT RmaxL rmsT rmsL pl1 pl2 res\n");
		fprintf(fp,"//   Manual       : 0= auto, 1= manually added\n");
		fprintf(fp,"//   Particle     : 1= muon, 2= charm, 3= electron, 4= e+e- pair\n");
		fprintf(fp,"//   Scanback     : 1= scanback, 0=not scanback\n");
		fprintf(fp,"//   Darkness     : 0= MIP, 1= BLACK, 2= GRAY\n");
		fprintf(fp,"//   OutOfBrick   : 1= pass through up, 2= edge out, 0 others\n");
		fprintf(fp,"//   Last plate   : 1= pass through down, n= edge out last plate number, 0 others\n");
		fprintf(fp,"//   Out_flag(res): 1= 1ry vertex track, 2= e+e- pair, 3= Low momentum (DTheta_RMS >20 mrad),\n");
		fprintf(fp,"//                  4= Scan Forth to be done, 5= Scan Forth  done\n");
		fprintf(fp,"// \n");
		fprintf(fp,"// Ipl x y z tx ty type irec ngrains\n");
		fprintf(fp,"//   Type         : 0 BT, 1 microtrack top, 2 microtrack bottom\n");
		fprintf(fp,"//   Irec         : 0 Automatic, 1 SB, 2 Manual\n");
		fprintf(fp,"// \n");
		fprintf(fp,"///////////////////////////////////////////// \n");
		fprintf(fp,"// ---------- Vertex Data ---------------- // \n");
		fprintf(fp,"///////////////////////////////////////////// \n");
	}
	// Vertices
	
	TObjArray *vertices = (gEDA->GetVertexSet())->GetVertices();
	
	printf("%d vertices\n", vertices->GetEntriesFast());
	EdbVertex *v_1ry=NULL;
	double z=1e9;
	for(int i=0;i<vertices->GetEntriesFast();i++){
		EdbVertex *v = (EdbVertex *) vertices->At(i);
		if(v->Z()<z) {
			z = v->Z();
			v_1ry=v;
		}
	}
	int vtxpl = gEDA->GetIPLZ(v_1ry->Z());
	printf("Vertex plate = %d. z=%7.1lf\n", vtxpl, gEDA->GetZ(vtxpl));
	if(eOutputFileMode==kBern) fprintf(fp, "// primary vertex plate = %d. dz = %f\n", vtxpl, gEDA->GetZ(vtxpl)-v_1ry->Z());
	
	if(eOutputFileMode==kBern) fprintf(fp, "// VID     X        Y        Z 1ry charm tau Ndw Nup OutOfBrick\n");
	
	for(int i=0; i<vertices->GetEntriesFast(); i++){
		EdbVertex *v = (EdbVertex *) vertices->At(i);
		v->SetID(i+1);
		int ISPRIMARY = v==v_1ry ? 1 : 0;
		int ISCHARM = 0;
		int ISTAU = 0;
		int Nup=0, Ndown=0;
		for(int j=0;j<v->N();j++){
			EdbTrackP *t = v->GetTrack(j);
			if(t->GetSegmentFirst()->Z()>v->Z()) Ndown++;
			if(t->GetSegmentLast()->Z()<v->Z()) Nup++;
		}
		int OutOfBrick = v->Z()<gEDA->GetZ(1)? 1:0; // 1: if vertex is outside of the brick. 0: others.
		
		if(eOutputFileMode==kBern){
			fprintf(fp,"%3d %8.1lf %8.1lf %8.1lf %3d %5d %3d %3d %3d %3d\n", 
					v->ID(), v->X(), v->Y(), v->Z(),ISPRIMARY,ISCHARM,ISTAU, Ndown, Nup, OutOfBrick);
		} else {
			fprintf(fp,"%d %8.1lf %8.1lf %8.1lf %3d %5d %3d %3d %3d %3d\n", 
					v->ID(), v->X(), v->Y(), v->Z(),ISPRIMARY,ISCHARM,ISTAU, Ndown, Nup, OutOfBrick);
		}
	}
	
	// Tracks + Segments
	int itrk=1;
	for(int i=0;i<gEDA->NTrackSets();i++){
		EdbEDATrackSet *set = gEDA->GetTrackSet(i);
		if(set->N()==0) continue;
		if(set->GetDraw()==0) continue;
		printf("TrackSet %s : %3d tracks\n", set->GetName(), set->N()); 
		if(eOutputFileMode==kBern){
			EdbID& id=set->GetID();
			fprintf(fp,"////////////////////////////////////////////////////////// \n");
			fprintf(fp,"// ---------- TrackSet : %s %6d 0 %d %4d ----------- //\n", set->GetName(), id.eBrick, id.eMajor, id.eMinor);
			fprintf(fp,"////////////////////////////////////////////////////////// \n");
		}
		for(int j=0;j<set->N();j++){
			EdbTrackP *t = set->GetTrack(j);
			EdbSegP *sf = t->GetSegmentFirst();
			EdbSegP *sl = t->GetSegmentLast();
			EdbVertex * v1 = NULL;
			EdbVertex * v2 = NULL;
			for(int k=0;k<vertices->GetEntriesFast();k++){
				EdbVertex *v = (EdbVertex *) vertices->At(k);
				for(int l=0;l<v->N();l++){
					EdbTrackP *tv = v->GetTrack(l);
					if(t==tv){
						if(NULL==v1) v1 = v;
						else {
							if(v1->Z()>v->Z()) { // v1 should 1ry.
								v2 = v1;
								v1 = v;
							}
							else v2 = v;
			}	}	}	}
			if(v1==NULL) v1=v_1ry; // all tracks are attached to 1ry vertex at this moment.
			
			// t  : pointer to the track
			// v1 : upstream vertex
			// v2 : downstream vertex
			// sf : first segment to be used for vertex calculation
			// fp : file pointer to print data.
			
			int ivtx1 = v1 ? v1->ID() : -1;
			int ivtx2 = v2 ? v2->ID() : -1;
			
			float ip1= v1? CalcIP(sf,v1) : -1.0; // v1->DistSeg(sf): -1.0;
			float ip2= v2? CalcIP(sl,v2) : -1.0; // v2->DistSeg(sl): -1.0;
			
			double p,pmin,pmax;
			CalcP(t, p, pmin, pmax); // pmin = pmin, pmax = pmax
			
			// Flags.... same as Feedback file.
			int MAN = 0;
			int SCANBACK=0, idscanback=-1;
			int MUON    =0;
			int DARKNESS=0; // no function to define.
			
			if(gEDA->GetTrackSet("MN")->FindTrack(t)) MAN=1;
			
			EdbTrackP * sbt = set==gEDA->GetTrackSet("SB")? t : gEDA->CheckScanback(t); // check if the track is a scanback track
			if(sbt!=NULL) idscanback = sbt->ID();
			if (idscanback!=-1) SCANBACK = 1;
			if (idscanback==gEDA->IdMuon()&&gEDA->IdMuon()!=-1) MUON =1;
			
			int OutOfBrick=0; // OutOfBrick: 1= pass through, 2=edge out, 0 others
			int LastPlate=0;  // Last plate: 1= pass though, n=edge out track last plate, 0 others 
			
			int nseg = t->N();
			
			// Small Kink search

			// remove Micro-tracks.
			double RmaxT=0.0, RmaxL=0.0, Rmax=0.0, dtmaxt=0.0, dtmaxl=0.0, rmstmax=-1.0, rmslmax=-1.0;
			int iplRmax1=0, iplRmax2=0;
			
			EdbTrackP *tt = new EdbTrackP;
			tt->Copy(*t);
			tt->SetCounters();
			
			for(int k=0;k<tt->N()&&tt->N()>1;k++){
				EdbSegP *s = tt->GetSegment(k);
				if(s->Side()!=0) tt->RemoveSegment(s); // if it's microtrack, remove.
			}
			
			int n = tt->N();
			
			if(n<=1) {
				// cannot calculate rms.
				RmaxT=RmaxL=-1.;
				iplRmax1=iplRmax2=-1;
			}
			else if( n<4 ){
				// if nseg <=3, cannot calculate kink, just calculate rms but not R.
				RmaxT=RmaxL=-1.;
				iplRmax1=iplRmax2=-1;
			}
			else {
				// Calculate rms and R.
				for(int k=0; k<n-1;k++){
					EdbSegP *ss1 = tt->GetSegment(k);
					EdbSegP *ss2 = tt->GetSegment(k+1);
					
					int ipl1 = ss1->Plate();
					int ipl2 = ss2->Plate();

					// if ss2 is on pl 57, skip.
					if(ipl2>=gEDA->GetLastPlate()) continue;
					// if ss1 is more downstream than 3 plates from vertex, skip. search only kink between 1-x,2-x,3-x.
					if(ipl1>vtxpl+2) continue;
					
					double dtt,dtl;
					CalcDTTransLongi(ss1, ss2, &dtt, &dtl);
					dtt=fabs(dtt);
					dtl=fabs(dtl);
					
					double r,rt,rl;
					int ndata;
					DTRMSTLGiven1Kink(tt, k, &r, &rt, &rl, &ndata);

					double RT =dtt/rt/sqrt((double)ipl2-ipl1);
					double RL =dtl/rl/sqrt((double)ipl2-ipl1);
					double R=sqrt(RT*RT+RL*RL);
					printf("ipl %2d-%2d dtt=%7.4lf dtl=%7.4lf rmst=%7.4lf rmsl=%7.4lf RT=%5.2lf RL=%5.2lf n=%d\n", 
							ipl1, ipl2, dtt, dtl, rt, rl, RT,RL,ndata);
					if(Rmax<R){
						Rmax  = R;
						RmaxT = RT;
						RmaxL = RL;
						iplRmax1=ipl1;
						iplRmax2=ipl2;
						dtmaxt = dtt;
						dtmaxl = dtl;
						rmstmax = rt;
						rmslmax = rl;
					}
				}
			}
			// nseg in 5 plate.
			int nin5pl = 0;
			for(int k=0;k<n;k++){
				EdbSegP *s = tt->GetSegment(k);
				if(s->Plate()>=vtxpl && s->Plate()<vtxpl+5) nin5pl++;
			}
			if( nin5pl<=2 ){
				// nseg <=2 in first 5 plate is out of criteria. just calculate RMS but not R.
				RmaxT=RmaxL=-1.;
				iplRmax1=iplRmax2=-1;
			}
			
			double rms, rmst, rmsl;
			if(n<=1) rms=rmst=rmsl=-1.;
			else DTRMSTL(tt, &rms, &rmst, &rmsl);
			
			delete tt;
			
			int result=1;
			
			const char *comment = set->GetComment(t);
			if(eOutputFileMode==kBern){
				// Bern format
				if(eOutputFileMode==kBern){
					fprintf(fp,"// tid v1 v2        x        y        z      tx      ty    ip1    ip2 ");
					fprintf(fp,"    p  pmin  pmax mn pa sb dk of ps ");
					fprintf(fp," n RmaxT RmaxL    rmsT    rmsL pl1 pl2 res\n");
				}
				fprintf(fp, "%6d %2d %2d %8.1f %8.1f %8.1f %7.4f %7.4f %6.1f %6.1f ",
					t->ID(), ivtx1, ivtx2, sf->X(), sf->Y(), sf->Z(), sf->TX(), sf->TY(), ip1, ip2);
				
			} else {
				// Official
				fprintf(fp, "%d %2d %2d %8.1f %8.1f %8.1f %7.4f %7.4f %6.1f %6.1f ",
					itrk++, ivtx1, ivtx2, sf->X(), sf->Y(), sf->Z(), sf->TX(), sf->TY(), ip1, ip2);
			}
			fprintf(fp, "%5.1lf %5.1lf %5.1lf %2d %2d %2d %2d %2d %2d ", 
				p,pmin,pmax, MAN, MUON, SCANBACK, DARKNESS, OutOfBrick, LastPlate);
			fprintf(fp, "%2d %5.2lf %5.2lf %7.4lf %7.4lf %3d %3d %3d\n",
				nseg, RmaxT, RmaxL, rmst, rmsl, iplRmax1, iplRmax2, result);

			// print on the shell.
			fprintf(stdout, "%5d %2d %2d %8.1f %8.1f %8.1f %7.4f %7.4f %6.1f %6.1f ",
				t->ID(), ivtx1, ivtx2, sf->X(), sf->Y(), sf->Z(), sf->TX(), sf->TY(), ip1, ip2);
			fprintf(stdout, "%5.1lf %5.1lf %5.1lf %2d %2d %2d %2d %2d %2d ", 
				p,pmin,pmax, MAN, MUON, SCANBACK, DARKNESS, OutOfBrick, LastPlate);
			fprintf(stdout, "%2d %5.2lf %5.2lf %7.4lf %7.4lf %3d %3d %3d\n",
				nseg, RmaxT, RmaxL, rmst, rmsl, iplRmax1, iplRmax2, result);
			
			// Check if the rms is smaller than angular resolution or not.
			// if it's smaller, then print a note in the feedback file.
			
			// minimum value of rms (angular resolution).
			double angle = sqrt(t->TX()*t->TX()+t->TY()*t->TY());
			double resT  = 0.0015*sqrt(2.0);
			double resL  = 0.0015*(1+angle*5)*sqrt(2.0);
			double RresT = RmaxT>0 ? dtmaxt/resT : -1.;
			double RresL = RmaxL>0 ? dtmaxl/resL : -1.;

			if( eOutputFileMode==kBern ){
				if( ( rmst<resT || rmsl<resL ) && rmst>0)
				     fprintf(fp, "//     rms for R = %7.4lf %7.4lf , note : rms is smaller than angular resolution.  "
							    "RresT RresL resT resL = %5.2lf %5.2lf %7.4lf %7.4lf\n", 
							    rmstmax, rmslmax, RresT, RresL, resT, resL);
				else fprintf(fp, "//     rms for R = %7.4lf %7.4lf                                                   "
							    "RresT RresL resT resL = %5.2lf %5.2lf %7.4lf %7.4lf\n", 
							    rmstmax, rmslmax, RresT, RresL, resT, resL);
				
				if( nin5pl<3 ) fprintf(fp, "//     nBT = %d. nBT in first five plate = %d.\n", n, nin5pl);
				fprintf(fp,      "//     track comment = \"%s\"\n", comment);
			}
			
			// Fill microtrack information
			
			int irec = 0;    // 0:TS, 1:PredScan, 2: Manual
			TString name = set->GetName();
			if( name.BeginsWith("TS") ) irec=0;
			if( name.BeginsWith("SB") ) irec=1;
			if( name.BeginsWith("SF") ) irec=1;
			if( name.BeginsWith("MN") ) irec=2;
			
			for(int k=0;k<t->N();k++){
				EdbSegP *s = t->GetSegment(k);

				int type = s->Side(); // 0:BT, 1:MTtop, 2:MTbottom
				EdbID idseg=s->ScanID();
				EdbEDATrackSet *setseg = gEDA->GetTrackSet( idseg);
				if(setseg){
					TString name = setseg->GetName();
					if( name.BeginsWith("TS") ) irec=0;
					if( name.BeginsWith("SB") ) irec=1;
					if( name.BeginsWith("SF") ) irec=1;
					if( name.BeginsWith("MN") ) irec=2;
				}
				
				if(eOutputFileMode==kBern){
					// Bern format
					if(k==0) fprintf(fp,"//       ipl        x        y        z      tx      ty type irec ngrains\n");
					fprintf(fp,"%12d %8.1f %8.1f %8.1f %7.4f %7.4f %4d %4d %4d\n",
						s->Plate(), s->X(), s->Y(), s->Z(), s->TX(), s->TY(), type, irec, (int) s->W());
				} else {
					// Official
					fprintf(fp,"%d %8.1f %8.1f %8.1f %7.4f %7.4f %4d %4d %4d\n",
						s->Plate(), s->X(), s->Y(), s->Z(), s->TX(), s->TY(), type, irec, (int) s->W());
				}
			}
			
		}
	}
	
	// Scanning area
	
	EdbEDATrackSet *set = gEDA->GetTrackSet("TS");
	EdbEDAAreaSet *aset = set->GetAreaSet();
	if(aset->N()!=0){
		fprintf(fp,"//////////////////////////////////////////\n");
		fprintf(fp,"// ------ Area information : %s ------- //\n", set->GetName());
		fprintf(fp,"//////////////////////////////////////////\n");
		for(int i=0;i<aset->N();i++){
			EdbEDAArea *a = aset->GetArea(i);
			fprintf(fp,"// area : TS %2d %lf %lf %lf %lf %lf\n", 
				a->Plate(), a->Xmin(), a->Xmax(),
				a->Ymin(), a->Ymax(), a->Z());
		}
	}
	
	fclose(fp);

	printf("file %s was written.\n", filename);

}

void EdbEDAIO::WriteFeedbackFile2008(EdbVertex *v, const char *filename){

// write feedback file. original code from Frederic.
// the flags are correct if you work in the correct directory.

	// extract vertex position and multiplicity
	// extract vertex tracks infos
	if(v==NULL) v = gEDA->GetSelectedVertex();
	if(filename==NULL){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Feedback file", "*.feedback", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		if(fi->fFilename!=NULL) filename = fi->fFilename;
	}
	
	if(NULL==filename||strlen(filename)==0) { printf("Cancel to write feedback\n"); return;}
	FILE *fp=fopen(filename,"wt");
	if(NULL==fp){
		printf("Couldn't open file : %s\n", filename);
		return;
	}


	// put first the primary vertex

	fprintf(fp,"%i %.1f %.1f %.1f %i %i %i \n",1,v->X(),v->Y(),v->Z(),1,0,0);

	// tracks from primary vertex
	for(int j=0; j<v->N(); j++){
		int MANUAL=0, MUON=0, SCANBACK=0, DARKNESS=0;
		EdbTrackP *t = v->GetTrack(j);
		EdbSegP *s = t->GetSegmentFirst();

		// check if the track is a scanback track
		int idscanback = -1;
		EdbTrackP *tsb = gEDA->CheckScanback(t);
		if(tsb) idscanback = tsb->ID();

		int idmuon     = IdMuon();
		if (idscanback!=-1) SCANBACK = 1;
		if (idscanback==idmuon&&idmuon!=-1) {
			MUON =1;
			printf("      Track: %d is muon\n", t->ID());
		}

		double p, pmin, pmax;
		CalcP(t, p, pmin, pmax); // pmin = pmin, pmax = pmax

		fprintf(fp,"%4i %i %i %8.1f %8.1f %8.1f %7.4f %7.4f %6.2f %3.1f %6.2lf %5.2lf %8.1lf %i %i %i %i \n",
		t->ID(),1,-1,s->X(),s->Y(),s->Z(),s->TX(), s->TY(), v->Impact(j),-1.,p,pmin,pmax,MANUAL,MUON,SCANBACK,DARKNESS);

	}
	printf("-- FeedbackFile : %s--\n", filename);
	fclose(fp);

#ifdef _WINDOWS
	char buf[256];
	sprintf(buf, "type %s",filename);
	for(int i=0;i<strlen(buf);i++) if(buf[i]=='/') buf[i]='\\';
	gSystem->Exec(buf);
#else
	gSystem->Exec(Form("cat %s",filename));
#endif
}

void EdbEDAIO::ReadFeedbackFile(const char *filename){
	// Read feedback file format (ver 2009 Oct).
	
	gEDA->Reset();
	
	
	char FlagPart[][10]={"", "mu", "charm", "e", "e-pair"};
	char FlagCS[][10]   ={"", "cs"};
	
	// if filename is not given, open file browser.
	if(filename==NULL){
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Feedback File", "*.feedback", "All files","*",0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve?gEve->GetMainWindow():0, kFDOpen, fi);
		filename = fi->fFilename;
	}
	
	if(NULL==filename) return;
	
	FILE *fp = fopen(filename,"rt");
	if(fp==NULL) {
		ErrorMessage(Form("File open error : %s . stop.\n", filename));
		return;
	}
	
	EdbPVRec *pvr = new EdbPVRec;
	EdbScanCond cond;
	EdbTrackP *t=0;
	EdbEDAText *tx = 0;
	
	char buf[256];
	
	char asetname[10];
	EdbEDAAreaSet *aset=NULL;
	
	while(fgets(buf,sizeof(buf),fp)){
		TString line=buf;
		printf("%s", line.Data());

		// track comment = ""
		if(line.Contains("track comment")){
			int istart = line.Index("\"");
			int iend   = line.Index("\"",1,istart+1, TString::kExact);
			if(istart!=kNPOS && iend!=kNPOS && istart+1!=iend){
				char buf[256];
				strncpy(buf, line.Data()+istart+1, iend-istart-1);
				buf[iend-istart-1]='\0';
				tx->AddText(buf);
			}
		}
		
		if(line.Contains("Area information :")){
			aset = new EdbEDAAreaSet;
			aset->SetDraw(kTRUE);
			sscanf( line.Data()+line.Index(":")+1, "%s", asetname);
			printf("Area infor %s\n", asetname);
		}
		if(aset){
			if(line.Contains(Form("area : %s", asetname))){
				int ipl;
				double xmin, xmax, ymin, ymax, z;
				sscanf(line.Data() + line.Index(Form("area : %s", asetname)) + strlen(Form("area : %s", asetname)),
					"%d %lf %lf %lf %lf %lf", &ipl, &xmin, &xmax, &ymin, &ymax, &z);
				aset->AddArea(ipl, xmin, xmax, ymin, ymax, z);
			}
		}
		
		// Remove comments
		line.ReplaceAll("\t"," ");
		int iposcmt = line.Index("//");
		if(iposcmt!=kNPOS) line.Remove(iposcmt, line.Length()-iposcmt);
		
		// Check number of columns.
		TObjArray *tokens = line.Tokenize(" ");
		int ntokens = tokens->GetEntries();
		delete tokens;
		
		
		if(ntokens==0) continue;
		
		else if( ntokens == 10 ){
			// Vertices
			float x,y,z; int id,isprimary,istau, ischarm;
			sscanf(line.Data(),"%d %f %f %f %d %d %d", &id, &x, &y, &z, &isprimary, &ischarm, &istau);
			EdbVertex *v = new EdbVertex();
			v->SetXYZ(x,y,z); v->SetID(id);
			v->SetFlag(isprimary);
			pvr->AddVertex(v);
			printf("Vertex %d %f %f %f\n",v->ID(), v->X(), v->Y(), v->Z());
		}

		else if( ntokens == 27 ){
			// Tracks
			float x,y,z,ax,ay, ip_upvtx, ip_downvtx,  p,pmin,pmax;
			int id_track, id_upvtx, id_downvtx,  manual, particle_id, scanback, darkness;
			int OutOfBrick, LastPlate, nseg, iplRmax1, iplRmax2, result;
			float RmaxT, RmaxL, rmst, rmsl;
			
			int n,nn;
			sscanf(line.Data(),         "%d %d %d %f %f %f %f %f %f%n", &id_track, &id_upvtx, &id_downvtx, &x, &y, &z, &ax, &ay, &ip_upvtx, &nn);
			sscanf(line.Data()+(n=nn),  "%f %f %f %f %d %d %d %d %d%n", &ip_downvtx, &p,&pmin,&pmax, &manual, &particle_id, &scanback, &darkness, &OutOfBrick, &nn);
			sscanf(line.Data()+(n+=nn), "%d %d %f %f %f %f %d %d %d",   &LastPlate, &nseg, &RmaxT, &RmaxL, &rmst, &rmsl, &iplRmax1, &iplRmax2, &result);
			
			// Create Track. "t" is defined out of main loop.
			t = new EdbTrackP;
			t->Set(id_track, x, y, ax, ay, 0, 0);
			t->SetZ(z);
			t->SetTrack(id_track);
			pvr->AddTrack(t);
			
			tx = new EdbEDAText(t->X(), t->Y(), t->Z(), "[");
			if(0<=particle_id && particle_id<=4) { if(tx->N()>1)tx->AddText(", "); tx->AddText(FlagPart[particle_id]);}
			if(scanback==1)                      { if(tx->N()>1)tx->AddText(", "); tx->AddText(FlagCS[scanback]);}
			tx->AddText("]");
			if(tx->N()==2) tx->SetText("");
			tx->SetReference(tx);
			gEDA->AddDrawObject(tx);
			
			// fill COV for vertexing
			t->SetErrors();
			cond.FillErrorsCov(t->TX(), t->TY(), t->COV());
			
			// associate to vertex.
			for(int i=0; i<pvr->Nvtx(); i++){
				EdbVertex *v = pvr->GetVertex(i);
				if(id_upvtx==v->ID()||id_downvtx==v->ID()){
					EdbVTA *vta = new EdbVTA(t, v);
					vta->SetZpos( t->Z()>v->Z() ? 1 : 0);
					vta->SetFlag(2);
					vta->SetImp( id_upvtx==v->ID()? ip_upvtx : ip_downvtx);
					v->AddVTA(vta);
				}
			}
		}
		
		else if( ntokens ==  9 ){
			// Segments
			int ipl, type, irec, ph;
			float x,y,z,ax,ay;
			sscanf(line.Data(),"%d %f %f %f %f %f %d %d %d", &ipl, &x, &y, &z, &ax, &ay, &type, &irec, &ph);
			
			EdbSegP *s = new EdbSegP(t->ID(),x,y,ax,ay,0,0);
			s->SetZ(z);
			s->SetPID(ipl);
			s->SetPlate(ipl);
			s->SetW(ph);
			s->SetTrack(t->ID());
			s->SetSide(type);

			// fill COV for vertexing
			s->SetErrors();
			cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
			
			// Add segment in PVRec and Track, keeping consistency of pointer in them.
			EdbSegP *s_in_pattern = pvr->AddSegment(*s);
			t->AddSegment( s_in_pattern);
			
			if(NULL!=tx){
				EdbSegP seg(*s);
				seg.PropagateTo(seg.Z()+1000);
				tx->SetXYZ(seg.X(),seg.Y(),seg.Z());
			}
		}
		
		
	}
	
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetPID(i);
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetSegmentsPID();
	for(int i=0;i<pvr->Ntracks(); i++) pvr->GetTrack(i)->SetCounters();
	
	printf("\nEdbPVRec summary. %d vertices, %d tracks.\n", pvr->Nvtx(), pvr->Ntracks());
	pvr->Print();
	
	fclose(fp);

	if(pvr==NULL) return;
	for(int i=0;i<pvr->Npatterns(); i++){
		EdbPattern *pat = pvr->GetPattern(i);
		pat->SetPID(gEDA->GetPID(pat->ID()));
	}

	gEDA->SetPVR(pvr);
	gEDA->GetTrackSet("TS")->AddTracksPVR(pvr);
	gEDA->GetVertexSet()->AddVertices(pvr->eVTX);
	
	if(aset){
		EdbEDATrackSet *set = gEDA->GetTrackSet(asetname);
		if(set) set->SetAreaSet(aset);
	}
	gEDA->Redraw();
	gEDA->UpdateScene();
}


EdbVertex * EdbEDAIO::ReadFeedbackFile2008(const char *filename){
	using namespace std;

	// if filename is not given, open file browser.
	if(filename==NULL){
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Feedback File", "*.feedback", "All files","*",0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDOpen, fi);
		filename = fi->fFilename;
	}
	
	if(NULL==filename) return NULL;
		
	// check the file.
	ifstream fin(filename);
	if(!fin) {
		printf("can not open file : %s\n", filename);
		return NULL;
	} printf("read feedback %s\n", filename);
	

	gEDA->ClearVertices();
	
	EdbEDATrackSet *set = gEDA->GetTrackSet("TS");
	set->ClearTracks();
	
	// read file line by line.
	char buf[512];
	while( fin && fin.getline(buf,512) ){
		string line(buf);
		
		// ignore comment
		if(line.find("//")!=line.npos) line[line.find("//")]='\0';
		
		// check the number of columnline.
		istringstream ss(line);
		string s2;
		int ncolumn;
		for(ncolumn=0;ss>>s2;ncolumn++){
			//cout <<s2 << endl;
		}
		
		// ncolumn == 7 --> vertex info.  ncolumn == 17 --> track info
		if(ncolumn==7){
			double x,y,z;
			int id,isprimary,istau, ischarm;
			sscanf(line.c_str(),"%d %lf %lf %lf %d %d %d", &id, &x, &y, &z, &isprimary, &ischarm, &istau);
			EdbVertex *v = new EdbVertex();
			v->SetXYZ(x,y,z);
			v->SetID(id);
			gEDA->AddVertex(v);
		}
		else if (ncolumn==17){
			double x,y,z,ax,ay, ip_upvtx, ip_downvtx;
			double p,pmin,pmax;
			int id_track, id_upvtx, id_downvtx;
			int manual, particle_id, scanback, darkness;
			
			istringstream ss(line);
			
			ss	>> id_track >> id_upvtx >> id_downvtx >> x >> y >> z >> ax >> ay >> ip_upvtx >> ip_downvtx 
				>> p >> pmin >> pmax >> manual >> particle_id >> scanback >> darkness;
			
			EdbSegP s(id_track, x, y, ax, ay);
			s.SetZ(z);
			
			EdbTrackP *t = set->FindTrack(&s);
			
			if(t) {
				printf("%s track %d\n", line.c_str(), t->ID());
				set->SetTrack(t);
				for(int i=0; i<gEDA->GetVertexSet()->N(); i++){
					EdbVertex *v = gEDA->GetVertexSet()->GetVertex(i);
					if(id_upvtx!=v->ID()) continue;
					EdbVTA *vta = new EdbVTA(t, v);
					vta->SetZpos(1);
					vta->SetFlag(2);
					vta->SetImp(EdbEDAUtil::CalcIP(&s,v));
					v->AddVTA(vta);
				}
			}
			else  printf("track not found. id = %d\n", id_track);
		}
	}
	gEDA->Redraw();
	if( gEDA->GetVertexSet()->N()) return gEDA->GetVertexSet()->GetVertex(0);
	else return NULL;
}

int EdbEDAIO::IdMuon(const char *filename){
	// code by GL
  char            buf[256];
  char            key[256];
  char m1[12], m2[10], m3[10],m4[10],m5[10],m6[10];
  static int val=-2;
  char mu[3];
  
  if(val>=0) return val;
  if(val==-1) return val;

  FILE *fp=fopen(filename,"r");
  if (fp==NULL)   {
    printf("ERROR open file: %s \n", filename);
    val = -1;
    return val;
  }
  
  // printf( "Read affine transformation from file: %s\n", fname );
  while( fgets(buf,256,fp)!=NULL ) {

    for( int i=0; i<(int)strlen(buf); i++ ) 
      if( buf[i]=='#' )  {
        buf[i]='0';                       // cut out comments starting from #
        break;
      }
    
   if( sscanf(buf,"%s",key)!=1 ) continue;
   if      ( !strcmp(key,"Prediction")   )
     {
       sscanf(buf+strlen(key),"%s %s %s %s %s %s ",m1,m2,m3,m4,m5,m6);
       mu[0] = m2[3];
       mu[1] = m2[4];
       mu[2] = m2[5];
       val = atol(mu);
       //       cout << m2[3] << " " << m2[4] << " " << m2[5]<< " " << val << endl;

     }
   if      ( !strcmp(key,"No")   )
     {
       val=-1;
     }

  }
  fclose(fp);	

  return val;
}

void EdbEDAIO::OpenTextEditor(const char *filename){
	TEveBrowser* browser = gEve->GetBrowser();
	browser->StartEmbedding(1);
	new TGTextEditor(filename, gClient->GetRoot());
	browser->StopEmbedding();
	
	unsigned int i,k=0;
	for(i=0;i<strlen(filename);i++){
		if(filename[i]=='\\'||filename[i]=='/') k=i+1;
	}
	browser->SetTabTitle(filename+k, 1);
}


EdbPVRec * EdbEDAIO::ReadFilteredText(const char *filename){
	int i,j;
	
	if(filename==NULL||strlen(filename)<1){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Filtered text", "*.txt", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDOpen, fi);
		if(fi->fFilename!=NULL) {
			filename = fi->fFilename;
		}
	}
	if(filename==NULL) return NULL;
	
	// check the file.
	ifstream fin(filename);
	if(!fin) {
		printf("can not open file : %s\n", filename);
		return NULL;
	} 
	else printf("Read Filtered Text File : %s\n", filename);
	

	TObjArray *vtxarr = new TObjArray;
	TObjArray *trkarr = new TObjArray;
	TObjArray *segarr = new TObjArray;	
	
	EdbPVRec *PVR = new EdbPVRec;
	
	EdbDataSet *dset = gEDA->GetDataSet();

	for(i=0;i<dset->N();i++){
		EdbDataPiece *piece = dset->GetPiece(i);
		EdbLayer* l = piece->GetLayer(0);
		EdbPattern *pat = new EdbPattern( 0.0, 0.0, l->Z());
		pat->SetPID(i);
		PVR->AddPattern(pat);
	}
	
	char buf[512];
	while( fin && fin.getline(buf,512) ){
		string line(buf);

		std::istringstream ssline(line);
		std::string column1;
		ssline >> column1;
		
		std::cout <<line <<std::endl;
		if(column1.find("//")!=column1.npos) continue;

		if(column1=="vertex"){
			int id;
			float x, y, z;
			ssline >> id >> x >> y >> z;
			EdbVertex *v = new EdbVertex;
			v->SetXYZ(x,y,z);
			v->SetID(id);
			vtxarr->Add(v);
			PVR->AddVertex(v);
		}
		
		if(column1=="track"){
			int id, flag, ipl1, ipl2, nseg, ivtx1, ivtx2;
			double tx,ty, ip_1ry, ip1, ip2;
			ssline >> id >> flag >> tx >> ty >> ipl1 >> ipl2 >> nseg >> ip_1ry >> ivtx1 >> ip1 >> ivtx2 >> ip2;
			EdbTrackP *t = new EdbTrackP;
			t->SetID(id);
			t->SetFlag(flag);
			for(i=0;i<vtxarr->GetEntriesFast();i++){
				EdbVertex *v = (EdbVertex *) vtxarr->At(i);
				if(v->ID()==ivtx1||v->ID()==ivtx2){
					EdbVTA *vta = new EdbVTA(t, v);
					vta->SetFlag(2);
					v->AddVTA(vta);
				}
			}
			trkarr->Add(t);
			PVR->AddTrack(t);
		}
		
		if(column1=="segment"){
			int itrk, ipl, iseg, sflag, tflag;
			float x, y, z, ax, ay, chi2, w;
			ssline >> itrk >> tflag >> ipl >> iseg >> x >> y >> z >> ax >> ay >> w >> chi2 >> sflag;
			
			EdbSegP *s = new EdbSegP (iseg, x, y, ax, ay, w, sflag);
			s->SetPID(gEDA->GetPID(ipl));
			s->SetZ(z);
			s->SetDZ(300);
			s->SetTrack(itrk);
			s->SetChi2(chi2);
			
			for(i=0;i<trkarr->GetEntriesFast();i++){
				EdbTrackP *t = (EdbTrackP *) trkarr->At(i);
				if(t->ID()==itrk&&t->Flag()==tflag) {
					t->AddSegment(s);
				}
			}
			
			EdbPattern *pat = PVR->GetPattern( gEDA->GetPID(ipl));
			pat->AddSegment(*s);
			segarr->Add(s);
			
		}
		
		if(column1=="area"){
			int ipl;
			double xmin,xmax,ymin,ymax,z;
			ssline >> ipl >> xmin >> xmax >> ymin >> ymax >> z;
			EdbEDAArea *a = new EdbEDAArea(ipl, xmin, xmax, ymin, ymax, z);
			gEDA->GetTrackSet("TS")->GetAreaSet()->AddArea(a);
		}
		
		
	}
	fin.close();

	printf("Vertices %4d\n", vtxarr->GetEntriesFast());
	printf("Tracks   %4d\n", trkarr->GetEntriesFast());
	printf("Segments %4d\n", segarr->GetEntriesFast());
	
	
	// Clear all track set
	for(i=0; i<gEDA->NTrackSets();i++){
		gEDA->GetTrackSet(i)->Clear();
	}
	
	for(i=0;i<trkarr->GetEntriesFast();i++){
		EdbTrackP *t = (EdbTrackP *) trkarr->At(i);
		t->SetCounters();
		EdbSegP *s = t->GetSegmentFirst();
		t->Set(t->ID(), s->X(), s->Y(), s->TX(), s->TY(), 0, t->Flag());
		t->SetZ(s->Z());
		
		if( 0 <= t->Flag()) gEDA->GetTrackSet("TS")->AddTrack(t);
		if( -100 < t->Flag() && t->Flag() <     0) gEDA->GetTrackSet("MN")->AddTrack(t);
		if( -200 < t->Flag() && t->Flag() <=  -99) gEDA->GetTrackSet("SB")->AddTrack(t);
		if(-1000 < t->Flag() && t->Flag() <= -200) gEDA->GetTrackSet("SF")->AddTrack(t);
	}
	
	EdbEDAVertexSet *vset = gEDA->GetVertexSet();
	gEDA->ClearSelectedVertex();
	vset->Clear();
	for(i=0;i<vtxarr->GetEntriesFast();i++){
		EdbVertex *v = (EdbVertex *) vtxarr->At(i);
		
		for(j=0;j<v->N();j++){
			EdbVTA    *vta = v->GetVTa(j);
			EdbTrackP *t   = vta->GetTrack();
			EdbSegP   *s   = t->GetSegmentFirst();
			vta->SetZpos( s->Z()>v->Z() ? 1 : 0 );
			vta->SetImp(CalcIP(s,v));
		}
		vset->AddVertex(v);
	}
	
	return PVR;
}

/*

void EdbEDAIO::PrintTrackFeedback(EdbTrackP *t, EdbVertex *v1, EdbVertex *v2, EdbSegP *s1, FILE *fp, bool print_segs){
	// t  : pointer to the track
	// v1 : upstream vertex
	// v2 : downstream vertex
	// s1 : first segment to be used for vertex calculation
	// fp : file pointer to print data.

	if(t==NULL){
		fprintf(fp,     "trkid v1 v2        x        y        z      tx      ty   ip1   ip2 ");
		fprintf(fp, "    p  pmin  pmax mn pa sb dk of ps ");
		fprintf(fp, " n RmaxT RmaxL   rmsT   rmsL pl1 pl2 res\n");
		
		if(print_segs){
			fprintf(fp, "        ipl        x        y        z      tx      ty type irec ngrains\n");
		}
		
		return;
	}
	int ivtx1 = v1 ? v1->ID() : -1;
	int ivtx2 = v2 ? v2->ID() : -1;

	if(s1==NULL) s1=t->GetSegmentFirst();
	EdbSegP *s2 = t->GetSegmentLast();

	float ip1= v1? v1->DistSeg(s1): -1.0;
	float ip2= v2? v1->DistSeg(s2): -1.0;
	
	double p,pmin,pmax;
	CalcP(t, p, pmin, pmax); // pmin = pmin, pmax = pmax

	// Flags.... same as Feedback file.
	int MAN = 0;
	int SCANBACK=0, idscanback=-1;
	int MUON    =0;
	int DARKNESS=0; // no function to define.
	if(gEDA->GetTrackSet("MN")->FindTrack(t)) MAN=1;
	
	EdbTrackP * sbt = gEDA->CheckScanback(t); // check if the track is a scanback track
	if(sbt!=NULL) idscanback = sbt->ID();
	if (idscanback!=-1) SCANBACK = 1;
	if (idscanback==gEDA->IdMuon()&&gEDA->IdMuon()!=-1) MUON =1;

	int OutOfBrick=0; // OutOfBrick: 1= pass through, 2=edge out, 0 others
	int LastPlate=0;  // Last plate: 1= pass though, n=edge out track last plate, 0 others 
	
	int nseg = t->N();
	double rms, rmst, rmsl;
	DTRMSTL1Kink(t, &rms, &rmst, &rmsl);
	
	double RmaxT, RmaxL;
	double Rmax=0.0;
	int iplRmax1, iplRmax2;
	if(nseg==1) {
		RmaxT=RmaxL=rmst=rmsl=0;
		iplRmax1=iplRmax2=s1->Plate();
	}
	else {
		for(int i=0; i<nseg-1;i++){
			EdbSegP *ss1 = t->GetSegment(i);
			EdbSegP *ss2 = t->GetSegment(i+1);
			double dtt,dtl;
			CalcDTTransLongi(ss1, ss2, &dtt, &dtl);
			dtt=fabs(dtt);
			dtl=fabs(dtl);
			double RT =dtt/rmst;
			double RL =dtl/rmsl;
			double R=sqrt(RT*RT+RL*RL);
			if(Rmax<R){
				Rmax  = R;
				RmaxT = RT;
				RmaxL = RL;
				iplRmax1=ss1->Plate();
				iplRmax2=ss2->Plate();
			}
		}
	}
	
	int result=1;
	
	fprintf(fp, "%5d %2d %2d %8.1f %8.1f %8.1f %7.4f %7.4f %5.1f %5.1f ",
		t->ID(), ivtx1, ivtx2, s1->X(), s1->Y(), s1->Z(), s1->TX(), s1->TY(), ip1, ip2);
	fprintf(fp, "%5.1lf %5.1lf %5.1lf %2d %2d %2d %2d %2d %2d ", 
		p,pmin,pmax, MAN, MUON, SCANBACK, DARKNESS, OutOfBrick, LastPlate);
	fprintf(fp, "%2d %5.2f %5.2f %6.4lf %6.4lf %3d %3d %3d\n",
		nseg, RmaxT, RmaxL, rmst, rmsl, iplRmax1, iplRmax2, result);
	
	if(print_segs){
		int type = 0; // 0:BT, 1:MTtop, 2:MTbottom
		int irec = 0;  // 0:TS, 1:PredScan, 2: Manual
		if( gEDA->GetTrackSet("TS")->FindTrack(t)) irec=0;
		if( gEDA->GetTrackSet("SB")->FindTrack(t)) irec=1;
		if( gEDA->GetTrackSet("SF")->FindTrack(t)) irec=1;
		if( gEDA->GetTrackSet("MN")->FindTrack(t)) irec=2;
		
		for(int i=0;i<t->N();i++){
			EdbSegP *s = t->GetSegment(i);
			fprintf(fp,"%11d %8.1f %8.1f %8.1f %7.4f %7.4f %4d %4d %4d\n",
				s->Plate(), s->X(), s->Y(), s->Z(), s->TX(), s->TY(), type, irec, s->W());
		}
	}
}

*/

void EdbEDAIO::WriteFilteredText(const char *filename){
	int i,j,k,l;
	FILE *fp;
	
	if(filename==NULL){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "filtered tracks", "*.txt", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		if(fi->fFilename!=NULL) {
			filename = fi->fFilename;
		}
	}
	if(filename==NULL) return;

	fp=fopen(filename,"wt");
	
	if(fp==NULL) {
		printf("File open error : %s\n", filename);
		return;
	}
	
	printf("Write Filtered Text File : %s\n", filename);

	// Vertices
	
	TObjArray *vertices = (gEDA->GetVertexSet())->GetVertices();

	EdbVertex *v_1ry=NULL;
	double z=1e9;
	for(i=0;i<vertices->GetEntriesFast();i++){
		EdbVertex *v = (EdbVertex *) vertices->At(i);
		if(v->Z()<z) {
			z = v->Z();
			v_1ry=v;
		}
	}

	fprintf(fp,"// vtx  id        X        Y        Z 1ry charm tau\n");
	for( i=0; i<vertices->GetEntriesFast(); i++){
		EdbVertex *v = (EdbVertex *) vertices->At(i);
		v->SetID(i);
		int ISPRIMARY = v==v_1ry ? 1 : 0;
		int ISCHARM = 0;
		int ISTAU = 0;
		fprintf(fp,"vertex %3d %8.1lf %8.1lf %8.1lf %3d %5d %3d\n", 
					v->ID(), v->X(), v->Y(), v->Z(),ISPRIMARY,ISCHARM,ISTAU);
	}
	
	printf("vertices %d\n", vertices->GetEntriesFast());
	
//	ID_Track ID_UpVtx ID_DownVtx X Y Z SX SY IP_UpVtx IP_DownVtx P Pmin Pmax Manual Mu/e Scanback Darknes

	fprintf(fp,"// trk    id flag      tx      ty pl1 pl2 nseg   ip1ry vtx1   ip1 vtx2   ip2      p   pmin       pmax man part sb dark\n");

	// Tracks + Segments
	
	for(i=0;i<gEDA->NTrackSets();i++){
		EdbEDATrackSet *set = gEDA->GetTrackSet(i);
		if(set->GetDraw()==0) continue;
		for(j=0;j<set->N();j++){
			EdbTrackP *t = set->GetTrack(j);
			EdbSegP *sf = t->GetSegmentFirst();
			EdbSegP *sl = t->GetSegmentLast();
			EdbVertex * vtx1st = NULL;
			EdbVertex * vtx2nd = NULL;
			for(k=0;k<vertices->GetEntriesFast();k++){
				EdbVertex *v = (EdbVertex *) vertices->At(k);
				for(l=0;l<v->N();l++){
					EdbTrackP *tv = v->GetTrack(l);
					if(t==tv){
						if(NULL==vtx1st) vtx1st = v;
						else {
							if(vtx1st->Z()>v->Z()) { // vtx1st should 1ry.
								vtx2nd = vtx1st;
								vtx1st = v;
							}
							else vtx2nd = v;
			}	}	}	}

			// Flags.... same as Feedback file.
			int MAN = 0;
			int SCANBACK=0, idscanback=-1;
			int MUON    =0;
			int DARKNESS=0; // no function to define.

			if(-99<=t->Flag()&&t->Flag()<0) MAN=1;
			
			EdbTrackP * sbt = gEDA->CheckScanback(t); // check if the track is a scanback track
			if(sbt!=NULL) idscanback = sbt->ID();
			if (idscanback!=-1) SCANBACK = 1;
			if (idscanback==IdMuon()&&IdMuon()!=-1) MUON =1;

			double p,pmin,pmax;
			CalcP(t, p, pmin, pmax); // pmin = pmin, pmax = pmax
			int ivtx1 = vtx1st ? vtx1st->ID() : -1;
			int ivtx2 = vtx2nd ? vtx2nd->ID() : -1;
			float ip1=-1.0, ip2=-1.0;
			if(vtx1st){
				for(k=0;k<vtx1st->N();k++){
					EdbVTA *vta = vtx1st->GetVTa(k);
					if(t==vta->GetTrack()) ip1 = vta->Imp();
				}
			}
			if(vtx2nd){
				for(k=0;k<vtx2nd->N();k++){
					EdbVTA *vta = vtx2nd->GetVTa(k);
					if(t==vta->GetTrack()) ip2 = vta->Imp();
				}
			}
			double ip_1ry = CalcIP(sf, v_1ry);
			
			fprintf(fp,"track %6d %4d %7.4f %7.4f %3d %3d %4d %7.1lf %4d %5.1f %4d %5.1f ",
				t->ID(), t->Flag(), t->TX(), t->TY(), sf->Plate(), sl->Plate(), t->N(), 
				ip_1ry, ivtx1, ip1, ivtx2, ip2);
			fprintf(fp,"%6.2lf %6.2lf %10.2lf %3d %4d %2d %4d\n", p,pmin,pmax, MAN, MUON, SCANBACK, DARKNESS);
		}
		printf("TrackSet : %s %d tracks\n", set->GetName(), set->N());
	}

	fprintf(fp,"// seg   trkid tflag ipl  segid        x        y        z      tx      ty  w   chi2 flag\n");
	
	for(i=0;i<gEDA->NTrackSets();i++){
		EdbEDATrackSet *set = gEDA->GetTrackSet(i);
		if(set->GetDraw()==0) continue;
		for(j=0;j<set->N();j++){
			EdbTrackP *t = set->GetTrack(j);
			for(k=0;k<t->N();k++){
				EdbSegP *s = t->GetSegment(k);
				fprintf(fp,"segment %6d %5d %3d %6d %8.1f %8.1f %8.1f %7.4f %7.4f %2d %6.4f %4d\n",
					t->ID(), t->Flag(), s->Plate(), s->ID(), s->X(),s->Y(),s->Z(),s->TX(),s->TY(),
					(int)s->W(),s->Chi2(),s->Flag());
			}
		}
	}
	
	
	// Scanning area
	fprintf(fp,"// area ipl xmin xmax ymin ymax z\n");
	EdbEDAAreaSet *aset = gEDA->GetTrackSet("TS")->GetAreaSet();
	for(i=0; i<aset->N(); i++){
		EdbEDAArea *a = aset->GetArea(i);
		fprintf(fp,"area %2d %8.1lf %8.1lf %8.1lf %8.1lf %8.1lf\n",
			a->Plate(), a->Xmin(), a->Xmax(), a->Ymin(), a->Ymax(), a->Z());
	}
	
	
	
	fclose(fp);
}





EdbTrackP * EdbEDASelection::GetTrack( TEveElement *e){
	if(IsTrack(e)) return (EdbTrackP*) e->GetUserData();
	if(IsSegment(e)){
		TEveElement *parent = *(e->BeginParents());
		if(IsTrack(parent)) return (EdbTrackP*) parent->GetUserData();
	}
	return NULL;
}


TObjArray * EdbEDASelection::SetSelected(TObjArray *selected)	{
	// fill eSelected by the current selected tracks.

	if(selected){ // fill by given object array.
		*eSelected = *selected;
		return eSelected;
	}
	ClearSelected();
	
	if(gEve->GetCurrentEvent()==NULL) return NULL; // if current event is not set, return.

	TEveSelection *sel = gEve->GetSelection();
	
	for(TEveElement::List_i it=sel->BeginChildren();it!=sel->EndChildren();it++){
		TEveElement *e = *it;
		if( IsSegment (e) || IsTrack(e) ){
			TObject *o = (TObject *)e->GetUserData();
			if(NULL==eSelected->FindObject(o)) eSelected->Add(o); // no double entry
			EdbTrackP *t = GetTrack(e);
			if(t!=NULL&&NULL==eSelectedTracks->FindObject(t)) eSelectedTracks->Add(t);
		}
		if( IsVertex (e)){
			// selection of vertex is managed in Added().
		}
	}
	return eSelected;
}


void EdbEDASelection::SelectAll(){
	TEveSelection *sel = gEve->GetSelection();
	sel->RemoveElements();
	TEveEventManager *ev = gEve->GetCurrentEvent();
	TEveElement::List_i it = ev->BeginChildren();
	for(int i=0;i<ev->NumChildren();i++, it++){
		TEveElement *e = *it;
		if (!IsTrack(e)) continue;
		TEveElement::List_i itseg = e->BeginChildren();
		for(int j=0;j<e->NumChildren();j++, itseg++){
			TEveElement *s = *itseg;
			if(e->FindChild(e->GetElementName())) continue;
			sel->AddElement(s);
		}
	}
}

void EdbEDASelection::SetSelection ( TObjArray * selected){
	if(selected) *eSelected = *selected;
	//if(GetTrackSet("TS")->N()==0) return;
	
	Disable();
	
	// TEveSelection clear
	TEveSelection *sel = gEve->GetSelection();
	sel->Clear();
	
	// loop for all elements in the current event.
	TEveEventManager *eEvent = gEve->GetCurrentEvent();
	if(eEvent==NULL) return;
	int nselected = eSelected->GetEntriesFast();

	for(TEveElement::List_i it = eEvent->BeginChildren(); it!=eEvent->EndChildren(); it++){
		TEveElement *e = *it;
		
		// only track or segment will be target to be selected.
		if(IsTrack(e)){
			// if this track is in the list of eSelected, add to TEveSelection
			for(int i=0; i<nselected;i++){
				if( eSelected->At(i) == e->GetUserData() ) sel->AddElement(e);
			}

			// loop over segments.
			if(e->NumChildren()==0) continue;
			for(TEveElement::List_i itj = e->BeginChildren(); itj!=e->EndChildren(); itj++){
				// if this track is in the list of eSelected, add to TEveSelection
				TEveElement *ee = *itj;
				for(int i=0;i<nselected;i++){
					if( eSelected->At(i) == ee->GetUserData() ) sel->AddElement(ee);
				}
			}
		}
		// coloring of selected vertex is done in EdbEDAVertexSet::DrawSingleVertex()
		/*
		if(IsVertex(e)){
			if(eVertex==(EdbVertex *) e->GetUserData()){
				TEvePointSet *ps = (TEvePointSet *) el;
				ps->SetMarkerColor(kRed);
			}
		}
		*/
	}
	Enable();
}

void EdbEDASelection::Added(TEveElement* el){
	// function which define the actions by mouse.
	// if the action is not by the mouse, return.
	if(eDisable) return; 
	
	printf("%s\n",el ? el->GetElementName() : "");
	gEve->SetStatusLine(el->GetElementName());
	if(el->GetUserData() ==NULL) return;

	// if a vertex was selected, the corresponding tracks will be selected.
	if(IsVertex(el)) {
		Disable();

		EdbVertex *v = (EdbVertex *) el->GetUserData();
		// to avoid double selection of the tracks
		eSel->RemoveElements();
		// add vertex to the selection, avoiding the loop in Added().
		eSel->AddElement(el);
		
		// add track into the selection.
		for(int i=0;i<v->N();i++){
			EdbTrackP *t = v->GetTrack(i);
			TEveElement * eltrack = GetTrackElement(t);
			if(NULL!=eltrack) {
				eSel->AddElement(eltrack);
				printf(" %s\n", eltrack->GetElementName());
			}
		}
		Enable();
		
		// change the color of the vertex point.
		if( eSelectedVertex!= v) {
			TEvePointSet *ps = GetVertexElement(eSelectedVertex);
			if(ps!=NULL) ps->SetMarkerColor(kOrange);
			eSelectedVertex=v;
			ps = (TEvePointSet *) el;
			ps->SetMarkerColor(kRed);
			gEve->Redraw3D();
		}

		// write feedback file
		gEDA->WriteFeedbackFile2008(v,"tmp.old.feedback");
	}

	// revise TObjArray * selected, *selected_tracks
	SetSelected();
	
	// revice extention
	gEDA->GetExtentionSet()->Draw(kFALSE);
	gEve->Redraw3D();
}

void EdbEDASelection::Removed(TEveElement* el){
	printf("Cancel %s\n", el ? el->GetElementName() : "");
	SetSelected();
}

void EdbEDASelection::Cleared(){
	SetSelected();
	if(eDisable) return; 
	printf("Cleared\n");
}
/*
EdbDataSet * EdbEDA::PVRec2DataSet(EdbPVRec *pvr){
	EdbDataSet *dset = new EdbDataSet();
	EdbDataPiece *piece=0;
	EdbDataPiece *pp=0;
	dset->Set0();
	
	for(int i=0; i<pvr->Npatterns(); i++){
		EdbPattern *pat = pvr->GetPattern(i);
		int plateID = pat->PID();
		int pieceID = 0;
		char filename[]="./*";
		int flag = 1;
		piece = new EdbDataPiece(plateID,pieceID,filename,flag);
		piece->MakeName();
		piece->GetMakeLayer(0)->SetZlayer(pat->Z(),0,0);
	}
	return dset;
}
*/
EdbDataProc * EdbEDA::PVRec2DataProc(EdbPVRec *pvr){
	FILE *fp;
	printf("write lnkdummy.def\n");
	fp=fopen("lnkdummy.def","wt");
	fprintf(fp, "INPUT_RUNS_LIST		lnkdummy.lst\n"
				"OUTPUT_DATA_DIR		./data/\n"
				"PARAMETERS_DIR		./pardummy/\n");
	fclose(fp);

	fp = fopen("lnkdummy.lst","wt");
	for(int i=0;i<pvr->Npatterns();i++){
		EdbPattern *pat = pvr->GetPattern(i);
		if(pat==NULL) continue;
		fprintf(fp, "%3d 1 ./data/* 1\n", pat->PID());
	}
	fclose(fp);
	
	printf("make pardummy directry\n");
	system("mkdir pardummy");
	printf("write parameters files\n");
	for(int i=0;i<pvr->Npatterns();i++){
		EdbPattern *pat = pvr->GetPattern(i);
		if(pat==NULL) continue;
		
		fp = fopen(Form("pardummy/%02d_001.par", pat->PID()),"wt");
		fprintf(fp,"ZLAYER 0 %.1lf 0. 0.\n", pat->Z());
		fprintf(fp,"AFFXY 0 1 0 0 1 0 0\n");
		fclose(fp);
	}
	
	EdbDataProc *proc = new EdbDataProc("lnkdummy.def");
	proc->SetPVR(pvr);
	
	system("rm lnkdummy.def");
	system("rm lnkdummy.lst");
	system("rm pardummy/*.*");
	system("rmdir pardummy");

	
	return proc;
}


#include<TGLUtil.h>
#include<TGLIncludes.h>

EdbEDAOverlay::EdbEDAOverlay():eDraw(1), eXa(80), eYa(30), eLength(120), eXb(50), eYb(20){}
void EdbEDAOverlay::SetXY(int pixel_x, int pixel_y){ if(gEDA->GetScaleZ()==1.0) eXb=pixel_x, eYb=pixel_y; else eXa=pixel_x, eYa=pixel_y;}
void EdbEDAOverlay::SetLength(int pixel_length){if(gEDA->GetScaleZ()==1.0) eLength=pixel_length; else eLength=pixel_length;}

void EdbEDAOverlay::Render(TGLRnrCtx& ctx){
	if(eDraw==kFALSE) return;
	if(gEDA->GetScaleZ()==1.0) RenderBar(ctx);
	else RenderAxes(ctx);
}

void EdbEDAOverlay::RenderAxes(TGLRnrCtx& ctx){
	// Render 3 scale axes.

	// view port to world.
	TGLCamera& cam = ctx.RefCamera();
	TGLLine3 l1 = cam.ViewportToWorld(eXa,eYa);
	TGLLine3 l2 = cam.ViewportToWorld(eXa+eLength,eYa);
	
	// find distance between 2 lines.
	const TGLVector3& glv1 = l1.Vector();
	const TGLVector3& glv12 = l2.Start()-l1.Start();
	TVector3 v1(glv1.X(), glv1.Y(), glv1.Z());
	v1.SetMag(1.);
	TVector3 v12(glv12.X(), glv12.Y(), glv12.Z());
	TVector3 vx = v1.Cross(v12);
	double d1 = vx.Mag();
	
	// find 
	double d;
	for(int i=0;i<10;i++){
		if( d1 < 2*TMath::Power(10,i) ) { d=1*TMath::Power(10,i); break; }
		if( d1 < 5*TMath::Power(10,i) ) { d=2*TMath::Power(10,i); break; }
		if( d1 < 10*TMath::Power(10,i) ) { d=5*TMath::Power(10,i); break; }
		d=d1;
	}
	//printf("%lf %lf\r", d, d1);
	
	
	// Draw Text
	TGLLine3 l3 = cam.ViewportToWorld(eXa+10,eYa+20);
	const TGLVertex3& glp31=l3.Start();
	const TGLVertex3& glp32=l3.End();
	TVector3 vec31(glp31.CArr());
	TVector3 vec32(glp32.CArr());
	TVector3 p31 = (vec31+vec32)*0.5;

	TString num=Form(" %d um",(int)d);
#if ROOT_VERSION_CODE>=ROOT_VERSION(5,26,00)
	TGLUtil::Color(kWhite);
	TGLFont font;
	ctx.RegisterFont(12, "arial", TGLFont::kBitmap, font);
	font.PreRender();
	font.Render(num, p31.X(), p31.Y(), p31.Z(), TGLFont::kLeft, TGLFont::kTop);
	font.PostRender();

	TGLUtil::LineWidth(1.0);
#else
	gEve->GetBrowser()->SetStatusText(num,0);
#endif

	// Light off
	TGLCapabilitySwitch lights_off(GL_LIGHTING, kFALSE); // need TGLIncludes.h
	
	// Draw Line
	const TGLVertex3& glp1=l1.Start();
	const TGLVertex3& glp2=l1.End();
	TVector3 vec1(glp1.CArr());
	TVector3 vec2(glp2.CArr());
	TVector3 p1 = (vec1+vec2)*0.5;
	TVector3 p2;

	TGLUtil::Color3f(1.,0.,0.);
	TGLVertex3 glpos(p1.X(), p1.Y(), p1.Z());
	TGLVertex3 glpos2(p1.X()+d, p1.Y(), p1.Z());
	TGLLine3 lx(glpos, glpos2);
	lx.Draw();

	TGLUtil::Color3f(0.,1.,0.);
	glpos2.Set(p1.X(), p1.Y()+d, p1.Z());
	TGLLine3 ly(glpos, glpos2);
	ly.Draw();

	TGLUtil::Color(kCyan);
	glpos2.Set(p1.X(), p1.Y(), p1.Z()+d*gEDA->GetScaleZ());
	TGLLine3 lz(glpos, glpos2);
	lz.Draw();
	
	
}

void EdbEDAOverlay::RenderBar(TGLRnrCtx& ctx){
	// Render a scale bar.

	// view port to world.
	TGLCamera& cam = ctx.RefCamera();
	TGLLine3 l1 = cam.ViewportToWorld(eXb,eYb);
	TGLLine3 l2 = cam.ViewportToWorld(eXb+eLength,eYb);
	
	// find distance between 2 lines.
	const TGLVector3& glv1 = l1.Vector();
	const TGLVector3& glv12 = l2.Start()-l1.Start();
	TVector3 v1(glv1.X(), glv1.Y(), glv1.Z());
	v1.SetMag(1.);
	TVector3 v12(glv12.X(), glv12.Y(), glv12.Z());
	TVector3 vx = v1.Cross(v12);
	double d1 = vx.Mag();
	
	// find 
	double d;
	for(int i=0;i<10;i++){
		if( d1 < 2*TMath::Power(10,i) ) { d=1*TMath::Power(10,i); break; }
		if( d1 < 5*TMath::Power(10,i) ) { d=2*TMath::Power(10,i); break; }
		if( d1 < 10*TMath::Power(10,i) ) { d=5*TMath::Power(10,i); break; }
		d=d1;
	}
	//printf("%lf %lf\r", d, d1);
	
	// position
	const TGLVertex3& p1=l1.Start();
	const TGLVertex3& p2=l2.End();
	TVector3 pp1(p1.X(), p1.Y(), p1.Z());
	TVector3 pp2(p2.X(), p2.Y(), p2.Z());
	TVector3 p=pp1*0.9+pp2*0.1-vx*0.05;
	TString num=Form("%d",(int)d);

	// Light off
	TGLCapabilitySwitch lights_off(GL_LIGHTING, kFALSE); // need TGLIncludes.h

#if ROOT_VERSION_CODE>=ROOT_VERSION(5,26,00)
	TGLUtil::Color(kWhite);
	TGLFont font;
	ctx.RegisterFont(12, "arial", TGLFont::kBitmap, font);
	font.Render(num, p.X(), p.Y(), p.Z(), TGLFont::kLeft, TGLFont::kTop);
	TGLUtil::LineWidth(2.0);
#else
	gEve->GetBrowser()->SetStatusText(num,0);
#endif

	
	// Draw Line
	const Float_t col[4]={1.0,0.0,0.0,0.0};
	TGLUtil::Color4fv(col);
	TGLVector3 vec12=(l2.End()-l1.Start())*(d/d1);
	TGLLine3 l(l1.Start(), vec12);
	l.Draw();
	
	// Draw Text
	
	
	
	// doesn't work.
	//TGLText tx;
	//tx.SetGLTextFont(41);
	//tx.SetTextColor(kWhite);
	//tx.SetTextSize(12);
	//tx.PaintGLText(p.X(), p.Y(), p.Z(),Form("%d", (int) d));

	// doesn't work.
	//	TGLUtil::DrawNumber(num, p); 

	// doesn't work.
	//	TGLFont f;
	//	ctx.RegisterFont(12, "arial", TGLFont::kBitmap, f);
	//	f.Render(num, p.X(), p.Y(), p.Z(), 0);

	// works but need OpenGL32.lib
	/*	const char *txt=Form("%d",(int)d);
		Float_t bb[6];
		TGLFont font;
		ctx.RegisterFont(12, "arial", TGLFont::kBitmap, font);
		font.BBox(txt, bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
		glPushMatrix();
		glTranslated(p.X(), p.Y(), p.Z());
		glRasterPos2i(0,0);
		font.Render(txt, );
		glPopMatrix();
	*/	

}


void EdbEDA::Init(void){
	// Initialization of TrackSets and VertexSet.
	gEDA = this;
	eScaleZ=1.0;
	eJapanese=kFALSE;
	
	// Setting of Path, Brick number.
	eScanProc = new EdbScanProc;
	eScanProc->eProcDirClient=FindProcDirClient();
	printf("%s\n", eScanProc->eProcDirClient.Data());
	int ibrick = FindBrickIDFromPath();
	EdbEDABrickData::SetBrick ( ibrick);

	//eScanSet  = new EdbScanSet();
	
	eDrawObjects      = new TObjArray;
	eTrackSets        = new TObjArray;
	
	EdbEDATrackSet *set;
	eTrackSets->Add(set = new EdbEDATrackSet("TS"));
	set->SetDrawAreaSet(kTRUE);
	
	eTrackSets->Add(set = new EdbEDATrackSet("SB"));
	set->SetTrackAttribute(kWhite);
	set->SetID(EdbID(ibrick,0,1,100));
	
	eTrackSets->Add(set = new EdbEDATrackSet("SF"));
	set->SetTrackAttribute(kGreen);
	set->SetID(EdbID(ibrick,0,1,200));
	
	eTrackSets->Add(set = new EdbEDATrackSet("BT"));
	set->SetInterpolation(kFALSE);   // disable to fill dashed line on holes.
	set->SetExtendMode(kExtendDown); // force extention down.

	eTrackSets->Add(set = new EdbEDATrackSet("MN"));
	set->SetTrackAttribute(kViolet);
	eVertexSet = new EdbEDAVertexSet();
	eExt = new EdbEDAExtentionSet();
	

	//AddDrawObject(new EdbEDAObject);
	
}


void EdbEDA::Run(){
	// Start Display.
	
	if(gEve==NULL){
		// TEve setting
#if ROOT_VERSION_CODE<ROOT_VERSION(5,22,00)
		TEveManager::Create();
		eViewer = gEve->GetGLViewer();

#else
		TEveManager::Create(kTRUE,"V");
		eViewer = gEve->GetDefaultGLViewer();
		gEve->GetWindowManager()->HideAllEveDecorations();
		
		gEve->GetDefaultViewer()->SetName("Viewer1");
		TGTab *tab=gEve->GetBrowser()->GetTabRight();
		tab->GetCurrentTab()->SetName("Viewer1");
		tab->Connect("Selected(Int_t)", "EdbEDA", this, "GetGLViewer()");

#endif 
		
		gEve->GetMainWindow()->Resize(1024,730); // resize for note book size 1280x800

		eViewer->SetCurrentCamera(TGLViewer::kCameraOrthoZOY); // Orthogonal Camera
		TGLOrthoCamera& cam = (TGLOrthoCamera &) eViewer->CurrentCamera();
		cam.SetEnableRotate(kTRUE);
		
		eViewer->SetClearColor(kBlack);
		SetOverlay();
		gEve->GetMainWindow()->SetWindowName(Form("OPERA Event Display and Analyser. Brick %d", EdbEDABrickData::GetBrick()));
		//gEve->GetMainWindow()->SetWindowName("OPERA Event Display and Analyser.");
		
		InitSelection(gEve->GetSelection());

		new EdbEDAListTab;
		eMainTab = new EdbEDAMainTab;

		eTrackSetTab = new EdbEDATrackSetTab;
		ePlotTab     = new EdbEDAPlotTab;
		eVertexTab   = new EdbEDAVertexTab;
		eDSTab       = new EdbEDADecaySearchTab;
		eOperationTab= new EdbEDAOperationTab;
		new EdbEDAPredTab;
		new EdbEDAShowerTab;
		new EdbEDAMCTab;
		// new EdbEDAShowerTab;
		//	gEve->GetBrowser()->GetTabBottom()->RemoveTab(0); // remove command tab
		gEve->GetBrowser()->GetTabBottom()->SetTab(0); // set TrackSetTab in front.
	}
	
	printf("GUI setting finished.\n");
	//if(eAreaSet->N()==0&&GetTrackSet("TS")->NBase()!=0) eAreaSet->SetAreas( GetTrackSet("TS")->GetTracksBase());
	Draw(kFALSE);
	if(gEve->GetCurrentEvent()!=NULL) gEve->Redraw3D(kTRUE);
}

void EdbEDA::Draw(int redraw){
	
	if(gEve->GetCurrentEvent()!=0) gEve->GetCurrentEvent()->RemoveElements();
	
	for(int i=0;i<NTrackSets();i++) GetTrackSet(i)->Draw(kFALSE);
	
	eVertexSet->Draw();
	
	for(int i=0;i<NDrawObjects();i++) if(GetDrawObject(i)) GetDrawObject(i)->Draw();
	
	SetSelection();
	
	SetSelected();
	eExt->Draw();
	
	if(redraw) gEve->Redraw3D();
	
	eMainTab->SetDrawCheckBox();
}

void EdbEDA::OpenFile(const char *filename, int datatype, TCut rcut){
	//  filename : filename for "LinkDef" or "*.set.root(ScanSet)" or "*.root(linked_track.root format)"
	//  datatype : only for "LinkDef". Data type for EdbDataProc::InitVolume(). put -1 for no InitVolume
	//            if datatype= 100 or 1000. linked tracks will be registred to "TS".
	//  rcut     : Track selection when read.
	
	// if filename is not given. open file browser.
	if(filename==NULL) {
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "LinkDef", "*.def",
		                            "Linkeda Tracks", "*.root",
		                            "ScanSet", "*.set.root",
		                            "Feedback file","*.feedback*",
		                            "Mxx file","*.all", 
		                            "All files","*",0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), 0, kFDOpen, fi);
		filename = fi->fFilename;
	}
	
	if(filename==NULL) return;
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
		printf("EDA Start with Scan Set : %s\n", s.Data());
		
		Init();
		TFile *f = TFile::Open(filename);
		eScanSet = (EdbScanSet *)f->Get("set");
		eID = *((EdbID *)eScanSet->eIDS.At(0));
		eID.ePlate = 0;
		
		TString s2;
		
		eScanProc->MakeFileName(s2, eID, "trk.root", kFALSE);
		int IsExist = gSystem->AccessPathName(s2.Data(), kReadPermission);
		if( IsExist==kFALSE) {
			GetTrackSet("TS")->SetScanSet(eScanSet);
			GetTrackSet("TS")->ReadTracksTree(eID,rcut);
			GetTrackSet("TS")->SetNsegCut(3);
			GetTrackSet("TS")->DoSelection();
		}
		else{
			GetTrackSet("SB")->ReadPredictionScan(eScanSet);
		}
	}
	else if(s.EndsWith(".def")){
		// Link def file
		printf("EDA Start with LinkDef : %s\n", s.Data());
		Init();
		EdbDataProc *dproc = new EdbDataProc(filename);
		eDataSet  = dproc->GetDataSet();
		
		if(datatype>=0) {
			dproc->InitVolume(datatype,rcut);
			SetDataType(datatype);
			GetTrackSet("TS")->ReadLinkedTracks(dproc);
			GetTrackSet("TS")->SetNsegCut(3);
			GetTrackSet("TS")->SetUpstreamPlate(GetPlateUp());
			GetTrackSet("TS")->DoSelection();
		}
	}
	else if(s.EndsWith(".root")){
		// Linked tracks file
		Init();
		printf("EDA Start with Linked_track : %s\n", s.Data());
		GetTrackSet("TS")->ReadLinkedTracksFile(filename,rcut);
		ePVR=GetTrackSet("TS")->GetPVRec();
		GetTrackSet("TS")->SetNsegCut(3);
		GetTrackSet("TS")->DoSelection();
	}
	else if(s.EndsWith(".feedback")){
		// Feedback
		Init();
		printf("EDA Start with Feedback : %s\n", s.Data());
		printf("Only for viewer!\n");
		ReadFeedbackFile(filename);
	}
	else if(s.EndsWith(".all")){
		// Mxx file (Nagoya style ascii)
		Init();
		SetJapanese();
		printf("EDA Start with Mxx file : %s\n", filename);
		EdbPVRec *pvr = EdbEDAUtil::ReadMxxPVR(filename);
		
		ePVR = pvr;
		GetTrackSet("TS")->AddTracksPVR(pvr);
		GetVertexSet()->AddVertices(pvr->eVTX);
	}
	
	
//		Run();
}


void EdbEDAText::Draw(Option_t *opt){
	if(GetDraw()==kFALSE) return;
	TEveText *tx = new TEveText();
	tx->SetElementName(eText.Data());
	tx->SetText(eText.Data());
	tx->PtrMainTrans()->SetPos(eX, eY, eZ*gEDA->GetScaleZ());
	
	tx->SetFontMode(TGLFont::kBitmap);
	tx->SetFontSize(eTextSize);
	tx->SetFontFile(eTextFont);
	tx->SetMainColor(eTextColor);
	gEve->AddElement(tx);
}

void EdbEDALine::Draw(Option_t *option){
	if(GetDraw()==kFALSE) return;
	TEveLine *l = new TEveLine;
	l->SetElementName(text.Data());
	l->SetNextPoint(X1,Y1,Z1*gEDA->GetScaleZ());
	l->SetNextPoint(X2,Y2,Z2*gEDA->GetScaleZ());
	l->SetLineColor(col);
	l->SetLineWidth(wid);
	gEve->AddElement(l);
}

