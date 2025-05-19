#include"EdbEDAUtil.h"


int EdbEDAUtil::IsVertex(TObject *o){
	if(o==NULL) return 0;
	if(strcmp(o->ClassName(),"EdbVertex")==0) return 1;
	return 0;
}


int EdbEDAUtil::IsVertex(TEveElement *e){
	if(e==NULL) return 0;
	TObject *o = (TObject *) e->GetUserData();
	return IsVertex(o);
}

int EdbEDAUtil::IsTrack(TObject *o){
	if(o==NULL) return 0;
	if(strcmp(o->ClassName(),"EdbTrackP")==0) return 1;
	return 0;
}

int EdbEDAUtil::IsTrack(TEveElement *e){
	if(e==NULL) return 0;
	TObject *o = (TObject *) e->GetUserData();
	return IsTrack(o);
}

int EdbEDAUtil::IsSegment(TObject *o){
	if(o==NULL) return 0;
	if(strcmp(o->ClassName(),"EdbSegP")==0) return 1;
	return 0;
}

int EdbEDAUtil::IsSegment(TEveElement *e){
	if(e==NULL) return 0;
	TObject *o = (TObject *) e->GetUserData();
	return IsSegment(o);
}

int EdbEDAUtil::IsIncludeCouples(EdbPattern *pat){
	// if pattern includes BT, return 1.
	// check all segments whether the seg->Track()==-1 or not.
	
	int flag=0;
	for(int j=0;j<pat->N();j++){
		EdbSegP *s = pat->GetSegment(j);
		if(s->Track()==-1){ // if only tracks, s->Track() should have track id.
			flag=1;
			break;
		}
	}
	return flag;
}

TEvePointSet * EdbEDAUtil::GetVertexElement ( EdbVertex * v){
	// Return the correspoinding TEvePointSet element to the EdbVertex.
	if (v==NULL) return NULL;
	TEveElement::List_i it = gEve->GetCurrentEvent()->BeginChildren();
	int n= gEve->GetCurrentEvent()->NumChildren();
	for(int i=0;i<n;i++, it++){
		TEvePointSet *e = (TEvePointSet *) *it;
		EdbVertex *vv = (EdbVertex *) e->GetUserData();
		if(v==vv) return e;
	}
	return NULL;
}

TEveCompound * EdbEDAUtil::GetTrackElement ( EdbTrackP * t){
	if(t==NULL) return NULL;
	TEveEventManager *ev = gEve->GetCurrentEvent();
	TEveElement::List_i it = ev->BeginChildren();
	int n= ev->NumChildren();
	for(int i=0; i<n; i++, it++){
		TEveElement *e = *it;
		EdbTrackP *tt = (EdbTrackP *) e->GetUserData();
		if(t==tt) return (TEveCompound *)e;
	}
	return NULL;
}




double EdbEDAUtil::CalcIP(EdbSegP *s, double x, double y, double z){
	// Calculate IP between a given segment and a given x,y,z.
	// return the IP value.
	
	double ax = s->TX();
	double ay = s->TY();
	double bx = s->X()-ax*s->Z();
	double by = s->Y()-ay*s->Z();

	double a;
	double r;
	double xx,yy,zz;
	
	a = (ax*(x-bx)+ay*(y-by)+1.*(z-0.))/(ax*ax+ay*ay+1.);
	xx = bx +ax*a;
	yy = by +ay*a;
	zz = 0. +1.*a;
	r = sqrt((xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z));
	
	return r;
}	

double EdbEDAUtil::CalcIP(EdbSegP *s, EdbVertex *v){
	// calculate IP between a given segment and a given vertex.
	// return the IP value.
	// this is used for IP cut.
	
	// if vertex is not given, use the selected vertex.
	// if(v==NULL) v=gEDA->GetSelectedVertex();
	
	if(v==NULL) return -1.;
	return CalcIP(s, v->X(), v->Y(), v->Z());
}


double EdbEDAUtil::CalcMinimumKinkAngle(TVector3 vertex, TVector3 daupos, TVector3 daumom){
	TVector3 parent = daupos-vertex;

	return daumom.Angle(parent);
}

double EdbEDAUtil::CalcMinimumKinkAngle( EdbVertex *v1ry, EdbSegP *tdaughter, int z_is_middle_of_base){
	// Minimum kink angle calculation in case that 1ry vertex is define && 1 daughter.
	
	EdbSegP td = *tdaughter;
	
	td.PrintNice();
	printf("%lf %lf %lf\n",v1ry->X(), v1ry->Y(), v1ry->Z());
	
	if(z_is_middle_of_base) td.PropagateTo( td.Z()-150);
	
	TVector3 vertex( v1ry->X(), v1ry->Y(), v1ry->Z());
	TVector3 daupos( td.X(), td.Y(), td.Z());
	TVector3 daumom( td.TX(), td.TY(), 1.0);
	
	return CalcMinimumKinkAngle( vertex, daupos, daumom);
}

double EdbEDAUtil::CalcMinimumKinkAngle( EdbSegP *t1ry, EdbSegP *tdaughter, int z_is_middle_of_base){
	// calculate minimum kink angle 
	// for the short decay with only 1 track from primary and 1 track from secondary vertex.
	// move virtual vertex in the lead, find minimimum kink angle.
	// if z_is_middle_of_base==1, move z position of tracks 150 microns upstream.
	
	double min = 1000000000000;
	EdbSegP t1=*t1ry;
	EdbSegP td=*tdaughter;
	
	if(z_is_middle_of_base){
		t1.PropagateTo( t1.Z()-150);
		td.PropagateTo( td.Z()-150);
	}
	
	float z1 = t1.Z();
	for(int iz=0; iz<1000; iz++) { // move virtual vertex point
		t1.PropagateTo(z1-iz);
		TVector3 vertex( t1.X(), t1.Y(), t1.Z());
		TVector3 daupos( td.X(), td.Y(), td.Z());
		TVector3 daumom( td.TX(), td.TY(), 1);
		double theta = CalcMinimumKinkAngle( vertex, daupos, daumom);
		if(min>theta) min=theta;
	}
	
	return min;
}

double EdbEDAUtil::CalcPtmin( TVector3 vertex, TVector3 daupos, TVector3 daumom){
	return daumom.Mag() * sin( CalcMinimumKinkAngle( vertex, daupos, daumom) );
}

double EdbEDAUtil::CalcPtmin( EdbVertex *v1ry, EdbSegP *tdaughter, int z_is_middle_of_base){
	// Ptmin calculation in case that 1ry vertex is define && 1 daughter.
	return tdaughter->P()*sin(CalcMinimumKinkAngle(v1ry, tdaughter, z_is_middle_of_base));
}

double EdbEDAUtil::CalcPtmin( EdbSegP *t1ry, EdbSegP *tdaughter, int z_is_middle_of_base){
	// Ptmin calculation in case that only 1ry track at 1ry vertex && 1 daughter.
	// move virtual vertex to find minimum kink angle. see also CalcMinimumKinkAngle()
	return tdaughter->P()*sin(CalcMinimumKinkAngle(t1ry, tdaughter, z_is_middle_of_base));
}

double EdbEDAUtil::CalcKinkAngle(EdbSegP *tparent, EdbSegP *tdaughter){
	// simple calcuration
	TVector3 parent(tparent->TX(), tparent->TY(), 1.0);
	TVector3 daughter(tdaughter->TX(), tdaughter->TY(), 1.0);
	
	return parent.Angle(daughter);
}

double EdbEDAUtil::CalcPt(EdbSegP *tparent, EdbSegP *tdaughter){
	
	return tdaughter->P()*sin( CalcKinkAngle(tparent, tdaughter));
	
}


double EdbEDAUtil::CalcEMCSelectron(EdbTrackP *t){
	
	if(t==NULL) return -99.0;
	
	int nseg = t->N();
	
	if(nseg>8) nseg=8;
	
	// trk 1
	double dtxsq=0., dtysq=0.;
	EdbSegP *s0 = t->GetSegmentFirst();
	
	for(int i=0;i<nseg-1;i++){
		EdbSegP *s1 = (EdbSegP *) t->GetSegment(i);
		EdbSegP *s2 = (EdbSegP *) t->GetSegment(i+1);
		
		double dz   = fabs( s2->Z()-s0->Z() ) / 1.3; // lead thickness
		double correction_factor = exp( -dz/5612 );
		double  dtx = (s1->TX()-s2->TX())*correction_factor;
		double  dty = (s1->TY()-s2->TY())*correction_factor;
		dtxsq += dtx*dtx;
		dtysq += dty*dty;
	}
	if(nseg>4){
		
		double dtrms = dtxsq/nseg + dtysq/nseg;
		
		if(dtrms<0.0) { return 99.;}
		else {
			dtrms = sqrt(dtrms/2.);
			return 0.0136/dtrms*sqrt(1/5.6)*(1+0.038*TMath::Log(1/5.6));
		}
	}
	return -99;
}



double EdbEDAUtil::CalcDmin(EdbSegP *seg1, EdbSegP *seg2, double *dminz){
// calculate minimum distance of 2 lines.
// use the data of (the selected object)->X(), Y(), Z(), TX(), TY().
// means, if the selected object == segment, use the data of the segment. or it == track, the use the fitted data.
// original code from Tomoko Ariga
// return deltaZ between Z(dmin) and Z(seg1). ==> if vertex is upstream, return plus value.
	double x1,y1,z1,ax1,ay1;
	double x2,y2,z2,ax2,ay2;
	double s1,s2,s1bunsi,s1bunbo,s2bunsi,s2bunbo;
	double p1x,p1y,p1z,p2x,p2y,p2z,p1p2;

	if(seg1->ID()==seg2->ID()&&seg1->PID()==seg2->PID()) return 0.0;
	
	x1 = seg1->X();
	y1 = seg1->Y();
	z1 = seg1->Z();
	ax1= seg1->TX();
	ay1= seg1->TY();

	x2 = seg2->X();
	y2 = seg2->Y();
	z2 = seg2->Z();
	ax2= seg2->TX();
	ay2= seg2->TY();

	s1bunsi=(ax2*ax2+ay2*ay2+1)*(ax1*(x2-x1)+ay1*(y2-y1)+z2-z1) - (ax1*ax2+ay1*ay2+1)*(ax2*(x2-x1)+ay2*(y2-y1)+z2-z1);
	s1bunbo=(ax1*ax1+ay1*ay1+1)*(ax2*ax2+ay2*ay2+1) - (ax1*ax2+ay1*ay2+1)*(ax1*ax2+ay1*ay2+1);
	s2bunsi=(ax1*ax2+ay1*ay2+1)*(ax1*(x2-x1)+ay1*(y2-y1)+z2-z1) - (ax1*ax1+ay1*ay1+1)*(ax2*(x2-x1)+ay2*(y2-y1)+z2-z1);
	s2bunbo=(ax1*ax1+ay1*ay1+1)*(ax2*ax2+ay2*ay2+1) - (ax1*ax2+ay1*ay2+1)*(ax1*ax2+ay1*ay2+1);
	s1=s1bunsi/s1bunbo;
	s2=s2bunsi/s2bunbo;
	p1x=x1+s1*ax1;
	p1y=y1+s1*ay1;
	p1z=z1+s1*1;
	p2x=x2+s2*ax2;
	p2y=y2+s2*ay2;
	p2z=z2+s2*1;
	p1p2=sqrt( (p1x-p2x)*(p1x-p2x)+(p1y-p2y)*(p1y-p2y)+(p1z-p2z)*(p1z-p2z) );

	if(dminz!=NULL) *dminz = z1-p1z;
	return p1p2;
}

EdbVertex * EdbEDAUtil::CalcVertex(TObjArray *segments){
	// calc vertex point with given segments. just topological calculation.
	// VTA is currently not set.
	// in case of Single-stop. make a vertex 650 micron upstream ob base.
	
	double xx,yy,zz,Det,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz;
	Ax=0.;Ay=0.;Az=0.;Bx=0.;By=0.;Bz=0.;Cx=0.;Cy=0.;Cz=0.;Dx=0.;Dy=0.;Dz=0.;

	if(segments->GetEntriesFast()==1){
		// in case of Single-stop. make a vertex 650 micron upstream ob base.
		EdbSegP *s = new EdbSegP(*((EdbSegP *) segments->At(0)));
		s->PropagateTo(s->Z()-650);
		xx=s->X();
		yy=s->Y();
		zz=s->Z();
		delete s;
	}
	else {
		for(int i=0;i<segments->GetEntriesFast();i++){
			EdbSegP *s = (EdbSegP *) segments->At(i);
			double ax = s->TX();
			double ay = s->TY();
			double az = 1.0;
			double x  = s->X();
			double y  = s->Y();
			double z  = s->Z();
			double a = ax*ax+ay*ay+az*az;
			double c = -ax*x-ay*y-az*z;
			double b = (ax*ax+ay*ay);
	//		double w = 1.0/a/a; // weight for small angle tracks.
			double w = 1.0; // no weight
			
			Ax+=2.0*w/b*( a*(ay*ay+az*az) );
			Bx+=2.0*w/b*( -a*ax*ay );
			Cx+=2.0*w/b*( -a*ax*az );
			Dx+=2.0*w/b*( -(a*x+c*ax)*(ax*ax-a)-(a*y+c*ay)*ax*ay-(a*z+c*az)*az*ax );

			Ay+=2.0*w/b*( -a*ay*ax );
			By+=2.0*w/b*( a*(az*az+ax*ax) );
			Cy+=2.0*w/b*( -a*ay*az );
			Dy+=2.0*w/b*( -(a*y+c*ay)*(ay*ay-a)-(a*z+c*az)*ay*az-(a*x+c*ax)*ax*ay );

			Az+=2.0*w/b*( -a*az*ax );
			Bz+=2.0*w/b*( -a*az*ay );
			Cz+=2.0*w/b*( a*b );
			Dz+=2.0*w/b*( -(a*z+c*az)*(az*az-a)-(a*x+c*ax)*az*ax-(a*y+c*ay)*ay*az );

		}

		Det=fabs( Ax*(By*Cz-Cy*Bz)-Bx*(Ay*Cz-Cy*Az)+Cx*(Ay*Bz-By*Az) );
		xx=( (By*Cz-Cy*Bz)*Dx-(Bx*Cz-Cx*Bz)*Dy+(Bx*Cy-Cx*By)*Dz)/Det;
		yy=(-(Ay*Cz-Cy*Az)*Dx+(Ax*Cz-Cx*Az)*Dy-(Ax*Cy-Cx*Ay)*Dz)/Det;
		zz=( (Ay*Bz-By*Az)*Dx-(Ax*Bz-Bx*Az)*Dy+(Ax*By-Bx*Ay)*Dz)/Det;
	}
	
	EdbVertex *v = new EdbVertex();
	v->SetXYZ(xx,yy,zz);
	
	return v;
}


void EdbEDAUtil::CalcPPartial(EdbTrackP *t0, EdbSegP *s1st, EdbSegP *slast, double& p, double& pmin, double& pmax, bool print){
	// momentum calculation partially using from s2 to slast.
	
	// copy the track. and create a new track.
	EdbTrackP *t = new EdbTrackP;
	t->Copy(*t0); 
	
	int i;
	// remove the segments before s1st.
	for(i=0;s1st->PID()!=t->GetSegmentFirst()->PID();i++) t->RemoveSegment( t->GetSegmentFirst());

	// remain the segments up to slast.
	for(i=0;slast->PID()!=t->GetSegment(i)->PID();i++) {}

	// remove the segment after the slast.
	i++;
	for(;t->GetSegment(i);) t->RemoveSegment( t->GetSegment(i));
	
	// calculate momentum of the new track.
	CalcP(t, p,pmin,pmax, print);
	
	delete t;
}


EdbMomentumEstimator *EdbEDAUtil::CalcP(EdbTrackP *t0, double& p, double& pmin, double& pmax, bool print){
	p = pmin = pmax = -99;
  if(t0==0) return NULL;
  EdbTrackP *t=CleanTrack(t0);  // remove microtrack and strange segment.
  t->SetCounters();
  
  	if(t->Npl()<t->N()){
		// in this case P estimation doesn't work.
		// this can happen if a track is formed from 2 different data set.
		// (PID definitions are defferent amond 2 data set. Npl is computed wrongly.)
		// put absolute plate number as PID. (temporary solution)
		for(int i=0;i<t->N();i++) t->GetSegment(i)->SetPID(t->GetSegment(i)->Plate());
	}
	
  if(t->N()<=2) return NULL;

  if(print) printf("CalcP : itrk %i ( %1.4f, %1.4f ) nBT= %i\n",t->ID(),t->TX(),t->TY(),t->N());

  EdbMomentumEstimator *tf= new EdbMomentumEstimator();
  tf->eMinEntr=3; 
  tf->PMS(*t);      

  if (sqrt(t->TX()*t->TX()+t->TY()*t->TY())<=0.2)  // use P3D
    {
      p=tf->eP;
      pmin=tf->ePmin;
      pmax=tf->ePmax;
	  if(print) printf("  ->P3D=%7.2f - %6.2f + %6.2f GeV 90%%C.L.\n", tf->eP, tf->ePmin, tf->ePmax);
	  if(print) printf("    PT =%7.2f - %6.2f + %6.2f GeV 90%%C.L.\n", tf->ePy, tf->ePYmin, tf->ePYmax);	   
    }
  if (sqrt(t->TX()*t->TX()+t->TY()*t->TY())>0.2)  // use PT
    {
      p=tf->ePy;
      pmin=tf->ePYmin;
      pmax=tf->ePYmax;
	  if(print) printf("    P3D=%7.2f - %6.2f + %6.2f GeV 90%%C.L.\n", tf->eP, tf->ePmin, tf->ePmax);
	  if(print) printf("  ->PT =%7.2f - %6.2f + %6.2f GeV 90%%C.L.\n", tf->ePy, tf->ePYmin, tf->ePYmax);	   
    }
  t0->SetP(p);
  
	if(print){
		double rms, rmst, rmsl;
		DTRMSTL(t, &rms, &rmst, &rmsl);
		printf("DTRMS (Sp,T,L) = ( %.4lf, %.4lf, %.4lf)\n", rms, rmst, rmsl);
	}

  

  if(t!=t0) delete t;
  return tf;
}


double EdbEDAUtil::CalcDistance(EdbSegP *s1,EdbSegP *s2, double z){
	// calculate distance
	double x1 = s1->X() + s1->TX()*(z-s1->Z());
	double y1 = s1->Y() + s1->TY()*(z-s1->Z());
	double x2 = s2->X() + s2->TX()*(z-s2->Z());
	double y2 = s2->Y() + s2->TY()*(z-s2->Z());
	return sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

double EdbEDAUtil::DTRMS(EdbTrackP *t){
	int nseg = t->N();
	if(nseg==1) return 0.0;
	double rms=0.0;
	for(int i=0;i<nseg-1;i++){
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);
		double dtx = s1->TX() - s2->TX();
		double dty = s1->TY() - s2->TY();
		rms+=dtx*dtx+dty*dty;
	}
	rms /= (nseg-1)*2;
	rms = sqrt(rms);
	return rms;
}

double EdbEDAUtil::DTRMSelectron(EdbTrackP *t){
	int nseg = t->N();
	if(nseg==1) return 0.0;
	double rms=0.0;
	EdbSegP *s0 = t->GetSegmentFirst();
	
	for(int i=0;i<nseg-1;i++){
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);
		double dtx  = s1->TX() - s2->TX();
		double dty  = s1->TY() - s2->TY();
		double dt2  = dtx*dtx+dty*dty;
		double dz   = fabs( s2->Z()-s0->Z() ) / 1.3; // lead thickness
		double correction_factor2 = exp( -2*dz/5612 );
		dt2 *= correction_factor2;
		rms+=dt2;
	}
	rms /= (nseg-1)*2;
	rms = sqrt(rms);
	return rms;
}


//	new TGMsgBox(gClient->GetRoot(),gEve->GetMainWindow(),"baka", "hoge", kMBIconAsterisk, kMBYes|kMBNo|kMBOk,&ret)
//		new TGMsgBox(gClient->GetRoot(),gEve->GetMainWindow(),"Confirmation", 
//			Form("%s already exist.\nDo you want to overwrite?",filename), kMBIconAsterisk, kMBYes|kMBNo,&ret);
//	kMBYes, kMBNo, kMBOk, kMBApply, kMBRetry, kMBIgnore, kMBCancel, kMBClose, kMBYesAll, kMBNoAll, kMBNewer, kMBAppend, kMBDismiss
//	kMBIconStop, kMBIconQuestion, kMBIconExclamation, kMBIconAsterisk

void EdbEDAUtil::ErrorMessage(const char *message){
	ErrorMessage("Error", message);
}
void EdbEDAUtil::ErrorMessage(const char *title, const char *message){
	printf("%s : %s\n", title, message);
	new TGMsgBox(gClient->GetRoot(), gEve? gEve->GetMainWindow() : 0,"Error", message, kMBIconAsterisk, kMBOk);
}

bool EdbEDAUtil::AskYesNo(const char *message){
	// Open a message dialog and ask you YES or No.
	// if YES, return kTRUE
	// if NO,  return kFALSE
	// else return kFALSE
	
	int retval=0;
	new TGMsgBox(gClient->GetRoot(), gEve? gEve->GetMainWindow() : 0,
	"EdbEDAUtil::AskYesNo()", message, kMBIconQuestion, kMBYes|kMBNo, &retval);
	
	if(retval==kMBYes) return kTRUE;
	if(retval==kMBNo) return kFALSE;
	else return kFALSE;
}

EdbTrackP * EdbEDAUtil::CleanTrack(EdbTrackP *t0){
	// remove fake-segment/Microtrack/pl57. return a pointer to newly created track.;
	// if Nseg<=4, no fake-segment will be removed.
	// if fake-segments exist, create a new track without the fake segments and return the pointer.
	// 1st segment will not be rejected as fake.
	// segment on pl57 will be rejected every time.
	
	EdbTrackP *t1 = new EdbTrackP;
	t1->Copy(*t0);
	// remove microtrack and pl 57.
	for(int i=0; i<t1->N()&&t1->N()>1; i++){
		EdbSegP *s = t1->GetSegment(i);
		if(s==NULL) break;
		if(s->Side()!=0) {
			printf("Microtrack removed pl%d\n", s->Plate());
			t1->RemoveSegment(s);
			i--;
		}
		if(s->Plate()>=57){
			printf("Segment on plate 57 removed.\n");
			t1->RemoveSegment(s);
			i--;
		}
	}
	
	t1->SetCounters();

	// if nseg<=4 return as it is.
	if(t1->N()<5) return t1;

	// calculate rms, rejecting each 1 of the segments.
	// DTRMS1Kink() give rms assuming 1 kink.
	double *rms = new double[t1->N()];
	for(int i=0; i<t1->N(); i++){
		EdbTrackP *t = new EdbTrackP;
		t->Copy(*t1);
		t->RemoveSegment( t->GetSegment(i));
		rms[i] = DTRMS1Kink(t);
		delete t;
	}
	
	// calculate the significance of strange segment.
	TObjArray *fakeseg = new TObjArray;
	for(int i=1; i<t1->N(); i++) {
		// think for the case to reject 1 segment.
		
		// calculate mean(rms) for the cases rejecting other segments.
		int k=0;
		double mean=0.0;
		for(int j=0; j<t1->N(); j++) if(j!=i) {mean+=rms[j]; k++;}
		mean/=k;
		
		// calculate error of rms (rms/sqrt((double)2*k))
		double rmserr = rms[i]/sqrt((double)2*k);
		
		// if mean-rms is bigger than 2 sigma, flag this segment to be rejected.
		int flag = mean-rms[i]>rmserr*2;
		EdbSegP *s = t1->GetSegment(i);
		if(flag)fakeseg->Add(s);
		
		if(flag)printf("CleanTrack Trkid %4d pid %2d id %6d rms %.4lf mean %.4lf rmserr %.4lf significance %4.1lf\n", t1->ID(), s->PID(), s->ID(), rms[i], mean, rmserr, (mean-rms[i])/rmserr);
	}
	delete rms;

	// if fake-segment exist, copy the track and remove the segment. return the pointer.
	if(fakeseg->GetEntriesFast()){
		EdbTrackP *t = new EdbTrackP;
		t->Copy(*t1);
		for(int i=0;i<fakeseg->GetEntriesFast();i++){
			t->RemoveSegment( (EdbSegP *) fakeseg->At(i));
		}
		delete t1;
		t1=t;
	}
	
	delete fakeseg;

	t1->SetCounters();
	
	// return the track microtrack are rejected.
	return t1;
}


double EdbEDAUtil::DTRMS1Kink(EdbTrackP *t, int *NKinkAngleUsed){
	// calculate rms of delta-theta assuming 1 kink exist(remove 1 biggest dtheta).
	// return rms(space), not projection. dtheta of 1 plate.
	// if nhale exists, the dtheta is calculated as dtheta=dtheta/sqrt(nhole+1)
	// rejection of 1 kink angle is for more than 4 segments.
	int n = t->N();
	double *DTX = new double [n-1];
	double *DTY = new double [n-1];
	
	int iMaxDT=0;       // isegment of biggest kink angle.
	double MaxDT=0.0; // biggest kink angle
	int ndata=0;      // number of combinations of dtx and dty.
	
	for(int i=0;i<n-1;i++){
		// calculation of kink angle.
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);
		
		double dtx = fabs(s1->TX()-s2->TX());
		double dty = fabs(s1->TY()-s2->TY());
		
		int dpid = s2->PID()-s1->PID();
		
		double dt = sqrt(dtx*dtx + dty*dty);
		if(MaxDT<dt) {
			// recoding of maximum kink angle.
			MaxDT = dt;
			iMaxDT = i;
		}

		dtx /= sqrt((double)dpid);
		dty /= sqrt((double)dpid);

		DTX[i] = dtx;
		DTY[i] = dty;
		ndata++;
	}
	
	// prepare data for rms calculation
	double rms=0.0;
	for(int i=0; i<ndata; i++){
		if(ndata>=2&&i==iMaxDT) {
			// for more than 4 segments.
			// remove a combination of dtx,dty with most biggest kink angle.
			continue;
			ndata--;
		}
		rms += DTX[i]*DTX[i];
		rms += DTY[i]*DTY[i];
	}
	delete DTX;
	delete DTY;
	
	rms /= 2*ndata;
	rms = sqrt(rms); // rms in projection
	rms *= sqrt(2.0);  // rms in space
	if(rms<0.0015) rms = 0.0015; // minimum value of rms (angular resolution).
	
//	printf("t->N()=%d n=%d rms = %lf\n", t->N(), 2*ndata, rms);
	
	if(NKinkAngleUsed) *NKinkAngleUsed = 2*ndata;
	
	return rms;
}

void EdbEDAUtil::CalcDTTransLongi(EdbSegP *s1, EdbSegP *s2, double *dtTransverse, double *dtLongitudinal){
	CalcDTTransLongi( s1->TX(), s1->TY(), s2->TX(), s2->TY(), dtTransverse, dtLongitudinal);
}

void EdbEDAUtil::CalcDTTransLongi(double tx1, double ty1, double tx2, double ty2, double *dtTransverse, double *dtLongitudinal){
	TVector3 v1,v2,v;
	v1.SetXYZ( tx1, ty1, -1. );
	v2.SetXYZ( tx2, ty2, -1. );
	
	float phi = v1.Phi();
	v1.RotateZ( -phi );
	v2.RotateZ( -phi );
	
	*dtLongitudinal = v2.X()-v1.X(); // longitudinal
	*dtTransverse   = v2.Y()-v1.Y(); // transverse. by definition v1->Y()=0.0
}

void EdbEDAUtil::CalcDXTransLongi(EdbSegP *s1, EdbSegP *s2, double *dxt, double *dxl){

	EdbSegP s11(*s1);
	s11.PropagateTo(s2->Z());
	
	TVector3 v1,vx1,vx2;
	v1.SetXYZ( s1->TX(), s2->TY(), -1. );
	float phi = v1.Phi();
	
	vx1.SetXYZ( s11.X()-s1->X(), s11.Y()-s1->Y(), -1.);
	vx2.SetXYZ( s2->X()-s1->X(), s2->Y()-s1->Y(), -1.);
	
	vx1.RotateZ( -phi );
	vx2.RotateZ( -phi );
	
	*dxl = vx2.X()-vx1.X(); // longitudinal
	*dxt = vx2.Y()-vx1.Y(); // transverse. by definition v1->Y()=0.0
	
}


double EdbEDAUtil::DTRMSTL(EdbTrackP *t, double *rmsspace, double *rmstransverse, double *rmslongitudinal, int *ndata){
	// calculate dtheta-RMS for each components (Spase angle, Transverse angle, Longitudinal angle)
	// normally RMS(Transverse) is smaller than RMS(Longitudinal). 
	// Transverse component doesn't have dependence on the angle.
	// ndata will be filled with number of data for transverse or longitudinal.
	// so that, actual ndata for rmsspace need to be multiplied by 2.
	int n = t->N()-1;
	
	double rmst = 0.0;
	double rmsl = 0.0;
	
	for(int i=0;i<n;i++){
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);
		
		// skip if s2 is on pl57
		if(s2->Plate()>=57){
			n--;
			continue;
		}
		
		double dtt,dtl;
		CalcDTTransLongi(s1,s2,&dtt, &dtl);
		
		float dplate = s2->Plate()-s1->Plate();
		dtt/=sqrt(dplate);
		dtl/=sqrt(dplate);
		
		rmst += dtt*dtt;
		rmsl += dtl*dtl;
	}
	
	if(n==0){
		rmst = -1.0;
		rmsl = -1.0;
	}
	else {
		rmst = sqrt(rmst/n); // n= t->N()-1
		rmsl = sqrt(rmsl/n);
	}
	
	// minimum value of rms (angular resolution).
	// double angle = sqrt(t->TX()*t->TX()+t->TY()*t->TY());
	// if(rmst<=0.0015*sqrt(2.0)) rmst=0.0015*sqrt(2.0); 
	// if(rmsl<=0.0015*(1+angle*5)*sqrt(2.0)) rmsl=0.0015*(1+angle*5)*sqrt(2.0);
	double rms = sqrt(rmsl*rmsl+rmst*rmst)/sqrt(2.0);

	// if rmst > rmsl, means to low nseg or scattering dominant.
	// to have more reliable value, use rms.
	// if(rmst>rmsl) rmst=rmsl=rms;
	
	*rmsspace = rms;
	*rmstransverse = rmst;
	*rmslongitudinal = rmsl;
	if(ndata) *ndata=n;

	return rmst;
}

double EdbEDAUtil::DTRMSTL1Kink(EdbTrackP *t0, double *rmsspace, double *rmstransverse, double *rmslongitudinal, int *NKinkAngleUsed){
	// calculate rms of delta-theta assuming 1 kink exist(remove 1 biggest dtheta).
	// return rms(transverse).
	// if nhale exists, the dtheta is calculated as dtheta=dtheta/sqrt(nhole+1).
	// rejection of 1 kink angle is for more than 4 segments.
	// if a segment is on pl57, this segment will not used.

	// remove Micro-tracks.
	EdbTrackP *t = new EdbTrackP (*t0);
	for(int i=0;i<t->N()&&t->N()>1;i++){
		EdbSegP *s = t->GetSegment(i);
		if(s->Side()!=0) t->RemoveSegment(s); // if it's microtrack, remove.
	}
	printf("microtrack remove nseg %d->%d\n", t0->N(), t->N());
	int nseg = t->N();
	
	int iMaxDT=0;       // isegment of biggest kink angle.
	double MaxDT=0.0;   // biggest kink angle
	
	for(int i=0;i<nseg-1;i++){
		// calculation of kink angle.
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);
		
		// skip if s2 is on pl57
		if(s2->Plate()>=57) continue;
		
		double dtl, dtt;
		CalcDTTransLongi(s1,s2, &dtt, &dtl);
		double dt = sqrt(dtt*dtt+dtl*dtl);

		if(MaxDT<dt) { // recoding of maximum kink angle.
			MaxDT = dt;
			iMaxDT = i;
		}
	}
	
	if(nseg<4)  DTRMSTL          (t,         rmsspace, rmstransverse, rmslongitudinal, NKinkAngleUsed);
	else        DTRMSTLGiven1Kink(t, iMaxDT, rmsspace, rmstransverse, rmslongitudinal, NKinkAngleUsed);
	
	delete t;
	
	return *rmstransverse;
}

double EdbEDAUtil::DTRMSTLGiven1Kink(EdbTrackP *t, int iKink, double *rmsspace, double *rmstransverse, double *rmslongitudinal, int *NKinkAngleUsed){
	// calculate rms of delta-theta assuming 1 "given" kink.
	// return rms(transverse).
	// if nhale exists, the dtheta is calculated as dtheta=dtheta/sqrt(nhole+1).
	// if a segment is on pl57, this segment will not used.

	int nseg = t->N();
	int ndata=0;
	double rmst=0.0, rmsl=0.0;
	for(int i=0;i<nseg-1;i++){
		// removal of the given delta-theta.
		if(i==iKink) continue;
		// calculation of kink angle.
		EdbSegP *s1 = t->GetSegment(i);
		EdbSegP *s2 = t->GetSegment(i+1);
		
		// skip if s2 is on pl57
		if(s2->Plate()>=57) continue;
		
		double dtl, dtt;
		CalcDTTransLongi(s1,s2, &dtt, &dtl);
		
		float dplate = s2->Plate()-s1->Plate();

		dtt /= sqrt(dplate); // if there is a hole, consider the effect of scatterning.
		dtl /= sqrt(dplate);
		
		rmst += dtt*dtt;
		rmsl += dtl*dtl;
		
		ndata++;
	}
	if(ndata==0) {
		*rmsspace=*rmstransverse=*rmslongitudinal=-1.;
		return -1;
	}
	
	if(ndata==0){
		rmst = -1.0;
		rmsl = -1.0;
	}
	else {
		rmst = sqrt(rmst/ndata); // n= t->N()-1
		rmsl = sqrt(rmsl/ndata);
	}
	
	double rms = sqrt(rmsl*rmsl+rmst*rmst)/sqrt(2.0);
	
	*rmsspace = rms;
	*rmstransverse = rmst;
	*rmslongitudinal = rmsl;
	if(NKinkAngleUsed) *NKinkAngleUsed = ndata;
	
	return rmst;
}

int EdbEDAUtil::InputNumberInteger(const char *message, int idefault){
	
	char buf[20];
	sprintf(buf,"%d",idefault);
	new TGInputDialog(gClient->GetRoot(), gEve?gEve->GetMainWindow():0, message, buf, buf);
	
	sscanf(buf,"%d",&idefault);
	return idefault;
	
}


double EdbEDAUtil::InputNumberReal(const char *message, double default_num, TGNumberFormat::EStyle es){
	// create a number entry window and return the number
	
	char buf[20];
	sprintf(buf,"%lf",default_num);
	new TGInputDialog(gClient->GetRoot(), gEve?gEve->GetMainWindow():0, message, buf, buf);
	
	sscanf(buf,"%lf",&default_num);
	return default_num;
}

//	EStyle :: kNESInteger, kNESRealOne, kNESRealTwo, kNESRealThree, kNESRealFour, kNESReal


int EdbEDAUtil::InputID(const char *message, EdbID &id){
	// create a number entry window and return the number
	// return 1 if successful, 0 if cancel.
	int ret=0;
	new EdbIDDialog(message, id, &ret);
	return ret;
}



char *EdbEDAUtil::FindProcDirClient(){
	// Find ProcDirClient from the current path.
	// "ONLINE" is the hard coded name for ProcDirClient.
	TString buf=gSystem->WorkingDirectory();
	
	int index = buf.Index("ONLINE");
	if(index==kNPOS) {
		printf("FindProcDirClient   : Not found \"ONLINE\" in the path. set eProcDirClient = \"..\"\n");
		return "..";
	}
	
	char *path = new char[256];
	strncpy(path, buf.Data(), index+6);
	path[index+6]='\0';
	printf("FindProcDirClient = %s\n", path);
	return path;
}

int EdbEDAUtil::FindBrickIDFromPath(){
	TString buf=gSystem->WorkingDirectory();
	int index = buf.Index("ONLINE");
	if(index==kNPOS) {
		printf("FindBrickIDFromPath : Not found \"ONLINE\" in the path. set eBrickID = 0\n");
		return 0;
	}
	int ibrick=0;
	sscanf(buf.Data()+index+8,"%d", &ibrick);
	printf("FindBrickIDFromPath : set eBrickID = %d\n", ibrick);
	return ibrick;
}


EdbPVRec * EdbEDAUtil::ReadFeedbackPVR(const char *filename){
	// Read feedback file format (ver 2009 Oct), and return EdbPVRec object.
	// eTracks, eVTX, patterns are filled.
	// stand alone function.
	
	// if filename is not given, open file browser.
	if(filename==NULL){
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Feedback File", "*.feedback", "All files","*",0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve?gEve->GetMainWindow():0, kFDOpen, fi);
		filename = fi->fFilename;
	}
	
	if(NULL==filename) {
		EdbEDAUtil::ErrorMessage("No filename. stop.");
		return NULL;
	}
	
	FILE *fp = fopen(filename,"rt");
	if(fp==NULL) {
		ErrorMessage(Form("File open error : %s . stop.\n", filename));
		return NULL;
	}
	
	EdbPVRec *pvr = new EdbPVRec;
	EdbScanCond cond;
	EdbTrackP *t=0;

	
	char buf[256];
	
	while(fgets(buf,sizeof(buf),fp)){
		TString line=buf;
		printf("%s", line.Data());

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
			t->SetFlag(particle_id);
			pvr->AddTrack(t);
			
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
		}
		else {
		    cout << "EdbEDAUtil::ReadFeedbackPVR Warning: Reading feedback file: items = " << ntokens << ". This is NOT specified!" << endl;
		}
	}
	
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetPID(i);
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetSegmentsPID();
	for(int i=0;i<pvr->Ntracks(); i++) pvr->GetTrack(i)->SetCounters();
	
	printf("\nEdbPVRec summary. %d vertices, %d tracks.\n", pvr->Nvtx(), pvr->Ntracks());
	pvr->Print();
	
	fclose(fp);
	return pvr;
}


EdbPVRec *EdbEDAUtil::ReadMxxPVR(const char *filename){

	if(filename==NULL){
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Mxx file", "*.all", 0,0};
		fi->fFileTypes=filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve?gEve->GetMainWindow():0, kFDOpen, fi);
		filename = fi->fFilename;
	}
	
	if(NULL==filename) {
		EdbEDAUtil::ErrorMessage("No filename. stop.");
		return NULL;
	}
	
	FILE *fp = fopen(filename,"rt");
	if(fp==NULL) {
		ErrorMessage(Form("File open error : %s . stop.\n", filename));
		return NULL;
	}
	
	int npid = 0, plid[100], pid[1000];
	double zlayer[100];

	EdbPVRec *pvr = new EdbPVRec;
	EdbScanCond cond;

	char buf[256];
	for(int i=0;i<100;i++) zlayer[i]=1.1E10;
	
	printf("reanding m-file : %s\n", filename);

	// Header
	// skip comment
	for(;NULL!=fgets(buf,sizeof(buf),fp);){
		if(buf[0]!='%') break;
	}
	
	// skip 1st line
	
	// nPlate
	fgets(buf,sizeof(buf),fp); 
	sscanf(buf,"%d", &npid);

	// pos list
	fgets(buf,sizeof(buf),fp); 
	int n=0,nn;
	for(int j=npid-1;j>=0;j--){
		sscanf(buf+n,"%d%n", &plid[j], &nn);
		plid[j]/=10;
		plid[j]=58-plid[j];
		pid[plid[j]]=j;
		n+=nn;
	}
	// Header end

	// track data
	for(int i=0;NULL!=fgets(buf,sizeof(buf),fp);i++){
		
		int itrk, nseg, pos1, pos2;
		sscanf(buf, "%d %d %d %d", &itrk, &nseg, &pos1, &pos2); // track headers
		
		EdbTrackP *t = new EdbTrackP();
		t->SetID(itrk);
		
		//printf("t %s", buf);
		EdbSegP *s=0;
		for(int j=0;j<nseg;j++){	// data of segmetns.
			int pos,isg,ph;
			double ax,ay,x,y,z;
			fgets(buf,sizeof(buf),fp);
			sscanf(buf,"%d %*d %d %d %lf %lf %lf %lf %lf",
						&pos, &isg, &ph, &ax, &ay, &x, &y, &z);
			pos/=10;
			pos=58-pos;
			if(zlayer[pid[pos]]>=1E10) zlayer[pid[pos]]=z;
			
			s = new EdbSegP(isg, x,y, ax,ay, ph/10000., 0);
			s->SetZ(z);
			s->SetPID(pid[pos]);
			s->SetPlate(pos);
			s->SetTrack(itrk);
			
			// fill COV for vertexing
			s->SetErrors();
			cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
			
			t->AddSegment(s);
			pvr->AddSegment(*s);
			//t->AddSegmentF( new EdbSegP(*((EdbSegP*)(s))) );
			
		}
		s = t->GetSegmentFirst();
		t->Set(itrk, s->X(), s->Y(), s->TX(), s->TY(), 0, 1);
		t->SetZ(s->Z());
		t->SetCounters();
		pvr->AddTrack(t);
	}

	
	for(int i=0;i<pvr->Ntracks();i++){
		EdbTrackP *t = pvr->GetTrack(i);
		// fill COV for vertexing
		t->SetErrors();
		cond.FillErrorsCov(t->TX(), t->TY(), t->COV());
	}
	
	
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetPID(i);
	for(int i=0;i<pvr->Npatterns(); i++) pvr->GetPattern(i)->SetSegmentsPID();

	pvr->Print();
	
	return pvr;
}

void EdbEDAUtil::WritePVRMxx(EdbPVRec *pvr, const char *filename){
	
	if(filename==NULL){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Mxx file", "*.all", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		if(fi->fFilename!=NULL) filename = fi->fFilename;
	}
	
	if(NULL==filename||strlen(filename)==0) {
		return;
	}
	FILE *fp=fopen(filename,"wt");
	if(NULL==fp){
		printf("Couldn't open file : %s\n", filename);
		return;
	}
	
	// Header
	
	// comment
	fprintf(fp,"%% Created by EDA\n");
	time_t tim;
	time(&tim);
	fprintf(fp,"%% on %s",ctime(&tim));
	
	// dummy module id
	fprintf(fp,"       0       0   3   0      0.0   0.0000\n");
	
	// number of plates
	fprintf(fp,"%d\n", pvr->Npatterns());
	
	// list of pos
	for(int i=pvr->Npatterns()-1;i>=0;i--){
		fprintf(fp,"%d ", (58-pvr->GetPattern(i)->Plate())*10+1);
	}
	fprintf(fp,"\n");
	
	// Tracks + Segments
	for(int i=0;i<pvr->Ntracks();i++){
		EdbTrackP *t = pvr->GetTrack(i);
		
		EdbSegP *s1 = t->GetSegmentFirst();
		EdbSegP *s2 = t->GetSegmentLast();
		
		fprintf(fp, "%3d %2d %3d %3d\n", t->ID(), t->N(), (58-s2->Plate())*10+1, (58-s1->Plate())*10+1);
		// write segment information
		for(int j=0;j<t->N();j++){
			EdbSegP *s = t->GetSegment(t->N()-j-1);
				fprintf(fp,"%3d %5d %7d %6d %7.4f %7.4f %8.1f %8.1f %8.1f\n",
					(58-s->Plate())*10+1, t->ID(), s->ID(), (int)(s->W()*10000), s->TX(), s->TY(), s->X(), s->Y(), s->Z());
			
		}
	}
	
	fclose(fp);
	
	printf("file %s was written.\n", filename);

}

void EdbEDAUtil::WriteTracksMxx(TObjArray *tracks, const char *filename){
	
	if(filename==NULL){ // File dialog
		TGFileInfo *fi=new TGFileInfo;
		fi->fIniDir    = StrDup(".");
		const char *filetypes[] = { "Mxx file", "*.all", "All files","*",0,0};
		fi->fFileTypes = filetypes;
		new TGFileDialog(gClient->GetRoot(), gEve->GetMainWindow(), kFDSave, fi);
		if(fi->fFilename!=NULL) filename = fi->fFilename;
	}
	
	if(NULL==filename||strlen(filename)==0) {
		return;
	}
	FILE *fp=fopen(filename,"wt");
	if(NULL==fp){
		printf("Couldn't open file : %s\n", filename);
		return;
	}
	
	// Header
	
	// comment
	fprintf(fp,"%% Created by EDA\n");
	time_t tim;
	time(&tim);
	fprintf(fp,"%% on %s",ctime(&tim));
	
	// dummy module id
	fprintf(fp,"       0       0   3   0      0.0   0.0000\n");
	
	
	// plates information
	std::vector<long > vec;
	std::vector<long >::iterator it;
	
	for(int i=0;i<tracks->GetEntriesFast();i++){
		EdbTrackP *t = (EdbTrackP *) tracks->At(i);
		for(int j=0;j<t->N();j++){
			EdbSegP *s = t->GetSegment(j);
			int ipl=(58-s->Plate())*10+1;
			it = find(vec.begin(), vec.end(), ipl);
			if(it==vec.end()) vec.push_back(  ipl);
		}
	}
	
	sort(vec.begin(), vec.end());
	
	// number of plates
	fprintf(fp,"%d\n", (int) vec.size());
	
	// list of pos
	for(int i=0;i<(int) vec.size();i++){
		fprintf(fp,"%d ", (int) vec[i]);
	}
	fprintf(fp,"\n");
	
	// Tracks + Segments
	for(int i=0;i<tracks->GetEntriesFast();i++){
		EdbTrackP *t = (EdbTrackP *) tracks->At(i);
		
		EdbSegP *s1 = t->GetSegmentFirst();
		EdbSegP *s2 = t->GetSegmentLast();
		
		fprintf(fp, "%3d %2d %3d %3d\n", t->ID(), t->N(), (58-s2->Plate())*10+1, (58-s1->Plate())*10+1);
		// write segment information
		for(int j=0;j<t->N();j++){
			EdbSegP *s = t->GetSegment(t->N()-j-1);
				fprintf(fp,"%3d %5d %7d %6d %7.4f %7.4f %8.1f %8.1f %8.1f\n",
					(58-s->Plate())*10+1, t->ID(), s->ID(), (int)(s->W()*10000), s->TX(), s->TY(), s->X(), s->Y(), s->Z());
			
		}
	}
	
	fclose(fp);
	
	printf("file %s was written.\n", filename);

}


void EdbEDAUtil::MakePVRFromTracksArray(TObjArray *tracks_or_segments, EdbPVRec& pvr){
	// convert tracks or segments array into pvr.
	// if pointers of segments will not be consistent.
	
	EdbScanCond cond;
	
	for(int i=0; i<tracks_or_segments->GetEntriesFast(); i++){
		TObject *o = tracks_or_segments->At(i);
		
		if(IsSegment(o)){
			EdbSegP *s = (EdbSegP *)o;
			
			// fill COV for vertexing
			if(s->SX()==0){
				s->SetErrors();
				cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
			}
			
			// Add segment in PVRec and Track, keeping consistency of pointer in them.
			EdbSegP *s_in_pattern = pvr.AddSegment(*s);
			
			EdbTrackP *t = new EdbTrackP(*s_in_pattern);
			t->AddSegment( s_in_pattern);
			pvr.AddTrack(t);
		}
		
		else if (IsTrack(o)){
			EdbTrackP *t = (EdbTrackP *) o;
			
			for(int j=0; j<t->N(); j++){
				EdbSegP *s = t->GetSegment(j);
				EdbSegP *s_in_pattern = pvr.AddSegment(*s);
			}
			
			pvr.AddTrack(t);
		}
	}
}

void EdbEDAUtil::FillTracksFromPatterns(EdbPVRec *pvr){
	// convert segments in patterns into Tracks in EdbPVRec object.
	// each track has 1 segment
	
	EdbScanCond cond;
	
	for(int i=0; i<pvr->Npatterns(); i++){
		EdbPattern *pat = pvr->GetPattern(i);
		
		for(int j=0; j<pat->N(); j++){
			EdbSegP *s = pat->GetSegment(j);
			// fill COV for vertexing
			if(s->SX()==0){
				s->SetErrors();
				cond.FillErrorsCov(s->TX(), s->TY(), s->COV());
			}
			EdbTrackP *t = new EdbTrackP(*s);
			t->AddSegment( s);
			pvr->AddTrack(t);
		}
	}
	
}

