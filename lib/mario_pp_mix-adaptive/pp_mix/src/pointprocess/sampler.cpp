#include "point2pattern.h"
#include "pointprocess.h"
#include <math.h>
#include <stdio.h>
#include "sampler.h"

void Sampler::GenerateDominatingPoisson(Point *headPoint, 
					long int *GeneratedPoints,
					long int *LivingPoints,
					long int *NoP,
					PointProcess *pp,
					gsl_rng *r)
{
  int i,l;
  double xtemp, ytemp, *marktemp, Xdim, Ydim;
  struct Point *TempPoint;
  char InWindow;
  *NoP = 0;
  //printf("Generating StraussProcess Poisson 1\n");
  Xdim = pp->Xmax-pp->Xmin;
  Ydim = pp->Ymax-pp->Ymin;
  *GeneratedPoints = gsl_ran_poisson(r, pp->TotalBirthRate);
  // printf("GeneratedPoints :%ld\n",*GeneratedPoints);
  *LivingPoints = *GeneratedPoints;
  //printf("Generating StraussProcess Poisson 2 %d\n",*GeneratedPoints);
  //FILE *out222;
  //out222 = fopen("dom.dat","w");

  // Generating a realisatin of the dominating Poisson process.
  // This is done by generating n points, where n is Poisson
  // distributed with mean TotalBirthRate. Of the generated
  // points we only keep the once with InWindow==1. The latter
  // is part of a thinning procedure.

  for (i=1; i<=*GeneratedPoints ; i++){
    //pp->NewEvent(&TempPoint->X,&TempPoint->Y,&InWindow,r);
    pp->NewEvent(&xtemp,&ytemp,&marktemp, &InWindow,r);
    if(InWindow==1){
      TempPoint = (struct Point *) malloc(sizeof *TempPoint);
      if(TempPoint == NULL){
	printf("Cannot allocate memory for TempPoint in Sampler::GenerateDominatingPoisson\n");
	exit(1);
      }
      TempPoint->X = xtemp;
      TempPoint->Y = ytemp;
      if(pp->MarkLength>0) TempPoint->Mark = marktemp;
      TempPoint->R = gsl_rng_uniform(r);
      TempPoint->No = 0;
      TempPoint->next = headPoint->next;
      headPoint->next = TempPoint;
      *NoP = *NoP + 1;
      //printf("mark: %f\n",*(TempPoint->Mark));
      /*
      fprintf(out222,"%f\t%f",TempPoint->X,TempPoint->Y);
      if(pp->MarkLength>0){
	for(l=0;l<pp->MarkLength;l++)
	  fprintf(out222,"\t%f",*(TempPoint->Mark+l));
      }
      fprintf(out222,"\n");
      */
    }
  }

  //fclose(out222);
  *LivingPoints = *NoP;
  /*
  FILE  *out;
  out = fopen("LivingPoints.dat","a");
  fprintf(out,"%ld\n",-(*LivingPoints));
  fclose(out);  
  printf("NoP: %ld\n",*NoP);
  */
  //exit(0);
};



void Sampler::Complement(Point2Pattern *Wpattern, Point2Pattern *Inputp2p, 
			 long int *xnop, PointProcess *pp, gsl_rng *r){

  long int NoOriginalWPoints, N, NoSurvivingWPoints, 
    dummy_int1, dummy_int2;
  double u, x, y, *mark, lnCondIntRatio;
  char InWindow;
  Point2 Point;
  Point2Pattern Xprocess(Inputp2p->Xmin,Inputp2p->Xmax,
			 Inputp2p->Ymin,Inputp2p->Ymax, 
			 Inputp2p->MaxXCell,Inputp2p->MaxYCell);

  Wpattern->Empty(); 
  struct Point *headOriginalW, *dummyOriginalW, *tempOriginalW;
  headOriginalW = (struct Point *) malloc(sizeof *headOriginalW);
  if(headOriginalW == NULL){
    printf("Cannot allocate memory for headOriginalW in Sampler::Complement\n");
    exit(1);
  }
  dummyOriginalW = (struct Point *) malloc(sizeof *dummyOriginalW);  
  if(dummyOriginalW == NULL ){
    printf("Cannot allocate memory for dummyOriginalW in Sampler::Complement\n\n");
    exit(1);
  }
  headOriginalW->next = dummyOriginalW; dummyOriginalW->next = dummyOriginalW;

  GenerateDominatingPoisson(headOriginalW, &dummy_int1,
			    &dummy_int2,
			    &NoOriginalWPoints,
			    pp,
			    r);
  printf("NoOriginalWPoints: %ld\n",NoOriginalWPoints);
  long int i;
  tempOriginalW = headOriginalW->next;
  i=0;
  while(tempOriginalW != dummyOriginalW){
    i++;
    tempOriginalW = tempOriginalW->next;
  }
  //printf("i: %ld\n",i);

  Inputp2p->NoP = Inputp2p->Count();
  Inputp2p->CopyTo(&Xprocess);
  pp->CalcBeta(&Xprocess);
  Xprocess.DumpToFile("sample3.dat");

  NoSurvivingWPoints = 0;
  //printf("NoOriginalWPoints: %ld\n",NoOriginalWPoints);
  printf("TotalBirthRate: %f\n",pp->TotalBirthRate);
  printf("Xprocess.NoP: %ld\n",Xprocess.NoP);
  //printf("NoP: %ld\n",Xprocess.Count());

  while(NoOriginalWPoints>0){
    //printf("%d\n",Xprocess.NoP);
    u = gsl_rng_uniform(r);
    if(u< ((double)Xprocess.NoP+pp->TotalBirthRate)/
       ((double)(Xprocess.NoP+NoOriginalWPoints)+pp->TotalBirthRate)){
      // Birth/death in X process
      u = gsl_rng_uniform(r);
      if(u<(((double) Xprocess.NoP)/
	    ((double) Xprocess.NoP + pp->TotalBirthRate))){
	// Death of point in X
	N = ((long int) floor(((double)Xprocess.NoP)*gsl_rng_uniform(r))) + 1;
	if(N>Xprocess.NoP) N = Xprocess.NoP;
	Xprocess.Remove(N);
      }else{
	// Potential birth of point in X
	pp->NewEvent(&x,&y,&mark, &InWindow,r);
	if(InWindow==1){
	  Point.X = x;
	  Point.Y = y;
	  pp->Beta(&Point);
	  lnCondIntRatio = pp->lnCondIntRatio(&Point,&Xprocess);
	  u = gsl_rng_uniform(r);
	  //printf("%f %f %f %f\n",Point.X,Point.Y,
	  //Point.Beta,exp(lnCondInt)/Point.Beta);
	  if(log(u)<lnCondIntRatio){ 
	    Xprocess.Add(&Point);
	  }
	}
      }
    }else{
      // Death in W process    
      if(headOriginalW->next == dummyOriginalW) printf("UPS!\n");
      tempOriginalW = headOriginalW->next;
      Point.X = headOriginalW->next->X;
      Point.Y = headOriginalW->next->Y;
      headOriginalW->next = tempOriginalW->next;
      free(tempOriginalW);

      pp->Beta(&Point);
      lnCondIntRatio = pp->lnCondIntRatio(&Point,&Xprocess);
      u = gsl_rng_uniform(r);
      if(log(u)>lnCondIntRatio){
	NoSurvivingWPoints++;
	Point.No = NoSurvivingWPoints;
	Wpattern->Add(&Point);
      }
      NoOriginalWPoints --;	
    }
    //printf("NoWPoints: %d\tNoOriginalWPoints: %d\n",
    //   NoWPoints,NoOriginalWPoints);
  }

/*
    Wprocess = poisson(beta);
    OriginalWpointStillAlive = NoP(Wprocess);
    
    u = runif
    
    
    

   */
  //printf("NoP Xprocess: %d\n",Xprocess.Count());
  *xnop = Xprocess.Count();

}

void Sampler::Complement(Point2Pattern *Wpattern, 
			 Point2Pattern *Inputp2p, 
			 Point2Pattern *ReturnXprocess, 
			 PointProcess *pp, gsl_rng *r){

  long int NoOriginalWPoints, N, NoSurvivingWPoints, 
    dummy_int1, dummy_int2;
  double u, x, y, *mark, lnCondIntRatio;
  char InWindow;
  Point2 Point;
  Point2Pattern Xprocess(Inputp2p->Xmin,Inputp2p->Xmax,
			 Inputp2p->Ymin,Inputp2p->Ymax, 
			 Inputp2p->MaxXCell,Inputp2p->MaxYCell);
  ReturnXprocess->Empty();
  Wpattern->Empty(); 
  struct Point *headOriginalW, *dummyOriginalW, *tempOriginalW;
  headOriginalW = (struct Point *) malloc(sizeof *headOriginalW);
  if(headOriginalW == NULL){
    printf("Cannot allocate memory for headOriginalW in Sampler::Complement\n");
    exit(1);
  }
  dummyOriginalW = (struct Point *) malloc(sizeof *dummyOriginalW);  
  if(dummyOriginalW == NULL){
    printf("Cannot allocate memory for dummyOriginalW in Sampler::Complement\n");
    exit(1);
  }
  headOriginalW->next = dummyOriginalW; dummyOriginalW->next = dummyOriginalW;

  GenerateDominatingPoisson(headOriginalW, &dummy_int1,
			    &dummy_int2,
			    &NoOriginalWPoints,
			    pp,
			    r);
  printf("NoOriginalWPoints: %ld\n",NoOriginalWPoints);
  long int i;
  tempOriginalW = headOriginalW->next;
  i=0;
  while(tempOriginalW != dummyOriginalW){
    i++;
    tempOriginalW = tempOriginalW->next;
  }
  printf("i: %ld\n",i);

  Inputp2p->NoP = Inputp2p->Count();
  Inputp2p->CopyTo(&Xprocess);
  pp->CalcBeta(&Xprocess);
  Xprocess.DumpToFile("sample3.dat");

  NoSurvivingWPoints = 0;
  //printf("NoOriginalWPoints: %ld\n",NoOriginalWPoints);
  printf("TotalBirthRate: %f\n",pp->TotalBirthRate);
  printf("Xprocess.NoP: %ld\n",Xprocess.NoP);
  //printf("NoP: %ld\n",Xprocess.Count());

  while(NoOriginalWPoints>0){
    //printf("%d\n",Xprocess.NoP);
    u = gsl_rng_uniform(r);
    if(u< ((double)Xprocess.NoP+pp->TotalBirthRate)/
       ((double)(Xprocess.NoP+NoOriginalWPoints)+pp->TotalBirthRate)){
      // Birth/death in X process
      u = gsl_rng_uniform(r);
      if(u<(((double) Xprocess.NoP)/
	    ((double) Xprocess.NoP + pp->TotalBirthRate))){
	// Death of point in X
	N = ((long int) floor(((double)Xprocess.NoP)*gsl_rng_uniform(r))) + 1;
	if(N>Xprocess.NoP) N = Xprocess.NoP;
	Xprocess.Remove(N);
      }else{
	// Potential birth of point in X
	pp->NewEvent(&x,&y,&mark, &InWindow,r);
	if(InWindow==1){
	  Point.X = x;
	  Point.Y = y;
	  pp->Beta(&Point);
	  lnCondIntRatio = pp->lnCondIntRatio(&Point,&Xprocess);
	  u = gsl_rng_uniform(r);
	  //printf("%f %f %f %f\n",Point.X,Point.Y,
	  //Point.Beta,exp(lnCondInt)/Point.Beta);
	  if(log(u)<lnCondIntRatio){ 
	    Xprocess.Add(&Point);
	  }
	}
      }
    }else{
      // Death in W process    
      if(headOriginalW->next == dummyOriginalW) printf("UPS!\n");
      tempOriginalW = headOriginalW->next;
      Point.X = headOriginalW->next->X;
      Point.Y = headOriginalW->next->Y;
      headOriginalW->next = tempOriginalW->next;
      free(tempOriginalW);

      pp->Beta(&Point);
      lnCondIntRatio = pp->lnCondIntRatio(&Point,&Xprocess);
      u = gsl_rng_uniform(r);
      if(log(u)>lnCondIntRatio){
	NoSurvivingWPoints++;
	Point.No = NoSurvivingWPoints;
	Wpattern->Add(&Point);
      }
      NoOriginalWPoints --;	
    }
    //printf("NoWPoints: %d\tNoOriginalWPoints: %d\n",
    //   NoWPoints,NoOriginalWPoints);
  }

/*
    Wprocess = poisson(beta);
    OriginalWpointStillAlive = NoP(Wprocess);
    
    u = runif
    
    
    

   */
  //printf("NoP Xprocess: %d\n",Xprocess.Count());
  //*xnop = Xprocess.Count();
  //Copy Xpross to ReturnXprocess
  if(ReturnXprocess != NULL){
    Xprocess.CopyTo(ReturnXprocess);
  }

}


void Sampler::Forward(long int TS, long int TT, int TX, int TY,
		      struct Point *Proposal, long int *DDD,
		      Point2Pattern *p2p, PointProcess *pp){

  long int XCell, YCell, iMark;
  double tmpR, TempGamma[2], TempI;
  double TempLnGamma[2], TempLnI;
  struct Point2 *TempCell, *TempCell2;

  /* Birth */
  if(TT==1){
    XCell = (int) floor((Proposal->X-p2p->Xmin)/p2p->XCellDim); 
    if(XCell>p2p->MaxXCell) XCell = p2p->MaxXCell;
    YCell = (int) floor((Proposal->Y-p2p->Ymin)/p2p->YCellDim);
    if(YCell>p2p->MaxYCell) YCell = p2p->MaxYCell;

    TempCell = (struct Point2 *) malloc(sizeof *TempCell);
    if(TempCell==NULL){
      printf("Cannot allocate memory for TempCell in Sampler::Forward\n");
      exit(1);
    }
    TempCell->No = Proposal->No;
    TempCell->X = Proposal->X;
    TempCell->Y = Proposal->Y;
    if(pp->MarkLength>0){
      TempCell->Mark = (double *) calloc(pp->MarkLength,sizeof(double));      
      if(TempCell->Mark==NULL){
	printf("Cannot allocate memory for Mark in Sampler::Forward\n");
	exit(1);
      }
      for(iMark=0;iMark<pp->MarkLength;iMark++){ 
      *(TempCell->Mark + iMark) = *(Proposal->Mark + iMark);
      //printf(":> %f\n",*(TempCell->Mark + iMark));
      }
    }


    tmpR = Proposal->R;
    TempCell->next = p2p->headCell[XCell][YCell]->next;
    p2p->headCell[XCell][YCell]->next = TempCell;
    TempCell->InLower[0]=0;
    TempCell->InLower[1]=0;

    TempGamma[0] = 1.0; TempGamma[1] = 1.0;    
    TempLnGamma[0] = 0.0; TempLnGamma[1] = 0.0;
    
    // Determine rectangular region of neighbours "measured in cells"
    long int rx,ry,fx,fy,tx,ty,xco,yco;
    rx = (int) floor(pp->InteractionRangeFct(TempCell)/p2p->XCellDim+1.0);
    ry = (int) floor(pp->InteractionRangeFct(TempCell)/p2p->YCellDim+1.0);

    // The rectangular region covers cell (fx,fy) to (tx,ty)
    if((XCell+rx)<=p2p->MaxXCell) tx=XCell+rx; else tx = p2p->MaxXCell;
    if((YCell+ry)<=p2p->MaxYCell) ty=YCell+ry; else ty = p2p->MaxYCell;
    if((XCell-rx)>=0) fx=XCell-rx; else fx = 0;
    if((YCell-ry)>=0) fy=YCell-ry; else fy = 0;
    //printf("InteractionRange: %f\n",pp->InteractionRange);
    //printf("rx: %ld ry: %ld\n",rx,ry);
    for(xco = fx; xco <= tx; xco++){
      for(yco = fy; yco <= ty; yco++){
	//printf("xco: %ld yco: %ld\n",xco,yco);
	TempCell2 = p2p->headCell[xco][yco]->next;
	while(TempCell2!=TempCell2->next){
	  TempI = pp->Interaction(TempCell,TempCell2);
	  if(TempCell2->InLower[0]==1) TempGamma[0] = TempGamma[0]*TempI;
	  if(TempCell2->InLower[1]==1) TempGamma[1] = TempGamma[1]*TempI;

	  TempLnI = log(pp->Interaction(TempCell,TempCell2));
	  if(TempCell2->InLower[0]==1) TempLnGamma[0] += TempLnI;
	  if(TempCell2->InLower[1]==1) TempLnGamma[1] += TempLnI;

	  TempCell2 = TempCell2->next; 
	  
	}
      }
    }
    //   TempCell->InLower[0],TempCell->InLower[1], tmpR);
    //if(tmpR <= TempGamma[1] ){ 
    if(log(tmpR) < TempLnGamma[1] ){ 
      TempCell->InLower[0]=1;
      p2p->UpperLiving[0] ++;
      pp->UpdateInteractionRange(TempCell);
    }
    //if(tmpR <= TempGamma[0] ){ 
    if(log(tmpR) < TempLnGamma[0] ){ 
      TempCell->InLower[1]=1;
      p2p->UpperLiving[1] ++;
    }
    //printf("> %f %f %d %d\n",TempGamma[0],TempGamma[1],
    //   TempCell->InLower[0],TempCell->InLower[1]);
  }
  /* Death */
  if(TT==0){
    TempCell=p2p->headCell[TX][TY];
    while(TempCell->next->No != *DDD){
      TempCell = TempCell->next;
      if(TempCell->next == TempCell){
	printf("hmm...\n"); 
	printf("TX: %d TY: %d\n",TX,TY);
	long ix,iy;
	for(ix=0;ix<=p2p->MaxXCell;ix++){
	  for(iy=0;iy<=p2p->MaxXCell;iy++){
	    printf("%ld %ld %ld\n",ix,iy,*DDD);
	    TempCell=p2p->headCell[ix][iy];
	    while(TempCell->next != TempCell){
	      //printf("%f %f %ld\n",TempCell->X,TempCell->Y,TempCell->No);
	      if(TempCell->No == *DDD){
		printf("TX: %d TY: %d\n",TX,TY);
		printf("!!! %lf %lf %ld\n",
		       TempCell->X,TempCell->Y,TempCell->No);
		printf("CellXDim %f YcellDim %f\n",
		       p2p->XCellDim,p2p->YCellDim);
		printf("%f %f\n",TempCell->X/p2p->XCellDim,
		       TempCell->Y/p2p->YCellDim);
	      }
	      TempCell = TempCell->next;
	    }
	  }
	}



	//P2P->Print(); 
	exit(0);}
    };
    if(*DDD!=TempCell->next->No) 
      printf("multi cell:  !!DDD:%ld TempUpper->No:%ld ",*DDD,TempCell->No);
    if(TempCell->next->InLower[0]==1) p2p->UpperLiving[0] --;
    if(TempCell->next->InLower[1]==1) p2p->UpperLiving[1] --;
    TempCell2 = TempCell->next;
    TempCell->next = TempCell2->next;
    if(pp->MarkLength>0) free(TempCell2->Mark);
    free(TempCell2);
    /* Common stuff */
    //KillCounter ++;
    *DDD = *DDD - 1;
  }
}


long int Sampler::BirthDeath(long int TimeStep,
			     struct Point *headLiving,
			     struct Point *headDeleted,
			     struct Point3 *headTransition,
			     Point2Pattern *p2p,
			     PointProcess *pp,
			     gsl_rng *r){
  float f1, f2, f3, f4;
  long int i,n;
  double xtemp,ytemp, *marktemp;
  char InWindow, Success;
  struct Point *TempPoint, *TempPoint2;
  struct Point3 *TempTransition;
  
  // Following line is there to 'shut up' gcc -Wall
  TempPoint2 = headLiving;

  f1 = LivingPoints; f2 = pp->TotalBirthRate; f3 = f2/(f1+f2);
  n = 0;
  Success = 0;

  //printf("LivingPoints: %d TotalBirthRate %f GeneratedPoints %d\n",
  // LivingPoints,PP->TotalBirthRate,GeneratedPoints);
  
  /* Birth */
  while(Success==0){
    f4 = gsl_rng_uniform(r);
    if(f4<f3){
      //printf("Ping 1 (BD)\n");
      pp->NewEvent(&xtemp, &ytemp, &marktemp, &InWindow, r);
      //printf("Ping 2 (BD)\n");
      if(InWindow==1){
	Success = 1;
	TempTransition = (struct Point3 *) malloc(sizeof *TempTransition);
	if(TempTransition == NULL){
	  printf("Cannot allocate memory for TempTransition in Sampler::BirthDeath (Birth)\n");
	  exit(1);
	}
	//printf("Ping 3 (BD)\n");
	TempTransition->Case = 0;
	LivingPoints ++;
	GeneratedPoints ++;
	TempPoint = (struct Point *) malloc(sizeof *TempPoint);
	if(TempPoint == NULL){
	  printf("Cannot allocate memory for TempPoint in Sampler::BirthDeath\n");
	  exit(1);
	}
	TempPoint->X = xtemp;
	TempPoint->Y = ytemp;
	TempPoint->Mark = marktemp;
	TempPoint->No = GeneratedPoints;
	TempPoint->R = gsl_rng_uniform(r);
	TempPoint->next = headLiving->next;
	headLiving->next = TempPoint;
	NoP ++;
	TempTransition->XCell = (int) floor((TempPoint->X-p2p->Xmin)/
					    p2p->XCellDim);
	TempTransition->YCell = (int) floor((TempPoint->Y-p2p->Ymin)/
					    p2p->YCellDim);
	if(TempTransition->XCell>p2p->MaxXCell){
	  TempTransition->XCell=p2p->MaxXCell;
	  // printf("X AvAvAv! %f\n",TempPoint->X);
	}
	if(TempTransition->YCell>p2p->MaxYCell){ 
	  TempTransition->YCell=p2p->MaxYCell;
	  // printf("Y AvAvAv! %f\n",TempPoint->Y);
	}
      TempTransition->next = headTransition->next;
      headTransition->next = TempTransition;
    }
  }
  /* Death */
  else{
    Success = 1;
    TempTransition = (struct Point3 *) malloc(sizeof *TempTransition);
    if(TempTransition == NULL){
      printf("Cannot allocate memory for TempTransition in Sampler::BirthDeath (Death)\n");
    }
    TempTransition->Case = 1;
    f1 = LivingPoints; f2 = f1*gsl_rng_uniform(r)+1.0;
    n = int(f2);
    if(n>LivingPoints){
      // printf("n=%ld!=%ld AvAvAv!\n",n,LivingPoints);
      n=LivingPoints;
    }
    TempPoint = headLiving;
    for(i=1; i<=n; i++){ 
      TempPoint2 = TempPoint;
      TempPoint = TempPoint->next;
      }
    TempPoint2->next = TempPoint->next;
    
    TempPoint->next = headDeleted->next;  
    headDeleted->next = TempPoint;

    LivingPoints --;
    NoP --;
    TempTransition->next = headTransition->next;
    headTransition->next = TempTransition;
  }
  }
  return(n);
}

void Sampler::Sim(Point2Pattern *p2p, PointProcess *pp, gsl_rng *r) {
  long int **summary;
  summary = (long int **) malloc(sizeof(long int*));
  *summary = (long int*) malloc(sizeof(long int));
  Sim(p2p,pp,summary,r);
  free(*summary);
}

void Sampler::Sim(Point2Pattern *p2p, PointProcess *pp, 
			     long int **Summary, gsl_rng *r) {

  //P2P = p2p;
  long int StartTime, EndTime, TimeStep, D0Time, D0Living;
  long int XCell, YCell, DDD, i, NoP, iMark;
  long int TimeSinceUpdateOfInteractionRange;
  FILE *out;
 
  //printf("Sim...\n");

  if((p2p->Xmax < pp->Xmax) || (p2p->Xmin > pp->Xmin) ||
     (p2p->Ymax < pp->Ymax) || (p2p->Ymin > pp->Ymin)){
    printf("point2pattern does not cover pointprocess!\n");
    printf("point2pattern: x range: %lf %lf  y range: %lf %lf\n",
	   p2p->Xmin,p2p->Xmax,p2p->Ymin,p2p->Ymax);
    printf("pointprocess:  x range: %lf %lf  y range: %lf %lf\n",
	   pp->Xmin,pp->Xmax,pp->Ymin,pp->Ymax);
    exit(0);
  }

  p2p->Empty();

  p2p->MarkLength = pp->MarkLength;
  /* Initialising linked listed for backward simulation */
  struct Point *headDeleted, *headLiving, *dummyDeleted, *dummyLiving;
  struct Point *TempPoint;
  headLiving = (struct Point *) malloc(sizeof *headLiving);
  if(headLiving == NULL){
    printf("Cannot allocate memory for headLiving in Sampler::Sim\n");
    exit(1);
  }
  dummyLiving = (struct Point *) malloc(sizeof *dummyLiving);  
  if(dummyLiving == NULL ){
    printf("Cannot allocate memory for dummyLiving in Sampler::Sim\n");
    exit(1);
  }
  headLiving->next = dummyLiving; dummyLiving->next = dummyLiving;

  headDeleted = (struct Point *) malloc(sizeof *headDeleted);
  if(headDeleted == NULL ){
    printf("Cannot allocate memory for headDeleted in Sampler::Sim\n");
    exit(1);
  }
  dummyDeleted = (struct Point *) malloc(sizeof *dummyDeleted);  
  if(dummyDeleted == NULL ){
    printf("Cannot allocate memory for dummyDeleted in Sampler::Sim\n");
    exit(1);
  }
  headDeleted->next = dummyDeleted; dummyDeleted->next = dummyDeleted;

  struct Point2 *TempCell2;

  struct Point3 *headTransition, *dummyTransition;
  headTransition = (struct Point3 *) malloc(sizeof *headTransition);
  if(headTransition == NULL){
    printf("Cannot allocate memory for headTransition in Sampler::Sim\n");
    exit(1);
  }
  dummyTransition = (struct Point3 *) malloc(sizeof *dummyTransition);
  if(dummyTransition == NULL){
    printf("Cannot allocate memory for dummyTransition in Sampler::Sim\n");
    exit(1);
  }
  headTransition->next = dummyTransition; 
  dummyTransition->next = dummyTransition;

  GenerateDominatingPoisson(headLiving, &GeneratedPoints,
			    &LivingPoints,
			    &NoP,
			    pp,
			    r); 
  *Summary = (long int *) realloc(*Summary,2*sizeof(long int));
  **Summary = GeneratedPoints;
  *(*Summary +1) = 0;
  //printf("D0 generated\n");
  StartTime=1;
  EndTime=1;

  TimeStep = 0; D0Time = 0;
  D0Living = LivingPoints;//GeneratedPoints;

  long int tmp, D0;
  
  do{
    tmp=BirthDeath(TimeStep, headLiving,
		   headDeleted,
		   headTransition,p2p,pp,r);
    if(tmp>0){ 
      if(tmp>(LivingPoints+1-D0Living)){
	D0Living --;
      }
    }
    //printf("D0Living: %d\n",D0Living);
    D0Time++;
  }while(D0Living>0);
  tmp=BirthDeath(TimeStep, headLiving,
		 headDeleted,
		 headTransition, p2p, pp, r); 
  StartTime=1; EndTime=D0Time+1; D0 = 0;
  //printf("D0Time: %ld\n",D0Time);
  do{	 
    if(D0==1){
      for(TimeStep=StartTime;TimeStep<=EndTime;TimeStep ++){
	tmp=BirthDeath(TimeStep, headLiving,
		       headDeleted,
		       headTransition,p2p, pp,r);      
      }
    }
    D0 = 1;
    p2p->Empty();

    *(*Summary+1) = *(*Summary+1) + 1;
    *Summary = (long int *) realloc(*Summary, 
				    2*(1+*(*Summary+1))*sizeof(long int));
    *((*Summary)+*(*Summary+1)*2) = GeneratedPoints;
    *((*Summary)+*(*Summary+1)*2+1) = LivingPoints;
    

    //printf("Yo!\n");

    /*
    headUpper->next = dummyUpper; dummyUpper->next = dummyUpper;
    for(XCell=0;XCell<=P2P->MaxXCell;XCell++){
      for(YCell=0;YCell<=P2P->MaxYCell;YCell++){
	headUpperCell[XCell][YCell]->next=dummyUpper;
      }
    }
    */
    
    p2p->UpperLiving[0] = LivingPoints;
    p2p->UpperLiving[1] = 0;

    /*
    out = fopen("LivingPoints.dat","a");
    fprintf(out,"%ld\n",LivingPoints);
    fclose(out);
    */

    // Copy SimPoints to Points
    //printf("Copy: MarkLength: %d\n",pp->MarkLength);
    p2p->NoP = 0;
    i=0;
    TempPoint = headLiving->next;
    while(TempPoint!=TempPoint->next){
      i++;
      TempCell2 = (struct Point2 *) malloc(sizeof *TempCell2);
      if(TempCell2 == NULL){
	printf("Cannot allocate memory for TempCell2 in Sampler::Sim\n");
	exit(1);
      }
      TempCell2->No = TempPoint->No;
      TempCell2->X = TempPoint->X;
      TempCell2->Y = TempPoint->Y;
      if(pp->MarkLength>0){
	TempCell2->Mark = (double *) calloc(pp->MarkLength,sizeof(double));
	if(TempCell2->Mark==NULL){
	  printf("Cannot allocate memory for Mark in Sampler::Sim\n");
	  exit(1);
	}
	for(iMark=0;iMark<pp->MarkLength;iMark++){ 
	  *(TempCell2->Mark + iMark) = *(TempPoint->Mark + iMark);	
	}
      }

      TempCell2->InLower[0] = 1;
      TempCell2->InLower[1] = 0;
      XCell = (int) floor((TempPoint->X-p2p->Xmin)/p2p->XCellDim);
      if(XCell>p2p->MaxXCell) XCell = p2p->MaxXCell;
      YCell = (int) floor((TempPoint->Y-p2p->Ymin)/p2p->YCellDim);
      if(YCell>p2p->MaxYCell) YCell = p2p->MaxYCell;     
      TempCell2->next = p2p->headCell[XCell][YCell]->next;
      
      p2p->headCell[XCell][YCell]->next = TempCell2;
      
      TempPoint = TempPoint->next;
      
    }

  struct Point3 *TempTransition;
  struct Point *Proposal;

  TempTransition = headTransition->next;
  Proposal = headDeleted->next;
  DDD = GeneratedPoints;
  
  pp->ResetInteractionRange();
  TimeSinceUpdateOfInteractionRange = 0;
  
  for(TimeStep=EndTime;TimeStep>=1;TimeStep--){
    Forward(TimeStep,TempTransition->Case,
	    TempTransition->XCell,TempTransition->YCell,
	    Proposal,&DDD, p2p, pp);
    if(TempTransition->Case == 1) Proposal = Proposal ->next;
    TempTransition = TempTransition->next;
    TimeSinceUpdateOfInteractionRange++;
    if(TimeSinceUpdateOfInteractionRange == 
       pp->InteractionRangeUpdateFrequency){ 
      pp->UpdateInteractionRange(p2p);
      TimeSinceUpdateOfInteractionRange = 0;
    }
  }

  /* Doubling strategy used!*/
  StartTime = EndTime+1;
  EndTime=EndTime*2;

  /*
  i=0;
  for(XCell=0;XCell<=p2p->MaxXCell;XCell++){
    for(YCell=0;YCell<=p2p->MaxYCell;YCell++){
      TempCell2 = p2p->headCell[XCell][YCell]->next;
      while(TempCell2->next!=TempCell2){
	i++;
	printf("%d: %f %f %d %d\n",i,TempCell2->X,TempCell2->Y,
	       TempCell2->InLower[0],TempCell2->InLower[1]);
	TempCell2 = TempCell2->next;
      }
      //headUpperCell[XCell][YCell]->next=dummyUpper;
    }
  }
  */
  // if(StartTime>10000)
    // printf("LL: %ld UL: %ld StartTime: %ld\n",
	  //  p2p->UpperLiving[0],p2p->UpperLiving[1],StartTime);

  }while(p2p->UpperLiving[0]!=p2p->UpperLiving[1]);
  p2p->Clean();
  i=0;
  struct Point *TempPoint2;
  TempPoint = headLiving->next;
  while(TempPoint!=TempPoint->next){
    TempPoint2 = TempPoint->next;
    i++;
    if(pp->MarkLength>0) free(TempPoint->Mark);
    free(TempPoint);
    TempPoint = TempPoint2;
  }
  free(TempPoint);
  free(headLiving);

  i = 0;
  TempPoint = headDeleted->next;
  while(TempPoint!=TempPoint->next){
    TempPoint2 = TempPoint->next;
    i++;
    if(pp->MarkLength>0) free(TempPoint->Mark);
    free(TempPoint);
    TempPoint = TempPoint2;
  }
  free(TempPoint);
  free(headDeleted);

  struct Point3 *TempTransition,*TempTransition2;

  i = 0;
  TempTransition = headTransition;
  TempTransition2 = headTransition->next;
  while(TempTransition!=TempTransition->next){
    i++;
    free(TempTransition);
    TempTransition = TempTransition2;
    TempTransition2 = TempTransition2->next;
  }
  free(TempTransition);
  //printf("%d ST: %d ET: %d\n",i,StartTime,EndTime);
  //scanf("%f",&f1);
};


long int Sampler::ForwardBD(double Time, Point2Pattern *p2p, PointProcess *pp,
	       gsl_rng *r){
  double RunTime, U, x,y,*mark,lnCIR;
  long int NoP, N, Transitions;
  char InWindow;
  Point2 Point;
  
  Transitions = 0;

  // Making sure the MarkLength matches for point pattern and
  // point process
  if(p2p->MarkLength != pp->MarkLength){
    printf("Mismatch in MarkLength!\n");
    p2p->Empty();
    p2p->MarkLength = pp->MarkLength;
  }

  RunTime = 0.0;
  NoP = p2p->Count();
  RunTime += gsl_ran_exponential(r,1.0/((double) NoP + pp->TotalBirthRate));
  while(RunTime < Time){
    U = gsl_rng_uniform(r);
    if(U < (((double)NoP)/((double) NoP + pp->TotalBirthRate))){
      // Death
      //printf("Death\n");
      N = gsl_rng_uniform_int(r,NoP)+1;
      p2p->Remove(N);
      NoP--;
      Transitions++;
    }else{
      // Birth
      //printf("Birth\n");
      //if(pp->MarkLength>0) 
      //mark = (double *) calloc(pp->MarkLength,sizeof(double));
      pp->NewEvent(&x,&y,&mark,&InWindow,r);
      if(InWindow==1){
	Point.X = x;
	Point.Y = y;
	Point.Mark = mark;
	pp->Beta(&Point); 
	lnCIR = pp->lnCondIntRatio(&Point,p2p);
	U = gsl_rng_uniform(r);
	if(log(U)<lnCIR){
	  p2p->Add(&Point);
	  NoP++;
	  Transitions++;
	}	
	//p2p->DumpToFile("test.dat");
      }	
      
      //if(pp->MarkLength>0) free(mark);     
    }
    RunTime += gsl_ran_exponential(r,1.0/((double) NoP + pp->TotalBirthRate));
  }
  return(Transitions);
}


void Sampler::ForwardMH(long int Steps, 
			Point2Pattern *p2p, PointProcess *pp,
			gsl_rng *r){
  long int Step, NoP, N, i;
  double ProposalX, ProposalY, *ProposalMark, logHR;
  char InWindow;
  struct Point2 TempPoint, *TempPointP;
  
  if(pp->MarkLength>0){ 
    TempPoint.Mark = (double *) calloc(pp->MarkLength,sizeof(double));
    if(TempPoint.Mark  == NULL){
      printf("Cannot allocate memory for TempPoint.Mark in Sampler::ForwardMH (Birth)\n");
      exit(1);
    }
  }
  NoP = p2p->Count();

  for(Step=0; Step<Steps;Step++){
    //printf("Step: %ld NoP: %ld\n",Step,NoP);
    if(gsl_rng_uniform(r)<0.5){
      // Birth proposal
      //printf("Birth proposal\n");
      pp->NewEvent(&ProposalX,&ProposalY,&ProposalMark,&InWindow,r);
      if(InWindow==1){
	TempPoint.X = ProposalX;
	TempPoint.Y = ProposalY;
	if(pp->MarkLength>0){
	  for(i=0;i<pp->MarkLength;i++) 
	    *(TempPoint.Mark + i) = *(ProposalMark + i);
	}
	pp->Beta(&TempPoint);
	//printf("ping! %lf %lf %lf\n",TempPoint.X,TempPoint.Y,*TempPoint.Mark);
	logHR = pp->lnCondInt(&TempPoint,p2p) - 
	  pp->logNewEventDensity(TempPoint.X,TempPoint.Y,TempPoint.Mark) -
	  log(NoP+1);
	//printf("pong!\n");
	/*
	printf("logHR: %lf %lf %lf\n",
	       logHR,pp->lnCondInt(&TempPoint,p2p),
	       pp->logNewEventDensity(TempPoint.X,TempPoint.Y,TempPoint.Mark));
	*/
	if(log(gsl_rng_uniform(r))<logHR){
	  p2p->Add(&TempPoint);
	  NoP++;
	}
      }
      free(ProposalMark);
    }else{
      // Death proposal
      //printf("Death proposal\n");
      if(NoP>0){
	N = gsl_rng_uniform_int(r,NoP)+1;
	TempPointP = p2p->Extract(N);
	if(TempPointP == TempPointP->next){ printf("URK!\n"); exit(0);}
	//if(TempPointP->next == TempPointP->next->next){ printf("KRU!\n"); exit(0);}
	logHR = -pp->lnCondInt(TempPointP,p2p)+ 
	  pp->logNewEventDensity(TempPointP->X,TempPointP->Y,TempPointP->Mark) +
	  log(NoP);
	if(log(gsl_rng_uniform(r))<logHR){
	  p2p->Remove(N);
	  NoP--;
	}      
      }
    }
  }
}
