#include <math.h>
#include <limits>
#include "point2pattern.h"
#include "pointprocess.h"
#include "strauss.h"

/****************
 StraussProcess
****************/

void StraussProcess::SetParameters(double b, double g, double Ri){
  beta = b; gamma = g; R = Ri; 
  InteractionRange = R;
  TotalBirthRate = beta*(this->Xmax-this->Xmin)*(this->Ymax-this->Ymin);
}

void StraussProcess::Print(){
  printf("StraussProcess with beta = %f, R = %f and gamma = %f\n",
	 beta,R,gamma);

}

StraussProcess::StraussProcess(double xmin, double xmax, 
			      double ymin, double ymax, 
			      double b, double g, double Ri) :
  PointProcess(xmin, xmax, ymin, ymax){
  beta = b; gamma = g; R = Ri; 
  InteractionRange = R;
  TotalBirthRate = beta*(xmax-xmin)*(ymax-ymin);
  }  

/*
double StraussProcess::Interaction(double r)
{
  double rtn;
  rtn = 1;
  if(r<R) rtn = gamma;
  return(rtn);
}
*/
double StraussProcess::Interaction(struct Point2 *Point1, struct Point2 *Point2)
{
  double rtn;
  if((sqrt(pow(Point1->X-Point2->X,2)+
	   pow(Point1->Y-Point2->Y,2)))<R){ rtn = gamma;}else{rtn = 1;}
  return(rtn);
}



void StraussProcess::NewEvent(double *x, double *y, 
			      double **mark, char *InWindow, gsl_rng *r)
{
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = gsl_rng_uniform(r)*Xdim+Xmin;
  *y = gsl_rng_uniform(r)*Ydim+Ymin;
  *InWindow = 1;
}


void StraussProcess::CalcBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
  long int i,j,k;
  k=0;
  printf("\nStrauss CalcBeta... %ld %ld\n",xsidepomm,ysidepomm);
  for(i=0; i<xsidepomm; i++){
    for(j=0; j<ysidepomm; j++){
      *(betapomm + i*ysidepomm + j) = this->beta;
      k++;
    }
  } 
}

void StraussProcess::CheckBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
  int dummyint;
  long int i,j,k;
  double d1;
  k=0;
  printf("\nStrauss CalcBeta... %ld %ld\n",xsidepomm,ysidepomm);
  for(i=0; i<xsidepomm; i++){
    for(j=0; j<ysidepomm; j++){
      if((fabs(*(betapomm + i*ysidepomm + j)- beta)>0.001) && (k==0)){
	printf("%f %f %f %ld %ld\n",fabs(*(betapomm + i*ysidepomm + j)- beta),
	       *(betapomm + i*ysidepomm + j),beta,i,j);
	k++;
	dummyint = scanf("%lf",&d1);
      }
    }
  } 
}

double StraussProcess::lnCondInt(struct Point2 *TempCell, 
				 Point2Pattern *p2p){
  double f1;
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx,k;
  double lnCI,dst;
  struct Point2 *TempCell2;
  
  xc = (int) floor((TempCell->X-p2p->Xmin)/p2p->XCellDim);
  if(xc>p2p->MaxXCell) xc = p2p->MaxXCell;
  yc = (int) floor((TempCell->Y-p2p->Ymin)/p2p->YCellDim); 
  if(yc>p2p->MaxYCell) yc = p2p->MaxYCell;
  
  //printf("x: %f y: %f  xc:%d  yc: %d\n",TempCell->X,TempCell->Y,
  //xc,yc);

  //dx = (Xmax-Xmin)/(double(p2p->MaxXCell+1));
  //dy = (Ymax-Ymin)/(double(p2p->MaxYCell+1));
  rx = (int) floor(this->InteractionRange/p2p->XCellDim+1.0);
  ry = (int) floor(this->InteractionRange/p2p->YCellDim+1.0);
  
  lnCI = log(TempCell->Beta);

  k = 0;
  
  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
  if((yc-ry)>=0) fy=yc-ry; else fy = 0;

  //printf("MCI! %ld %ld %ld %ld\n",fx,tx,fy,ty);

  for(xco = fx; xco <= tx; xco++){
    for(yco = fy; yco <= ty; yco++){
      //printf("%d %d\n",xco,yco);
      TempCell2 = p2p->headCell[xco][yco]->next;
      while(TempCell2!=TempCell2->next){
	if(TempCell2 != TempCell){
	  k++;
	  //dst = sqrt(pow(TempCell->X-TempCell2->X,2)+
	  //     pow(TempCell->Y-TempCell2->Y,2));
	  //printf("%f %f\n",dst,log(Interaction(dst)));
	  //lnCI += log(Interaction(dst));
	  lnCI += log(Interaction(TempCell,TempCell2));
	}
	TempCell2 = TempCell2->next; 
	
      }
    }
  }
  return(lnCI);
}

double StraussProcess::lnCondIntRatio(struct Point2 *TempCell, 
				 Point2Pattern *p2p){
  double f1;
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx,k;
  double lnCI,dst;
  struct Point2 *TempCell2;
  
  xc = (int) floor((TempCell->X-p2p->Xmin)/p2p->XCellDim);
  if(xc>p2p->MaxXCell) xc = p2p->MaxXCell;
  yc = (int) floor((TempCell->Y-p2p->Ymin)/p2p->YCellDim); 
  if(yc>p2p->MaxYCell) yc = p2p->MaxYCell;
  
  //printf("x: %f y: %f  xc:%d  yc: %d\n",TempCell->X,TempCell->Y,
  //xc,yc);

  //dx = (Xmax-Xmin)/(double(p2p->MaxXCell+1));
  //dy = (Ymax-Ymin)/(double(p2p->MaxYCell+1));
  rx = (int) floor(this->InteractionRange/p2p->XCellDim+1.0);
  ry = (int) floor(this->InteractionRange/p2p->YCellDim+1.0);
  
  //printf("InteractionRange: %lf\n",this->InteractionRange);

  lnCI = 0;

  k = 0;
  
  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
  if((yc-ry)>=0) fy=yc-ry; else fy = 0;

  //printf("MCI! %ld %ld %ld %ld\n",fx,tx,fy,ty);

  for(xco = fx; xco <= tx; xco++){
    for(yco = fy; yco <= ty; yco++){
      //printf("%d %d\n",xco,yco);
      TempCell2 = p2p->headCell[xco][yco]->next;
      while(TempCell2!=TempCell2->next){
	if(TempCell2 != TempCell){
	  k++;
	  //dst = sqrt(pow(TempCell->X-TempCell2->X,2)+
	  //     pow(TempCell->Y-TempCell2->Y,2));
	  //printf("%f %f\n",dst,log(Interaction(dst)));
	  //lnCI += log(Interaction(dst));
	  lnCI += log(Interaction(TempCell,TempCell2));
	}
	TempCell2 = TempCell2->next; 
	
      }
    }
  }
  return(lnCI);
}


void StraussProcess::Beta(struct Point2 *TempCell){
  TempCell->Beta = beta;
}

void StraussProcess::CalcBeta(Point2Pattern *p2p){
  long int xco,yco;
  struct Point2 *TempMother;

  for(xco = 0; xco <= p2p->MaxXCell; xco++){
    for(yco = 0; yco <= p2p->MaxYCell; yco++){
      TempMother = p2p->headCell[xco][yco]->next;
      while(TempMother!=TempMother->next){
	TempMother->Beta = this->beta;
	TempMother = TempMother->next;
      }
    }
  }
}

long int StraussProcess::Pairs(Point2Pattern *p2p){
  //long int xc,yc,rx,ry,tx,ty,fx,xco,yco,i;
  struct Point2 *Point, *Point2;
  long int Pairs;

  Pairs = 0;
  Point = p2p->First();
  while(Point!=NULL){
    Point2 = p2p->Next(Point,R,Point);
    while(Point2!=NULL){
      if(sqrt((Point->X-Point2->X)*(Point->X-Point2->X)+
	      (Point->Y-Point2->Y)*(Point->Y-Point2->Y))<R)
	Pairs++;
      Point2 = p2p->Next(Point2,R,Point);
    }
    Point = p2p->Next(Point);
  }

  /*
  i=0;
  // Same cell
  for(xc = 0; xc <= p2p->MaxXCell; xc++){
    for(yc = 0; yc <= p2p->MaxYCell; yc++){
      Point1 = p2p->headCell[xc][yc]->next;
      while(Point1 != Point1->next){
	Point2 = Point1->next;
	while(Point2!=Point2->next){
	  if( sqrt((Point1->X-Point2->X)*(Point1->X-Point2->X) + 
		   (Point1->Y-Point2->Y)*(Point1->Y-Point2->Y))
	      < R ) i++;
	  Point2 = Point2->next;
	}
	Point1 = Point1->next;
      }
    }
  }
  // Different cells
  for(xc = 0; xc <= p2p->MaxXCell; xc++){
    for(yc = 0; yc <= p2p->MaxYCell; yc++){

      rx = (int) floor(R/p2p->XCellDim+1.0);
      ry = (int) floor(R/p2p->YCellDim+1.0);
      if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
      if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
      if((xc-rx)>=0) fx=xc-rx; else fx = 0;

      Point1 = p2p->headCell[xc][yc]->next;
      while(Point1 != Point1->next){
	for(xco = fx; xco <= tx; xco++){
	  for(yco = yc; yco <= ty; yco++){
	    if(!((xco<=xc) && (yco==yc))){
	      Point2 = p2p->headCell[xco][yco]->next;
	      while(Point2!=Point2->next){
		if( sqrt((Point1->X-Point2->X)*(Point1->X-Point2->X) + 
			 (Point1->Y-Point2->Y)*(Point1->Y-Point2->Y))
		    < R ) i++;		
		Point2 = Point2->next;
	      }
	    }
	  }
	}
	Point1 = Point1->next;       
      }
    }
  } 
  */
  return(Pairs);
}



/**********************
  StraussProcessInhom
**********************/

void StraussProcessInhom::SummaryToFile(const char *FileName){
  FILE *out;

  out = fopen(FileName,"w");

  fprintf(out,"Xmin:\t%f\n",Xmin);
  fprintf(out,"Xmax:\t%f\n",Xmax);
  fprintf(out,"Ymin:\t%f\n",Ymin);
  fprintf(out,"Ymax:\t%f\n",Ymax);

  fprintf(out,"InteractionRange:\t%f\n",InteractionRange);

  fprintf(out,"beta:\t%f\n",beta);
  fprintf(out,"R:\t%f\n",R);
  fprintf(out,"gamma:\t%f\n",gamma);
  
  fprintf(out,"intercept:\t%f\n",intercept);
  fprintf(out,"slopex:\t%f\n",slopex);
  fprintf(out,"slopex2:\t%f\n",slopex2);
  fprintf(out,"slopey:\t%f\n",slopey);
  fprintf(out,"slopey2:\t%f\n",slopey2);
  fprintf(out,"slopexy:\t%f\n",slopexy);

  fprintf(out,"betaMax:\t%f\n",betaMax);


  fclose(out);
};

void StraussProcessInhom::SetParameters(double b, double g, double Ri,
					double ntrcpt, 
					double slpx, double slpx2,
					double slpy, double slpy2, 
					double slpxy){
  //change?
  
  beta = b; gamma = g; R = Ri; 
  InteractionRange = R;
  TotalBirthRate = beta*(this->Xmax-this->Xmin)*(this->Ymax-this->Ymin);
}

void StraussProcessInhom::SetParameters(double b, double g, double Ri){
  beta = b; gamma = g; R = Ri; 
  InteractionRange = R;
  TotalBirthRate = beta*(this->Xmax-this->Xmin)*(this->Ymax-this->Ymin);
}


void StraussProcessInhom::Print(){
  //change?
  printf("StraussProcessInhom with beta = %f, R = %f and gamma = %f \n",
	 beta,R,gamma);
  printf("slopex = %f slopex2 = %f slopey = %f slopey2 = %f\n",
	 slopex,slopex2,slopey,slopey2);
  printf("Xmin = %f Xmax = %f Ymin = %f Ymax = %f\n",
	 Xmin,Xmax,Ymin,Ymax);


}

StraussProcessInhom::StraussProcessInhom(
		 double xmin, double xmax, double ymin, double ymax, 
		 double b, double g, double Ri,
		 double ntrcpt, double slpx, double slpx2,
		 double slpy, double slpy2, double slpxy) :
  PointProcess(xmin, xmax, ymin, ymax){
  double xmode, ymode, tmp;
  char modeinside, fail;


  beta = b; gamma = g; R = Ri; 
  InteractionRange = R;

  intercept = ntrcpt;
  slopex = slpx; slopex2 = slpx2;
  slopey = slpy; slopey2 = slpy2;
  slopexy = slpxy;

  printf("slopex: %f slopex2: %f\n",slopex,slopex2);
  printf("slopey: %f slopey2: %f\n",slopey,slopey2);
  printf("slopexy: %f\n",slopexy);
  betaMin = std::numeric_limits<double>::max();
  betaMax = -betaMin;
  printf("betaMax: %f\n",betaMax);
  printf("betaMin: %f\n",betaMin);

  modeinside = 0;
  fail = 0;
  tmp = 4*slopex2*slopey2-slopexy*slopexy;
  printf("tmp: %f\n",tmp);
  if(tmp!=0){
    xmode = (slopey*slopexy - 2*slopey2*slopex)/tmp;
    ymode = (slopex*slopexy - 2*slopex2*slopey)/tmp;
    if(xmode>=this->Xmin & xmode<=this->Xmax & 
       ymode>=this->Ymin & ymode <= this->Ymax){ 
      modeinside = 1;
      tmp = Beta(xmode,ymode);
      printf("tmp: %f\n",tmp);
      if(tmp < betaMin) betaMin = tmp;
      if(tmp > betaMax) betaMax = tmp;      
      if(tmp<0) fail = 1;
    }
  }

  printf("betaMax: %f\n",betaMax);
  printf("betaMin: %f\n",betaMin);


  if(modeinside==0 & fail==0){
    // Bottom
    tmp = Beta(-(1/2)*(slopex+slopexy*this->Ymin)/slopex2,this->Ymin);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      
    // Top
    tmp = Beta(-(1/2)*(slopex+slopexy*this->Ymax)/slopex2,this->Ymax);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      
    // Bottom
    tmp = Beta(this->Xmin,-(1/2)*(slopey+slopexy*this->Xmin)/slopey2);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      
    // Top
    tmp = Beta(this->Xmax,-(1/2)*(slopey+slopexy*this->Xmax)/slopey2);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      

    // Top left
    tmp = Beta(this->Xmin,this->Ymax);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      
    // Top right
    tmp = Beta(this->Xmax,this->Ymax);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      
    // Bottom right
    tmp = Beta(this->Xmax,this->Ymin);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      
    // Bottom left
    tmp = Beta(this->Xmin,this->Ymin);
    printf("tmp: %f\n",tmp);
    if(tmp < betaMin) betaMin = tmp;
    if(tmp > betaMax) betaMax = tmp;      

  }
  //if(fail==1){
  //  printf("Invalid parameters StraussProcessInhom!\n");
  //  exit(0);
  //}

  printf("betaMax: %f\n",betaMax);
  printf("betaMin: %f\n",betaMin);

  TotalBirthRate = (xmax-xmin)*(ymax-ymin)*betaMax;
  printf("TotalBirthRate: %f\n",TotalBirthRate);
}  

/*
double StraussProcessInhom::Interaction(double r)
{
  //change?
  double rtn;
  rtn = 1;
  if(r<R) rtn = gamma;
  return(rtn);
}
*/

double StraussProcessInhom::Interaction(struct Point2 *Point1, struct Point2 *Point2)
{
  //change?
  double rtn;
  if((sqrt(pow(Point1->X-Point2->X,2)+
	   pow(Point1->Y-Point2->Y,2)))<R){ rtn = gamma;}else{rtn = 1;}
  return(rtn);
}


void StraussProcessInhom::NewEvent(double *x, double *y, double **mark,
				   char *InWindow, gsl_rng *r)
{
  //change?
  double Xdim, Ydim;
  Xdim = Xmax-Xmin;
  Ydim = Ymax-Ymin;
  *x = gsl_rng_uniform(r)*Xdim+Xmin;
  *y = gsl_rng_uniform(r)*Ydim+Ymin;
  if(gsl_rng_uniform(r)<= (Beta(*x,*y)/betaMax)){ 
    *InWindow = 1;
  }else{
    *InWindow = 0;
  }
}


void StraussProcessInhom::CalcBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
  //change?
  long int i,j,k;
  k=0;
  printf("\nStrauss CalcBeta... %ld %ld\n",xsidepomm,ysidepomm);
  for(i=0; i<xsidepomm; i++){
    for(j=0; j<ysidepomm; j++){
      *(betapomm + i*ysidepomm + j) = this->beta;
      k++;
    }
  } 
}

void StraussProcessInhom::CheckBeta(long int xsidepomm, long int ysidepomm, 
		   double *betapomm){ 
  //change?
  int dummyint;
  long int i,j,k;
  double d1;
  k=0;
  printf("\nStrauss CalcBeta... %ld %ld\n",xsidepomm,ysidepomm);
  for(i=0; i<xsidepomm; i++){
    for(j=0; j<ysidepomm; j++){
      if((fabs(*(betapomm + i*ysidepomm + j)- beta)>0.001) && (k==0)){
	printf("%f %f %f %ld %ld\n",fabs(*(betapomm + i*ysidepomm + j)- beta),
	       *(betapomm + i*ysidepomm + j),beta,i,j);
	k++;
	dummyint = scanf("%lf",&d1);
      }
    }
  } 
}

double StraussProcessInhom::lnCondInt(struct Point2 *TempCell, 
				 Point2Pattern *p2p){
  //change?
  double f1;
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx,k;
  double lnCI,dst;
  struct Point2 *TempCell2;
  
  f1 = (TempCell->X-p2p->Xmin)/p2p->XCellDim;  xc = int(f1);
  if(xc>p2p->MaxXCell) xc = p2p->MaxXCell;
  f1 = (TempCell->Y-p2p->Ymin)/p2p->YCellDim;  yc = int(f1);
  if(yc>p2p->MaxYCell) yc = p2p->MaxYCell;

  //printf("x: %f y: %f  xc:%d  yc: %d\n",TempCell->X,TempCell->Y,
  // xc,yc);

  rx = int(this->InteractionRange/p2p->XCellDim+1.0);
  ry = int(this->InteractionRange/p2p->YCellDim+1.0);

  
  
  lnCI = log(TempCell->Beta);

  k = 0;
  
  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
  if((yc-ry)>=0) fy=yc-ry; else fy = 0;

  //printf("MCI! %ld %ld %ld %ld\n",fx,tx,fy,ty);

  for(xco = fx; xco <= tx; xco++){
    for(yco = fy; yco <= ty; yco++){
      TempCell2 = p2p->headCell[xco][yco]->next;
      while(TempCell2!=TempCell2->next){
	if(TempCell2 != TempCell){
	  k++;
	  lnCI += log(Interaction(TempCell,TempCell2));
	  
	}
	TempCell2 = TempCell2->next; 
	
      }
    }
  }
  return(lnCI);
}

double StraussProcessInhom::lnCondIntRatio(struct Point2 *TempCell, 
				 Point2Pattern *p2p){
  //change?
  double f1;
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx,k;
  double lnCI,dst;
  struct Point2 *TempCell2;
  
  f1 = (TempCell->X-p2p->Xmin)/p2p->XCellDim;  xc = int(f1);
  if(xc>p2p->MaxXCell) xc = p2p->MaxXCell;
  f1 = (TempCell->Y-p2p->Ymin)/p2p->YCellDim;  yc = int(f1);
  if(yc>p2p->MaxYCell) yc = p2p->MaxYCell;

  //printf("x: %f y: %f  xc:%d  yc: %d\n",TempCell->X,TempCell->Y,
  // xc,yc);

  rx = int(this->InteractionRange/p2p->XCellDim+1.0);
  ry = int(this->InteractionRange/p2p->YCellDim+1.0);

  
  
  lnCI = 0;

  k = 0;
  
  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
  if((yc-ry)>=0) fy=yc-ry; else fy = 0;

  //printf("MCI! %ld %ld %ld %ld\n",fx,tx,fy,ty);

  for(xco = fx; xco <= tx; xco++){
    for(yco = fy; yco <= ty; yco++){
      TempCell2 = p2p->headCell[xco][yco]->next;
      while(TempCell2!=TempCell2->next){
	if(TempCell2 != TempCell){
	  k++;
	  lnCI += log(Interaction(TempCell,TempCell2));
	  
	}
	TempCell2 = TempCell2->next; 
	
      }
    }
  }
  return(lnCI);
}


 
double StraussProcessInhom::Beta(double x, double y){
  return(beta*exp(intercept + slopex*x + slopex2*x*x + 
	     slopey*y + slopey2*y*y + slopexy*x*y));
};

void StraussProcessInhom::Beta(struct Point2 *TempCell){
  //change?
  TempCell->Beta = Beta(TempCell->X,TempCell->Y);
    /*
    beta*exp(intercept + 
			    slopex*TempCell->X + 
			    slopex2*TempCell->X*TempCell->X + 
			    slopey*TempCell->Y + 
			    slopey2*TempCell->Y*TempCell->Y + 
			    slopexy*TempCell->X*TempCell->Y);
    */
}

void StraussProcessInhom::CalcBeta(Point2Pattern *p2p){
  //change?
  long int xco,yco;
  struct Point2 *TempMother;

  for(xco = 0; xco <= p2p->MaxXCell; xco++){
    for(yco = 0; yco <= p2p->MaxYCell; yco++){
      TempMother = p2p->headCell[xco][yco]->next;
      while(TempMother!=TempMother->next){
	TempMother->Beta = beta*exp(intercept + 
			    slopex*TempMother->X + 
			    slopex2*TempMother->X*TempMother->X + 
			    slopey*TempMother->Y + 
			    slopey2*TempMother->Y*TempMother->Y + 
			    slopexy*TempMother->X*TempMother->Y);
	TempMother = TempMother->next;
      }
    }
  }
}
