#include <math.h>
#include <stdio.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "point2pattern.h"
#include "pointprocess.h"

double PointProcess::lnDens(Point2Pattern *p2p){  
  long int xco,yco,xc,yc,fx,tx,fy,ty,ry,rx;
  double lnDens,dst;
  struct Point2 *TempCell, *TempCell2;

  rx = int(InteractionRange/p2p->XCellDim+1.0);
  ry = int(InteractionRange/p2p->YCellDim+1.0);
  
  //printf("1:%f 2:%f 3:%d 4:%d 5:%f 6:%f\n",dx,dy,rx,ry,
  // this->InteractionRange,InteractionRange);
  //printf("mx:%d my:%d\n",p2p->MaxXCell,p2p->MaxYCell);

  lnDens = 0;

  //printf("lnDens: %f (0)\n",lnDens);
  
  for(xc = 0; xc <= p2p->MaxXCell; xc++){
    for(yc = 0; yc <= p2p->MaxYCell; yc++){
      //if(xc==1) printf("%d %d\n",xc,yc);
      TempCell = p2p->headCell[xc][yc]->next;
      while(TempCell != TempCell->next){
	lnDens += log(TempCell->Beta);
	//printf("lnDens: %f (1) %d %d %d %d Beta %f\n",lnDens,xc,yc,
	//       p2p->MaxXCell,p2p->MaxYCell,TempCell->Beta);
	//if(lnDens<(-100000)){printf("%f",lnDens); scanf("%f",&f1);}
	if(InteractionRange>0){
	  if((xc+rx)<=p2p->MaxXCell) tx=xc+rx; else tx = p2p->MaxXCell;
	  if((yc+ry)<=p2p->MaxYCell) ty=yc+ry; else ty = p2p->MaxYCell;
	  if((xc-rx)>=0) fx=xc-rx; else fx = 0;
	  if((yc-ry)>=0) fy=yc-ry; else fy = 0;
	  for(xco = fx; xco <= tx; xco++){
	    for(yco = fy; yco <= ty; yco++){
	      //if(xc==1) printf("%d %d %d %d %d %d\n",xco,yco,fx,tx,fy,ty);
	      TempCell2 = p2p->headCell[xco][yco]->next;
	      while(TempCell2!=TempCell2->next){
		if(TempCell2 != TempCell){
		  dst = sqrt(pow(TempCell->X-TempCell2->X,2)+
			     pow(TempCell->Y-TempCell2->Y,2));
		  lnDens += log(Interaction(TempCell,TempCell2));
		}
		TempCell2 = TempCell2->next; 
	      }
	    }
	  }
	  //printf("lnDens: %f\n",lnDens);
	}
	TempCell = TempCell->next;
      }
    }
  }
  return(lnDens);

}

void PointProcess::GenerateDominatingPoisson(Point2Pattern *p2p, gsl_rng *r){
  p2p->Empty();
  long int NoDP,i;
  double xtemp, ytemp,*marktemp;
  Point2 TempPoint;
  char InWindow;
  
  // Set mark length of the resulting point pattern 
  // to be the same as that of the dominating Poisson process
  p2p->MarkLength = MarkLength;

  NoDP = gsl_ran_poisson(r, TotalBirthRate);  
  p2p->NoP=0;
  for (i=1; i<=NoDP ; i++){
    NewEvent(&xtemp, &ytemp, &marktemp, &InWindow,r);
    if(InWindow==1){
      TempPoint.X = xtemp;
      TempPoint.Y = ytemp;
      TempPoint.No = 0;
      if(MarkLength>0) TempPoint.Mark = marktemp;
      p2p->Add(&TempPoint);
      p2p->NoP = p2p->NoP + 1;
    }else{
      if(MarkLength>0) free(marktemp);
    }
  }
}

void PointProcess::ResetInteractionRange(void){
  // printf("BLA!\n");
};

void PointProcess::UpdateInteractionRange(Point2Pattern *p2p){
  //
};

void PointProcess::UpdateInteractionRange(struct Point2 *Point){
};

double PointProcess::InteractionRangeFct(struct Point2 *Point){
  return(InteractionRange);
};
