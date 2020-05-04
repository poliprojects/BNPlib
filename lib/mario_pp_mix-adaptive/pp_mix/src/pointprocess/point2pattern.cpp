#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream> 
#include <fstream> 
#include <string>
#include <sstream>
#include <algorithm>
#include "point2pattern.h"
using namespace std; 

Point2Pattern::Point2Pattern(double xmin, double xmax,
				  double ymin, double ymax, 
				  long int mxc, long int myc){
    long int i,j;
    MarkLength = 0;
    UpperLiving[0] = 0;
    UpperLiving[1] = 0;
    Xmin = xmin; Xmax = xmax;
    Ymin = ymin; Ymax = ymax;
    if(Xmin>Xmax){ printf("Xmin>Xmax\n"); exit(1);}
    if(Ymin>Ymax){ printf("Ymin>Ymax\n"); exit(1);}
    NoP = 0;
    dummyCell = (struct Point2 *) malloc(sizeof *dummyCell);
    if(dummyCell == NULL){
      printf("Cannot allocate memory for dummyCell in Point2Pattern::Point2Pattern\n");
      exit(1);
    }
    dummyCell->next = dummyCell;
    dummyCell->No = 0;
    MaxXCell = mxc; MaxYCell = myc;
    if(MaxXCell>9) MaxXCell = 9;
    if(MaxYCell>9) MaxYCell = 9;
    for(i=0;i<=MaxXCell;i++){
      for(j=0;j<=MaxYCell;j++){
	headCell[i][j] = 
	  (struct Point2 *) malloc(sizeof *headCell[i][j]);
	if(headCell[i][j] == NULL){
	  printf("Cannot allocate memory for headCell in Point2Pattern::Point2Pattern\n");
	  exit(1);
	}
	headCell[i][j]->next=dummyCell;
      }
    }
    XCellDim = (Xmax-Xmin)/((double)(MaxXCell+1));
    YCellDim = (Ymax-Ymin)/((double)(MaxYCell+1));
  };

void Point2Pattern::Diagnostics(){
  printf("MaxXCell: %ld  MaxYCell: %ld\n",MaxXCell,MaxYCell);
  printf("NoP: %ld\n",NoP);
  printf("XCellDim: %f  YCellDim: %f\n",XCellDim,YCellDim);
  printf("Xmin: %f  Xmax: %f\n",Xmin,Ymax);
  printf("Ymin: %f  Ymax: %f\n",Ymin,Ymax);
}

void Point2Pattern::CopyTo(Point2Pattern *p2p){
  long int xc,yc,i;
  Point2 *Point1, *Point2;
  p2p->Empty();
  p2p->MarkLength = this->MarkLength;
  for(xc=0;xc<=MaxXCell;xc++){
    for(yc=0;yc<=MaxYCell;yc++){
      Point1 = headCell[xc][yc];
      while(Point1->next!=Point1->next->next){
	Point2 = (struct Point2 *) malloc(sizeof *Point2);
	if(Point2 == NULL){
	  printf("Cannot allocate memory for Point2 in Sampler::CopyTo\n");
	  exit(1);
	}
	Point2->X = Point1->next->X;
	Point2->Y = Point1->next->Y;
	Point2->Beta = Point1->next->Beta;
	Point2->TempBeta = Point1->next->TempBeta;
	Point2->next = p2p->headCell[xc][yc]->next;
	p2p->headCell[xc][yc]->next = Point2;
	if(this->MarkLength>0){
	  Point2->Mark = (double *) calloc(MarkLength,sizeof(double));
	  if(Point2->Mark == NULL){
	    printf("Cannot allocate memory for mark in Point2Pattern::CopyTo\n");
	    exit(1);
	  }
	  for(i=0;i<MarkLength;i++){
	    *(Point2->Mark+i) = *(Point1->next->Mark+i);   
	  }
	}
	Point1 = Point1->next;
      }      
    }    
  }
  p2p->NoP = NoP;
}

void Point2Pattern::Add(Point2 *Point){
  long int XCell, YCell, i;
  double *mark;
  Point2 *Point2;
  Point2 = (struct Point2 *) malloc(sizeof *Point2);
  if(Point2 == NULL){
    printf("Cannot allocate memory for Point2 in Point2Pattern::Add\n");
    exit(1);
  }
  Point2->X = Point->X;
  Point2->Y = Point->Y;
  Point2->No = Point->No;
  Point2->Beta = Point->Beta;
  if(MarkLength>0){
    mark = (double *) calloc(MarkLength,sizeof(double));
    if(mark == NULL){
      printf("Cannot allocate memory for mark in Point2Pattern::Add\n");
      exit(1);
    }    
    for(i=0;i<MarkLength;i++){
      *(mark+i) = *(Point->Mark+i);   
    }
    Point2->Mark = mark;
  }

  XCell = //(long int) 
    (floor((Point->X-this->Xmin)/this->XCellDim));  
  if(XCell>this->MaxXCell) XCell = this->MaxXCell;
  YCell = //(long int) 
    (floor((Point->Y-this->Ymin)/this->YCellDim));  
  if(YCell>this->MaxYCell) YCell = this->MaxYCell;
  
  Point2->next = this->headCell[XCell][YCell]->next;
  this->headCell[XCell][YCell]->next = Point2;

  NoP++;
}

void Point2Pattern::Remove(long int N){
  // Remove Nth point in point pattern
  long int tmpN,xc,yc;
  Point2 *Point1, *Point2;

  tmpN = 0;
  for(xc=0;xc<=MaxXCell;xc++){
    for(yc=0;yc<=MaxYCell;yc++){
      //printf("N: %d tmpN: %d\n",N,tmpN);
      Point1 = headCell[xc][yc];
      while(Point1->next!=Point1->next->next){
	tmpN++;
	if(tmpN==N){
	  // Remove Nth point
	  Point2 = Point1->next;
	  Point1->next = Point2->next;
	  if(MarkLength>0) free(Point2->Mark);
	  free(Point2);
	}else{
	  Point1 = Point1->next;
	}
	if(tmpN==N) break;
      }
      if(tmpN==N) break;
    }
    if(tmpN==N) break;
  }
  NoP--;
}

struct Point2* Point2Pattern::Extract(long int N){
  // Remove Nth point in point pattern
  long int tmpN,xc,yc;
  Point2 *Point1, *Point2;

  tmpN = 0;
  for(xc=0;xc<=MaxXCell;xc++){
    for(yc=0;yc<=MaxYCell;yc++){
      //printf("N: %d tmpN: %d\n",N,tmpN);
      Point1 = headCell[xc][yc];
      while(Point1->next!=Point1->next->next){
	tmpN++;
	if(tmpN==N){
	  // Remove Nth point
	  Point2 = Point1->next;
	}else{
	  Point1 = Point1->next;
	}
	if(tmpN==N) break;
      }
      if(tmpN==N) break;
    }
    if(tmpN==N) break;
  }
  return(Point2);
}


void Point2Pattern::Print(){
  long int i,j,k;
  k =0;
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	k++;
	printf("%f %f %ld %ld %ld=%d %ld=%d UL0 %d UL1 %d %f\n",
	       TempCell->X,TempCell->Y,k,
	       TempCell->No,
	       i,int(TempCell->X/XCellDim),
	       j,int(TempCell->Y/YCellDim),
	       TempCell->InLower[0],TempCell->InLower[1],
	       TempCell->Beta);
	TempCell = TempCell->next;
      }
    }
  }
  printf("Printed %ld points.\n",k);
}

long int Point2Pattern::Count(){
  long int i,j,k;
  k =0;
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	k++;
	TempCell = TempCell->next;
      }
    }
  }
  //printf("Printed %d points.\n",k);
  return(k);
}


void Point2Pattern::Empty(){
  struct Point2 *TempCell, *TempCell2;
  long int i,j;
  
  //k=0;

  for(i=0; i<=this->MaxXCell; i++){
    for(j=0; j<=this->MaxYCell; j++){
      TempCell = headCell[i][j]->next;
      while(TempCell!=TempCell->next){	
	//k++; printf("%d %d %d\n",i,j,k);
	TempCell2 = TempCell->next;
	if(MarkLength>0) free(TempCell->Mark);
	free(TempCell);
	TempCell = TempCell2;
      }
      headCell[i][j]->next = dummyCell;
    }
  }
}

void Point2Pattern::Clean(){
  struct Point2 *TempCell, *TempCell2;
  long int i,j;
  
  for(i=0; i<=MaxXCell; i++){
    for(j=0; j<=MaxYCell; j++){
      TempCell = headCell[i][j];
      TempCell2 = headCell[i][j]->next;
      while(TempCell2!=TempCell2->next){
	TempCell2->No = 0;
	if(TempCell2->InLower[0]==0){
	  TempCell->next = TempCell2->next;
	  if(MarkLength>0) free(TempCell2->Mark);
	  free(TempCell2);
	  TempCell2 = TempCell->next;
	}
	else{
	  TempCell2 = TempCell2->next;
	  TempCell = TempCell->next;
	}
      }
    }
  }
}

void Point2Pattern::AppendToFile(const char *FileName){
  FILE *out;
  long int i,j,iMark;
  out = fopen(FileName,"a");

  // Produce a row of -1 to separate patterns
  fprintf(out,"%f\t%f\t", -1.0,-1.0);
  for(iMark=0;iMark<MarkLength;iMark++)
    fprintf(out,"%lf\t",-1.0);
  fprintf(out,"%d\n",-1);

  // Dump pattern to disc
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	fprintf(out,"%f\t%f\t",
	       TempCell->X,TempCell->Y);
	for(iMark=0;iMark<MarkLength;iMark++)
	  fprintf(out,"%lf\t",*(TempCell->Mark + iMark));
	fprintf(out,"%ld\n",TempCell->No);
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}


void Point2Pattern::AppendToFile(const char *FileName, 
				 double xmin, double xmax,
				 double ymin,double ymax){
  FILE *out;
  long int i,j;
  out = fopen(FileName,"a");
  fprintf(out,"%f\t%f\t%d\n",
	  -1.0,-1.0,-1);
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	if((TempCell->X>=xmin) & (TempCell->X<=xmax) &
	   (TempCell->Y>=ymin) & (TempCell->Y<=ymax))
	  fprintf(out,"%f\t%f\t%ld\n",
		  TempCell->X,TempCell->Y,TempCell->No);
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}



void Point2Pattern::DumpToFile(const char *FileName){
  FILE *out;
  long int i,j,iMark;
  long int k;
  out = fopen(FileName,"w");
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%ld %ld %ld:\n",i,j,MarkLength);
      TempCell = headCell[i][j]->next;
      k=0;
      while(TempCell->next != TempCell){
        fprintf(out,"%f\t%f\t",
              TempCell->X,TempCell->Y);
        for(iMark=0;iMark<MarkLength;iMark++){
          fprintf(out,"%lf\t",*(TempCell->Mark + iMark));
        }
        fprintf(out,"%ld\n",TempCell->No);
        TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}

Eigen::MatrixXd Point2Pattern::to_eigen_mat() {
  std::vector<Eigen::VectorXd> out;
  long int i, j, iMark;
  long int k;
  struct Point2 *TempCell;
  double currx, curry;
  for (i = 0; i <= MaxXCell; i++)
  {
    for (j = 0; j <= MaxYCell; j++)
    {
      TempCell = headCell[i][j]->next;
      k = 0;
      while (TempCell->next != TempCell)
      {
        currx = TempCell->X;
        curry = TempCell->Y;
        Eigen::VectorXd currpoint(2);
        currpoint << currx, curry;
        out.push_back(currpoint);
        for (iMark = 0; iMark < MarkLength; iMark++) {
          // fprintf(out, "%lf\t", *(TempCell->Mark + iMark));
        }
        // fprintf(out, "%ld\n", TempCell->No);
        TempCell = TempCell->next;
      }
    }
  }
  // fclose(out);
  int nrows = out.size();
  int ncols = 2;

  MatrixXd out2(nrows, ncols);
  for (int i = 0; i < nrows; i++)
    out2.row(i) = out[i].transpose();

  return out2;
}

void Point2Pattern::DumpToFileNoNo(const char *FileName){
  FILE *out;
  long int i,j,iMark;
  long int k;
  out = fopen(FileName,"w");
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%ld %ld %ld:\n",i,j,MarkLength);
      TempCell = headCell[i][j]->next;
      k=0;
      while(TempCell->next != TempCell){
	fprintf(out,"%f\t%f\t",
	       TempCell->X,TempCell->Y);
	for(iMark=0;iMark<MarkLength;iMark++){
	  fprintf(out,"%lf\t",*(TempCell->Mark + iMark));
	}
	fprintf(out,"\n");
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}



void Point2Pattern::DumpToFileBeta(const char *FileName){
  FILE *out;
  long int i,j,iMark;
  out = fopen(FileName,"w");
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	fprintf(out,"%f\t%f\t",
	       TempCell->X,TempCell->Y);
	for(iMark=0;iMark<MarkLength;iMark++)
	  fprintf(out,"%lf\t",*(TempCell->Mark + iMark));
	fprintf(out,"%ld\t",TempCell->No);
	fprintf(out,"%lf\n",TempCell->Beta);
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}



void Point2Pattern::DumpToFile(const char *FileName,
			       double xmin, double xmax,
			       double ymin, double ymax){
  FILE *out;
  long int i,j;
  out = fopen(FileName,"w");
  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	if((TempCell->X>=xmin) & (TempCell->X<=xmax) &
	   (TempCell->Y>=ymin) & (TempCell->Y<=ymax))
	  fprintf(out,"%f\t%f\t%ld\n",
		  TempCell->X,TempCell->Y,TempCell->No);
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}


void Point2Pattern::DumpToFile(const char *FileName,
			       const char *WriteMode,
			       double xmin, double xmax,
			       double ymin, double ymax){
  FILE *out;
  long int i,j;
  out = fopen(FileName,WriteMode);
  if(WriteMode=="a")
    fprintf(out,"%f\t%f\t%d\n",0.0,0.0,-1);

  struct Point2 *TempCell;
  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      //printf("%d %d:\n",i,j);
      TempCell = headCell[i][j]->next;
      while(TempCell->next != TempCell){
	if((TempCell->X>=xmin) & (TempCell->X<=xmax) &
	   (TempCell->Y>=ymin) & (TempCell->Y<=ymax))
	  fprintf(out,"%f\t%f\t%ld\n",
		  TempCell->X,TempCell->Y,TempCell->No);
	TempCell = TempCell->next;
      }
    }
  }
  fclose(out);
}



void Point2Pattern::ReadFromFile(const char *FileName){
  ReadFromFile(FileName,0);
    /*
  FILE *out;
  int dummyint;
  long int k,XCell,YCell;
  float xs,ys, f1;
  out = fopen(FileName,"r");
  struct Point2 *TempCell;
  k=0;
  while(feof(out)==0){
    k++;
    dummyint = fscanf(out,"%f%f\n",&xs,&ys);
    //printf("%f %f\n",xs,ys);
    TempCell = (struct Point2 *) malloc(sizeof *TempCell);
    if(TempCell==NULL) printf("NULL\n");
    TempCell->No = k;
    TempCell->X = xs;
    TempCell->Y = ys;
    TempCell->InLower[0] = 1;
    TempCell->InLower[1] = 1;

    f1 = (xs-Xmin)/XCellDim;  XCell = int(f1);
    if(XCell>MaxXCell) XCell = MaxXCell;
    f1 = (ys-Ymin)/YCellDim;  YCell = int(f1);
    if(YCell>MaxYCell) YCell = MaxYCell;

    TempCell->next = headCell[XCell][YCell]->next;
    headCell[XCell][YCell]->next = TempCell;

  }
  fclose(out);
  printf("%ld points loaded.\n",k);
    */
}

void Point2Pattern::ReadFromFile(const char *FileName,
				 long int marklength){
  FILE *out;
  int dummyint;
  long int k,XCell,YCell, marki;
  double xs,ys, f1, tempd;
  double *TempMark;
  MarkLength = marklength;
  out = fopen(FileName,"r");
  if(out==NULL){
    printf("No file named %s\n",FileName);
    exit(1);
  }
  struct Point2 *TempCell;
  k=0;
  while(feof(out)==0){
    dummyint = fscanf(out,"%lf%lf",&xs,&ys);
    if(dummyint==2){
    k++;
    //printf("%lf %lf",xs,ys);
    if((xs<Xmin) || (xs>Xmax) || (ys<Ymin)|| (ys>Ymax)){
      printf("Point in file outside pattern window!\n");
      printf("Point:  x = %lf  y = %lf\n",xs,ys);
      printf("Window: x = %lf .. %lf\n",Xmin,Xmax);
      printf("        y = %lf .. %lf\n",Ymin,Ymax);
      exit(1);
    }
    
    TempCell = (struct Point2 *) malloc(sizeof *TempCell);
    if(TempCell==NULL){ 
      printf("Cannot allocate Point2 in Point2Pattern::ReadFromFile\n");
      exit(1);
    }
    TempCell->No = k;
    TempCell->X = xs;
    TempCell->Y = ys;
    TempCell->InLower[0] = 1;
    TempCell->InLower[1] = 1;
    if(MarkLength>0){
      TempMark = (double *) calloc(MarkLength,sizeof(double));
      if(TempMark==NULL){
	printf("Cannot allocate memory for mark in Point2Pattern::ReadFromFile\n");
	exit(1);
      }
      TempCell->Mark = TempMark;
      for(marki = 0; marki<MarkLength;marki++){
	dummyint = fscanf(out,"%lf",&tempd);
	*(TempCell->Mark+marki) = tempd;
	//printf(" %lf",tempd); 
      }
    }
    //printf("\n");
    //f1 = (xs-Xmin)/XCellDim;  XCell = int(f1);
    XCell = floor((TempCell->X-this->Xmin)/this->XCellDim);  
    if(XCell>MaxXCell) XCell = MaxXCell;
    //f1 = (ys-Ymin)/YCellDim;  YCell = int(f1);
    YCell = floor((TempCell->Y-this->Ymin)/this->YCellDim);  
    if(YCell>MaxYCell) YCell = MaxYCell;

    //printf("Hertil %ld %ld\n",XCell,YCell);

    TempCell->next = headCell[XCell][YCell]->next;
    headCell[XCell][YCell]->next = TempCell;
    }
  }
  fclose(out);
  printf("%ld points loaded.\n",k);

}

void Point2Pattern::ReadFromFileExperimental(const char *FileName,
					     long int marklength){
  FILE *out;
  int dummyint;
  long int k,XCell,YCell, marki, RowElement, LineNumber;
  double xs,ys, f1, tempd;
  string line;
  double value;
  double *TempMark;
  struct Point2 *TempPoint;

  // Check if input file exists
  ifstream file(FileName);
  if(file.fail()){
    printf("No such file!\n");
    exit(1);
  }

  // Set mark length of point pattern to that specified by input
  MarkLength = marklength;
  LineNumber = 0;

  while ( getline(file, line) ){
    istringstream iss(line);
    LineNumber++;
    // Allocate memory for point
    TempPoint = (struct Point2 *) malloc(sizeof *TempPoint);
    RowElement=0;
    // If marked point allocate memory for mark vector
    if(MarkLength>0){
      TempMark = (double *) calloc(MarkLength,sizeof(double));
      if(TempMark==NULL){
	printf("Cannot allocate memory for mark in Point2Pattern::ReadFromFileExperimental\n");
	exit(1);
      }
      TempPoint->Mark = TempMark;
    }
    while(iss >> value){      
      RowElement++;
      if(RowElement==1) xs = value;
      if(RowElement==2) ys = value;      
      if(RowElement>2) *(TempPoint->Mark+RowElement-3) = value;
    }
    //printf("(x,y)=(%lf,%lf)\n",xs,ys);
    if((RowElement-2)<MarkLength){
      printf("Not enough columns in file\n");
      exit(1);
    }
    //printf("RowElement: %ld\n",RowElement);

    // Check if point is inside window
    if((xs<Xmin) || (xs>Xmax) || (ys<Ymin)|| (ys>Ymax)){
      printf("Point in file outside pattern window!\n");
      printf("Point:  x = %lf  y = %lf\n",xs,ys);
      printf("Window: x = %lf .. %lf\n",Xmin,Xmax);
      printf("        y = %lf .. %lf\n",Ymin,Ymax);
      exit(1);
    }
    TempPoint->No = LineNumber;
    TempPoint->X = xs;
    TempPoint->Y = ys;
    TempPoint->InLower[0] = 1;
    TempPoint->InLower[1] = 1;

    // Place point into subcell structure
    f1 = (xs-Xmin)/XCellDim;  XCell = int(f1);
    if(XCell>MaxXCell) XCell = MaxXCell;
    f1 = (ys-Ymin)/YCellDim;  YCell = int(f1);
    if(YCell>MaxYCell) YCell = MaxYCell;

    TempPoint->next = headCell[XCell][YCell]->next;
    headCell[XCell][YCell]->next = TempPoint;

  }

  printf("%ld points loaded.\n",LineNumber);

}



void Point2Pattern::Truncate(double xlower, double xupper,
			     double ylower, double yupper){

  long int i,j;
  struct Point2 *Point1, *Point3;

  for(i=0;i<=MaxXCell;i++){
    for(j=0;j<=MaxYCell;j++){
      Point1 = headCell[i][j];
      while(Point1->next->next != Point1->next){
	Point3 = Point1->next;
	if((Point3->X<xlower) || (Point3->X>xupper) ||
	   (Point3->Y<ylower) || (Point3->Y>yupper)){
	  // Point to be truncated
	  Point1->next = Point3->next;
	  if(MarkLength>0) free(Point3->Mark);
	  free(Point3);
	}else{
	  Point1 = Point1->next;
	}
      }
    }
  }
}


struct Point2* Point2Pattern::FirstSloppy(double R, struct Point2 *CentrePoint){
  // Unlike the First function FirstSloppy does not take into account
  // if the found point is within distance R of CentrePoint. The idea
  // is to save calculating the distance if this is being calculated
  // anyway by the code that calls the function. 
  long int px,py,fx,fy,tx,ty;
  double dx,dy;
  struct Point2 *ReturnPoint;

  fx = max((long int) floor((CentrePoint->X-R-Xmin) / XCellDim),(long int) 0);
  fy = max((long int) floor((CentrePoint->Y-R-Ymin) / YCellDim),(long int) 0);

  tx = min((long int) floor((CentrePoint->X+R-Xmin) / XCellDim), MaxXCell);
  ty = min((long int) floor((CentrePoint->Y+R-Ymin) / YCellDim), MaxYCell);
  
  px = fx; py = fy;

  ReturnPoint = this->headCell[px][py]->next;

  while((ReturnPoint == CentrePoint) || (ReturnPoint == ReturnPoint->next)){
    if(ReturnPoint == CentrePoint) ReturnPoint = ReturnPoint->next;
    ReturnPoint = ReturnPoint->next;
    if(ReturnPoint->next == ReturnPoint){
      px++;
      if(px>tx){ px = fx; py++;}
      if(py>ty) return(NULL);
      ReturnPoint = this->headCell[px][py]->next;
    }
  }
  return(ReturnPoint);
}

struct Point2* Point2Pattern::NextSloppy(struct Point2 *OldPoint, 
				   double R, struct Point2 *CentrePoint){
  struct Point2 *ReturnPoint;
  long int px,py,fx,fy,tx,ty;
  ReturnPoint = OldPoint->next;
  if((ReturnPoint != CentrePoint) && (ReturnPoint!=ReturnPoint->next)){
    return(ReturnPoint);
  }
  px = floor((OldPoint->X-Xmin) / XCellDim);
  py = floor((OldPoint->Y-Ymin) / YCellDim);

  fx = max((long int) floor((CentrePoint->X-R-Xmin) / XCellDim),(long int) 0);

  tx = min((long int) floor((CentrePoint->X+R-Xmin) / XCellDim), MaxXCell);
  ty = min((long int) floor((CentrePoint->Y+R-Ymin) / YCellDim), MaxYCell);

  while((ReturnPoint == CentrePoint) || (ReturnPoint == ReturnPoint->next)){
    if(ReturnPoint == CentrePoint) ReturnPoint = ReturnPoint->next;
    ReturnPoint = ReturnPoint->next;
    if(ReturnPoint->next == ReturnPoint){
      px++;
      if(px>tx){ px = fx; py++;}
      if(py>ty) return(NULL);
      ReturnPoint = this->headCell[px][py]->next;
    }
  }
  return(ReturnPoint);
}


struct Point2* Point2Pattern::First(double R, struct Point2 *CentrePoint){
  // This function finds the "first" point within distance R from CentrePoint
  // which is also different from CentrePoint. If no such point is found
  // the function return NULL. The search for the "first" point happens
  // in a rectangular region which contains a cricle centred at CentrePoint
  // with radius R. The remainig points are found  using the
  // Next function repeatedly. 
  // This way detailed knowledge about the variable structure of the
  // Point2Pattern type is not needed.

  long int px,py,fx,fy,tx,ty;
  double dx,dy;
  double DistSqr, Rsqr;
  Rsqr = R*R;
  struct Point2 *ReturnPoint;
  
  fx = max((long int) floor((CentrePoint->X-R-Xmin) / XCellDim),(long int) 0);
  fy = max((long int) floor((CentrePoint->Y-R-Ymin) / YCellDim),(long int) 0);

  tx = min((long int) floor((CentrePoint->X+R-Xmin) / XCellDim), MaxXCell);
  ty = min((long int) floor((CentrePoint->Y+R-Ymin) / YCellDim), MaxYCell);
  
  px = fx; py = fy;

  ReturnPoint = this->headCell[px][py]->next;

  DistSqr = (ReturnPoint->X-CentrePoint->X)*(ReturnPoint->X-CentrePoint->X)
    + (ReturnPoint->Y-CentrePoint->Y)*(ReturnPoint->Y-CentrePoint->Y);

  while((ReturnPoint == CentrePoint) || (ReturnPoint == ReturnPoint->next)
	|| (DistSqr>Rsqr)){
    if(ReturnPoint == CentrePoint) ReturnPoint = ReturnPoint->next;
    ReturnPoint = ReturnPoint->next;
    if(ReturnPoint->next == ReturnPoint){
      px++;
      if(px>tx){ px = fx; py++;}
      if(py>ty)	return(NULL);
      ReturnPoint = this->headCell[px][py]->next;
    }
    DistSqr = (ReturnPoint->X-CentrePoint->X)*(ReturnPoint->X-CentrePoint->X)
      + (ReturnPoint->Y-CentrePoint->Y)*(ReturnPoint->Y-CentrePoint->Y);
  }
  return(ReturnPoint);
}

struct Point2* Point2Pattern::Next(struct Point2 *OldPoint, 
				   double R, struct Point2 *CentrePoint){
  struct Point2 *ReturnPoint;
  long int px,py,fx,fy,tx,ty;
  double DistSqr, Rsqr;
  Rsqr = R*R;
  ReturnPoint = OldPoint->next;
  if((ReturnPoint != CentrePoint) && (ReturnPoint!=ReturnPoint->next)){
    DistSqr = (ReturnPoint->X-CentrePoint->X)*(ReturnPoint->X-CentrePoint->X)
      + (ReturnPoint->Y-CentrePoint->Y)*(ReturnPoint->Y-CentrePoint->Y);
    if(DistSqr<Rsqr) return(ReturnPoint);
  }
  px = floor((OldPoint->X-Xmin) / XCellDim);
  py = floor((OldPoint->Y-Ymin) / YCellDim);

  fx = max((long int) floor((CentrePoint->X-R-Xmin) / XCellDim),(long int) 0);

  tx = min((long int) floor((CentrePoint->X+R-Xmin) / XCellDim), MaxXCell);
  ty = min((long int) floor((CentrePoint->Y+R-Ymin) / YCellDim), MaxYCell);
  
  while((ReturnPoint == CentrePoint) || (ReturnPoint == ReturnPoint->next)
    	|| (DistSqr>Rsqr)){
    if(ReturnPoint == CentrePoint) ReturnPoint = ReturnPoint->next;
    ReturnPoint = ReturnPoint->next;
    if(ReturnPoint->next == ReturnPoint){
      px++;
      if(px>tx){ px = fx; py++;}
      if(py>ty) return(NULL);
      ReturnPoint = this->headCell[px][py]->next;
    }
    DistSqr = (ReturnPoint->X-CentrePoint->X)*(ReturnPoint->X-CentrePoint->X)
      + (ReturnPoint->Y-CentrePoint->Y)*(ReturnPoint->Y-CentrePoint->Y);
  }
  return(ReturnPoint);
}

struct Point2* Point2Pattern::First(){
  // This function finds the "first" point in the point pattern.
  // The remainig points are found  using the
  // Next function repeatedly. 
  // This way detailed knowledge about the variable structure of the
  // Point2Pattern type is not needed.

  long int px,py;
  struct Point2 *ReturnPoint;
  ReturnPoint = headCell[0][0]->next;
  px = 0; py = 0;
  while(ReturnPoint == ReturnPoint->next){
    px++;
    if(px > this->MaxXCell){ px = 0; py++;}
    if(py > this->MaxYCell){
      //printf("No first point!\n");
      return(NULL);
    }
    ReturnPoint = this->headCell[px][py]->next;
  }
  return(ReturnPoint);
}

struct Point2* Point2Pattern::Next(struct Point2 *OldPoint){
  struct Point2 *ReturnPoint;
  long int px,py;
  ReturnPoint = OldPoint->next;
  px = floor((OldPoint->X-Xmin) / XCellDim);
  py = floor((OldPoint->Y-Ymin) / YCellDim);

  if(ReturnPoint!= ReturnPoint->next) return(ReturnPoint);

  px = floor((OldPoint->X-Xmin) / XCellDim);
  py = floor((OldPoint->Y-Ymin) / YCellDim);

  while(ReturnPoint == ReturnPoint->next){
    px++;
    if(px > this->MaxXCell){ px = 0; py++;}
    if(py > this->MaxYCell) return(NULL);
    ReturnPoint = this->headCell[px][py]->next;
  }
  return(ReturnPoint);
}

