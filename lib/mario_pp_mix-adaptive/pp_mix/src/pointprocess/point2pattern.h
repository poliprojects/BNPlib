#ifndef GUARD_point2pattern_h
#define GUARD_point2pattern_h

#include <Eigen/Dense>
#include "../utils.hpp"

struct Point2{ long int No; double X; double Y; double *Mark;
  int InLower[2]; 
  double Beta; double TempBeta; struct Point2 *next; }; 

class Point2Pattern {
public:
  long int UpperLiving[2];
  long int MaxXCell, MaxYCell, NoP, MarkLength;
  double XCellDim, YCellDim, Xmin, Xmax, Ymin, Ymax;
  struct Point2 *headCell[10][10],*dummyCell;
  char DirX[10], DirY[10]; // what is this? this is extra to what existed in point2pattern.c class-definition
 
  Point2Pattern(double xmin, double xmax,
		double ymin, double ymax, 
		long int mxc, long int myc);
  ~Point2Pattern(){}
  void Print();
  long int Count();
  void Empty();
  void Clean();
  void Add(Point2 *Point);
  void CopyTo(Point2Pattern *p2p);
  void Remove(long int N);
  struct Point2* Extract(long int N);
  void DumpToFile(const char *FileName);
  Eigen::MatrixXd to_eigen_mat();
  void DumpToFileNoNo(const char *FileName);
  void DumpToFileBeta(const char *FileName);
  void DumpToFile(const char *FileName, 
		  double xmin, double xmax,
		  double ymin,double ymax);
  void DumpToFile(const char *FileName, 
		  const char *WriteMode,
		  double xmin, double xmax,
		  double ymin,double ymax);
  void AppendToFile(const char *FileName);
  void AppendToFile(const char *FileName, 
		    double xmin, double xmax,
		    double ymin,double ymax);
  void ReadFromFile(const char *FileName);
  void ReadFromFile(const char *FileName, long int MarkLength);
  void ReadFromFileExperimental(const char *FileName, long int MarkLength);
  void Diagnostics();
  void Truncate(double, double, double, double);
  struct Point2* FirstSloppy(double R, Point2 *CentrePoint);
  struct Point2* NextSloppy(Point2 *OldPoint, double R, Point2 *CentrePoint);
  struct Point2* First(double R, Point2 *CentrePoint);
  struct Point2* Next(Point2 *OldPoint, double R, Point2 *CentrePoint);
  struct Point2* First();
  struct Point2* Next(Point2 *OldPoint);
};
#endif 
