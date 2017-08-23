
#ifndef SUGEOM_H
#define SUGEOM_H

#include <math.h>

#include "su.h"
#include "segy.h"


char textbuffer[120];      /* buffer for SPS records read and process */

int verbose;               /* flag to get advisory messages */

/* Structure for point information*/
struct  PointInfo {

   float  Line;     /*  1- Line name*/
   float  Point;    /*  2- Point number (station, PT, etc)*/
   int    StaCor;   /*  5- Static correction */
   float  PDepth;   /*  6- Point depth in relation to surface (ex buried SP)*/
   double X;        /* 10- X coordinate*/
   double Y;        /* 11- Y coordinate*/
   float  Elev;     /* 12- Surface elevation at point*/

};

/* Structure for relation (aka register) information*/
struct  RegInfo {

   int   Num;       /*  2- Field Record Number*/
   int   Inc;       /*  3- Field Record Increment*/
   float SLine;     /*  5- Source Line*/
   float SPoint;    /*  6- Source Point*/
   int   FChan;     /*  8- Start Channel*/
   int   TChan;     /*  9- To Channel*/
   int   IncChan;   /* 10- Channel Increment for multicomponent data*/
   float RLine;     /* 11- Receiver Line*/
   float FRecv;     /* 12- From Receiver Station*/
   float TRecv;     /* 13- To receiver Station*/
   int   RInc;      /* 14- Receiver Index*/

};


/* function proto.*/
int countRec(FILE*,  char );
int getPoints(FILE*, char ,int , struct  PointInfo*);
int getFiles(FILE*,  int ,struct  RegInfo*);
int getSPSfield(char *, int , int , void*);
int GetPointIndex(struct PointInfo* ,int , float ,float );

segy tr;


#endif   /* SUGEOM_H*/

