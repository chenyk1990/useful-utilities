
/*********************************************************************
 *                                                                   *
 *  int getSPSfield(char *Line, const int rev, \                    *
 *                int const idx, void *addr)                         *
 *                                                                   *
 *  Author:  Fernando Roxo <petro@roxo.org>                          *
 *                                                                   *
 *  Description:                                                     *
 *      This function takes as input one single line record from an  *
 *    SPS file and extracts one single field from it.  The only      *
 *    supported records are point records 'R' or 'S' and the relation*
 *    record 'X'.                                                    *
 *                                                                   *
 *  Parameters:                                                      *
 *    char *Line :  line record to be processed.  First character    *
 *                  must be the label of one of supported records.   *
 *    const int rev :  SPS version being processed.  Supported:     *
 *                  0 - SPS revision 0.                              *
 *                  2 - SPS revision 2.1                             *
 *    int const idx  :  Index of the field to be extracted.  This is *
 *                  the item number as described by standard.        *
 *                  The label is item 0 (zero) and so one.           *
 *    void *addr     :  This is a pointer for the recovered item.    *
 *                  The item type must be according to the format    *
 *                  used by function, that differs from standard     *
 *                  for some fields.                                 *
 *                                                                   * 
 *                                                                   *
 *********************************************************************/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

struct linepars {
   int st;
   int len;
   char *fmt;
   char *name;
};

struct linepars *pparse;
struct linepars *xparse;

//  Version 0 record description:
struct linepars pparse0[] = {
      {  0, 1, "%c","RecTyp"   },  //  0 Label (R or S)
      {  1,16, "%f","LName"    },  //  1 Line name       (string in std)
      { 17, 8, "%f","PNum"     },  //  2 Point number    (string in std)
      { 25, 1, "%d","PInd"     },  //  3 Point index
      { 26, 2, "%d","PCode"    },  //  4 Point code      (string in std)
      { 28, 4, "%d","StaCor"   },  //  5 Static Correction
      { 32, 4, "%f","PDepth"   },  //  6 Point Depth
      { 36, 4, "%d","Datum"    },  //  7 Seismic datum
      { 40, 2, "%d","UpHole"   },  //  8 Uphole Time
      { 42, 4, "%f","WDept"    },  //  9 Water Depth
      { 46, 9, "%lf","MapX"    },  // 10 Map Grid Easting
      { 55,10, "%lf","MapY"    },  // 11 Map Grid Northing
      { 65, 6, "%f","Elev"     },  // 12 Surface Elevation
      { 71, 3, "%d","YDay"     },  // 13 Day of Year
      { 74, 6, "%d","Time"     },  // 14 Time of day     (2I2 in std)
};


struct linepars xparse0[] = {
      {  0, 1, "%c","RecTyp"   },  //  0 Label X
      {  1, 6, "%s","TapeNo"   },  //  1 Tape Number     (String in std)
      {  7, 4, "%d","RecNum"   },  //  2 Field Record Number
      { 11, 1, "%d","RecInc"   },  //  3 Field Record Increment
      { 12, 1, "%c","InstCode" },  //  4 Instrument Code (string in std)
      { 13,16, "%f","SLName"   },  //  5 Source Line Name(string in std)
      { 29, 8, "%f","PNum"     },  //  6 Point number    (string in std)
      { 37, 1, "%d","PIndex"   },  //  7 Point Index
      { 38, 4, "%d","FromChan" },  //  8 From Channel
      { 42, 4, "%d","ToChan"   },  //  9 To Channel
      { 46, 1, "%d","ChanInc"  },  // 10 Channel Increment
      { 47,16, "%f","RLName"   },  // 11 Receiver Line name(string in std)
      { 63, 8, "%f","FromRecv" },  // 12 From Receiver   (string in std)
      { 71, 8, "%f","ToRecv"   },  // 13 To Receiver     (string in std)
      { 79, 1, "%d","RenIx"    }   // 14 Receiver Index

};

//  Version 2.1 record description:
struct linepars pparse2[] = {
      {  0, 1, "%c","RecTyp"   },  //  0 Label (R or S)
      {  1,10, "%f","LName"    },  //  1 Line name       (float in std)
      { 11,10, "%f","PNum"     },  //  2 Point number    (float in std)
      { 23, 1, "%d","PInd"     },  //  3 Point index
      { 24, 2, "%d","PCode"    },  //  4 Point code      (string in std)
      { 26, 4, "%d","StaCor"   },  //  5 Static Correction
      { 30, 4, "%f","PDepth"   },  //  6 Point Depth
      { 34, 4, "%d","Datum"    },  //  7 Seismic datum
      { 38, 2, "%d","UpHole"   },  //  8 Uphole Time
      { 40, 6, "%f","WDept"    },  //  9 Water Depth
      { 46, 9, "%lf","MapX"    },  // 10 Map Grid Easting
      { 55,10, "%lf","MapY"    },  // 11 Map Grid Northing
      { 65, 6, "%f","Elev"     },  // 12 Surface Elevation
      { 71, 3, "%d","YDay"     },  // 13 Day of Year
      { 74, 6, "%d","Time"     }   // 14 Time of day     (3I2 in std)

};

struct linepars xparse2[] = {
      {  0, 1, "%c","RecTyp"   }, //  0 Label X
      {  1, 6, "%s","TapeNo"   }, //  1 Tape Number     (Sting in std)
      {  7, 8, "%d","RecNum"   }, //  2 Field Record Number
      { 15, 1, "%d","RecInc"   }, //  3 Field Record Increment
      { 16, 1, "%c","InstCode" }, //  4 Instrument Code (string in std)
      { 17,10, "%f","SLName"   }, //  5 Source Line Name(float in std)
      { 27,10, "%f","PNum"     }, //  6 Point number    (float in std)
      { 37, 1, "%d","PIndex"   }, //  7 Point Index
      { 38, 5, "%d","FromChan" }, //  8 From Channel
      { 43, 5, "%d","ToChan"   }, //  9 To Channel
      { 48, 1, "%d","ChanInc"  }, // 10 Channel Increment
      { 49,10, "%f","RLName"   }, // 11 Receiver Line name(float in std)
      { 59,10, "%f","FromRecv" }, // 12 From Receiver   (float in std)
      { 69,10, "%f","ToRecv"   }, // 13 To Receiver     (float in std)
      { 79, 1, "%d","RenIx"    }  // 14 Receiver Index

};

int getSPSfield(char *Line, int rev, int idx, void *addr)
{

   int ilast, Ini,nconv;
   char tmp;

   if (rev != 0 && rev != 2) {
         // throw a tantrum
         fprintf(stderr,"****  Unsupported SPS format %d!!\n\n",rev);
         return 1;
   }
  /*  Define the structure to use */
   if ( Line[0] == 'R' || Line[0] == 'S' ) {/* this is a point record */
      if ( rev == 0 ) {
         pparse = pparse0;
      } else if ( rev == 2 ) {
         pparse = pparse2;
      }
   }
   else if ( Line[0] == 'X' ) {             /*this is a relation record */
      if ( rev == 0 ) {
         pparse = xparse0;
      } else if ( rev == 2 ) {
         pparse = xparse2;
      }
   }
   else  {                                  /* what is this? */
      /* throw a tantrum */
      fprintf(stderr,"****  Unsupported SPS record %c!!\n\n",Line[0]);
      return 1;
   }

   Ini = pparse[idx].st;
   ilast = Ini + pparse[idx].len;

   if ( strlen(Line) < Ini )
   {   /* No field left to be read */
      /* Just return with target untouched */
      return -1;  /*  Return a -1 like sscanf would if failed. */
   } else if ( strlen(Line) < ilast) {
      ilast = strlen(Line);  /*  Just reduce the parse length */
   }

   // Save the starting character of next field
   tmp = Line[ilast];

   Line[ilast] = '\0';

   /* This option will produce a zero if conversion worked, nonzero otherwise. */
   /* nconv = 1-sscanf(&Line[Ini], pparse[idx].fmt, addr); */
   nconv = sscanf(&Line[Ini], pparse[idx].fmt, addr);
   /* Debug */
   /* fprintf (stderr,"%s: %s '%s' = %d\n",pparse[idx].name,pparse[idx].fmt, &Line[Ini],nconv);  */

   /* recover the starting character of next field */
   Line[ilast] = tmp;

   return nconv;
}


