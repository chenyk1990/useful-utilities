/*=========================================================
  calculate 2D-interpolation coeff.
  the value at point (x,y) in the squre = Sum coef[i]*z[i] i=0,3
    0-------4-------1
    |       |       |
    |       |   x   |
    |       |       |
    5-------8-------6
    |       |       |
    |       |       |
    |       |       |
    2-------7-------3

  algorithm: divide the squre by 4 and compute the contributions
  for each grid to the center of the sub-square where the point
  is located. then doing the same to the sub-square, ....

  Author	Lupei Zhu
  Revision
	12/04/2000	Initial coding
=========================================================*/
void intp2D(float x, float y,	/* In: location of the point [0-1]*/
		float coef[],	/* Out: coefficients */
		int n		/* In: number of divisions */
	)
{
   int i, j;
   float a[4][4];
   void intp2Dsub(float,float,float a[4][4],int);
   for(i=0; i<4; i++) {
      for(j=0; j<4; j++) {
	 a[i][j] = 0;
      }
      a[i][i] = 1.;
   } 
   intp2Dsub(x,y,a,n);
   for(i=0; i<4; i++) {
      coef[i] = 0.;
      for(j=0; j<4; j++) coef[i] += 0.25*a[j][i];
   } 
   return;
}

/*=====================================================
recursive part, a[i][j] is the contribution from grid j
of the original square to the node i of the sub-square
=======================================================*/
void intp2Dsub(float x, float y, float a[4][4], int n) {
   int i, j;
   void tsfm(int, float a[4][4]);
   i = x>0.5;
   j = y>0.5;
   tsfm(2*i+j,a);
   x = 2*x-i;
   y = 2*y-j;
   n--;
   if (n==0) return;
   intp2Dsub(x, y, a, n);
}

void tsfm(int subbox, float a[4][4]) {
   int   i, j, k;
   static int  indx[4][4] = {{5,8,2,7},{0,4,5,8},{8,6,7,3},{4,1,8,6}};
   static float u[9][4] = { {4,0,0,0},{0,4,0,0},{0,0,4,0},{0,0,0,4},
			    {2,2,0,0},{2,0,2,0},{0,2,0,2},{0,0,2,2},{1,1,1,1} };
   float tmp[4][4];
   for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
	 tmp[i][j] = 0.;
	 for(k=0;k<4;k++) {
	   tmp[i][j] += u[indx[subbox][i]][k]*a[k][j];
	 }
      }
   }
   for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
	 a[i][j]=0.25*tmp[i][j];
      }
   }
   return;
}
