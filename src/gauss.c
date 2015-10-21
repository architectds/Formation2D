/*Gauss Elimination */
# include <stdio.h>
//# include <conio.h>
# include <math.h>
# include <stdlib.h>
# define MAX 3
int main(){
  FILE * infile = fopen("input/gauss.in","r");
  int n,i,j, err;
  fscanf(infile,"%d",&n);
  fprintf(stderr, " Matrix of size %d read from input.\n", n);
  double b[n], x[n], temp;
  double ** a = (double**) malloc( n * sizeof(double*));
  for (i = 0; i < n ; i ++)
    a[i] = (double*) malloc( n* sizeof(double));
  for (i = 0; i < n ; i ++){
    for ( j = 0; j < n ; j ++){
      fscanf(infile, "%lf", &a[i][j]);
      //fscanf(infile, "%lf", &temp);
      //fprintf(stderr, "%4.2f\t", temp);
      fprintf(stderr, "%4.2f\t", a[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "Right Hand Side Vector is :\n");
  for (i = 0 ; i < n; i ++){
    fscanf(infile, "%lf", &b[i]);
    fprintf(stderr, "%6.4f\t", b[i]);
  }
  fprintf(stderr, "\n");
  err = gauss(a,b,n  );

  if (err == 1) fprintf(stderr, " Gaussian Elimination Failed!\n");
  if (err == 0){
    fprintf(stderr, " Gaussian Elimination Succeeded!\n");
    for (i = 0 ; i < n; i ++)
      fprintf(stderr, "%6.4f\t", b[i]);
    fprintf(stderr, "\n");
  }
  return(0);
}



int gauss(double ** a, double * b, int n)
{
  /* solving ax = b type linear equation using gaussian elimination */
  int i,j,k,row, nrows = n;
  double temp,ainv,sum=0, ratio;
  double ** A = a;
  double * B = b;
  double * X = (double*) malloc(n*sizeof(double));
  for ( i = 1; i<nrows; i++) {
    ainv = 1.0 / A[i-1][i-1];
    for ( row = i; row<nrows; row++) {
      ratio = A[row][i-1] * ainv;
      B[row] = B[row] - ratio * B[i-1];

      for ( j = 0; j<nrows; j++) {
        A[row][j] = A[row][j] - ratio * A[i-1][j];
      }
    }
  }
  for ( i = 0; i < nrows; i ++) if (A[i][i] == 0) return(1);
  X[nrows-1] = B[nrows-1] / A[nrows-1][nrows-1];
  for ( i = nrows-2; i >= 0 ; i--) {
    //int i = nrows - 2 - ii;
    ainv = 1.0 / A[i][i];
    X[i] = B[i] * ainv;
    for ( j=i+1; j < nrows; j++) {
      X[i] = X[i] - A[i][j] * X[j] * ainv;
    }
  }
  for (i = 0 ; i < nrows; i++)
    b[i] = X[i];
  return(0);
}
