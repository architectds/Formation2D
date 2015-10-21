/******************************************************************************
 * PIHM-RT is a finite volume based, reactive transport module that operate
 * on top of the hydrological processes described by PIHM. PIHM-RT track the 
 * transportation and reaction in a given watershed. PIHM-RT uses operator
 * splitting technique to couple transport and reaction. 
 *
 * PIHM-RT requires two additional input files: 
 *     a. chemical condition file:     projectname.chem
 *     b. index of initial conditions: projectname.cini
 *
 *
 *
 * If you have any questions, concerns, suggestions, please contact me at 
 * the following address:
 *
 *     Developer: Chen Bao <baochen.d.s@gmail.com>
 *     Version  : 0.2
 *     Date     : Feb, 2014   
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>



#define UNIT_C 1440
#define ZERO   1E-20
#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80
#define TIME_MIN  1E-5
#define EPS       0.05
#define INFTYSMALL  1E-6
#define RTdepth 5.0
#define MIN(a,b) (((a)<(b))? (a):(b))
#define MAX(a,b) (((a)>(b))? (a):(b))

/* Functions declarations and usage */
int  keymatch(const char *, const char * , double * , char ** );
static double timer();
int  realcheck(const char* words);
void ConditionAssign(int condition, char* str, int * index);
/* Fucntion declarations finished   */

// Timer

static double timer() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}


int realcheck(const char* words){

  int flg = 1, i;
  if ((( words[0] < 58)&&(words[0] > 47))||(words[0] == 46)||(words[0]==45)||(words[0]==43)){
    for ( i = 0 ; i < strlen(words); i ++)
      if ( (words[i] >57 || words[i] < 43) && (words[i]!=69) && ( words[i]!=101) &&(words[i]!=10) &&(words[i]!=13))
	flg = 0;
  }
  else flg = 0;
  return(flg);
}


int keymatch(const char * line, const char * keyword, double * value, char ** strval){
  /* A very general and convinient way of reading datafile and input file */
  /* find keyword in line, assign the value after keyword to value array if there is any */
  /* store both numbers and strings in order for later use, buffer required */
  /* if is keyword not found return 0. If comments, return 2. Otherwise return 1 */
  int i;

  for(i = 0; i <WORDS_LINE; i++) value[i] = 0.0;

  if ((line[0] =='!')|| (line[0] =='#')){

    return(2);
    /* assign a special flag for comments */
  }

  int j, k, line_width, word_width = WORD_WIDTH, quoteflg = 0;
  int words_line = WORDS_LINE;
  int keyfoundflag = 0;

  char ** words;
  words = (char **) malloc ( WORDS_LINE * sizeof(char* ));
  
  for (i = 0 ; i < WORDS_LINE ; i++){
    words[i] = (char*) malloc ( WORD_WIDTH * sizeof(char));
    memset( words[i], 0, WORD_WIDTH);
  }
  i = j = k = 0 ;

  /* Partition the line into words */
  while ( i < strlen(line)){
    if(line[i]!=39){
      while (line[i]!=9 && line[i]!=0 && line[i]!=10 && line[i]!=32 && line[i]!=13){
	words[k][j++] = line[i++];
	if ( line[i]==9 || line[i]==32 || line[i] ==13){
	  k ++;
	  j =0;
	}
      }
    }
    else{
      words[k][j++]= line[i++];
      while(line[i]!=39){
	words[k][j++] = line[i++];
      }
      words[k++][j] = line[i++];
      j = 0;
    }
    i ++;
  }

  words_line = k + 1;

  for ( i = 0 ; i < words_line; i ++)
    if (strcmp(words[i], keyword) == 0) 
      keyfoundflag = 1;

  j = k = 0;
  for ( i = 0; i < words_line; i ++){
    //    fprintf(stderr, "word#%d=%s, length=%d\n" , i, words[i], strlen(words[i]));
    strcpy(strval[k++],words[i]);
    //    if ((( words[i][0] < 58)&&(words[i][0] > 47))||(words[i][0] == 46)||(words[i][0]==45)||(words[i][0]==43))
    if (realcheck(words[i]) == 1)
      value[j++] = atof(words[i]);
  }
  
  for (i = 0 ; i < WORDS_LINE ; i++)
    free(words[i]) ;
  free(words);
   return(keyfoundflag);

}

void ConditionAssign(int condition, char* str, int * index){
  /* This subroutine takes in input strings and output an index array that record the conditions each blocks assigned to */
  /* input strings could use separators like - and , */

  int i, j , k , l, length = strlen(str);
  char ** words = (char **) malloc( length* sizeof(char*));
  for ( i = 0; i < length; i++){
    words[i]= (char* ) malloc ( WORD_WIDTH * sizeof(char));
    memset( words[i], 0, WORD_WIDTH);
  }
  char * tmpstr    = (char*) malloc(WORD_WIDTH * sizeof(char));

  char * separator = (char*) malloc( length * sizeof(char));
  int  * value     = (int *) malloc( length * sizeof(char));



  i = j = k = l = 0;
  while (i < length){
    while (str[i]!= 0 && str[i] !=45 && str[i] != 44 ){
      //      fprintf(stderr, "%s %c\n",words[k], str[i]);
      words[k][j++] = str[i++];
      if ( str[i]==45 || str[i]==44){
        k ++;
        j =0;
	separator[l++] = str[i];
      }
    }
    i++;
  }
  
  for ( i = 0; i <= k; i ++){
    strcpy(tmpstr, words[i]);
    value[i] = atoi(tmpstr);
    //    fprintf(stderr, " tset %d %s\n", value[i], words[i]);
  }
  /* 
  for ( i = 0; i < length; i ++)
    fprintf(stderr, "%s\t%c\n",words[i],separator[i]);

  for ( i = 0; i < length; i ++)
    fprintf(stderr, "%d\n",value[i]);
 
  fprintf(stderr, "condition = %d\n", condition);
  */
  for ( i = 0; i < l ; i ++){
    if( separator[i] == ',' || separator[i] == 0){
      index[value[i]] = condition;
      fprintf(stderr, " CS %d: %d\n", value[i], index[value[i]]);
    }
    if( separator[i] == '-')
      for ( j = value[i]; j <= value[i+1]; j ++){
	index[j] = condition;
	//	fprintf(stderr, " CS %d: %d\n", j, index[j]);

      }
  }
  
  for (i = 0 ; i < length ; i++)
    free(words[i]) ;
  free(words);
  free(separator);
  free(value);
  free(tmpstr);
  
}
