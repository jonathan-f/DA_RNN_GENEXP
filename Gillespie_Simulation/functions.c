///////////////////////////////////////////////
//
//
// 5/2019
//
//   Title: Functions of Gillespie algorithm
//
//	Author: Michele
//
//
/////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "functions.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

double pi = 3.1415926;

double Hill_activation(double gene, double K, double k, double h, double treshold){
/*
this function returns the hill function for an activation function of hill kind plus treshold
*/

return k*(pow(gene,h)/(pow(gene,h)+ pow(K,h))) + treshold;

}

double Hill_repression(double gene, double K, double k, double h, double treshold){
/*
this function returns the hill function for repression function of hill kind plus treshold
*/

return k*(1-pow(gene,h)/(pow(gene,h)+ pow(K,h))) + treshold;

}

// select the stochiometry of the system:
// stochio[i][j] is the stochiometric coefficient of specie j in the reaction i

void fill_stochio( int **stochio, int n, int m){

int i;
//general stochiometric matrix for a gene regulatory netwotk. From 0 to n you have the gene epression, n - 2n gene stop exressing,
// 2n - 3n protein translation, 3n - 4n protein dilution
  for(i=0;i<n;i++){

  // expression reactions gene from 0 to n
stochio[i][i]=1;

  // stop expresison of genes from n to 2 n
  stochio[i+n][i]=-1;

  // protein translation
  stochio[i+2*n][i]=1;

  //protein Dilution

  stochio[i+3*n][i]=-1;



  }
}

int step_function(int x){
  if(x>0)return 1;
  if(x<=0)return 0;
}

void fill_rates( double *rate, double k, double mu, int n){

// simpler way to fill the reaction rates

int i;
  for(i=0;i<n;i++){

//ecpression reaction rates
rate[i]= k;

// stop gene expression reactions
rate[i+n] = 0.1*mu;

// protein translation  reaction rates
rate[n+n+i]= 2*k;

// protein Dilution reaction rates
rate[2*n+n+i]= 0.01*mu;



  }

}

void fill_h( double *h,  int m){
int i;
// select the hill coefficients array to be fed into the Activation-Repression reactions

  for(i=0;i<m;i++){
h[i]= 2;

  }
h[0]=10;
}

void fill_treshold( double *treshold, int m){
int i;
// select the trechold array to feed the Activation-repression functions of the gene expression
  for(i =0;i<m;i++){
  treshold[i]=0;
  }

}

void selectK(int *K, int m){
int i;

for(i=0;i<m;i++){
K[i]=50;
}
		}

    // RANDOM generators of parameters
void fill_rates_RND( double *rate, double k, double mu, int n, int n1){

    // think a way to fill the rate in a smart way

        gsl_rng * r;
        const gsl_rng_type * T;
        T = gsl_rng_mt19937;
        r = gsl_rng_alloc(T);
        //setting the seed to the rng
        gsl_rng_set(r, 3*time(0));
        srand48(time(0));

    int i;
    double x, sig, k1;
    sig = k/10.;
      for(i=0;i<n1;i++){
    //ecpression reaction rates from 0 to n1

    k1 = k+ gsl_ran_gaussian(r, sig);

    rate[i]= k;

    //translation reactions from n1 to n1+n

    rate[i+n] = 100*k;

    // gene Dilution reaction rates from n1 + n to n + 2*n1
    rate[n1+n+i]= mu;
    // Protein dilution rates from 2* n1 + n to 2* n1 + 2* n
    if(i<n)rate[2*n1+n+i]= mu;

    // active degradation genes from 2* n1 + 2*n to 3* n1 + 2* n
    rate[2*n+2*n1+i] = 1*k;

      }

    }

void fill_h_RND( double *h,  int m){
    int i;
    // select the hill coefficients array to be fed into the Activation-Repression reactions

      for(i=0;i<m;i++){
    h[i]= 2;

      }
    h[0]=10;
    }

void fill_treshold_RND( double *treshold, int m){
    int i;
    // select the trechold array to feed the Activation-repression functions of the gene expression
      for(i =0;i<m;i++){
      treshold[i]=0;
      }

    }

void selectK_RND(int *K, int m){
    int i,x;

    gsl_rng * r;
    const gsl_rng_type * T;
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);
    //setting the seed to the rng
    gsl_rng_set(r, 3*time(0));
    srand48(time(0));

    for(i=0;i<m;i++){
      x=((int)lrand48()%100);

      K[i]=x;
    }
    		}

void propensityJij(int *Kpg, double *pro, int *gene,  int *protein,  double *rate, double *hpg, double *treshold, int **Jij_PG, int n){

    // Expresison DNA reactions from 0 to n

    int i,j;
    // protein - gene network: n genes (gene promoters) up to n Transcription factors (protein involved)
    for(i=0;i<n;i++){
      pro[i]=0;
      //if(i>3){
      for(j=0;j<n;j++){
        // random expression mechanism
      if(Jij_PG[i][j]==1){
        if(pro[i]!=0)pro[i]*= Hill_activation(protein[j], Kpg[i], 1, hpg[i], treshold[i]);
        if(pro[i]==0)pro[i]= Hill_activation(protein[j], Kpg[i], 1, hpg[i], treshold[i]);
    }
      if(Jij_PG[i][j]==-1){
        if(pro[i]!=0)pro[i]*= Hill_repression(protein[j], Kpg[i], 1, hpg[i], treshold[i]);
        if(pro[i]==0)pro[i]= Hill_repression(protein[j], Kpg[i], 1, hpg[i], treshold[i]);
    }
    }

      pro[i]=pro[i]*rate[i];

      // gene Dilution from n+ n1 to n + 2*n1
      pro[i+n]= rate[i+n]*gene[i];

}

// Translation reactions from n1 to n+n1

// Tranlastion reaction gene --> Protein from n1 to n1+n


for(i=0;i<n;i++){

pro[i+2*n]=rate[i+2*n]*gene[i];


}


//Protein Dilution

for(i=0;i<n;i++){

pro[i+3*n]=rate[i+3*n]*protein[i];


}


// expression dynamic of the core clock
      pro[0]= Hill_repression(protein[3], Kpg[0], rate[0], hpg[0], treshold[0]);
      pro[1]= Hill_activation(protein[0], Kpg[0], rate[1], hpg[1], treshold[0]);
      pro[2]= Hill_activation(protein[1], Kpg[0], rate[2], hpg[2], treshold[0]);
      pro[3]= Hill_activation(protein[2], Kpg[0], rate[3], hpg[3], treshold[0]);

    }

void interactionMatrix_RND(int **Jij, int n, int n1, double reps, double acts){
// Matrix n1 x n

  //reps percetnage of repression links
  //acts percetnage of activation links

      gsl_rng * r;
      const gsl_rng_type * T;
      T = gsl_rng_mt19937;
      r = gsl_rng_alloc(T);
      //setting the seed to the rng
      gsl_rng_set(r, 3*time(0));
      srand48(time(0));

double x;

int i, j;
for(i=0;i<n1;i++){
  for(j=0;j<n;j++){
    x =  gsl_rng_uniform(r);
    if(x<reps)Jij[i][j] = - 1;
    if(x>reps && x<reps+acts)Jij[i][j] = 1;
    }
  }

}

void interactionMatrix_RNDexp(int **Jij, int n, int n1, double reps, double acts){

//TO BE EDITED

//reps percetnage of repression links
//acts percetnage of activation links

      gsl_rng * r;
      const gsl_rng_type * T;
      T = gsl_rng_mt19937;
      r = gsl_rng_alloc(T);
      //setting the seed to the rng
      gsl_rng_set(r, 3*time(0));
      srand48(time(0));

double x;

int i, j;
for(i=0;i<n1;i++){
  for(j=0;j<n;j++){
    x =  gsl_rng_uniform(r);
    if(x<reps)Jij[i][j] = - 1;
    if(x>reps && x<reps+acts)Jij[i][j] = 1;
  }
}

}

void interactionMatrix_fromFIle(int **Jij, int n,int n1, char filename[85]){

  /*
  Jij: matrix to be read n1 x n
  filename: input file of ints, n x n1 numbers that will be inserted in Jij

  */

    FILE *file;

    file=fopen(filename,"r");
    if (file==NULL)
        {
            printf("no such file.");
        }
    double x;

int i, j;
for(i=0;i<n1;i++){
  for(j=0;j<n;j++){
    x =  fscanf(file,"%d",&Jij[i][j]);

  }
}

}

void saveMatrix_int(int **matrix, char filename[65], int m, int n){
    /*
     Save the matrix n*m of int on a file called filename
     */
// matrix Jij[i][j] i: 0-n j: 0-m
    int i , j;
    FILE *file;

    file=fopen(filename,"w");

    for (i=0; i<n; i++) {

        for(j=0;j<m;j++){
        fprintf(file,"%d\t",matrix[i][j]);
        }
    fprintf(file,"\n");

    }

    fclose(file);

}

void saveMatrix_double(double **matrix, char filename[65], int m, int n){
    /*
     Save the matrix n*m of int on a file called filename
     */

    int i , j;
    FILE *file;

    file=fopen(filename,"w");

    for(i=0;i<m;i++){
        for (j=0; j<n; j++) {
            fprintf(file,"%lf\t",matrix[i][j]);
        }
        fprintf(file,"\n");

    }

    fclose(file);

}

void compute_atot(double *atot, double *pro, int m){
	        int i;
	        (*atot)=0;
        for(i=0;i<m;i++){
            (*atot)= (*atot)+pro[i];
        }
	}

void select_reaction(int *i, double *pro, double atot){

	double x,z;
        x=((double)lrand48()/RAND_MAX)*atot;
        // could be possible also use other algorithms for extracting the reaction, like a cumulative algorithm
        z=0;
        (*i)=0;
        do{
            z+=pro[(*i)];
            (*i)++;

        }while(x>z);

        (*i)=(*i)-1;
}

void upgrade_species( int *gene, int *protein, int **stochio, int i, int n){


//Upgrade specie for the Gillespie System
            int j;

					for(j=0;j<n;j++){
							gene[j]=Upgradegene(stochio[i][j], gene[j]);
              protein[j]= protein[j] + stochio[i][j];
						}
						// in just one step I am upgrading everything!!
}

void check_Overloading(int S, double *t_delay){

	int j;
	int k=0;
 for(j=0;j<S;j++){
			if(t_delay[j]==0){
			j=S;
			//printf("NO Overloading\n");
		k++;

		}
	}
if(k==0) printf("OVERLOADING!!!!!!\n");

}

void print_out(double atot, double t_inc, double t,int i, int t_step){

      if(t_step%100000 ==0){
        printf("%lf\t time inc: %lf\t time: %lf \t fired reaction = %d\t time step = %d \n",atot,t_inc,t,i,t_step-1);
		//printf("atp=%lf\n",atp);
  }

}

int Upgradegene(int x, int y){

  int a;

  if(x>0){
    a =1;
  }

  if(x==0){
    a =y;
  }


  if(x<0){
    a =0;
  }

return a;

}

void printPropensityFormulas(int **Jij_PG,int **Jij_PP, int **Jij_PR, int **Jij_RR, int n, int n1){

    // Expresison reactions from 0 to n1

    int i,j;
    // protein - gene network: n1 genes (gene promoters) up to n Transcription factors (protein involved)
    for(i=0;i<n1;i++){
printf("Reaction: %d \n",i);      //if(i>3){
      for(j=0;j<n;j++){
        // random expression mechanism
      if(Jij_PG[i][j]==1){
        printf("p[%d]/(p[%d]+Kp[%d]) ",j,j,i);//Hill_activation(protein[j], Kpg[i], 1, hpg[i], treshold[i]);

    }
      if(Jij_PG[i][j]==-1){
        printf("Kp[%d]/(p[%d]+Kp[%d]) ",i,i,j);//Hill_repression(protein[j], Kpg[i], 1, hpg[i], treshold[i]);
    }
    }

      printf("\n");
}



for(i=0;i<n1;i++){
  printf("Reaction: %d \n",i+n1);
  if(i<n)printf("rate[%d] r[%d] ",i+n1,i);
  for(j=0;j<n1;j++){
      // random binding mechanism protein gene, it can enhance or depress the rate of translation. Generally repress.
      if(j<n){
    if(Jij_PR[i][j]==1){
      printf("pr[%d]/(pr[%d]+Kpr[%d]) ",j,j,i+n1);
      }
    if(Jij_PR[i][j]==-1){
      printf("Kpr[%d]/(pr[%d]+Kpr[%d]) ",i+n1,i+n1,j);\
      }
}
      // random binding mechanism gene gene, it can enhance or repress the rate of translation. Generally repress, like microgene effects

    if(Jij_RR[i][j]==1){
      printf("r[%d]/(r[%d]+Kr[%d]) ",j,j,i+n1);
    }
    if(Jij_RR[i][j]==-1){
      printf("Kr[%d]/(r[%d]+Kr[%d]) ",i+n1,i+n1,j);
    }
  }

 printf("\n");

}
    for(i=0;i<n1;i++){
      printf("Reaction: %d \n",i+2*n1+n);


          if(i<n){

            // Protein Dilution
            printf("rate[%d]*p[%d]\n ",i+2*n1+n,i);
          }
          // gene Dilution
          printf("Reaction: %d \n",i+n1+n);

          printf("rate[%d]*r[%d] \n",i+n+n1,i);


          //active degradation - constant rate for the gene
          printf("Reaction: %d \n",i+2*n1+2*n);     

          printf("rate[%d]*StepFunction[ r[%d]] \n ",i+2*n1+2*n,i);


}
    }
