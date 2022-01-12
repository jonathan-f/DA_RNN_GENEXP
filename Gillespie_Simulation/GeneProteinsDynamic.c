///////////////////////////////////////////////
//
//
//		5/2021
//
//   Title:  Gillespie Algorithm on a Gene regulatory network comptued over a Random interaction matrix
//	Author: Michele Monti
//
//
/////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "functions.h"
#include <stdbool.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

int main(){

//////////////////////////////////////////////////////////////////////////////////

    // setting variables of the system

    gsl_rng * r;
    const gsl_rng_type * T;
    T = gsl_rng_mt19937;
    //initialising instance of rng
    r = gsl_rng_alloc(T);
    //setting the seed to the rng
    gsl_rng_set(r, 3*time(0));
    srand48(time(0));

    int i,j,t_max, t_step, count, nr, na;
    double k,dt, x,y, atot, t_out, acts, reps;
		double t,t_inc,mu;
    double tgrad ;
    double tt=0.;

    int checkMatrixSaveLoad=0;    // if 1 will load the saved RND matrix

	x=0;
	y=0;

	int n , m;

	n= 10;								// number of Protein expressed




  int *gene;
  gene = (int*)calloc(n, sizeof(int*));

  int *protein;
  protein = (int*)calloc(n, sizeof(int*));

				// number of reactions

  m=4*n;


// Stochio is the stochiometric matrix of the system. stochio[i][j] means how the gene-protein j changes if the reaction i is fired
// important to notice that stochio[i][j] is a n+n1 * m matrix where the first n1 lines represent the stochiometry of the gene and
// and the second n lines represent the stochiometry of the Proteins

int **stochio;
stochio = (int **)calloc(m, sizeof(int*));
for(i = 0; i < m; i++){
    stochio[i] = (int *)calloc(n, sizeof(int));
}

// Jij_PG interaction matrix between the protein (n) and the gene-promoter (n1)
int **Jij_PG;
Jij_PG = (int **)calloc(n, sizeof(int*));
for(i = 0; i < n; i++){
    Jij_PG[i] = (int *)calloc(n, sizeof(int));
}



    // pro(i) is an array of the propensity function of each reaction(i)
    // rate(i) is the array of the rate of the reaction(i)
    // hpg hpr hrr are the array of the hill coefficients
    // treshold array of tresholds for the hill functinos
    // gene array where at each position is encoded the protein concentrations at time t
    // K array of Kinetic contants of the Hill functions
    // isto array of the histogra of the reacitons

    double *pro;
    pro = (double*)calloc(m, sizeof(double*));
		double *rate;
    rate= (double*)calloc(m, sizeof(double*));

// think one second if it is better to introduce multiple h array of the hill coefficients of the different kind of reactions.
// hpp protein-protein binding
// ppg protein gene binding
// hpr proteion gene binding

    double *hpg;
		hpg= (double*)calloc(m, sizeof(double*));


    double *h;
    h= (double*)calloc(m, sizeof(double*));
		double *treshold;
		treshold= (double*)calloc(m, sizeof(double*));

	// K is the treshold of the gene expression reactions

    int *Kpg;
    Kpg = (int*)calloc(m, sizeof(int*));



    // isto array to fill with the frequency of having a reaction through th simulation
    int *isto;
    isto = (int *)calloc(m, sizeof(int*));

    // gene <--- array of number of proteins for each genes

    FILE *fileP;

    FILE *fileR;

    FILE *istog;

  char filenameP[85];
  char filenameR[85];
	char istoname[85];

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initialization of the system:

    //initial conditions
for (i=0; i< n;i++) {
gene[i]=1;
protein[i]=2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //set the algorithm and output precision parameters parameters

// final integration time
  t_max=1000;

// counter of time steps initialized to 0
  t_step=0;

// initialized orinting out time step
	dt=0;

// starting time
  t=0.;

// frquency on how you update the extegenel input
  tgrad = 0.1;

// frequency on which you print out
	t_out=0.5;

// counting variable
	count=0;

//initialization of the input variable starting time
  tt=0;

// here you can start to loop over the biochemical parameters

// Biochemical Parameters: na, nr #Activator and Repressor on average per genes

  // nr=(int)n/2;
  // na=(int)n/2;
  // reps = (double)nr/n;
  // acts = (double)na/n;

    nr=3;
    na=3;
    reps = (double)nr/(n);
    acts = (double)na/(n);

// mu dilution rate, k base rate of expression

  mu = .1;
	k=.1;

// Output file names
    sprintf(filenameP,"Protein_N_%d_nr%d_na%d.txt", n,nr,na);

    sprintf(filenameR,"gene_N_%d_nr%d_na%d.txt", n,nr,na);

  	sprintf(istoname,"Histo_reactions_N_%d_nr%d_na%d.txt", n, nr,na);

      fileP=fopen(filenameP,"w");

      fileR=fopen(filenameR,"w");

      istog=fopen(istoname, "w");

  // set the stochiometry of the system and relative rates

  interactionMatrix_RND(Jij_PG, n,n, reps, acts);


  fill_stochio( stochio, n,m);

	  // compute the propensity function

		//fill_rates_RND(rate, k, mu,n,n);

fill_rates(rate, k, mu, n);


		fill_h(h, m);

		fill_treshold(treshold,m);

    selectK(Kpg,m);


    propensityJij(Kpg, pro, gene,  protein, rate, hpg, treshold, Jij_PG,n);



    compute_atot(&atot,pro,m);


        if(n<11){

    for(i=0;i<m;i++){
      if(i==0)printf("gene Expression Reactions \n");
      if(i==n)printf("protein translation Reactions \n");
      if(i==n+n )printf("genes stop Expression Reactions \n");
      if(i==n+2*n)printf("Proteins dilution Reactions \n");

      printf("pro[%d] = %lf \t rate[%d] = %lf \n",i,pro[i],i,rate[i]);
//printf("random Hill = %lf \n rate[%d]= %lf",Hill_activation(gene[i%n], K[i], 1, h[i], treshold[i]), i,rate[i] );
    }





  printf("J_PG matrix\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++){
  printf("%d\t",Jij_PG[i][j]);
      }
      printf("\n");
    }




    printf("Stochio matrix\n");
        for(i=0;i<m;i++){
          for(j=0;j<n;j++){
    printf("%d\t",stochio[i][j]);
        }
        printf("\n");
      }


      }

// Save the Jij_PG matriix of interaction protein genes
      char filematrix1[65];

      sprintf(filematrix1,"InteractionMatrix_JijPG.txt");

      saveMatrix_int(Jij_PG, filematrix1, n, n);



if(checkMatrixSaveLoad == 1){
    printf("input matrix\n");

        for(i=0;i<n;i++){
          for(j=0;j<n;j++){
    printf("%d\t",Jij_PG[i][j]);
        }
        printf("\n");
      }
    interactionMatrix_fromFIle(Jij_PG,  n,n,  filematrix1);

    printf("Loaded matrix\n");
        for(i=0;i<n;i++){
          for(j=0;j<n;j++){
    printf("%d\t",Jij_PG[i][j]);
        }
        printf("\n");
      }

    }



    //printing on the data file all the initial conditions



printf("Architecture of the GRN:\n Number of Genes : %d \n  Number of Average reprssors and Activators : %d , %d \n",n,na, nr );
printf("Probability of having an Activating or Repressing node in the Random Matrixes : %lf , %lf \n",acts, reps );
printf("Total Simulation time: %d\n",t_max);

    printf("start!!\n");

// time integration of Gillespie algorithm
    do{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        t_step++;

		//compute the level of phosphorylation that you need for the propensity function

		//propensity2(K, pro, gene,  rate, h, treshold);
    //propensityJij(K, pro, gene,  rate, h, treshold, Jij, n);
    propensityJij(Kpg,  pro, gene,  protein, rate, hpg,  treshold, Jij_PG,  n);

        // computation of the total propensity function

    compute_atot(&atot, pro, m);

        // select time

    t_inc=  gsl_ran_exponential (r, 1/atot);

// compute the extegenel input signal
//
    /*
			if(tt+tgrad< t +t_inc) {

				tt=tt+tgrad;
				t = tt;
// extegenel input behaviour here
			}
			//else{
      */
		// Select reaction

        select_reaction(&i,pro,atot);

		// Upgrade genes and proteins


		upgrade_species(gene, protein, stochio,  i,  n);
    for(j=0;j<n;j++){

      if(gene[j]<0){
        printf("gene[%d] <0 was fired %d\n",j,i);
        break;
}
    }
			t=t+t_inc;

			isto[i]++;

		//}

  // printing data in the output file selecting the precision of the graph (basically I have more than 10^6 point so I can print out just 1 per 10)

  print_out(atot,t_inc,t,i,t_step);


  // print on the files with the precision that yoiu have to set on t_step%precision, or print each dt, better if you want to make autocorrelation,
  // in fact in that cae you need to have a fixed time step in your array of data


// print on file
if(t>dt){
  //  printf("%d\n",step_function(gene[1]));
        //printing the protein concentration behaviour
    do{

      fprintf(fileP,"%lf\t",dt);
      fprintf(fileR,"%lf\t",dt);
        for(j=0;j<n+n;j++){
            if(j<n)fprintf(fileR,"%d\t",gene[j]);

            if(j>=n)fprintf(fileP,"%d\t",protein[j-n]);

        }
        fprintf(fileP,"\n");
        fprintf(fileR,"\n");

            dt=dt+t_out;
        }while(t>dt);
}

		}while(t<t_max);

    printf("finish!!\n total time steps: %d\n time %lf \n",t_step,t);

    for(j=0;j<m;j++){

        fprintf(istog,"%lf\n",(double)isto[j]/t_step);

    }

    fclose(istog);

    fclose(fileP);

    fclose(fileR);

//here you need to end the loops over the parameters


    return 0;

    free(h);
    free(hpg);


		free(treshold);
    free(Kpg);

    free(pro);
    free(rate);
    free(stochio);


    free(Jij_PG);

    free(gene);
    free(isto);

}
