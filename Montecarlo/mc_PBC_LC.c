#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <string.h>


/* -------------------------------------------------------------------------- */
/* PROGRAM CONSTANTS */
/* -------------------------------------------------------------------------- */

// DEPRECATED CONSTANTS DEFINITION
// #define NP 1000
// #define LBOX 1.0
// #define DELTAR LBOX/50.0
// #define T 1.0
// #define NSTEP 100.0
// #define NSAMP 1000.0
// #define ENERGY 1.0
// #define SIGMA 0.05
// #define CUTOFF 2.5*SIGMA

static const float LBOX = 1.0;
static const float LCELL = 1.0/4;
static const int NP = 1000;
static const int NCELL = 64;
static const float DELTAR = 1.0/50.0;
static const float T = 1.0;
static const float NSTEP = 100.0;
static const float NSAMP = 1000.0;
static const float EPSILON = 1.0;
static const float SIGMA = 0.1;
static const float CUTOFF = 0.1*2.5;

/* -------------------------------------------------------------------------- */
/* FUNCTIONS */
/* -------------------------------------------------------------------------- */

float start_pos();
float get_dist( float pos1[3], float pos2[3] );
int find_cell_atom_i( float pos[3] );
float calc_energy_particle( float r[3], float pos[NP][3], int list[NP], int head[NCELL] );
float calc_LJpot( float r[3], float pos[3] );
float prob( float energy );
float step_correction ( float deltar, float ratio );
void write_pos_to_file( float pos[NP][3], char filename[20] );
float apply_pbc( float r );

/* -------------------------------------------------------------------------- */
/* MAIN BODY */
/* -------------------------------------------------------------------------- */

main(){

  float deltaR = DELTAR;
  float init_pos[NP][3], pos[NP][3], new_pos[NP][3];
  float rold[3], rnew[3];
  float chi;
  float enerold, enernew, ratio;
  float pold, pnew;
  float naccept, ntry;
  int neighbours_old[NP-1], neighbours_new[NP-1];
  int nsamp = (int)NSAMP;
  int posit;
  int list[NP];
  int head[NCELL];
  int icel;
  int counting = 0;
  float tot_en;

  srand(time(NULL)); /* Change seed to ensure new results every run */

  /* CALL INITIAL CONDITIONS SUBROUTINE --> READ FROM FILE (T, LBOX, NSAMP, NSTEP, SIGMA[Particle radius]?) */

  /*   /|
      / |
        |
        |
       _|_  PLACE PARTICLES --> INITIALIZE STARTING POSITIONS; INITIALIZE NEIGHBOUR LISTS */

  for (int part = 0; part < NP; part++){
    for (int dim = 0; dim < 3; dim++){
      init_pos[part][dim] = start_pos();
      pos[part][dim] = init_pos[part][dim];
      new_pos[part][dim] = init_pos[part][dim];
    }
  }

  write_pos_to_file(init_pos,"first_pos.dat");

  /* for (int part = 0; part < NP; part++){
    printf("%f\t%f\t%f\n",pos[part][0],pos[part][1],pos[part][2]);
  } Works perfectly */

  for (int i = 1; i <= NP; i++){ /* Initialize list vector */
    list[i] = 0;
  }

  for (int i = 1; i <= NCELL; i++){ /* Initialize head vector */
    head[i] = 0;
  }

  for (int i = 1; i <= NP; i++){
    icel = find_cell_atom_i(pos[i-1]);
    // printf("Atom %i is placed in cell %i\n",i,icel); It works but *TAKE CARE OF INDEXES
    list[i] = head[icel];
    head[icel] = i;
  }

  /*  ______
            |
       _____|
      |
      |______   DO MONTECARLO LOOP --> CALL ENERGY SUBROUTINE ON NSAMP STEP */

  for (int iter = 0; iter < NP*NSTEP; iter++){ /* MC start */

    /* _____
            |
        ____|
            |
       _____|  SELECT RANDOM PARTICLE AND RANDOM STEP */

    posit = random();
    posit = posit % NP; /* Random position selected for movement */
    /* printf("\nPARTICLE %i\n",posit); Check randomly-selected particle */

    /* Random step is done in each dimension */
    for (int dim = 0; dim < 3; dim++){
      chi = rand();
      chi = 2*(chi/RAND_MAX - 0.5);

      rold[dim] = pos[posit][dim];
      rnew[dim] = rold[dim] + chi*deltaR;

      new_pos[posit][dim] = rnew[dim];
    }

    /* printf("%f\t%f\t%f\t\t\t%f\t%f\t%f\n",rold[0],rold[1],rold[2],rnew[0],rnew[1],rnew[2]); Comparison between rold-rnew --> WORKS FINE */
    /* float dist = get_dist(rold,rnew);
    printf("%f\n",dist); */


 /*  |
     |     |
     |_____|
           |
           |  CALCULATE ENERGY OF NEIGHBOURHOOD AND PROBABILITIES */




    /* Lennard-Jones potential --> energy*((sigma/r)^12 - sigma/r)^6) */

    enerold = calc_energy_particle( rold, pos, list, head );
    enernew = calc_energy_particle( rnew, new_pos, list, head);

    pold = prob(enerold);
    pnew = prob(enernew);


    /* CHECK PROBABILITY OF TRANSITION */
    if (pnew/pold > 1){
      ratio = 1;
    } else {
      ratio = pnew/pold;
    }



    /*   _____
        |
        |_____
              |
         _____|  CHECK PROBABILITIES, TAKE/NOT STEPS, AND UPDATE SYSTEM */


    chi = rand();
    chi = chi/RAND_MAX;
    if (chi < ratio){
      for (int dim = 0; dim < 3; dim++){
        rold[dim] = rnew[dim];

        /* APPLY PBC; otherwise, errors may occur when representing 3d scatterplots of particles positions */
        rold[dim] = apply_pbc(rold[dim]);
      }
      naccept++;
      /* printf("TRANSITION ACCEPTED WITH P = %f\n",ratio);
      printf("ENERGY OLD: %f \t ENERGY NEW: %f\n",enerold,enernew);
      printf("PROB. OLD: %f \t PROB. NEW: %f\n",pold,pnew); */
    } else {
      for (int dim = 0; dim < 3; dim++){
        rold[dim] = rold[dim];
      }
      /* printf("TRANSITION REJECTED WITH P = %f\n",ratio);
      printf("ENERGY OLD: %f \t ENERGY NEW: %f\n",enerold,enernew);
      printf("PROB. OLD: %f \t PROB. NEW: %f\n",pold,pnew); */
    }
    ntry++;

    /* Need to update particle position list --> done below */
    for (int dim = 0; dim < 3; dim++){
      pos[posit][dim] = rold[dim];
    }


    /* Save trajectories each NSAMP times */
    if (iter > 1 && counting < (int)NSAMP){
      if ( iter % (int)NSAMP == 0){
        deltaR = step_correction(deltaR,naccept/ntry); /* Control step size */

        printf("%3.0f %% completed\n",100*((float)iter+1)/(NP*NSTEP)); /* Percentage until completed MC run */
        printf("\tTOTAL ENERGY OF THE SYSTEM: %f\n",enernew/(NP*2));

        char filename[20];
        sprintf(filename, "%d", iter);
        strcat(filename,"step_pos.dat");
        write_pos_to_file(pos,filename);
        counting++; /* Save trajectories into file */
      }
    }

  } /* MC end */

  /* Write result to file */
  write_pos_to_file(pos,"last_pos.dat");

  printf("\n---------------------------------\nACCEPTANCE RATIO: %f\n\n",naccept/ntry);
}


/* ---------------------------------------------------------------------------- */
/* SUBROUTINES */
/* ---------------------------------------------------------------------------- */

float start_pos() /* Subroutine to generate a random number between -LBOX/2 and LBOX/2 */
{
  float chi = random();
  chi = LBOX*(chi/RAND_MAX - 0.5);
  return chi;
}


float get_dist( float pos1[3], float pos2[3] ) /* Subroutine to get euclidean distance between two particles */
{
  float dist, xdif, ydif, zdif;

  xdif = (pos2[0] - pos1[0]);
  ydif = (pos2[1] - pos1[1]);
  zdif = (pos2[2] - pos1[2]);

  /* Need to apply minimum image convention, so as to keep PBC  */
  xdif = apply_pbc(xdif);
  ydif = apply_pbc(ydif);
  zdif = apply_pbc(zdif);

  dist = sqrt(pow(xdif,2) + pow(ydif,2) + pow(zdif,2));
  /* printf("Dist: %f, x: %f, y: %f, z: %f\n",dist,xdif,ydif,zdif); */
  return dist;
}


int find_cell_atom_i( float pos[3] )
{
  int mx, my, mz;

  mx = (int)(LBOX / LCELL);
  my = (int)(LBOX / LCELL);
  mz = (int)(LBOX / LCELL);
  /* printf("Cell size: %f %f %f\n",mx,my,mz); */

  int icelx, icely, icelz, icel;

  icelx = 1 + (int)((0.5 + pos[0]/LBOX)*mx);
  icely = 1 + (int)((0.5 + pos[1]/LBOX)*my);
  icelz = 1 + (int)((0.5 + pos[2]/LBOX)*mz);

  icel = icelx + (icely-1)*mx + (icelz-1)*my*mx;

  return icel;
}


float calc_energy_particle( float r[3], float pos[NP][3], int list[NP], int head[NCELL] ) // IT WORKS!! FINALLY :,)
{
  int mx, my, mz;
  float energy = 0;

  mx = (int)(LBOX / LCELL);
  my = (int)(LBOX / LCELL);
  mz = (int)(LBOX / LCELL);

  int icelx, icely, icelz, icel, i, jcel, j;
  int jx, jy, jz;

  icelx = 1 + (int)((0.5 + r[0]/LBOX)*mx);
  icely = 1 + (int)((0.5 + r[1]/LBOX)*my);
  icelz = 1 + (int)((0.5 + r[2]/LBOX)*mz);
  icel = icelx + (icely-1)*mx + (icelz-1)*my*mx;
  i = head[icel];

  for (int jcelx = icelx - 1; jcelx < icelx + 2; jcelx++){
    for (int jcely = icely - 1; jcely < icely + 2; jcely++){
      for (int jcelz = icelz - 1; jcelz < icelz + 2; jcelz++){

        jx = jcelx;
        jy = jcely;
        jz = jcelz;

        // Apply PBCs by left boundaries
        if(jx == 0){
          jx = mx;
        }
        if(jy == 0){
          jy = my;
        }
        if(jz == 0){
          jz = mz;
        }

        // Apply PBCs by right boundaries
        if(jx == mx+1){
          jx = 1;
        }
        if(jy == my+1){
          jy = 1;
        }
        if(jz == mz+1){
          jz = 1;
        }

        jcel = jx + (jy-1)*mx + (jz-1)*my*mx;
        j = head[jcel];

        if (j != 0 && j != i){ // Avoid comparing a particle with itself
          while (j != 0){
            if (get_dist(pos[i-1],pos[j-1]) < CUTOFF){
              energy = energy + calc_LJpot(pos[i-1],pos[j-1]);
            }
            j = list[j];
          }
        }

      }
    }
  }

  return energy;
}


float calc_LJpot( float r[3], float pos[3] ) /* Subroutine to get Lennard-Jones potential between two interacting particles */
{
  float LJpot = 4 * EPSILON * ( pow(SIGMA/get_dist(r,pos),12) - pow(SIGMA/get_dist(r,pos),6) );
  return LJpot;
}


float prob( float energy ) /* Function to compute probability of a given state with a given energy according to a Boltzmann distribution */
{
  float p = exp(-energy/T);
  return p;
}


float step_correction ( float deltar, float ratio )
{
  if (ratio < 0.5){
    deltar = deltar/1.01;
  }else{
    deltar = deltar*1.01;
  }
  return deltar;
}


void write_pos_to_file( float pos[NP][3], char filename[20] )
{
  FILE *fp;
  fp = fopen(filename,"w+");
  for (int p = 0; p < NP; p++){
    /* Apply PBC --> Fold back to primary box if particle leaves it
    for (int dim = 0; dim < 3; dim++){
      pos[p][dim] = pos[p][dim] - (int)((pos[p][dim]/LBOX)+0.5)*LBOX;
    } */
    fprintf(fp, "%f\t%f\t%f\n", pos[p][0], pos[p][1], pos[p][2] );
  }
  fclose(fp);
}


float apply_pbc( float r )
{
  r-=(int)( ( (r/LBOX<0)?-0.5:0.5) + r)*LBOX;
  // r-=(int)( ( (r<0)?-0.5:0.5 ) + r); /* TRuncates correctly negative float values */
  /* r = r - (int)((r/LBOX)+0.5)*LBOX; */
  return r;
}
