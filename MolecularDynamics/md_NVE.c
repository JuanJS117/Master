#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

/* -------------------------------------------------------------------------- */
/* MOLECULAR DYNAMICS CONSTANTS */
/* -------------------------------------------------------------------------- */

static const double SIGMA = 1.0;              // Particle diameter
static const double EPSILON = 1.0;            // Energy unit
static const double MASS = 1.0;               // Particle mass
//static const double LBOX = 10.0;            // Box size --> Always set LBOX to be 10 times SIGMA (LBOX = 10*SIGMA) if working with NP = 500
static const double LBOX = 7.368063;          // Box size --> Always set LBOX to be around 7.36 times SIGMA (LBOX = 7.36*SIGMA) if working with NP = 200
//static const double VOLBOX = 10.0*10.0*10.0;// Box volume --> If working with NP = 500 and density = 0.5
static const double VOLBOX = 400.0;           // Box volume --> If working with NP = 200 and density = 0.5
//static const int NP = 500;                  // Number of particles --> Always set NP = 500, in order to reach a density of 0.5 with LBOX = 10*SIGMA
static const int NP = 200;                    // Number of particles --> Always set NP = 200, in order to reach a density of 0.5 with LBOX = 7.36*SIGMA
static const double T = 2.0;                  // Temperature --> To be reached by Langevin thermostat (This value ensures a T = 2.0)
static const double dt = 0.003;               // Time step (fs)
static const double NSTEP = 10000.0;           // Number of steps --> For a simulation of 1 ps, use NSTEP = 1000000
static const double NSAMP = 1000.0;            // Number of steps after which we store an image of the system --> Usually NSAMP = NSTEP/100
static const double CUTOFF = 1.0*2.5;         // Lennard-Jones potential cutoff --> Always set it to be 2.5*SIGMA
static const double UCUT = -0.01631689;       // Lennard-Jones potential at the cutoff distance --> Needed for truncated Lennard-Jones potential
static const double FCUT = -0.01559979;       // Lennard-Jones force at the cutoff distance --> Needed for truncated Lennard-Jones force
static const double PI = 3.14159265359;

/* -------------------------------------------------------------------------- */
/* FUNCTIONS */
/* -------------------------------------------------------------------------- */

void print_welcome();                                                           // Friendly message at the start
double norm( double pos[3] );                                                   // Calc the euclidean norm of a vector
double Xdif( double pos1[3], double pos2[3] );                                  // Calc the x difference between two points
double Ydif( double pos1[3], double pos2[3] );                                  // Calc the y difference between two points
double Zdif( double pos1[3], double pos2[3] );                                  // Calc the z difference between two points
double get_dist( double pos1[3], double pos2[3] );                              // Calc the distance between 2 particles
double calc_LJpot( double r[3], double pos[3] );                                // Calc the Lennard-Jones potential between 2 particles
double calc_LJforce( double r[3], double pos[3] );                              // Calc the Lennard-Jones force between 2 particles
double get_force( double pos[NP][3], int i );                                   // Get all the internal forces acting on a particle
double apply_pbc( double r );                                                   // Apply Periodic Boundary Conditions to a point in space
void save_trajectory( double pos[NP][3], double vel[NP][3], double force[NP][3], double ener, double temp, double P, char filename[20] );   // Save trajectory in file
double randn();                                                                 // Subroutine to obtain normally distributed random numbers between -1 and 1


/* -------------------------------------------------------------------------- */
/* MAIN BODY */
/* -------------------------------------------------------------------------- */


void main()
{

  double pos[NP][3];
  double vel[NP][3];
  double v_half[NP][3];
  double force[NP][3];
  double force_next[NP][3];
  double F;
  double TotEn, KinEn, PotEn;
  double MomX, MomY, MomZ;
  double Temp;
  double P;

  srandom(time(NULL));
  // srandom(123456);


  /* 1 --> SYSTEM INITIALIZATION */

  // 1.1 --> Initialize positions


  // // (OPTION 1) Initial positions are extracted from file (MonteCarlo equilibrium configuration at desired T)
  // FILE *ifp;
  // ifp = fopen("initial.txt","r");
  // double x,y,z;
  // int count = 0;
  // while (fscanf(ifp, "%lf\t%lf\t%lf", &x, &y, &z) == 3){
  //   pos[count][0] = x;
  //   pos[count][1] = y;
  //   pos[count][2] = z;
  //   count = count + 1;
  //   // printf("%lf\t%lf\t%lf\n",x,y,z);
  // }
  // fclose(ifp);

  // (OPTION 2) Initialize equally-spaced lattice of particles
  double np_side = floor(pow(NP,(1.0/3.0)) + 1.0);
  int np_count = 0;
  //printf("%f\n",np_side);
  //printf("The %i particles will be placed at the following coordinates:\n\n",NP);
  for (int x = 0; x < np_side; x++){
    for (int y = 0; y < np_side; y++){
      for (int z = 0; z < np_side; z++){
        if (np_count < NP){
          pos[np_count][0] = -(LBOX/2.0) + (x+0.5)*(LBOX/(np_side));
          pos[np_count][1] = -(LBOX/2.0) + (y+0.5)*(LBOX/(np_side));
          pos[np_count][2] = -(LBOX/2.0) + (z+0.5)*(LBOX/(np_side));
          //printf("\t%f\t%f\t%f\n",pos[np_count][0],pos[np_count][1],pos[np_count][2]);
          np_count = np_count + 1;
        } else {
          break;
        }
      }
    }
  }


  // 1.2 --> Initialize velocities

  // printf("\nVELOCITIES:\n");
  for (int i = 0; i < NP; i++){
    for (int dim = 0; dim < 3; dim++){
      vel[i][dim] = sqrt(T/5.0)*randn();
      v_half[i][dim] = 0.0;
    }
    // printf("%lf\t%lf\t%lf\n",vel[i][0],vel[i][1],vel[i][2]); // Works fine, uncomment to check
  }

  // 1.3 --> Initialize forces

  // First of all, we set all forces to zero
  for (int i = 0; i < NP; i++){
    force[i][0] = 0.0;
    force[i][1] = 0.0;
    force[i][2] = 0.0;
  }

  // Then, we get the forces acting over each particle

  for (int i = 0; i < NP - 1; i++){
    for (int j = i+1; j < NP; j++){
      F = calc_LJforce(pos[i],pos[j]);

      force[i][0] = force[i][0] + F*Xdif(pos[i],pos[j]);
      force[i][1] = force[i][1] + F*Ydif(pos[i],pos[j]);
      force[i][2] = force[i][2] + F*Zdif(pos[i],pos[j]);

      force[j][0] = force[j][0] - F*Xdif(pos[i],pos[j]);
      force[j][1] = force[j][1] - F*Ydif(pos[i],pos[j]);
      force[j][2] = force[j][2] - F*Zdif(pos[i],pos[j]);
    }
  }


  /* 2 --> MOLECULAR DYNAMICS SIMULATION */

  printf("Simulation_Percentage\tIteration\tTime\tTotEn\tPotEn\tKinEn\tTemperature\tMomentumX\tMomentumY\tMomentumZ\n");

  for (int iter = 0; iter < NSTEP; iter++){

    // 2.1 --> First Verlet step: update v_half and pos

    //printf("\nITERATION %i\n\n",iter+1);
    //printf("\tPOSITIONS:\n");
    for (int i = 0; i < NP; i++){
      for (int dim = 0; dim < 3; dim++){
        v_half[i][dim] = vel[i][dim] + 0.5*force[i][dim]*dt;
        pos[i][dim] = pos[i][dim] + v_half[i][dim]*dt;
        pos[i][dim] = apply_pbc(pos[i][dim]);
      }
      //printf("\t%lf\t%lf\t%lf\n",pos[i][0],pos[i][1],pos[i][2]);
    }

    // 2.2 --> Second Verlet step: retrieve forces

    for (int i = 0; i < NP; i++){
      force[i][0] = 0.0;
      force[i][1] = 0.0;
      force[i][2] = 0.0;
    }

    for (int i = 0; i < NP - 1; i++){
      for (int j = i + 1; j < NP; j++){
        F = calc_LJforce(pos[i],pos[j]);

        force[i][0] = force[i][0] + F*Xdif(pos[i],pos[j]);
        force[i][1] = force[i][1] + F*Ydif(pos[i],pos[j]);
        force[i][2] = force[i][2] + F*Zdif(pos[i],pos[j]);

        force[j][0] = force[j][0] - F*Xdif(pos[i],pos[j]);
        force[j][1] = force[j][1] - F*Ydif(pos[i],pos[j]);
        force[j][2] = force[j][2] - F*Zdif(pos[i],pos[j]);
      }
    }

    // 2.3 --> Third Verlet step: update velocities

    //printf("\n\t\t\tVELOCITIES:\n");
    for (int i = 0; i < NP; i++){
      for (int dim = 0; dim < 3; dim++){
        vel[i][dim] = v_half[i][dim] + 0.5*force[i][dim]*dt;
      }
      //printf("\t\t\t%lf\t%lf\t%lf\n",vel[i][0],vel[i][1],vel[i][2]);
    }

    // 2.4 --> Store trajectories

    double print_iter = floor(NSTEP/NSAMP);
    if (iter % (int)print_iter == 0){

      TotEn = 0.0;
      KinEn = 0.0;
      PotEn = 0.0;

      MomX = 0.0;
      MomY = 0.0;
      MomZ = 0.0;

      Temp = 0.0;
      //P = 0.0;

      for (int i = 0; i < NP; i++){
        KinEn = KinEn + 0.5*(vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );

        MomX = MomX + vel[i][0];
        MomY = MomY + vel[i][1];
        MomZ = MomZ + vel[i][2];

        Temp = Temp + (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] );
        //P = P - (1.0/3.0)*norm(pos[i])*norm(force[i]);

        for (int j = 0; j < NP; j++){
          if (i != j){
            PotEn = PotEn + 0.5*calc_LJpot(pos[j],pos[i]);
          } else {
            PotEn = PotEn + 0.0;
          }
        }
      }

      PotEn = PotEn/NP;
      KinEn = KinEn/NP;
      TotEn = TotEn + PotEn + KinEn;

      Temp = Temp/NP;
      //P = P + Temp*NP;
      //P = P / VOLBOX;

      printf("%.2lf%%\t%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%.14lf\t%.14lf\t%.14lf\n",(iter/NSTEP)*100.0,iter,dt*(double)iter,TotEn,PotEn,KinEn,Temp,MomX,MomY,MomZ);

      // Uncomment code below only if a verbose output is desired
      // printf("Simulation status: %f %%\t...\t%1.12f fs\n", (100.0*iter/NSTEP),dt*iter);
      // printf("\tTotal energy: %f\tPotential energy: %f\tKinetic energy: %f\n",TotEn,PotEn,KinEn);
      // printf("\t\tTotal momentum: %f\t%f\t%f\n",MomX,MomY,MomZ);
      // printf("\t\t\tTemperature of the system: %f\n",T);
      // printf("\t\t\tTemperature: %f\tPressure: %f\n",Temp,P);
      // printf("\n--------------------------------------------------------\n");

      char filename[20];
      sprintf(filename, "%d", iter);
      strcat(filename,"_traj.dat");
      save_trajectory(pos,vel,force,TotEn,Temp,P,filename);
    }

  }

}


/* -------------------------------------------------------------------------- */
/* FUNCTIONS CODE */
/* -------------------------------------------------------------------------- */



double randn() /* Subroutine to obtain normally distributed random numbers
                  by implementing the Box-Muller polar algorithm */
{
  double x1, x2, w, y1, y2;


  for (int i = 0; i < 2; i++){
    while (w >= 1.0){
      x1 = random();
      x2 = random();
      x1 = 2.0*(x1/RAND_MAX - 0.5);
      x2 = 2.0*(x2/RAND_MAX - 0.5);
      w = x1*x1 + x2*x2;
      // printf("%f\t%f\n",x1,x2);
    }

    w = sqrt((-2.0*log(w))/w);
    y1 = x1*w;
    y2 = x2*w;

    // printf("%f\t%f\t\t%f\t\t%f\t%f\n",x1,x2,w,y1,y2);

    if (i != 0){ // The first iteration always goes wrong, so we discard it
      return y1;
    }
  }
}



void print_welcome() /* Subroutine to print a friendly message at the start */
{

  printf("\n--------------------------------------------------------\n");
  printf("\tMOLECULAR DYNAMICS SIMULATION -- LENNARD-JONES POTENTIAL\n");
  printf("\t\t\t\t\t\tBy Juan Jimenez Sanchez\n\t\t\t\t\t(juan.jimenezs@estudiante.uam.es)\n\n");
  printf("\tThis simulation is carried out keeping NVT constant (Langevin dynamics) under the following conditions:\n");
  printf("\tNumber of particles: %i;\tBox size: %.1f sigma;\tDensity: %f;\tParticle mass: %.1f\n",NP,LBOX,NP/VOLBOX,MASS);
  printf("\tTime step: %f;\tSimulation steps: %.1f\n",dt,NSTEP);
  printf("\n--------------------------------------------------------\n\n");

}



double norm( double pos[3] ) /* Subroutine to get the euclidean norm of a vector */
{
  double norm;

  norm = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);

  return norm;
}



double Xdif( double pos1[3], double pos2[3] ) /* Subroutine to get the difference between 2 points in their x-coord. Useful when getting forces */
{
  double xdif;

  xdif = (pos2[0] - pos1[0]);
  /* Need to apply minimum image convention, so as to keep PBC  */
  xdif = apply_pbc(xdif);

  return xdif;
}



double Ydif( double pos1[3], double pos2[3] ) /* Subroutine to get the difference between 2 points in their y-coord. Useful when getting forces */
{
  double ydif;

  ydif = (pos2[1] - pos1[1]);
  /* Need to apply minimum image convention, so as to keep PBC  */
  ydif = apply_pbc(ydif);

  return ydif;
}



double Zdif( double pos1[3], double pos2[3] ) /* Subroutine to get the difference between 2 points in their z-coord. Useful when getting forces */
{
  double zdif;

  zdif = (pos2[2] - pos1[2]);
  /* Need to apply minimum image convention, so as to keep PBC  */
  zdif = apply_pbc(zdif);

  return zdif;
}



double get_dist( double pos1[3], double pos2[3] ) /* Subroutine to get euclidean distance between two particles */
{
  double dist, xdif, ydif, zdif;

  xdif = Xdif(pos1,pos2);
  ydif = Ydif(pos1,pos2);
  zdif = Zdif(pos1,pos2);

  // dist = sqrt(pow(xdif,2) + pow(ydif,2) + pow(zdif,2));
  dist = (xdif*xdif + ydif*ydif + zdif*zdif);
  // printf("Dist: %lf, x: %lf, y: %lf, z: %lf\n",dist,xdif,ydif,zdif);
  return dist;
}




double calc_LJpot( double r[3], double pos[3] ) /* Subroutine to get Lennard-Jones potential between two interacting particles */
{
  double LJpot, dist;
  dist = get_dist(r,pos);

  double dist2 = dist*dist;
  double dist6 = dist2*dist2*dist2;
  double dist12 = dist6*dist6;

  if (dist < CUTOFF){
    //LJpot = 4.0 * EPSILON * ( pow(SIGMA/dist,12) - pow(SIGMA/dist,6) );
    LJpot = 4.0 * EPSILON * ( 1.0/dist12 - 1.0/dist6 ) - UCUT;
  } else {
    LJpot = 0.0;
  }
  return LJpot;
}



double calc_LJforce( double r[3], double pos[3] ) /* Subroutine to get Lennard-Jones potential force between two interacting particles */
{
  double LJforce, dist;
  dist = get_dist(r,pos);

  // As SIGMA is equal to 1, to save time, we will use 1.0 instead of the variable
  double dist2 = dist*dist;
  double dist6 = dist2*dist2*dist2;
  double dist12 = dist6*dist6;

  if (dist <= CUTOFF){
    //LJforce = - 24.0 * EPSILON * ( 2.0*(pow(SIGMA,12)/pow(dist,14)) - (pow(SIGMA,6)/pow(dist,8)) );
    // CLARIFICATION --> Minus sign before LJforce is required, as it is the way to grant that we are calculating the force exerted by other particle over OUR particle
    // Otherwise, we get the force OUR particle does over another particle, thus leading to opposite Newton laws, progressive approachment, and collisions between particles.
    LJforce = -24.0 * EPSILON * ( 2.0/dist12 - 1.0/dist6 ) / dist2 - FCUT;
  } else {
    LJforce = 0.0;
  }
  return LJforce;
}


double apply_pbc( double r ) /* Subroutine to apply Periodic Boundary Conditions to a given set of coordinates */
{
  // A lot of crap in this comments, just previous versions until working one was achieved 
  //r-=(int)( ( (r/LBOX<0)?-0.5:0.5) + r)*LBOX;
  //printf("%f\t",r);
  r = r - floor(r/LBOX + 0.5)*LBOX;
  //printf("%f\n",r);

  // r-=(int)( ( (r<0)?-0.5:0.5 ) + r); /* TRuncates correctly negative double values */
  /* r = r - (int)((r/LBOX)+0.5)*LBOX; */
  return r;

  // Uncomment this little piece of code to check that PBCs work correctly
  /* for (double r = -10.0; r < 10.1; r = r + 0.2){
    double r_pbc = apply_pbc(r);
    printf("%lf\t%lf\n",r,r_pbc);
  } */
}



void save_trajectory( double pos[NP][3], double vel[NP][3], double force[NP][3], double ener, double temp, double P, char filename[20] ) /* Subroutine to store simulation data into a file */
{
  FILE *fp;
  fp = fopen(filename,"w+");
  fprintf(fp, "PosX\tPosY\tPosZ\tVelX\tVelY\tVelZ\tForceX\tForceY\tForceZ\tEner\tT\tP\n");
  for (int i = 0; i < NP; i++){
    /* Apply PBC --> Fold back to primary box if particle leaves it
    for (int dim = 0; dim < 3; dim++){
      pos[p][dim] = pos[p][dim] - (int)((pos[p][dim]/LBOX)+0.5)*LBOX;
    } */
    fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", pos[i][0], pos[i][1], pos[i][2], vel[i][0], vel[i][1], vel[i][2], force[i][0], force[i][1], force[i][2], ener, temp, P);
  }
  fclose(fp);
}