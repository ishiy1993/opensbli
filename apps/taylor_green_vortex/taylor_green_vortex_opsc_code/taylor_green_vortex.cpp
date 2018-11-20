#include <stdlib.h>
#include <string.h>
#include <math.h>
// Global constants in the equations are
int nx2;
int nx1;
double deltai2;
double rinv8;
double rinv9;
double Minf;
double rinv1;
double rinv4;
double Pr;
double rinv12;
double rinv13;
double rinv10;
double rinv11;
double rinv14;
double rknew[3];
double rinv15;
double rc5;
double rc6;
double rc0;
double rc2;
double rc3;
int nx0;
double deltai1;
double deltai0;
double Re;
double deltat;
double gama;
double rkold[3];
double rinv7;
// OPS header file
#define OPS_3D
#include "ops_seq.h"
#include "ops_hdf5.h"
#include "taylor_green_vortex_block_0_kernel.h"

// main program start
int main (int argc, char **argv) 
{

   gama = 1.40000000000000;
   Pr = 0.710000000000000;
   nx2 = 64;
   nx0 = 64;
   Re = 1600;
   deltat = 0.00338500000000000;
   nx1 = 64;
   Minf = 0.100000000000000;
   rc5 = 5.0/2.0;
   rc6 = 4.0/3.0;
   rc0 = 1.0/2.0;
   rc2 = 1.0/12.0;
   rc3 = 2.0/3.0;
   rkold[0] = 1.0/4.0;
   rkold[1] = 3.0/20.0;
   rkold[2] = 3.0/5.0;
   rknew[0] = 2.0/3.0;
   rknew[1] = 5.0/12.0;
   rknew[2] = 3.0/5.0;
   rinv13 = 1.0/(gama - 1);
   rinv11 = 1.0/Re;
   rinv14 = pow(Minf, -2);
   rinv15 = 1.0/(gama*pow(Minf, 2));
   rinv12 = 1.0/Pr;
   deltai2 = 0.03125*M_PI;
   deltai1 = 0.03125*M_PI;
   deltai0 = 0.03125*M_PI;
   rinv10 = pow(deltai0, -2);
   rinv8 = 1.0/deltai0;
   rinv9 = pow(deltai2, -2);
   rinv1 = 1.0/deltai2;
   rinv4 = pow(deltai1, -2);
   rinv7 = 1.0/deltai1;

   // Initializing OPS 
   ops_init(argc,argv,1);

   ops_decl_const("nx2" , 1, "int", &nx2);
   ops_decl_const("nx1" , 1, "int", &nx1);
   ops_decl_const("deltai2" , 1, "double", &deltai2);
   ops_decl_const("rinv8" , 1, "double", &rinv8);
   ops_decl_const("rinv9" , 1, "double", &rinv9);
   ops_decl_const("Minf" , 1, "double", &Minf);
   ops_decl_const("rinv1" , 1, "double", &rinv1);
   ops_decl_const("rinv4" , 1, "double", &rinv4);
   ops_decl_const("Pr" , 1, "double", &Pr);
   ops_decl_const("rinv12" , 1, "double", &rinv12);
   ops_decl_const("rinv13" , 1, "double", &rinv13);
   ops_decl_const("rinv10" , 1, "double", &rinv10);
   ops_decl_const("rinv11" , 1, "double", &rinv11);
   ops_decl_const("rinv14" , 1, "double", &rinv14);
   ops_decl_const("rinv15" , 1, "double", &rinv15);
   ops_decl_const("rc5" , 1, "double", &rc5);
   ops_decl_const("rc6" , 1, "double", &rc6);
   ops_decl_const("rc0" , 1, "double", &rc0);
   ops_decl_const("rc2" , 1, "double", &rc2);
   ops_decl_const("rc3" , 1, "double", &rc3);
   ops_decl_const("nx0" , 1, "int", &nx0);
   ops_decl_const("deltai1" , 1, "double", &deltai1);
   ops_decl_const("deltai0" , 1, "double", &deltai0);
   ops_decl_const("Re" , 1, "double", &Re);
   ops_decl_const("deltat" , 1, "double", &deltat);
   ops_decl_const("gama" , 1, "double", &gama);
   ops_decl_const("rinv7" , 1, "double", &rinv7);

   // Defining block in OPS Format
   ops_block taylor_green_vortex_block;

   // Initialising block in OPS Format
   taylor_green_vortex_block = ops_decl_block(3, "taylor_green_vortex_block");

   // Define dataset
   ops_dat p;
   ops_dat rhou1_old;
   ops_dat wk20;
   ops_dat wk21;
   ops_dat wk28;
   ops_dat wk61;
   ops_dat wk0;
   ops_dat wk15;
   ops_dat wk35;
   ops_dat wk12;
   ops_dat wk34;
   ops_dat wk30;
   ops_dat wk39;
   ops_dat wk44;
   ops_dat wk40;
   ops_dat wk45;
   ops_dat rhou2_old;
   ops_dat wk66;
   ops_dat wk3;
   ops_dat wk43;
   ops_dat wk2;
   ops_dat wk50;
   ops_dat wk14;
   ops_dat wk26;
   ops_dat wk64;
   ops_dat wk49;
   ops_dat wk27;
   ops_dat wk65;
   ops_dat wk5;
   ops_dat wk23;
   ops_dat wk9;
   ops_dat rhou1;
   ops_dat wk51;
   ops_dat rhou0;
   ops_dat T;
   ops_dat wk56;
   ops_dat wk13;
   ops_dat wk36;
   ops_dat wk47;
   ops_dat wk29;
   ops_dat wk19;
   ops_dat wk55;
   ops_dat wk18;
   ops_dat rhoE_old;
   ops_dat rhoE;
   ops_dat wk8;
   ops_dat wk37;
   ops_dat wk42;
   ops_dat wk10;
   ops_dat rho_old;
   ops_dat wk46;
   ops_dat wk62;
   ops_dat u2;
   ops_dat wk48;
   ops_dat wk25;
   ops_dat wk63;
   ops_dat wk7;
   ops_dat wk1;
   ops_dat wk54;
   ops_dat wk33;
   ops_dat wk6;
   ops_dat rho;
   ops_dat wk11;
   ops_dat rhou2;
   ops_dat wk32;
   ops_dat wk59;
   ops_dat wk38;
   ops_dat wk31;
   ops_dat u1;
   ops_dat wk60;
   ops_dat u0;
   ops_dat wk67;
   ops_dat wk22;
   ops_dat wk24;
   ops_dat wk53;
   ops_dat wk4;
   ops_dat wk57;
   ops_dat wk41;
   ops_dat wk52;
   ops_dat wk17;
   ops_dat wk16;
   ops_dat rhou0_old;
   ops_dat wk58;

   // Initialise/allocate OPS dataset.
   int halo_p[] = {2, 2, 2};
   int halo_m[] = {-2, -2, -2};
   int size[] = {nx0, nx1, nx2};
   int base[] = {0, 0, 0};
   double* val = NULL;
   p = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "p");
   rhou1_old = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhou1_old");
   wk20 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk20");
   wk21 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk21");
   wk28 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk28");
   wk61 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk61");
   wk0 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk0");
   wk15 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk15");
   wk35 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk35");
   wk12 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk12");
   wk34 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk34");
   wk30 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk30");
   wk39 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk39");
   wk44 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk44");
   wk40 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk40");
   wk45 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk45");
   rhou2_old = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhou2_old");
   wk66 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk66");
   wk3 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk3");
   wk43 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk43");
   wk2 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk2");
   wk50 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk50");
   wk14 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk14");
   wk26 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk26");
   wk64 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk64");
   wk49 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk49");
   wk27 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk27");
   wk65 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk65");
   wk5 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk5");
   wk23 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk23");
   wk9 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk9");
   rhou1 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhou1");
   wk51 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk51");
   rhou0 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhou0");
   T = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "T");
   wk56 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk56");
   wk13 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk13");
   wk36 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk36");
   wk47 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk47");
   wk29 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk29");
   wk19 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk19");
   wk55 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk55");
   wk18 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk18");
   rhoE_old = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhoE_old");
   rhoE = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhoE");
   wk8 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk8");
   wk37 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk37");
   wk42 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk42");
   wk10 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk10");
   rho_old = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rho_old");
   wk46 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk46");
   wk62 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk62");
   u2 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "u2");
   wk48 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk48");
   wk25 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk25");
   wk63 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk63");
   wk7 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk7");
   wk1 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk1");
   wk54 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk54");
   wk33 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk33");
   wk6 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk6");
   rho = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rho");
   wk11 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk11");
   rhou2 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhou2");
   wk32 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk32");
   wk59 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk59");
   wk38 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk38");
   wk31 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk31");
   u1 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "u1");
   wk60 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk60");
   u0 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "u0");
   wk67 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk67");
   wk22 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk22");
   wk24 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk24");
   wk53 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk53");
   wk4 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk4");
   wk57 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk57");
   wk41 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk41");
   wk52 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk52");
   wk17 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk17");
   wk16 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk16");
   rhou0_old = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "rhou0_old");
   wk58 = ops_decl_dat(taylor_green_vortex_block, 1, size, base, halo_m, halo_p, val, "double", "wk58");

   // Declare all the stencils used 
   int stencil3_temp[] = {0,-2,0,0,-1,0,0,1,0,0,2,0};
   ops_stencil stencil3 = ops_decl_stencil(3,4,stencil3_temp,"0,-2,0,0,-1,0,0,1,0,0,2,0");
   int stencil6_temp[] = {-2,0,0,-1,0,0,0,0,0,1,0,0,2,0,0};
   ops_stencil stencil6 = ops_decl_stencil(3,5,stencil6_temp,"-2,0,0,-1,0,0,0,0,0,1,0,0,2,0,0");
   int stencil2_temp[] = {0,-2,0,0,-1,0,0,0,0,0,1,0,0,2,0};
   ops_stencil stencil2 = ops_decl_stencil(3,5,stencil2_temp,"0,-2,0,0,-1,0,0,0,0,0,1,0,0,2,0");
   int stencil5_temp[] = {0,0,-2,0,0,-1,0,0,0,0,0,1,0,0,2};
   ops_stencil stencil5 = ops_decl_stencil(3,5,stencil5_temp,"0,0,-2,0,0,-1,0,0,0,0,0,1,0,0,2");
   int stencil1_temp[] = {0,0,-2,0,0,-1,0,0,1,0,0,2};
   ops_stencil stencil1 = ops_decl_stencil(3,4,stencil1_temp,"0,0,-2,0,0,-1,0,0,1,0,0,2");
   int stencil0_temp[] = {0,0,0};
   ops_stencil stencil0 = ops_decl_stencil(3,1,stencil0_temp,"0,0,0");
   int stencil4_temp[] = {-2,0,0,-1,0,0,1,0,0,2,0,0};
   ops_stencil stencil4 = ops_decl_stencil(3,4,stencil4_temp,"-2,0,0,-1,0,0,1,0,0,2,0,0");

   ops_reduction enstrophy = ops_decl_reduction_handle(sizeof(double), "double", "reduction_enstrophy");
   ops_reduction rhomean = ops_decl_reduction_handle(sizeof(double), "double", "reduction_rhomean");
   ops_reduction ke = ops_decl_reduction_handle(sizeof(double), "double", "reduction_ke");

   // Boundary condition exchange code
   ops_halo_group halo_exchange0 ;
   {
      int halo_iter[] = {2, nx1 + 4, nx2 + 4};
      int from_base[] = {0, -2, -2};
      int to_base[] = {nx0, -2, -2};
      int dir[] = {1, 2, 3};
      ops_halo halo0 = ops_decl_halo(rho, rho, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo1 = ops_decl_halo(rhou0, rhou0, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo2 = ops_decl_halo(rhou1, rhou1, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo3 = ops_decl_halo(rhou2, rhou2, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo4 = ops_decl_halo(rhoE, rhoE, halo_iter, from_base, to_base, dir, dir);
      ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
      halo_exchange0 = ops_decl_halo_group(5,grp);
   }
   // Boundary condition exchange code
   ops_halo_group halo_exchange1 ;
   {
      int halo_iter[] = {2, nx1 + 4, nx2 + 4};
      int from_base[] = {nx0 - 2, -2, -2};
      int to_base[] = {-2, -2, -2};
      int dir[] = {1, 2, 3};
      ops_halo halo0 = ops_decl_halo(rho, rho, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo1 = ops_decl_halo(rhou0, rhou0, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo2 = ops_decl_halo(rhou1, rhou1, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo3 = ops_decl_halo(rhou2, rhou2, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo4 = ops_decl_halo(rhoE, rhoE, halo_iter, from_base, to_base, dir, dir);
      ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
      halo_exchange1 = ops_decl_halo_group(5,grp);
   }
   // Boundary condition exchange code
   ops_halo_group halo_exchange2 ;
   {
      int halo_iter[] = {nx0 + 4, 2, nx2 + 4};
      int from_base[] = {-2, 0, -2};
      int to_base[] = {-2, nx1, -2};
      int dir[] = {1, 2, 3};
      ops_halo halo0 = ops_decl_halo(rho, rho, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo1 = ops_decl_halo(rhou0, rhou0, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo2 = ops_decl_halo(rhou1, rhou1, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo3 = ops_decl_halo(rhou2, rhou2, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo4 = ops_decl_halo(rhoE, rhoE, halo_iter, from_base, to_base, dir, dir);
      ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
      halo_exchange2 = ops_decl_halo_group(5,grp);
   }
   // Boundary condition exchange code
   ops_halo_group halo_exchange3 ;
   {
      int halo_iter[] = {nx0 + 4, 2, nx2 + 4};
      int from_base[] = {-2, nx1 - 2, -2};
      int to_base[] = {-2, -2, -2};
      int dir[] = {1, 2, 3};
      ops_halo halo0 = ops_decl_halo(rho, rho, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo1 = ops_decl_halo(rhou0, rhou0, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo2 = ops_decl_halo(rhou1, rhou1, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo3 = ops_decl_halo(rhou2, rhou2, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo4 = ops_decl_halo(rhoE, rhoE, halo_iter, from_base, to_base, dir, dir);
      ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
      halo_exchange3 = ops_decl_halo_group(5,grp);
   }
   // Boundary condition exchange code
   ops_halo_group halo_exchange4 ;
   {
      int halo_iter[] = {nx0 + 4, nx1 + 4, 2};
      int from_base[] = {-2, -2, 0};
      int to_base[] = {-2, -2, nx2};
      int dir[] = {1, 2, 3};
      ops_halo halo0 = ops_decl_halo(rho, rho, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo1 = ops_decl_halo(rhou0, rhou0, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo2 = ops_decl_halo(rhou1, rhou1, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo3 = ops_decl_halo(rhou2, rhou2, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo4 = ops_decl_halo(rhoE, rhoE, halo_iter, from_base, to_base, dir, dir);
      ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
      halo_exchange4 = ops_decl_halo_group(5,grp);
   }
   // Boundary condition exchange code
   ops_halo_group halo_exchange5 ;
   {
      int halo_iter[] = {nx0 + 4, nx1 + 4, 2};
      int from_base[] = {-2, -2, nx2 - 2};
      int to_base[] = {-2, -2, -2};
      int dir[] = {1, 2, 3};
      ops_halo halo0 = ops_decl_halo(rho, rho, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo1 = ops_decl_halo(rhou0, rhou0, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo2 = ops_decl_halo(rhou1, rhou1, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo3 = ops_decl_halo(rhou2, rhou2, halo_iter, from_base, to_base, dir, dir);
      ops_halo halo4 = ops_decl_halo(rhoE, rhoE, halo_iter, from_base, to_base, dir, dir);
      ops_halo grp[] = {halo0,halo1,halo2,halo3,halo4};
      halo_exchange5 = ops_decl_halo_group(5,grp);
   }

   // Init OPS partition
   ops_partition("");

   int iter_range88[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
   ops_par_loop(taylor_green_vortex_block0_88_kernel, "Initialisation", taylor_green_vortex_block, 3, iter_range88,
   ops_arg_dat(rhou1, 1, stencil0, "double", OPS_WRITE),
   ops_arg_dat(rhoE, 1, stencil0, "double", OPS_WRITE),
   ops_arg_dat(rho, 1, stencil0, "double", OPS_WRITE),
   ops_arg_dat(rhou2, 1, stencil0, "double", OPS_WRITE),
   ops_arg_dat(rhou0, 1, stencil0, "double", OPS_WRITE),
   ops_arg_idx());



   // Boundary condition exchange calls
   ops_halo_transfer(halo_exchange0);
   // Boundary condition exchange calls
   ops_halo_transfer(halo_exchange1);
   // Boundary condition exchange calls
   ops_halo_transfer(halo_exchange2);
   // Boundary condition exchange calls
   ops_halo_transfer(halo_exchange3);
   // Boundary condition exchange calls
   ops_halo_transfer(halo_exchange4);
   // Boundary condition exchange calls
   ops_halo_transfer(halo_exchange5);

   double cpu_start, elapsed_start;
   ops_timers(&cpu_start, &elapsed_start);

   for (int iteration=0; iteration<5909; iteration++){


      int iter_range87[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
      ops_par_loop(taylor_green_vortex_block0_87_kernel, "Save equations", taylor_green_vortex_block, 3, iter_range87,
      ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
      ops_arg_dat(rhoE, 1, stencil0, "double", OPS_READ),
      ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
      ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
      ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
      ops_arg_dat(rhou1_old, 1, stencil0, "double", OPS_WRITE),
      ops_arg_dat(rhou2_old, 1, stencil0, "double", OPS_WRITE),
      ops_arg_dat(rhou0_old, 1, stencil0, "double", OPS_WRITE),
      ops_arg_dat(rho_old, 1, stencil0, "double", OPS_WRITE),
      ops_arg_dat(rhoE_old, 1, stencil0, "double", OPS_WRITE));



      for (int stage=0; stage<3; stage++){


         int iter_range0[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_0_kernel, "Grouped Formula Evaluation", taylor_green_vortex_block, 3, iter_range0,
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u0, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(u1, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(u2, 1, stencil0, "double", OPS_WRITE));


         int iter_range1[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_1_kernel, "Non-Grouped Formula Evaluation", taylor_green_vortex_block, 3, iter_range1,
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhoE, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(p, 1, stencil0, "double", OPS_WRITE));


         int iter_range2[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_2_kernel, "Non-Grouped Formula Evaluation", taylor_green_vortex_block, 3, iter_range2,
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(p, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(T, 1, stencil0, "double", OPS_WRITE));


         int iter_range3[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_3_kernel, "D(rhoE x2)", taylor_green_vortex_block, 3, iter_range3,
         ops_arg_dat(rhoE, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk0, 1, stencil0, "double", OPS_WRITE));


         int iter_range4[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_4_kernel, "D(u0 x1 x1)", taylor_green_vortex_block, 3, iter_range4,
         ops_arg_dat(u0, 1, stencil2, "double", OPS_READ),
         ops_arg_dat(wk1, 1, stencil0, "double", OPS_WRITE));


         int iter_range5[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_5_kernel, "D(rhou1 x2)", taylor_green_vortex_block, 3, iter_range5,
         ops_arg_dat(rhou1, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk2, 1, stencil0, "double", OPS_WRITE));


         int iter_range6[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_6_kernel, "D(rhou0 x1)", taylor_green_vortex_block, 3, iter_range6,
         ops_arg_dat(rhou0, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk3, 1, stencil0, "double", OPS_WRITE));


         int iter_range7[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_7_kernel, "D(u0 x2)", taylor_green_vortex_block, 3, iter_range7,
         ops_arg_dat(u0, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk4, 1, stencil0, "double", OPS_WRITE));


         int iter_range8[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_8_kernel, "D(rhou1 x0)", taylor_green_vortex_block, 3, iter_range8,
         ops_arg_dat(rhou1, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk5, 1, stencil0, "double", OPS_WRITE));


         int iter_range9[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_9_kernel, "D(rhou0 x0)", taylor_green_vortex_block, 3, iter_range9,
         ops_arg_dat(rhou0, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk6, 1, stencil0, "double", OPS_WRITE));


         int iter_range10[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_10_kernel, "rho*u1[x0 x1 x2 t]", taylor_green_vortex_block, 3, iter_range10,
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range11[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_11_kernel, "D(rho*u1[x0 x1 x2 t] x1)", taylor_green_vortex_block, 3, iter_range11,
         ops_arg_dat(wk63, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk7, 1, stencil0, "double", OPS_WRITE));


         int iter_range12[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_12_kernel, "D(u0 x2 x2)", taylor_green_vortex_block, 3, iter_range12,
         ops_arg_dat(u0, 1, stencil5, "double", OPS_READ),
         ops_arg_dat(wk8, 1, stencil0, "double", OPS_WRITE));


         int iter_range13[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_13_kernel, "D(u2 x2 x2)", taylor_green_vortex_block, 3, iter_range13,
         ops_arg_dat(u2, 1, stencil5, "double", OPS_READ),
         ops_arg_dat(wk9, 1, stencil0, "double", OPS_WRITE));


         int iter_range14[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_14_kernel, "D(rho x1)", taylor_green_vortex_block, 3, iter_range14,
         ops_arg_dat(rho, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk10, 1, stencil0, "double", OPS_WRITE));


         int iter_range15[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_15_kernel, "D(u1[x0 x1 x2 t] x2 x2)", taylor_green_vortex_block, 3, iter_range15,
         ops_arg_dat(u1, 1, stencil5, "double", OPS_READ),
         ops_arg_dat(wk11, 1, stencil0, "double", OPS_WRITE));


         int iter_range16[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_16_kernel, "p*u1", taylor_green_vortex_block, 3, iter_range16,
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(p, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range17[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_17_kernel, "D(p*u1 x1)", taylor_green_vortex_block, 3, iter_range17,
         ops_arg_dat(wk63, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk12, 1, stencil0, "double", OPS_WRITE));


         int iter_range18[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_18_kernel, "D(u2 x0 x0)", taylor_green_vortex_block, 3, iter_range18,
         ops_arg_dat(u2, 1, stencil6, "double", OPS_READ),
         ops_arg_dat(wk13, 1, stencil0, "double", OPS_WRITE));


         int iter_range19[] = {0, nx0, 0, nx1, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_19_kernel, "D(u2 x1)", taylor_green_vortex_block, 3, iter_range19,
         ops_arg_dat(u2, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk14, 1, stencil0, "double", OPS_WRITE));


         int iter_range20[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_20_kernel, "rhou1*u0", taylor_green_vortex_block, 3, iter_range20,
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range21[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_21_kernel, "D(rhou1*u0 x0)", taylor_green_vortex_block, 3, iter_range21,
         ops_arg_dat(wk63, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk15, 1, stencil0, "double", OPS_WRITE));


         int iter_range22[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_22_kernel, "D(u0 x1)", taylor_green_vortex_block, 3, iter_range22,
         ops_arg_dat(u0, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk16, 1, stencil0, "double", OPS_WRITE));


         int iter_range23[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_23_kernel, "D(p x2)", taylor_green_vortex_block, 3, iter_range23,
         ops_arg_dat(p, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk17, 1, stencil0, "double", OPS_WRITE));


         int iter_range24[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_24_kernel, "D(rhou2 x0)", taylor_green_vortex_block, 3, iter_range24,
         ops_arg_dat(rhou2, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk18, 1, stencil0, "double", OPS_WRITE));


         int iter_range25[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_25_kernel, "D(rho x0)", taylor_green_vortex_block, 3, iter_range25,
         ops_arg_dat(rho, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk19, 1, stencil0, "double", OPS_WRITE));


         int iter_range26[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_26_kernel, "D(rhoE x1)", taylor_green_vortex_block, 3, iter_range26,
         ops_arg_dat(rhoE, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk20, 1, stencil0, "double", OPS_WRITE));


         int iter_range27[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_27_kernel, "rho*u0", taylor_green_vortex_block, 3, iter_range27,
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range28[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_28_kernel, "D(rho*u0 x0)", taylor_green_vortex_block, 3, iter_range28,
         ops_arg_dat(wk63, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk21, 1, stencil0, "double", OPS_WRITE));


         int iter_range29[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_29_kernel, "D(rhou1 x1)", taylor_green_vortex_block, 3, iter_range29,
         ops_arg_dat(rhou1, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk22, 1, stencil0, "double", OPS_WRITE));


         int iter_range30[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_30_kernel, "D(u1[x0 x1 x2 t] x0 x0)", taylor_green_vortex_block, 3, iter_range30,
         ops_arg_dat(u1, 1, stencil6, "double", OPS_READ),
         ops_arg_dat(wk23, 1, stencil0, "double", OPS_WRITE));


         int iter_range31[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_31_kernel, "rho*u2", taylor_green_vortex_block, 3, iter_range31,
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range32[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_32_kernel, "D(rho*u2 x2)", taylor_green_vortex_block, 3, iter_range32,
         ops_arg_dat(wk63, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk24, 1, stencil0, "double", OPS_WRITE));


         int iter_range33[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_33_kernel, "D(u0 x0 x0)", taylor_green_vortex_block, 3, iter_range33,
         ops_arg_dat(u0, 1, stencil6, "double", OPS_READ),
         ops_arg_dat(wk25, 1, stencil0, "double", OPS_WRITE));


         int iter_range34[] = {0, nx0, 0, nx1, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_34_kernel, "D(u1[x0 x1 x2 t] x1)", taylor_green_vortex_block, 3, iter_range34,
         ops_arg_dat(u1, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk26, 1, stencil0, "double", OPS_WRITE));


         int iter_range35[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_35_kernel, "D(p x0)", taylor_green_vortex_block, 3, iter_range35,
         ops_arg_dat(p, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk27, 1, stencil0, "double", OPS_WRITE));


         int iter_range36[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_36_kernel, "D(T x1 x1)", taylor_green_vortex_block, 3, iter_range36,
         ops_arg_dat(T, 1, stencil2, "double", OPS_READ),
         ops_arg_dat(wk28, 1, stencil0, "double", OPS_WRITE));


         int iter_range37[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_37_kernel, "rhou0*u1[x0 x1 x2 t]", taylor_green_vortex_block, 3, iter_range37,
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range38[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_38_kernel, "D(rhou0*u1[x0 x1 x2 t] x1)", taylor_green_vortex_block, 3, iter_range38,
         ops_arg_dat(wk63, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk29, 1, stencil0, "double", OPS_WRITE));


         int iter_range39[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_39_kernel, "rhou1*u2", taylor_green_vortex_block, 3, iter_range39,
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range40[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_40_kernel, "D(rhou1*u2 x2)", taylor_green_vortex_block, 3, iter_range40,
         ops_arg_dat(wk63, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk30, 1, stencil0, "double", OPS_WRITE));


         int iter_range41[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_41_kernel, "D(u2 x1 x1)", taylor_green_vortex_block, 3, iter_range41,
         ops_arg_dat(u2, 1, stencil2, "double", OPS_READ),
         ops_arg_dat(wk31, 1, stencil0, "double", OPS_WRITE));


         int iter_range42[] = {0, nx0, 0, nx1, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_42_kernel, "D(u2 x0)", taylor_green_vortex_block, 3, iter_range42,
         ops_arg_dat(u2, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk32, 1, stencil0, "double", OPS_WRITE));


         int iter_range43[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_43_kernel, "p*u2", taylor_green_vortex_block, 3, iter_range43,
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(p, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range44[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_44_kernel, "D(p*u2 x2)", taylor_green_vortex_block, 3, iter_range44,
         ops_arg_dat(wk63, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk33, 1, stencil0, "double", OPS_WRITE));


         int iter_range45[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_45_kernel, "D(rhoE x0)", taylor_green_vortex_block, 3, iter_range45,
         ops_arg_dat(rhoE, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk34, 1, stencil0, "double", OPS_WRITE));


         int iter_range46[] = {0, nx0, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_46_kernel, "D(u0 x0)", taylor_green_vortex_block, 3, iter_range46,
         ops_arg_dat(u0, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk35, 1, stencil0, "double", OPS_WRITE));


         int iter_range47[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_47_kernel, "D(rhou2 x2)", taylor_green_vortex_block, 3, iter_range47,
         ops_arg_dat(rhou2, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk36, 1, stencil0, "double", OPS_WRITE));


         int iter_range48[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_48_kernel, "rhou2*u2", taylor_green_vortex_block, 3, iter_range48,
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range49[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_49_kernel, "D(rhou2*u2 x2)", taylor_green_vortex_block, 3, iter_range49,
         ops_arg_dat(wk63, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk37, 1, stencil0, "double", OPS_WRITE));


         int iter_range50[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_50_kernel, "D(u1[x0 x1 x2 t] x1 x1)", taylor_green_vortex_block, 3, iter_range50,
         ops_arg_dat(u1, 1, stencil2, "double", OPS_READ),
         ops_arg_dat(wk38, 1, stencil0, "double", OPS_WRITE));


         int iter_range51[] = {0, nx0, -2, nx1 + 2, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_51_kernel, "D(u1[x0 x1 x2 t] x0)", taylor_green_vortex_block, 3, iter_range51,
         ops_arg_dat(u1, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk39, 1, stencil0, "double", OPS_WRITE));


         int iter_range52[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_52_kernel, "D(p x1)", taylor_green_vortex_block, 3, iter_range52,
         ops_arg_dat(p, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk40, 1, stencil0, "double", OPS_WRITE));


         int iter_range53[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_53_kernel, "rhou2*u1", taylor_green_vortex_block, 3, iter_range53,
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range54[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_54_kernel, "D(rhou2*u1 x1)", taylor_green_vortex_block, 3, iter_range54,
         ops_arg_dat(wk63, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk41, 1, stencil0, "double", OPS_WRITE));


         int iter_range55[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_55_kernel, "rhoE*u1[x0 x1 x2 t]", taylor_green_vortex_block, 3, iter_range55,
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhoE, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range56[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_56_kernel, "D(rhoE*u1[x0 x1 x2 t] x1)", taylor_green_vortex_block, 3, iter_range56,
         ops_arg_dat(wk63, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk42, 1, stencil0, "double", OPS_WRITE));


         int iter_range57[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_57_kernel, "rhou0*u0", taylor_green_vortex_block, 3, iter_range57,
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range58[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_58_kernel, "D(rhou0*u0 x0)", taylor_green_vortex_block, 3, iter_range58,
         ops_arg_dat(wk63, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk43, 1, stencil0, "double", OPS_WRITE));


         int iter_range59[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_59_kernel, "D(T x2 x2)", taylor_green_vortex_block, 3, iter_range59,
         ops_arg_dat(T, 1, stencil5, "double", OPS_READ),
         ops_arg_dat(wk44, 1, stencil0, "double", OPS_WRITE));


         int iter_range60[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_60_kernel, "D(rho x2)", taylor_green_vortex_block, 3, iter_range60,
         ops_arg_dat(rho, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk45, 1, stencil0, "double", OPS_WRITE));


         int iter_range61[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_61_kernel, "rhoE*u2", taylor_green_vortex_block, 3, iter_range61,
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhoE, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range62[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_62_kernel, "D(rhoE*u2 x2)", taylor_green_vortex_block, 3, iter_range62,
         ops_arg_dat(wk63, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk46, 1, stencil0, "double", OPS_WRITE));


         int iter_range63[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_63_kernel, "D(u2 x2)", taylor_green_vortex_block, 3, iter_range63,
         ops_arg_dat(u2, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk47, 1, stencil0, "double", OPS_WRITE));


         int iter_range64[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_64_kernel, "D(T x0 x0)", taylor_green_vortex_block, 3, iter_range64,
         ops_arg_dat(T, 1, stencil6, "double", OPS_READ),
         ops_arg_dat(wk48, 1, stencil0, "double", OPS_WRITE));


         int iter_range65[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_65_kernel, "rhou2[x0 x1 x2 t]*u0", taylor_green_vortex_block, 3, iter_range65,
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range66[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_66_kernel, "D(rhou2[x0 x1 x2 t]*u0 x0)", taylor_green_vortex_block, 3, iter_range66,
         ops_arg_dat(wk63, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk49, 1, stencil0, "double", OPS_WRITE));


         int iter_range67[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_67_kernel, "D(u1[x0 x1 x2 t] x2)", taylor_green_vortex_block, 3, iter_range67,
         ops_arg_dat(u1, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk50, 1, stencil0, "double", OPS_WRITE));


         int iter_range68[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_68_kernel, "rhoE*u0", taylor_green_vortex_block, 3, iter_range68,
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhoE, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range69[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_69_kernel, "D(rhoE*u0 x0)", taylor_green_vortex_block, 3, iter_range69,
         ops_arg_dat(wk63, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk51, 1, stencil0, "double", OPS_WRITE));


         int iter_range70[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_70_kernel, "rhou0*u2[x0 x1 x2 t]", taylor_green_vortex_block, 3, iter_range70,
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range71[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_71_kernel, "D(rhou0*u2[x0 x1 x2 t] x2)", taylor_green_vortex_block, 3, iter_range71,
         ops_arg_dat(wk63, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk52, 1, stencil0, "double", OPS_WRITE));


         int iter_range72[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_72_kernel, "p[x0 x1 x2 t]*u0", taylor_green_vortex_block, 3, iter_range72,
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(p, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range73[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_73_kernel, "D(p[x0 x1 x2 t]*u0 x0)", taylor_green_vortex_block, 3, iter_range73,
         ops_arg_dat(wk63, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk53, 1, stencil0, "double", OPS_WRITE));


         int iter_range74[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_74_kernel, "rhou1*u1[x0 x1 x2 t]", taylor_green_vortex_block, 3, iter_range74,
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range75[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_75_kernel, "D(rhou1*u1[x0 x1 x2 t] x1)", taylor_green_vortex_block, 3, iter_range75,
         ops_arg_dat(wk63, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk54, 1, stencil0, "double", OPS_WRITE));


         int iter_range76[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_76_kernel, "D(rhou2 x1)", taylor_green_vortex_block, 3, iter_range76,
         ops_arg_dat(rhou2, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk55, 1, stencil0, "double", OPS_WRITE));


         int iter_range77[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_77_kernel, "D(rhou0 x2)", taylor_green_vortex_block, 3, iter_range77,
         ops_arg_dat(rhou0, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk56, 1, stencil0, "double", OPS_WRITE));


         int iter_range78[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_78_kernel, "D(u0 x0 x1)", taylor_green_vortex_block, 3, iter_range78,
         ops_arg_dat(wk35, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk57, 1, stencil0, "double", OPS_WRITE));


         int iter_range79[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_79_kernel, "D(u2 x1 x2)", taylor_green_vortex_block, 3, iter_range79,
         ops_arg_dat(wk14, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk58, 1, stencil0, "double", OPS_WRITE));


         int iter_range80[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_80_kernel, "D(u1[x0 x1 x2 t] x0 x1)", taylor_green_vortex_block, 3, iter_range80,
         ops_arg_dat(wk39, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk59, 1, stencil0, "double", OPS_WRITE));


         int iter_range81[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_81_kernel, "D(u1[x0 x1 x2 t] x1 x2)", taylor_green_vortex_block, 3, iter_range81,
         ops_arg_dat(wk26, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk60, 1, stencil0, "double", OPS_WRITE));


         int iter_range82[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_82_kernel, "D(u0 x0 x2)", taylor_green_vortex_block, 3, iter_range82,
         ops_arg_dat(wk35, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk61, 1, stencil0, "double", OPS_WRITE));


         int iter_range83[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_83_kernel, "D(u2 x0 x2)", taylor_green_vortex_block, 3, iter_range83,
         ops_arg_dat(wk32, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk62, 1, stencil0, "double", OPS_WRITE));


         int iter_range84[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_84_kernel, "Residual of equation", taylor_green_vortex_block, 3, iter_range84,
         ops_arg_dat(wk58, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk20, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk47, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk21, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk28, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk61, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk52, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk29, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk55, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk19, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk59, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk35, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk18, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk11, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk12, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk31, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk8, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk37, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk34, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk42, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk10, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk30, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk39, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk44, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk40, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk15, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk46, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk45, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk62, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk48, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk25, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk3, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk7, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk54, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk33, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk6, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk50, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk32, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk38, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk14, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhoE, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk60, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk26, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk43, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk49, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk22, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk24, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk27, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk5, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk23, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk9, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk36, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk4, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk57, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk41, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk51, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk17, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk56, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk13, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk53, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk16, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk65, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(wk64, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(wk67, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(wk66, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_WRITE));


         int iter_range85[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_85_kernel, "RK new (subloop) update", taylor_green_vortex_block, 3, iter_range85,
         ops_arg_dat(wk64, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk67, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou1_old, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk65, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhoE_old, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2_old, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk66, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0_old, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho_old, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(rhoE, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_WRITE),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_WRITE),
         ops_arg_gbl(&rknew[stage], 1, "double", OPS_READ));


         int iter_range86[] = {-2, nx0 + 2, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_86_kernel, "RK old update", taylor_green_vortex_block, 3, iter_range86,
         ops_arg_dat(wk64, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk67, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk65, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk66, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk63, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou1_old, 1, stencil0, "double", OPS_RW),
         ops_arg_dat(rhou2_old, 1, stencil0, "double", OPS_RW),
         ops_arg_dat(rhou0_old, 1, stencil0, "double", OPS_RW),
         ops_arg_dat(rho_old, 1, stencil0, "double", OPS_RW),
         ops_arg_dat(rhoE_old, 1, stencil0, "double", OPS_RW),
         ops_arg_gbl(&rkold[stage], 1, "double", OPS_READ));



         // Boundary condition exchange calls
         ops_halo_transfer(halo_exchange0);
         // Boundary condition exchange calls
         ops_halo_transfer(halo_exchange1);
         // Boundary condition exchange calls
         ops_halo_transfer(halo_exchange2);
         // Boundary condition exchange calls
         ops_halo_transfer(halo_exchange3);
         // Boundary condition exchange calls
         ops_halo_transfer(halo_exchange4);
         // Boundary condition exchange calls
         ops_halo_transfer(halo_exchange5);

      }



      if(fmod(1.0*iteration + 1.0,739.000000000000) == 0)
      {
         char buf[100];
         sprintf(buf,"taylor_green_vortex_%d.h5",iteration);
         ops_fetch_block_hdf5_file(taylor_green_vortex_block, buf);
         ops_fetch_dat_hdf5_file(rho, buf);
         ops_fetch_dat_hdf5_file(rhou0, buf);
         ops_fetch_dat_hdf5_file(rhou1, buf);
         ops_fetch_dat_hdf5_file(rhou2, buf);
         ops_fetch_dat_hdf5_file(rhoE, buf);
      }if(fmod(iteration,100) == 0)
      {
         int iter_range89[] = {-2, nx0 + 2, -2, nx1 + 2, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_89_kernel, "Non-Grouped Formula Evaluation", taylor_green_vortex_block, 3, iter_range89,
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u2, 1, stencil0, "double", OPS_WRITE));


         int iter_range90[] = {0, nx0, -2, nx1 + 2, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_90_kernel, "Non-Grouped Formula Evaluation", taylor_green_vortex_block, 3, iter_range90,
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rhou0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u0, 1, stencil0, "double", OPS_WRITE));


         int iter_range91[] = {-2, nx0 + 2, 0, nx1, -2, nx2 + 2};
         ops_par_loop(taylor_green_vortex_block0_91_kernel, "Non-Grouped Formula Evaluation", taylor_green_vortex_block, 3, iter_range91,
         ops_arg_dat(rhou1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u1, 1, stencil0, "double", OPS_WRITE));


         int iter_range92[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_92_kernel, "D(u1[x0 x1 x2 t] x0)", taylor_green_vortex_block, 3, iter_range92,
         ops_arg_dat(u1, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk0, 1, stencil0, "double", OPS_WRITE));


         int iter_range93[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_93_kernel, "D(u0 x1)", taylor_green_vortex_block, 3, iter_range93,
         ops_arg_dat(u0, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk1, 1, stencil0, "double", OPS_WRITE));


         int iter_range94[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_94_kernel, "D(u2 x1)", taylor_green_vortex_block, 3, iter_range94,
         ops_arg_dat(u2, 1, stencil3, "double", OPS_READ),
         ops_arg_dat(wk2, 1, stencil0, "double", OPS_WRITE));


         int iter_range95[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_95_kernel, "D(u0 x2)", taylor_green_vortex_block, 3, iter_range95,
         ops_arg_dat(u0, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk3, 1, stencil0, "double", OPS_WRITE));


         int iter_range96[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_96_kernel, "D(u2 x0)", taylor_green_vortex_block, 3, iter_range96,
         ops_arg_dat(u2, 1, stencil4, "double", OPS_READ),
         ops_arg_dat(wk4, 1, stencil0, "double", OPS_WRITE));


         int iter_range97[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_97_kernel, "D(u1[x0 x1 x2 t] x2)", taylor_green_vortex_block, 3, iter_range97,
         ops_arg_dat(u1, 1, stencil1, "double", OPS_READ),
         ops_arg_dat(wk5, 1, stencil0, "double", OPS_WRITE));


         int iter_range98[] = {0, nx0, 0, nx1, 0, nx2};
         ops_par_loop(taylor_green_vortex_block0_98_kernel, "Reduction equations", taylor_green_vortex_block, 3, iter_range98,
         ops_arg_dat(wk3, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u0, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk5, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(rho, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk4, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u1, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(u2, 1, stencil0, "double", OPS_READ),
         ops_arg_dat(wk0, 1, stencil0, "double", OPS_READ),
         ops_arg_reduce(ke, 1, "double", OPS_INC),
         ops_arg_reduce(enstrophy, 1, "double", OPS_INC),
         ops_arg_reduce(rhomean, 1, "double", OPS_INC));


         double ke_reduction = 0.0; 
         ops_reduction_result(ke, &ke_reduction);
         double enstrophy_reduction = 0.0; 
         ops_reduction_result(enstrophy, &enstrophy_reduction);
         double rhomean_reduction = 0.0; 
         ops_reduction_result(rhomean, &rhomean_reduction);
         ops_printf("%g, %g, %g, %g\n", (iteration + 1)*deltat, ke_reduction, enstrophy_reduction, rhomean_reduction);
      }

   }

   double cpu_end, elapsed_end;
   ops_timers(&cpu_end, &elapsed_end);

   ops_printf("\nTimings are:\n");
   ops_printf("-----------------------------------------\n");
   ops_printf("Total Wall time %lf\n",elapsed_end-elapsed_start);

   ops_fetch_block_hdf5_file(taylor_green_vortex_block, "taylor_green_vortex_5909.h5");
   ops_fetch_dat_hdf5_file(rho, "taylor_green_vortex_5909.h5");
   ops_fetch_dat_hdf5_file(rhou0, "taylor_green_vortex_5909.h5");
   ops_fetch_dat_hdf5_file(rhou1, "taylor_green_vortex_5909.h5");
   ops_fetch_dat_hdf5_file(rhou2, "taylor_green_vortex_5909.h5");
   ops_fetch_dat_hdf5_file(rhoE, "taylor_green_vortex_5909.h5");

   // Exit OPS 
   ops_exit();

}
