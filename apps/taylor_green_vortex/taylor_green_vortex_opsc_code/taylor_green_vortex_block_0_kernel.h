#ifndef block_0_KERNEL_H
#define block_0_KERNEL_H

void taylor_green_vortex_block0_0_kernel(const double *rhou1 , const double *rho , const double *rhou2 , const double *rhou0 , double *u0 , double *u1 , double *u2)
{
u0[OPS_ACC4(0,0,0)] = rhou0[OPS_ACC3(0,0,0)]/rho[OPS_ACC1(0,0,0)];
u2[OPS_ACC6(0,0,0)] = rhou2[OPS_ACC2(0,0,0)]/rho[OPS_ACC1(0,0,0)];
u1[OPS_ACC5(0,0,0)] = rhou1[OPS_ACC0(0,0,0)]/rho[OPS_ACC1(0,0,0)];
}


void taylor_green_vortex_block0_1_kernel(const double *u0 , const double *u1 , const double *u2 , const double *rho , const double *rhoE , double *p)
{
p[OPS_ACC5(0,0,0)] = (gama - 1)*(-rc0*(pow(u0[OPS_ACC0(0,0,0)], 2) + pow(u1[OPS_ACC1(0,0,0)], 2) + pow(u2[OPS_ACC2(0,0,0)], 2))*rho[OPS_ACC3(0,0,0)] + rhoE[OPS_ACC4(0,0,0)]);
}


void taylor_green_vortex_block0_2_kernel(const double *rho , const double *p , double *T)
{
T[OPS_ACC2(0,0,0)] = gama*p[OPS_ACC1(0,0,0)]*pow(Minf, 2)/rho[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_3_kernel(const double *rhoE , double *wk0)
{
wk0[OPS_ACC1(0,0,0)] = rinv1*((rc2)*rhoE[OPS_ACC0(0,0,-2)] - rc3*rhoE[OPS_ACC0(0,0,-1)] + (rc3)*rhoE[OPS_ACC0(0,0,1)] - rc2*rhoE[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_4_kernel(const double *u0 , double *wk1)
{
wk1[OPS_ACC1(0,0,0)] = rinv4*(-rc5*u0[OPS_ACC0(0,0,0)] - rc2*u0[OPS_ACC0(0,-2,0)] + (rc6)*u0[OPS_ACC0(0,-1,0)] + (rc6)*u0[OPS_ACC0(0,1,0)] - rc2*u0[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_5_kernel(const double *rhou1 , double *wk2)
{
wk2[OPS_ACC1(0,0,0)] = rinv1*((rc2)*rhou1[OPS_ACC0(0,0,-2)] - rc3*rhou1[OPS_ACC0(0,0,-1)] + (rc3)*rhou1[OPS_ACC0(0,0,1)] - rc2*rhou1[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_6_kernel(const double *rhou0 , double *wk3)
{
wk3[OPS_ACC1(0,0,0)] = rinv7*((rc2)*rhou0[OPS_ACC0(0,-2,0)] - rc3*rhou0[OPS_ACC0(0,-1,0)] + (rc3)*rhou0[OPS_ACC0(0,1,0)] - rc2*rhou0[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_7_kernel(const double *u0 , double *wk4)
{
wk4[OPS_ACC1(0,0,0)] = rinv1*((rc2)*u0[OPS_ACC0(0,0,-2)] - rc3*u0[OPS_ACC0(0,0,-1)] + (rc3)*u0[OPS_ACC0(0,0,1)] - rc2*u0[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_8_kernel(const double *rhou1 , double *wk5)
{
wk5[OPS_ACC1(0,0,0)] = rinv8*((rc2)*rhou1[OPS_ACC0(-2,0,0)] - rc3*rhou1[OPS_ACC0(-1,0,0)] + (rc3)*rhou1[OPS_ACC0(1,0,0)] - rc2*rhou1[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_9_kernel(const double *rhou0 , double *wk6)
{
wk6[OPS_ACC1(0,0,0)] = rinv8*((rc2)*rhou0[OPS_ACC0(-2,0,0)] - rc3*rhou0[OPS_ACC0(-1,0,0)] + (rc3)*rhou0[OPS_ACC0(1,0,0)] - rc2*rhou0[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_10_kernel(const double *u1 , const double *rho , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rho[OPS_ACC1(0,0,0)]*u1[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_11_kernel(const double *wk63 , double *wk7)
{
wk7[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk63[OPS_ACC0(0,-2,0)] - rc3*wk63[OPS_ACC0(0,-1,0)] + (rc3)*wk63[OPS_ACC0(0,1,0)] - rc2*wk63[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_12_kernel(const double *u0 , double *wk8)
{
wk8[OPS_ACC1(0,0,0)] = rinv9*(-rc5*u0[OPS_ACC0(0,0,0)] - rc2*u0[OPS_ACC0(0,0,-2)] + (rc6)*u0[OPS_ACC0(0,0,-1)] + (rc6)*u0[OPS_ACC0(0,0,1)] - rc2*u0[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_13_kernel(const double *u2 , double *wk9)
{
wk9[OPS_ACC1(0,0,0)] = rinv9*(-rc5*u2[OPS_ACC0(0,0,0)] - rc2*u2[OPS_ACC0(0,0,-2)] + (rc6)*u2[OPS_ACC0(0,0,-1)] + (rc6)*u2[OPS_ACC0(0,0,1)] - rc2*u2[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_14_kernel(const double *rho , double *wk10)
{
wk10[OPS_ACC1(0,0,0)] = rinv7*((rc2)*rho[OPS_ACC0(0,-2,0)] - rc3*rho[OPS_ACC0(0,-1,0)] + (rc3)*rho[OPS_ACC0(0,1,0)] - rc2*rho[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_15_kernel(const double *u1 , double *wk11)
{
wk11[OPS_ACC1(0,0,0)] = rinv9*(-rc5*u1[OPS_ACC0(0,0,0)] - rc2*u1[OPS_ACC0(0,0,-2)] + (rc6)*u1[OPS_ACC0(0,0,-1)] + (rc6)*u1[OPS_ACC0(0,0,1)] - rc2*u1[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_16_kernel(const double *u1 , const double *p , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = p[OPS_ACC1(0,0,0)]*u1[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_17_kernel(const double *wk63 , double *wk12)
{
wk12[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk63[OPS_ACC0(0,-2,0)] - rc3*wk63[OPS_ACC0(0,-1,0)] + (rc3)*wk63[OPS_ACC0(0,1,0)] - rc2*wk63[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_18_kernel(const double *u2 , double *wk13)
{
wk13[OPS_ACC1(0,0,0)] = rinv10*(-rc5*u2[OPS_ACC0(0,0,0)] - rc2*u2[OPS_ACC0(-2,0,0)] + (rc6)*u2[OPS_ACC0(-1,0,0)] + (rc6)*u2[OPS_ACC0(1,0,0)] - rc2*u2[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_19_kernel(const double *u2 , double *wk14)
{
wk14[OPS_ACC1(0,0,0)] = rinv7*((rc2)*u2[OPS_ACC0(0,-2,0)] - rc3*u2[OPS_ACC0(0,-1,0)] + (rc3)*u2[OPS_ACC0(0,1,0)] - rc2*u2[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_20_kernel(const double *rhou1 , const double *u0 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou1[OPS_ACC0(0,0,0)]*u0[OPS_ACC1(0,0,0)];
}


void taylor_green_vortex_block0_21_kernel(const double *wk63 , double *wk15)
{
wk15[OPS_ACC1(0,0,0)] = rinv8*((rc2)*wk63[OPS_ACC0(-2,0,0)] - rc3*wk63[OPS_ACC0(-1,0,0)] + (rc3)*wk63[OPS_ACC0(1,0,0)] - rc2*wk63[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_22_kernel(const double *u0 , double *wk16)
{
wk16[OPS_ACC1(0,0,0)] = rinv7*((rc2)*u0[OPS_ACC0(0,-2,0)] - rc3*u0[OPS_ACC0(0,-1,0)] + (rc3)*u0[OPS_ACC0(0,1,0)] - rc2*u0[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_23_kernel(const double *p , double *wk17)
{
wk17[OPS_ACC1(0,0,0)] = rinv1*((rc2)*p[OPS_ACC0(0,0,-2)] - rc3*p[OPS_ACC0(0,0,-1)] + (rc3)*p[OPS_ACC0(0,0,1)] - rc2*p[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_24_kernel(const double *rhou2 , double *wk18)
{
wk18[OPS_ACC1(0,0,0)] = rinv8*((rc2)*rhou2[OPS_ACC0(-2,0,0)] - rc3*rhou2[OPS_ACC0(-1,0,0)] + (rc3)*rhou2[OPS_ACC0(1,0,0)] - rc2*rhou2[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_25_kernel(const double *rho , double *wk19)
{
wk19[OPS_ACC1(0,0,0)] = rinv8*((rc2)*rho[OPS_ACC0(-2,0,0)] - rc3*rho[OPS_ACC0(-1,0,0)] + (rc3)*rho[OPS_ACC0(1,0,0)] - rc2*rho[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_26_kernel(const double *rhoE , double *wk20)
{
wk20[OPS_ACC1(0,0,0)] = rinv7*((rc2)*rhoE[OPS_ACC0(0,-2,0)] - rc3*rhoE[OPS_ACC0(0,-1,0)] + (rc3)*rhoE[OPS_ACC0(0,1,0)] - rc2*rhoE[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_27_kernel(const double *u0 , const double *rho , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rho[OPS_ACC1(0,0,0)]*u0[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_28_kernel(const double *wk63 , double *wk21)
{
wk21[OPS_ACC1(0,0,0)] = rinv8*((rc2)*wk63[OPS_ACC0(-2,0,0)] - rc3*wk63[OPS_ACC0(-1,0,0)] + (rc3)*wk63[OPS_ACC0(1,0,0)] - rc2*wk63[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_29_kernel(const double *rhou1 , double *wk22)
{
wk22[OPS_ACC1(0,0,0)] = rinv7*((rc2)*rhou1[OPS_ACC0(0,-2,0)] - rc3*rhou1[OPS_ACC0(0,-1,0)] + (rc3)*rhou1[OPS_ACC0(0,1,0)] - rc2*rhou1[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_30_kernel(const double *u1 , double *wk23)
{
wk23[OPS_ACC1(0,0,0)] = rinv10*(-rc5*u1[OPS_ACC0(0,0,0)] - rc2*u1[OPS_ACC0(-2,0,0)] + (rc6)*u1[OPS_ACC0(-1,0,0)] + (rc6)*u1[OPS_ACC0(1,0,0)] - rc2*u1[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_31_kernel(const double *u2 , const double *rho , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rho[OPS_ACC1(0,0,0)]*u2[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_32_kernel(const double *wk63 , double *wk24)
{
wk24[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk63[OPS_ACC0(0,0,-2)] - rc3*wk63[OPS_ACC0(0,0,-1)] + (rc3)*wk63[OPS_ACC0(0,0,1)] - rc2*wk63[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_33_kernel(const double *u0 , double *wk25)
{
wk25[OPS_ACC1(0,0,0)] = rinv10*(-rc5*u0[OPS_ACC0(0,0,0)] - rc2*u0[OPS_ACC0(-2,0,0)] + (rc6)*u0[OPS_ACC0(-1,0,0)] + (rc6)*u0[OPS_ACC0(1,0,0)] - rc2*u0[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_34_kernel(const double *u1 , double *wk26)
{
wk26[OPS_ACC1(0,0,0)] = rinv7*((rc2)*u1[OPS_ACC0(0,-2,0)] - rc3*u1[OPS_ACC0(0,-1,0)] + (rc3)*u1[OPS_ACC0(0,1,0)] - rc2*u1[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_35_kernel(const double *p , double *wk27)
{
wk27[OPS_ACC1(0,0,0)] = rinv8*((rc2)*p[OPS_ACC0(-2,0,0)] - rc3*p[OPS_ACC0(-1,0,0)] + (rc3)*p[OPS_ACC0(1,0,0)] - rc2*p[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_36_kernel(const double *T , double *wk28)
{
wk28[OPS_ACC1(0,0,0)] = rinv4*(-rc5*T[OPS_ACC0(0,0,0)] - rc2*T[OPS_ACC0(0,-2,0)] + (rc6)*T[OPS_ACC0(0,-1,0)] + (rc6)*T[OPS_ACC0(0,1,0)] - rc2*T[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_37_kernel(const double *u1 , const double *rhou0 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou0[OPS_ACC1(0,0,0)]*u1[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_38_kernel(const double *wk63 , double *wk29)
{
wk29[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk63[OPS_ACC0(0,-2,0)] - rc3*wk63[OPS_ACC0(0,-1,0)] + (rc3)*wk63[OPS_ACC0(0,1,0)] - rc2*wk63[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_39_kernel(const double *rhou1 , const double *u2 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou1[OPS_ACC0(0,0,0)]*u2[OPS_ACC1(0,0,0)];
}


void taylor_green_vortex_block0_40_kernel(const double *wk63 , double *wk30)
{
wk30[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk63[OPS_ACC0(0,0,-2)] - rc3*wk63[OPS_ACC0(0,0,-1)] + (rc3)*wk63[OPS_ACC0(0,0,1)] - rc2*wk63[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_41_kernel(const double *u2 , double *wk31)
{
wk31[OPS_ACC1(0,0,0)] = rinv4*(-rc5*u2[OPS_ACC0(0,0,0)] - rc2*u2[OPS_ACC0(0,-2,0)] + (rc6)*u2[OPS_ACC0(0,-1,0)] + (rc6)*u2[OPS_ACC0(0,1,0)] - rc2*u2[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_42_kernel(const double *u2 , double *wk32)
{
wk32[OPS_ACC1(0,0,0)] = rinv8*((rc2)*u2[OPS_ACC0(-2,0,0)] - rc3*u2[OPS_ACC0(-1,0,0)] + (rc3)*u2[OPS_ACC0(1,0,0)] - rc2*u2[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_43_kernel(const double *u2 , const double *p , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = p[OPS_ACC1(0,0,0)]*u2[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_44_kernel(const double *wk63 , double *wk33)
{
wk33[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk63[OPS_ACC0(0,0,-2)] - rc3*wk63[OPS_ACC0(0,0,-1)] + (rc3)*wk63[OPS_ACC0(0,0,1)] - rc2*wk63[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_45_kernel(const double *rhoE , double *wk34)
{
wk34[OPS_ACC1(0,0,0)] = rinv8*((rc2)*rhoE[OPS_ACC0(-2,0,0)] - rc3*rhoE[OPS_ACC0(-1,0,0)] + (rc3)*rhoE[OPS_ACC0(1,0,0)] - rc2*rhoE[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_46_kernel(const double *u0 , double *wk35)
{
wk35[OPS_ACC1(0,0,0)] = rinv8*((rc2)*u0[OPS_ACC0(-2,0,0)] - rc3*u0[OPS_ACC0(-1,0,0)] + (rc3)*u0[OPS_ACC0(1,0,0)] - rc2*u0[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_47_kernel(const double *rhou2 , double *wk36)
{
wk36[OPS_ACC1(0,0,0)] = rinv1*((rc2)*rhou2[OPS_ACC0(0,0,-2)] - rc3*rhou2[OPS_ACC0(0,0,-1)] + (rc3)*rhou2[OPS_ACC0(0,0,1)] - rc2*rhou2[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_48_kernel(const double *u2 , const double *rhou2 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou2[OPS_ACC1(0,0,0)]*u2[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_49_kernel(const double *wk63 , double *wk37)
{
wk37[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk63[OPS_ACC0(0,0,-2)] - rc3*wk63[OPS_ACC0(0,0,-1)] + (rc3)*wk63[OPS_ACC0(0,0,1)] - rc2*wk63[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_50_kernel(const double *u1 , double *wk38)
{
wk38[OPS_ACC1(0,0,0)] = rinv4*(-rc5*u1[OPS_ACC0(0,0,0)] - rc2*u1[OPS_ACC0(0,-2,0)] + (rc6)*u1[OPS_ACC0(0,-1,0)] + (rc6)*u1[OPS_ACC0(0,1,0)] - rc2*u1[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_51_kernel(const double *u1 , double *wk39)
{
wk39[OPS_ACC1(0,0,0)] = rinv8*((rc2)*u1[OPS_ACC0(-2,0,0)] - rc3*u1[OPS_ACC0(-1,0,0)] + (rc3)*u1[OPS_ACC0(1,0,0)] - rc2*u1[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_52_kernel(const double *p , double *wk40)
{
wk40[OPS_ACC1(0,0,0)] = rinv7*((rc2)*p[OPS_ACC0(0,-2,0)] - rc3*p[OPS_ACC0(0,-1,0)] + (rc3)*p[OPS_ACC0(0,1,0)] - rc2*p[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_53_kernel(const double *u1 , const double *rhou2 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou2[OPS_ACC1(0,0,0)]*u1[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_54_kernel(const double *wk63 , double *wk41)
{
wk41[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk63[OPS_ACC0(0,-2,0)] - rc3*wk63[OPS_ACC0(0,-1,0)] + (rc3)*wk63[OPS_ACC0(0,1,0)] - rc2*wk63[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_55_kernel(const double *u1 , const double *rhoE , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhoE[OPS_ACC1(0,0,0)]*u1[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_56_kernel(const double *wk63 , double *wk42)
{
wk42[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk63[OPS_ACC0(0,-2,0)] - rc3*wk63[OPS_ACC0(0,-1,0)] + (rc3)*wk63[OPS_ACC0(0,1,0)] - rc2*wk63[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_57_kernel(const double *u0 , const double *rhou0 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou0[OPS_ACC1(0,0,0)]*u0[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_58_kernel(const double *wk63 , double *wk43)
{
wk43[OPS_ACC1(0,0,0)] = rinv8*((rc2)*wk63[OPS_ACC0(-2,0,0)] - rc3*wk63[OPS_ACC0(-1,0,0)] + (rc3)*wk63[OPS_ACC0(1,0,0)] - rc2*wk63[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_59_kernel(const double *T , double *wk44)
{
wk44[OPS_ACC1(0,0,0)] = rinv9*(-rc5*T[OPS_ACC0(0,0,0)] - rc2*T[OPS_ACC0(0,0,-2)] + (rc6)*T[OPS_ACC0(0,0,-1)] + (rc6)*T[OPS_ACC0(0,0,1)] - rc2*T[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_60_kernel(const double *rho , double *wk45)
{
wk45[OPS_ACC1(0,0,0)] = rinv1*((rc2)*rho[OPS_ACC0(0,0,-2)] - rc3*rho[OPS_ACC0(0,0,-1)] + (rc3)*rho[OPS_ACC0(0,0,1)] - rc2*rho[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_61_kernel(const double *u2 , const double *rhoE , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhoE[OPS_ACC1(0,0,0)]*u2[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_62_kernel(const double *wk63 , double *wk46)
{
wk46[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk63[OPS_ACC0(0,0,-2)] - rc3*wk63[OPS_ACC0(0,0,-1)] + (rc3)*wk63[OPS_ACC0(0,0,1)] - rc2*wk63[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_63_kernel(const double *u2 , double *wk47)
{
wk47[OPS_ACC1(0,0,0)] = rinv1*((rc2)*u2[OPS_ACC0(0,0,-2)] - rc3*u2[OPS_ACC0(0,0,-1)] + (rc3)*u2[OPS_ACC0(0,0,1)] - rc2*u2[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_64_kernel(const double *T , double *wk48)
{
wk48[OPS_ACC1(0,0,0)] = rinv10*(-rc5*T[OPS_ACC0(0,0,0)] - rc2*T[OPS_ACC0(-2,0,0)] + (rc6)*T[OPS_ACC0(-1,0,0)] + (rc6)*T[OPS_ACC0(1,0,0)] - rc2*T[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_65_kernel(const double *u0 , const double *rhou2 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou2[OPS_ACC1(0,0,0)]*u0[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_66_kernel(const double *wk63 , double *wk49)
{
wk49[OPS_ACC1(0,0,0)] = rinv8*((rc2)*wk63[OPS_ACC0(-2,0,0)] - rc3*wk63[OPS_ACC0(-1,0,0)] + (rc3)*wk63[OPS_ACC0(1,0,0)] - rc2*wk63[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_67_kernel(const double *u1 , double *wk50)
{
wk50[OPS_ACC1(0,0,0)] = rinv1*((rc2)*u1[OPS_ACC0(0,0,-2)] - rc3*u1[OPS_ACC0(0,0,-1)] + (rc3)*u1[OPS_ACC0(0,0,1)] - rc2*u1[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_68_kernel(const double *u0 , const double *rhoE , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhoE[OPS_ACC1(0,0,0)]*u0[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_69_kernel(const double *wk63 , double *wk51)
{
wk51[OPS_ACC1(0,0,0)] = rinv8*((rc2)*wk63[OPS_ACC0(-2,0,0)] - rc3*wk63[OPS_ACC0(-1,0,0)] + (rc3)*wk63[OPS_ACC0(1,0,0)] - rc2*wk63[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_70_kernel(const double *u2 , const double *rhou0 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou0[OPS_ACC1(0,0,0)]*u2[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_71_kernel(const double *wk63 , double *wk52)
{
wk52[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk63[OPS_ACC0(0,0,-2)] - rc3*wk63[OPS_ACC0(0,0,-1)] + (rc3)*wk63[OPS_ACC0(0,0,1)] - rc2*wk63[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_72_kernel(const double *u0 , const double *p , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = p[OPS_ACC1(0,0,0)]*u0[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_73_kernel(const double *wk63 , double *wk53)
{
wk53[OPS_ACC1(0,0,0)] = rinv8*((rc2)*wk63[OPS_ACC0(-2,0,0)] - rc3*wk63[OPS_ACC0(-1,0,0)] + (rc3)*wk63[OPS_ACC0(1,0,0)] - rc2*wk63[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_74_kernel(const double *rhou1 , const double *u1 , double *wk63)
{
wk63[OPS_ACC2(0,0,0)] = rhou1[OPS_ACC0(0,0,0)]*u1[OPS_ACC1(0,0,0)];
}


void taylor_green_vortex_block0_75_kernel(const double *wk63 , double *wk54)
{
wk54[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk63[OPS_ACC0(0,-2,0)] - rc3*wk63[OPS_ACC0(0,-1,0)] + (rc3)*wk63[OPS_ACC0(0,1,0)] - rc2*wk63[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_76_kernel(const double *rhou2 , double *wk55)
{
wk55[OPS_ACC1(0,0,0)] = rinv7*((rc2)*rhou2[OPS_ACC0(0,-2,0)] - rc3*rhou2[OPS_ACC0(0,-1,0)] + (rc3)*rhou2[OPS_ACC0(0,1,0)] - rc2*rhou2[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_77_kernel(const double *rhou0 , double *wk56)
{
wk56[OPS_ACC1(0,0,0)] = rinv1*((rc2)*rhou0[OPS_ACC0(0,0,-2)] - rc3*rhou0[OPS_ACC0(0,0,-1)] + (rc3)*rhou0[OPS_ACC0(0,0,1)] - rc2*rhou0[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_78_kernel(const double *wk35 , double *wk57)
{
wk57[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk35[OPS_ACC0(0,-2,0)] - rc3*wk35[OPS_ACC0(0,-1,0)] + (rc3)*wk35[OPS_ACC0(0,1,0)] - rc2*wk35[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_79_kernel(const double *wk14 , double *wk58)
{
wk58[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk14[OPS_ACC0(0,0,-2)] - rc3*wk14[OPS_ACC0(0,0,-1)] + (rc3)*wk14[OPS_ACC0(0,0,1)] - rc2*wk14[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_80_kernel(const double *wk39 , double *wk59)
{
wk59[OPS_ACC1(0,0,0)] = rinv7*((rc2)*wk39[OPS_ACC0(0,-2,0)] - rc3*wk39[OPS_ACC0(0,-1,0)] + (rc3)*wk39[OPS_ACC0(0,1,0)] - rc2*wk39[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_81_kernel(const double *wk26 , double *wk60)
{
wk60[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk26[OPS_ACC0(0,0,-2)] - rc3*wk26[OPS_ACC0(0,0,-1)] + (rc3)*wk26[OPS_ACC0(0,0,1)] - rc2*wk26[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_82_kernel(const double *wk35 , double *wk61)
{
wk61[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk35[OPS_ACC0(0,0,-2)] - rc3*wk35[OPS_ACC0(0,0,-1)] + (rc3)*wk35[OPS_ACC0(0,0,1)] - rc2*wk35[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_83_kernel(const double *wk32 , double *wk62)
{
wk62[OPS_ACC1(0,0,0)] = rinv1*((rc2)*wk32[OPS_ACC0(0,0,-2)] - rc3*wk32[OPS_ACC0(0,0,-1)] + (rc3)*wk32[OPS_ACC0(0,0,1)] - rc2*wk32[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_84_kernel(const double *wk58 , const double *wk20 , const double *wk47 , const double *wk21 , const double *wk28 , const double *wk61 , const double *wk52 , const double *wk29 , const double *wk55 , const double *wk19 , const double *wk59 , const double *wk35 , const double *wk18 , const double *wk11 , const double *wk12 , const double *wk31 , const double *wk8 , const double *wk37 , const double *wk34 , const double *wk42 , const double *wk10 , const double *wk30 , const double *wk39 , const double *wk0 , const double *wk44 , const double *u0 , const double *wk40 , const double *wk15 , const double *wk46 , const double *wk45 , const double *wk62 , const double *u2 , const double *wk48 , const double *wk25 , const double *wk3 , const double *wk7 , const double *wk1 , const double *wk54 , const double *wk33 , const double *wk6 , const double *rho , const double *rhou2 , const double *wk50 , const double *wk32 , const double *wk38 , const double *wk14 , const double *rhoE , const double *u1 , const double *wk60 , const double *wk26 , const double *wk43 , const double *wk49 , const double *wk22 , const double *wk24 , const double *wk27 , const double *wk5 , const double *wk23 , const double *wk9 , const double *wk36 , const double *rhou1 , const double *wk4 , const double *wk2 , const double *wk57 , const double *wk41 , const double *wk51 , const double *rhou0 , const double *wk17 , const double *wk56 , const double *wk13 , const double *wk53 , const double *wk16 , double *wk65 , double *wk64 , double *wk67 , double *wk66 , double *wk63)
{
wk63[OPS_ACC75(0,0,0)] = -0.5*(wk26[OPS_ACC49(0,0,0)] + wk35[OPS_ACC11(0,0,0)] + wk47[OPS_ACC2(0,0,0)])*rho[OPS_ACC40(0,0,0)] - 0.5*wk10[OPS_ACC20(0,0,0)]*u1[OPS_ACC47(0,0,0)] - 0.5*wk19[OPS_ACC9(0,0,0)]*u0[OPS_ACC25(0,0,0)] - 0.5*wk21[OPS_ACC3(0,0,0)] - 0.5*wk24[OPS_ACC53(0,0,0)] - 0.5*wk45[OPS_ACC29(0,0,0)]*u2[OPS_ACC31(0,0,0)] - 0.5*wk7[OPS_ACC35(0,0,0)];
wk64[OPS_ACC72(0,0,0)] = 1.0*rinv11*(wk1[OPS_ACC36(0,0,0)] + wk59[OPS_ACC10(0,0,0)]) + 1.0*rinv11*(wk62[OPS_ACC30(0,0,0)] + wk8[OPS_ACC16(0,0,0)]) + 1.0*rinv11*((rc6)*wk25[OPS_ACC33(0,0,0)] - rc3*wk59[OPS_ACC10(0,0,0)] - rc3*wk62[OPS_ACC30(0,0,0)]) - 0.5*(wk26[OPS_ACC49(0,0,0)] + wk35[OPS_ACC11(0,0,0)] + wk47[OPS_ACC2(0,0,0)])*rhou0[OPS_ACC65(0,0,0)] - wk27[OPS_ACC54(0,0,0)] - 0.5*wk29[OPS_ACC7(0,0,0)] - 0.5*wk3[OPS_ACC34(0,0,0)]*u1[OPS_ACC47(0,0,0)] - 0.5*wk43[OPS_ACC50(0,0,0)] - 0.5*wk52[OPS_ACC6(0,0,0)] - 0.5*wk56[OPS_ACC67(0,0,0)]*u2[OPS_ACC31(0,0,0)] - 0.5*wk6[OPS_ACC39(0,0,0)]*u0[OPS_ACC25(0,0,0)];
wk65[OPS_ACC71(0,0,0)] = 1.0*rinv11*(wk11[OPS_ACC13(0,0,0)] + wk58[OPS_ACC0(0,0,0)]) + 1.0*rinv11*(wk23[OPS_ACC56(0,0,0)] + wk57[OPS_ACC62(0,0,0)]) + 1.0*rinv11*((rc6)*wk38[OPS_ACC44(0,0,0)] - rc3*wk57[OPS_ACC62(0,0,0)] - rc3*wk58[OPS_ACC0(0,0,0)]) - 0.5*(wk26[OPS_ACC49(0,0,0)] + wk35[OPS_ACC11(0,0,0)] + wk47[OPS_ACC2(0,0,0)])*rhou1[OPS_ACC59(0,0,0)] - 0.5*wk15[OPS_ACC27(0,0,0)] - 0.5*wk2[OPS_ACC61(0,0,0)]*u2[OPS_ACC31(0,0,0)] - 0.5*wk22[OPS_ACC52(0,0,0)]*u1[OPS_ACC47(0,0,0)] - 0.5*wk30[OPS_ACC21(0,0,0)] - wk40[OPS_ACC26(0,0,0)] - 0.5*wk5[OPS_ACC55(0,0,0)]*u0[OPS_ACC25(0,0,0)] - 0.5*wk54[OPS_ACC37(0,0,0)];
wk66[OPS_ACC74(0,0,0)] = 1.0*rinv11*(wk13[OPS_ACC68(0,0,0)] + wk61[OPS_ACC5(0,0,0)]) + 1.0*rinv11*(wk31[OPS_ACC15(0,0,0)] + wk60[OPS_ACC48(0,0,0)]) + 1.0*rinv11*(-rc3*wk60[OPS_ACC48(0,0,0)] - rc3*wk61[OPS_ACC5(0,0,0)] + (rc6)*wk9[OPS_ACC57(0,0,0)]) - 0.5*(wk26[OPS_ACC49(0,0,0)] + wk35[OPS_ACC11(0,0,0)] + wk47[OPS_ACC2(0,0,0)])*rhou2[OPS_ACC41(0,0,0)] - wk17[OPS_ACC66(0,0,0)] - 0.5*wk18[OPS_ACC12(0,0,0)]*u0[OPS_ACC25(0,0,0)] - 0.5*wk36[OPS_ACC58(0,0,0)]*u2[OPS_ACC31(0,0,0)] - 0.5*wk37[OPS_ACC17(0,0,0)] - 0.5*wk41[OPS_ACC63(0,0,0)] - 0.5*wk49[OPS_ACC51(0,0,0)] - 0.5*wk55[OPS_ACC8(0,0,0)]*u1[OPS_ACC47(0,0,0)];
wk67[OPS_ACC73(0,0,0)] = 1.0*rinv11*rinv12*rinv13*rinv14*wk28[OPS_ACC4(0,0,0)] + 1.0*rinv11*rinv12*rinv13*rinv14*wk44[OPS_ACC24(0,0,0)] + 1.0*rinv11*rinv12*rinv13*rinv14*wk48[OPS_ACC32(0,0,0)] + 1.0*rinv11*(wk1[OPS_ACC36(0,0,0)] + wk59[OPS_ACC10(0,0,0)])*u0[OPS_ACC25(0,0,0)] + 1.0*rinv11*(wk11[OPS_ACC13(0,0,0)] + wk58[OPS_ACC0(0,0,0)])*u1[OPS_ACC47(0,0,0)] + 1.0*rinv11*(wk13[OPS_ACC68(0,0,0)] + wk61[OPS_ACC5(0,0,0)])*u2[OPS_ACC31(0,0,0)] + 1.0*rinv11*(wk14[OPS_ACC45(0,0,0)] + wk50[OPS_ACC42(0,0,0)])*wk14[OPS_ACC45(0,0,0)] + 1.0*rinv11*(wk14[OPS_ACC45(0,0,0)] + wk50[OPS_ACC42(0,0,0)])*wk50[OPS_ACC42(0,0,0)] + 1.0*rinv11*(wk16[OPS_ACC70(0,0,0)] + wk39[OPS_ACC22(0,0,0)])*wk16[OPS_ACC70(0,0,0)] + 1.0*rinv11*(wk16[OPS_ACC70(0,0,0)] + wk39[OPS_ACC22(0,0,0)])*wk39[OPS_ACC22(0,0,0)] + 1.0*rinv11*(wk23[OPS_ACC56(0,0,0)] + wk57[OPS_ACC62(0,0,0)])*u1[OPS_ACC47(0,0,0)] + 1.0*rinv11*(wk31[OPS_ACC15(0,0,0)] + wk60[OPS_ACC48(0,0,0)])*u2[OPS_ACC31(0,0,0)] + 1.0*rinv11*(wk32[OPS_ACC43(0,0,0)] + wk4[OPS_ACC60(0,0,0)])*wk32[OPS_ACC43(0,0,0)] + 1.0*rinv11*(wk32[OPS_ACC43(0,0,0)] + wk4[OPS_ACC60(0,0,0)])*wk4[OPS_ACC60(0,0,0)] + 1.0*rinv11*(wk62[OPS_ACC30(0,0,0)] + wk8[OPS_ACC16(0,0,0)])*u0[OPS_ACC25(0,0,0)] + 1.0*rinv11*((rc6)*wk25[OPS_ACC33(0,0,0)] - rc3*wk59[OPS_ACC10(0,0,0)] - rc3*wk62[OPS_ACC30(0,0,0)])*u0[OPS_ACC25(0,0,0)] + 1.0*rinv11*(-rc3*wk26[OPS_ACC49(0,0,0)] - rc3*wk35[OPS_ACC11(0,0,0)] + (rc6)*wk47[OPS_ACC2(0,0,0)])*wk47[OPS_ACC2(0,0,0)] + 1.0*rinv11*(-rc3*wk26[OPS_ACC49(0,0,0)] + (rc6)*wk35[OPS_ACC11(0,0,0)] - rc3*wk47[OPS_ACC2(0,0,0)])*wk35[OPS_ACC11(0,0,0)] + 1.0*rinv11*((rc6)*wk26[OPS_ACC49(0,0,0)] - rc3*wk35[OPS_ACC11(0,0,0)] - rc3*wk47[OPS_ACC2(0,0,0)])*wk26[OPS_ACC49(0,0,0)] + 1.0*rinv11*((rc6)*wk38[OPS_ACC44(0,0,0)] - rc3*wk57[OPS_ACC62(0,0,0)] - rc3*wk58[OPS_ACC0(0,0,0)])*u1[OPS_ACC47(0,0,0)] + 1.0*rinv11*(-rc3*wk60[OPS_ACC48(0,0,0)] - rc3*wk61[OPS_ACC5(0,0,0)] + (rc6)*wk9[OPS_ACC57(0,0,0)])*u2[OPS_ACC31(0,0,0)] - 0.5*(wk26[OPS_ACC49(0,0,0)] + wk35[OPS_ACC11(0,0,0)] + wk47[OPS_ACC2(0,0,0)])*rhoE[OPS_ACC46(0,0,0)] - 0.5*wk0[OPS_ACC23(0,0,0)]*u2[OPS_ACC31(0,0,0)] - wk12[OPS_ACC14(0,0,0)] - 0.5*wk20[OPS_ACC1(0,0,0)]*u1[OPS_ACC47(0,0,0)] - wk33[OPS_ACC38(0,0,0)] - 0.5*wk34[OPS_ACC18(0,0,0)]*u0[OPS_ACC25(0,0,0)] - 0.5*wk42[OPS_ACC19(0,0,0)] - 0.5*wk46[OPS_ACC28(0,0,0)] - 0.5*wk51[OPS_ACC64(0,0,0)] - wk53[OPS_ACC69(0,0,0)];
}


void taylor_green_vortex_block0_85_kernel(const double *wk64 , const double *wk67 , const double *rhou1_old , const double *wk65 , const double *rhoE_old , const double *rhou2_old , const double *wk66 , const double *rhou0_old , const double *rho_old , const double *wk63 , double *rhou1 , double *rhoE , double *rho , double *rhou2 , double *rhou0 , const double *rknew)
{
rho[OPS_ACC12(0,0,0)] = deltat*rknew[0]*wk63[OPS_ACC9(0,0,0)] + rho_old[OPS_ACC8(0,0,0)];
rhou0[OPS_ACC14(0,0,0)] = deltat*rknew[0]*wk64[OPS_ACC0(0,0,0)] + rhou0_old[OPS_ACC7(0,0,0)];
rhou1[OPS_ACC10(0,0,0)] = deltat*rknew[0]*wk65[OPS_ACC3(0,0,0)] + rhou1_old[OPS_ACC2(0,0,0)];
rhou2[OPS_ACC13(0,0,0)] = deltat*rknew[0]*wk66[OPS_ACC6(0,0,0)] + rhou2_old[OPS_ACC5(0,0,0)];
rhoE[OPS_ACC11(0,0,0)] = deltat*rknew[0]*wk67[OPS_ACC1(0,0,0)] + rhoE_old[OPS_ACC4(0,0,0)];
}


void taylor_green_vortex_block0_86_kernel(const double *wk64 , const double *wk67 , const double *wk65 , const double *wk66 , const double *wk63 , double *rhou1_old , double *rhou2_old , double *rhou0_old , double *rho_old , double *rhoE_old , const double *rkold)
{
rho_old[OPS_ACC8(0,0,0)] = deltat*rkold[0]*wk63[OPS_ACC4(0,0,0)] + rho_old[OPS_ACC8(0,0,0)];
rhou0_old[OPS_ACC7(0,0,0)] = deltat*rkold[0]*wk64[OPS_ACC0(0,0,0)] + rhou0_old[OPS_ACC7(0,0,0)];
rhou1_old[OPS_ACC5(0,0,0)] = deltat*rkold[0]*wk65[OPS_ACC2(0,0,0)] + rhou1_old[OPS_ACC5(0,0,0)];
rhou2_old[OPS_ACC6(0,0,0)] = deltat*rkold[0]*wk66[OPS_ACC3(0,0,0)] + rhou2_old[OPS_ACC6(0,0,0)];
rhoE_old[OPS_ACC9(0,0,0)] = deltat*rkold[0]*wk67[OPS_ACC1(0,0,0)] + rhoE_old[OPS_ACC9(0,0,0)];
}


void taylor_green_vortex_block0_87_kernel(const double *rhou1 , const double *rhoE , const double *rho , const double *rhou2 , const double *rhou0 , double *rhou1_old , double *rhou2_old , double *rhou0_old , double *rho_old , double *rhoE_old)
{
rho_old[OPS_ACC8(0,0,0)] = rho[OPS_ACC2(0,0,0)];
rhou0_old[OPS_ACC7(0,0,0)] = rhou0[OPS_ACC4(0,0,0)];
rhou1_old[OPS_ACC5(0,0,0)] = rhou1[OPS_ACC0(0,0,0)];
rhou2_old[OPS_ACC6(0,0,0)] = rhou2[OPS_ACC3(0,0,0)];
rhoE_old[OPS_ACC9(0,0,0)] = rhoE[OPS_ACC1(0,0,0)];
}


void taylor_green_vortex_block0_88_kernel(double *rhou1 , double *rhoE , double *rho , double *rhou2 , double *rhou0 , const int *idx)
{
double x = deltai0*idx[0];
double y = deltai1*idx[1];
double z = deltai2*idx[2];
double u = sin(x)*cos(y)*cos(z);
double v = -cos(x)*sin(y)*cos(z);
double w = 0.0;
double p = 1.0*rinv15 + 0.0625*(cos(2.0*x) + cos(2.0*y))*(2.0 + cos(2.0*z));
double r = gama*pow(Minf, 2)*p;
rho[OPS_ACC2(0,0,0)] = r;
rhou0[OPS_ACC4(0,0,0)] = r*u;
rhou1[OPS_ACC0(0,0,0)] = r*v;
rhou2[OPS_ACC3(0,0,0)] = 0.0;
rhoE[OPS_ACC1(0,0,0)] = rinv13*p + 0.5*r*(pow(u, 2) + pow(v, 2) + pow(w, 2));
}


void taylor_green_vortex_block0_89_kernel(const double *rho , const double *rhou2 , double *u2)
{
u2[OPS_ACC2(0,0,0)] = rhou2[OPS_ACC1(0,0,0)]/rho[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_90_kernel(const double *rho , const double *rhou0 , double *u0)
{
u0[OPS_ACC2(0,0,0)] = rhou0[OPS_ACC1(0,0,0)]/rho[OPS_ACC0(0,0,0)];
}


void taylor_green_vortex_block0_91_kernel(const double *rhou1 , const double *rho , double *u1)
{
u1[OPS_ACC2(0,0,0)] = rhou1[OPS_ACC0(0,0,0)]/rho[OPS_ACC1(0,0,0)];
}


void taylor_green_vortex_block0_92_kernel(const double *u1 , double *wk0)
{
wk0[OPS_ACC1(0,0,0)] = rinv8*((rc2)*u1[OPS_ACC0(-2,0,0)] - rc3*u1[OPS_ACC0(-1,0,0)] + (rc3)*u1[OPS_ACC0(1,0,0)] - rc2*u1[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_93_kernel(const double *u0 , double *wk1)
{
wk1[OPS_ACC1(0,0,0)] = rinv7*((rc2)*u0[OPS_ACC0(0,-2,0)] - rc3*u0[OPS_ACC0(0,-1,0)] + (rc3)*u0[OPS_ACC0(0,1,0)] - rc2*u0[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_94_kernel(const double *u2 , double *wk2)
{
wk2[OPS_ACC1(0,0,0)] = rinv7*((rc2)*u2[OPS_ACC0(0,-2,0)] - rc3*u2[OPS_ACC0(0,-1,0)] + (rc3)*u2[OPS_ACC0(0,1,0)] - rc2*u2[OPS_ACC0(0,2,0)]);
}


void taylor_green_vortex_block0_95_kernel(const double *u0 , double *wk3)
{
wk3[OPS_ACC1(0,0,0)] = rinv1*((rc2)*u0[OPS_ACC0(0,0,-2)] - rc3*u0[OPS_ACC0(0,0,-1)] + (rc3)*u0[OPS_ACC0(0,0,1)] - rc2*u0[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_96_kernel(const double *u2 , double *wk4)
{
wk4[OPS_ACC1(0,0,0)] = rinv8*((rc2)*u2[OPS_ACC0(-2,0,0)] - rc3*u2[OPS_ACC0(-1,0,0)] + (rc3)*u2[OPS_ACC0(1,0,0)] - rc2*u2[OPS_ACC0(2,0,0)]);
}


void taylor_green_vortex_block0_97_kernel(const double *u1 , double *wk5)
{
wk5[OPS_ACC1(0,0,0)] = rinv1*((rc2)*u1[OPS_ACC0(0,0,-2)] - rc3*u1[OPS_ACC0(0,0,-1)] + (rc3)*u1[OPS_ACC0(0,0,1)] - rc2*u1[OPS_ACC0(0,0,2)]);
}


void taylor_green_vortex_block0_98_kernel(const double *wk3 , const double *u0 , const double *wk1 , const double *wk2 , const double *wk5 , const double *rho , const double *wk4 , const double *u1 , const double *u2 , const double *wk0 , double *ke , double *enstrophy , double *rhomean)
{
*ke = *ke + (rc0)*(pow(u0[OPS_ACC1(0,0,0)], 2) + pow(u1[OPS_ACC7(0,0,0)], 2) + pow(u2[OPS_ACC8(0,0,0)], 2))*rho[OPS_ACC5(0,0,0)];
*enstrophy = *enstrophy + (rc0)*(pow(wk0[OPS_ACC9(0,0,0)] - wk1[OPS_ACC2(0,0,0)], 2) + pow(wk2[OPS_ACC3(0,0,0)] - wk5[OPS_ACC4(0,0,0)], 2) + pow(wk3[OPS_ACC0(0,0,0)] - wk4[OPS_ACC6(0,0,0)], 2))*rho[OPS_ACC5(0,0,0)];
*rhomean = *rhomean + rho[OPS_ACC5(0,0,0)];
}


#endif