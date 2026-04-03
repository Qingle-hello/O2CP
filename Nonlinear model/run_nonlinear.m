clc;clear;close all;

addpath('Solver');



fprintf('\n 1\n\n');
cL = 1; J1 = 20;
Ex2_parareal(@ SDIRK2, J1,cL);
Ex2_parareal(@ OCP, J1,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J1/2,cL);

cL = 1; J2 = 50;
Ex2_parareal(@ SDIRK2, J2,cL);
Ex2_parareal(@ OCP, J2,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J2/2,cL);

%%

fprintf('\n 2\n\n');
cL = 5; J1 = 20;
Ex2_parareal(@ SDIRK2, J1,cL);
Ex2_parareal(@ OCP, J1,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J1/2,cL);

cL = 5; J2 = 50;
Ex2_parareal(@ SDIRK2, J2,cL);
Ex2_parareal(@ OCP, J2,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J2/2,cL);

%%

fprintf('\n 3\n\n');
cL = 10; J1 = 20;
Ex2_parareal(@ SDIRK2, J1,cL);
Ex2_parareal(@ OCP, J1,cL); % not stable
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J1/2,cL);

cL = 10; J2 = 50;
Ex2_parareal(@ SDIRK2, J2,cL);
Ex2_parareal(@ OCP, J2,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J2/2,cL);

%

fprintf('\n 4\n\n');
cL = 15; J1 = 20;
Ex2_parareal(@ SDIRK2, J1,cL);
Ex2_parareal(@ OCP, J1,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J1/2,cL);

cL = 15; J2 = 50;
Ex2_parareal(@ SDIRK2, J2,cL);
Ex2_parareal(@ OCP, J2,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J2/2,cL);

%

fprintf('\n 5\n\n');
cL = 20; J1 = 20;
Ex2_parareal(@ SDIRK2, J1,cL);
Ex2_parareal(@ OCP, J1,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J1/2,cL);

cL = 20; J2 = 50;
Ex2_parareal(@ SDIRK2, J2,cL);
Ex2_parareal(@ OCP, J2,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J2/2,cL);

%

fprintf('\n 6\n\n');
cL = 40; J1 = 20;
Ex2_parareal(@ SDIRK2, J1,cL);
Ex2_parareal(@ OCP, J1,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J1/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J1/2,cL);

cL = 40; J2 = 50;
Ex2_parareal(@ SDIRK2, J2,cL);
Ex2_parareal(@ OCP, J2,cL);
Ex2_mparareal(@ BDF2_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver,@ LIIIC2_solver,J2/2,cL);
Ex2_mparareal(@ MOCP1_solver_Extplt,@ LIIIC2_solver,J2/2,cL);
