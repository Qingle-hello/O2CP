clc;clear;close all;

addpath('Solver');


fprintf('1\n\n');
Ex1_parareal_nonsmooth(@ SDIRK2, 20);
Ex1_parareal_nonsmooth(@ OCP, 20);
Ex1_mparareal_nonsmooth(@ BDF2_solver,@ LIIIC2_solver,20/2);
Ex1_mparareal_nonsmooth(@ MOCP1_solver,@ LIIIC2_solver,20/2);

fprintf('2\n\n');
Ex1_parareal_nonsmooth(@ SDIRK2, 50);
Ex1_parareal_nonsmooth(@ OCP, 50);
Ex1_mparareal_nonsmooth(@ BDF2_solver,@ LIIIC2_solver,50/2);
Ex1_mparareal_nonsmooth(@ MOCP1_solver,@ LIIIC2_solver,50/2);

fprintf('3\n\n');
Ex1_parareal_smooth(@ SDIRK2, 20);
Ex1_parareal_smooth(@ OCP, 20);
Ex1_mparareal_smooth(@ BDF2_solver,@ LIIIC2_solver,20/2);
Ex1_mparareal_smooth(@ MOCP1_solver,@ LIIIC2_solver,20/2);

fprintf('4\n\n');
Ex1_parareal_smooth(@ SDIRK2, 50);
Ex1_parareal_smooth(@ OCP, 50);
Ex1_mparareal_smooth(@ BDF2_solver,@ LIIIC2_solver,50/2);
Ex1_mparareal_smooth(@ MOCP1_solver,@ LIIIC2_solver,50/2);

%%
