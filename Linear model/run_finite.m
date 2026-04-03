clc;clear;close all;

addpath('Solver');


fprintf('1\n\n');
Ex1_parareal_finite(@ SDIRK2, 20);
Ex1_parareal_finite(@ OCP, 20);
Ex1_mparareal_finite(@ BDF2_solver,@ LIIIC2_solver,20/2);
Ex1_mparareal_finite(@ MOCP1_solver,@ LIIIC2_solver,20/2);

fprintf('2\n\n');
Ex1_parareal_finite(@ SDIRK2, 50);
Ex1_parareal_finite(@ OCP, 50);
Ex1_mparareal_finite(@ BDF2_solver,@ LIIIC2_solver,50/2);
Ex1_mparareal_finite(@ MOCP1_solver,@ LIIIC2_solver,50/2);
