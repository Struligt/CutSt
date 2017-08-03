%module loadparams just loads the params for function to be used. 
% for example, can clear workspace, then load params , then run the desired
% matlab function eg loadparams; dancut(ws, W,bs);

% calibration against Bertsimas MIT unbounded knapsack example
% ws = [6,4,3,1]'; %length of each piece [2,4,5]
% W = 25; %stock length (eg 4,20 m lumber) 12
% p = [11, 7, 5, 1]';
% bs = [100,100,100,100]'; %unbounded

% building climbing wall lengths
% bs = 4*ones(1,5); %quantity demanded for each length
% ws = [2100,2174,903,1806,557]; %length of each piece [2,4,5]
% W = 4200; %stock length (eg 4,20 m lumber) 12


% building climbing wall lengths, used 19/05-17
%delayed column generation produces a singular matrix for some reason
% ie knapsack solver finds "best" pattern that 
% ws = [2190,70,1829,2000]; %length of each piece [2,4,5]
% W = 4200; %stock length (eg 4,20 m lumber) 12
% bs = [6,6,6,2]; %demand for each cut'
% Cs = W.*[14.33]/1000;
% 
%used 19-24/05-17, for testing dancut2.m 
% ws = [2190,70,1829,2000]; %length of each piece [2,4,5]
% bs = [6,6,6,2]; %demand for each cut'
% W = 4800; %stock length (eg 4,20 m lumber) 12
% Cs = W.*[15.74]/1000;
 
% % building climbing wall lengths, used 19/05-17, for testing dancut2.m
ws = [2190,70,1829,2000]; %length of each piece [2,4,5]
ws = [70,1829,2000,2190];
bs = [6,6,2,6]; %demand for each cut'
% this ratio of values gives "mixing" of stock lengths chosen in final soln :
W = [4200, 4800]; %stock length (eg 4,20 m lumber) 12
Cs = W.*[14.33 15.74*0.7966169]/1000;
% actual available lengths & prices :
% C14 och G04 virke: f0rsta ar byggregel, andra ar formvirke
W = [2500, 3600, 4200, 4800]; 
Cs = W.*[24.20 15.26 15.27 16.79]/1000;

% C24 hog hallfasthet virke :
W = [3000 3600 4200 4800 5400];
Cs = W.*[20.98 19.08 19.08 20.98 20.97]/1000;


% ws = [2190,70,1829,2000]; %length of each piece [2,4,5]
% W = [3000, 3600, 4200, 4800]; %stock length (eg 4,20 m lumber) 12
% Cs = W.*[15.75 14.33 14.33 15.74]/1000;
% bs = [6,6,6,2]; %demand for each cut'