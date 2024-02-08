%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA108
% Project Title: Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

%% Problem Settings

CostFunction=@Ackley;   % Cost Function

nVar=10;                % Number of Unknown (Decision) Variables

VarMin=-10;             % Lower Bound of Decision Variables
VarMax= 10;             % Upper Bound of Decision Variables

%% CMA-ES Settings

% Maximum Number of Iterations
MaxIt=300;

BestCost = solve(CostFunction, nVar, MaxIt, VarMax, VarMin);

%% Display Results

figure;
% plot(BestCost, 'LineWidth', 2);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

