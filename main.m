clc;
clear;
close all;

%% Laith Abualigah, Enhanced Aquila Optimizer for Global Optimization and Data Clustering, 
%% Sintafic Report, Springer, 2025
%% Problem Definition% D=normr(data); iris 3  vowel 10
dat = xlsread('iris.xlsx');
[k]=3;
D=dat(:,1:(size(dat,2)-1));
%% 
N=5;    
M_Iter=100;
LB=1;  UB=k;
Dim=length(D);


[Best_FF,Best_P,con7]=LOBLAO(N,M_Iter,LB,UB,Dim,k,D);

 
 T=1:M_Iter;
semilogy(T,con7,'LineWidth', 2)
title('Iris Dataset')
xlabel('No. of Iterations');
ylabel('Best fitness value');
legend('LOBLAO')