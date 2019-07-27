clear;close all;clc;format short g;

alfaset=[0.1];
for i=1:length(alfaset)
    %  set of parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the paramaters of this model%%%%%%%%%%%%%%%%%
    n=500;   %the number of agent number
    N=50;    %the number of star number      
    L=N;    %the range of the star's position
    agentwidth=8;
    alfa=alfaset(i);
    timestep=10;           %the total timesteps
    agenttime=2;    %the number of agents' evoltion
    Thermo=L*0;
    
    autostar(n,N,L,agentwidth,alfa,timestep,agenttime,Thermo);
    
end
