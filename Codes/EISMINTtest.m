
function EISMINTtest

% test on modelling EISMINT
% Kori v0.9

%----------------------------------------
clear;
close all;

% Fixed margin experiment

ctr.imax=31; % number of grid points in y-direction
ctr.jmax=31; % number of grid points in x-direction
ctr.delta=50.e3; % grid size (m)
ctr.nsteps=2001; % number of iterations (total time is nsteps*dt)
ctr.dt=25; % time step (years)
ctr.snapshot=10;
ctr.plotH=1;

Li=(ctr.imax-1)*ctr.delta/1e3;
Lj=(ctr.jmax-1)*ctr.delta/1e3;
[X,Y] = meshgrid(0:ctr.delta/1e3:Li,0:ctr.delta/1e3:Lj);
dist=max(abs(X-Lj/2.),abs(Y-Li/2.));
Mb=zeros(ctr.imax,ctr.jmax)+0.3;
Ts=239.+8e-8*dist.^3.-273.15;
save('EISMINT','Mb','Ts');

KoriModel('EISMINT','EISMINT01',ctr);


% Moving margin experiment

ctr.MbType=4;
ctr.TsType=2;

KoriModel('EISMINT','EISMINT02',ctr);

end

