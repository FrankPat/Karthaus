
function GreenlandRun

% test on modelling Greenland ice sheet 10km
% f.ETISh v.1.8

%----------------------------------------
clear;
close all;

% maxNumCompThreads(2);

% Initialization SIA

ctr.inverse=1;
ctr.imax=281;
ctr.jmax=151;
ctr.delta=10.e3;
ctr.nsteps=2001;
ctr.dt=20;

ctr.Tcalc=2;
ctr.Tinit=1; % SET 1 FOR INITIALIZATION !
ctr.Asin=zeros(ctr.imax,ctr.jmax)+3e-9; %3e-9
ctr.Ao=5.0e-17; % 5.0e-17
ctr.m=3;
ctr.calving=4;

KoriModel('Green10km','GRinit10a',ctr);


% forcing

ctr.inverse=0;
ctr.Tinit=0;
ctr.MbType=2;
ctr.TsType=1;
ctr.PDDcalc=1;
ctr.nsteps=501;
ctr.plotH=1;
fc.DeltaT=zeros(ctr.nsteps,1)+5;

KoriModel('GRinit10a','GRrun10a',ctr,fc);



end

