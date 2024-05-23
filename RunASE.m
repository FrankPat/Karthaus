function RunASE

% test on modelling PIG and Thwaites at 3km
% assembling input files
% Kori-ULB

clear;
close all;

%% First SIA inversion (As)
ctr.imax=318;
ctr.jmax=270;
ctr.kmax=21;
ctr.delta=3.e3;
ctr.nsteps=8001;
ctr.dt=5;
ctr.Ao=5.0e-17; % 5.0e-17
ctr.m=3;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+3e-9;
ctr.calving=4;
ctr.Tinit=1;
ctr.Tcalc=2;
ctr.basin=1;
ctr.inverse=1;

% KoriModel('ASE3km','ASEinta1',ctr);


%% Second SSA inversion (As, MeltInv) - short run

ctr.inverse=2;
ctr.meltfunc=1;
ctr.GroundedMelt=1; % better for basins!
ctr.shelf=1;
ctr.SSA=2;
ctr.nsteps=101;
ctr.dt=0.01;
ctr.Tinit=0; % SET 1 FOR INITIALIZATION !
ctr.Tinv=10;

KoriModel('ASEinta1','ASEintb1',ctr);

% Second SSA inversion (As, MeltInv) - long run
ctr.nsteps=10001;
ctr.dt=0.2;

% KoriModel('ASEintb1','ASEintc1',ctr);

%% Control run with optimized melt (MeltInv)
ctr.inverse=0;
ctr.nsteps=1001;
ctr.meltfunc=11; % use MeltInv as melt
ctr.calving=2;

% KoriModel('ASEintc1','ASEruna1',ctr);

%% Forcing run
ctr.meltfunc=3;
ctr.gammaT=1e-4;

KoriModel('ASEruna1','ASErunb1',ctr); % start from control run
% KoriModel('ASEintc1','ASErunc1',ctr); % start from optimized run

% ctr.timeslice=1;
ctr.nsteps=251;
ctr.dt=0.1;
ctr.damage=1;
ctr.calving=5;
% KoriModel('ASEintc','ASErunf1',ctr);

% ctr.GeoidCalc=1;
% ctr.BedAdj=1;
% KoriModel('ASEintc','ASErune',ctr);

% ctr.damage=0;
% ctr.meltfunc=8;
% fc.butfac=zeros(ctr.nsteps,1);
% KoriModel('ASEintc','ASErung',ctr,fc);


end


