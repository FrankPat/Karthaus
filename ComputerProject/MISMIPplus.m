function MISMIPplus

% MISMIP+ experiment with Kori-ULB
% Details in Cornford et al (2020)

clear;
close all;

% Initial ice sheet creation

ctr.delta=2e3;
ctr.imax=23; % need number of cells + 2
ctr.jmax=352; % need number of cells + 1
ctr.m=3;
ctr.dt=1;
ctr.nsteps=8001;
ctr.mismip=1; % indicates that BCs need to be applied
ctr.SSA=1;
ctr.shelf=1;
ctr.Ao=2e-17; % 4.0e-17
ctr.shelftune=1;
ctr.Asin=zeros(ctr.imax,ctr.jmax)+1e-6; % 1e-5

Li=(ctr.imax-2)*ctr.delta;
Lj=(ctr.jmax-2)*ctr.delta;
[X,Y]=meshgrid(-ctr.delta:ctr.delta:Lj,-ctr.delta:ctr.delta:Li);
B=Bx(X)+By(Y,4.0e3,5.0e2,24.0e3);
B=max(B,-720.0);
% periodic BCs
B(1,:)=B(3,:);
B(ctr.imax,:)=B(ctr.imax-2,:);
H=zeros(ctr.imax,ctr.jmax)+10;
Mb=zeros(ctr.imax,ctr.jmax)+0.3;
Ts=zeros(ctr.imax,ctr.jmax)-5.0;
save('MISMIPplusin','B','H','Mb','Ts');

%KoriModel('MISMIPplusin','mismiptestin_2km',ctr);


% Refining grounding line position on upward bed slope: Ice0

ctr.Ao=4e-17;
ctr.nsteps=8001;
KoriModelAll('mismiptestin_2km','mismiptest1_2km',ctr);

% Ice1r experiment for 100 year with melting

% ctr.nsteps=101;
% ctr.meltfunc=10;
% ctr.meltfac=1;
% ctr.BetaIter=ctr.nsteps;
% KoriModel('mismiptest1_2km','mismiptest1r_2km',ctr);
% 
% % Ice1rr experiment for 200 year with melting
% 
% ctr.nsteps=201;
% KoriModel('mismiptest1_2km','mismiptest1rr_2km',ctr);
% 
% % Ice1ra experiment starting from Ice1r without melt
% 
% ctr.nsteps=101;
% ctr.meltfunc=0;
% KoriModel('mismiptest1r_2km','mismiptest1ra_2km',ctr);
% 
% % Plot results
% 
% load mismiptest1_2km_toto;
% Li=(ctr.imax-2)*ctr.delta/1e3;
% Lj=(ctr.jmax-2)*ctr.delta/1e3;
% [Xo,Yo] = meshgrid(-ctr.delta/1e3:ctr.delta/1e3:Lj,-ctr.delta/1e3:ctr.delta/1e3:Li);
% figure;
% contour(Xo,-Yo,MASK,[1],'k','linewidth',2);
% hold on;
% load mismiptest1r_2km;
% contour(Xo,-Yo,MASK,[1],'r','linewidth',2);
% load mismiptest1rr_2km;
% contour(Xo,-Yo,MASK,[1],'g','linewidth',2);
% load mismiptest1ra_2km;
% contour(Xo,-Yo,MASK,[1],'b','linewidth',2);
% xlim([300 550]);
% xlabel('x (km)');
% ylabel('y (km)');
% grid on;
% legend('Ice0','Ice1r','Ice1rr','Ice1ra');


end


function [Bxx]=Bx(x)
    B0=-150;
    B2=-728.8;
    B4=343.91;
    B6=-50.57;
    xx=x./300.0e3;
    xx2=xx.^2;
    xx4=xx2.^2;
    xx6=xx4.*xx2;
    Bxx=B0+B2*xx2+B4*xx4+B6*xx6;
end

function [Byy]=By(y,fc,dc,wc)
    Byy=dc./(1.0+exp(-2.0*(y-wc)/fc)) + dc./(1.0+exp(2.0*(y+wc)./fc));
end
