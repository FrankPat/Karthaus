% Tests on grounding-line resolution

clear all;
close all;

L=500;
time_end=5000;
dx=.05;
A=1e-24;
timestep=20.;
maxx=round(L/dx);
grlpos=zeros(3,1);
h=zeros(maxx,1)+10;
[grlpos(1),h]=FlowlineSSAnew(maxx,dx,timestep,A,h,time_end);
A=0.2e-24;
[grlpos(2),h]=FlowlineSSAnew(maxx,dx,timestep,A,h,time_end);
A=1e-24;
[grlpos(3),h]=FlowlineSSAnew(maxx,dx,timestep,A,h,time_end);

%%
max=20;
time_end=5000;
L=500;
dx=linspace(3.5,9.5,max)';
dx=exp(dx)/1000;
maxx=round(L./dx);
A=1e-24;
timestep=20.;

grlpos=zeros(length(maxx),1);
for i=1:length(maxx)
    h=zeros(maxx(i),1)+10;
    [grlpos(i),h]=FlowlineSSAnew(maxx(i),dx(i),timestep,A,h,time_end);
end

save toto;

load toto

figure;
semilogy(grlpos,dx,'o-');
ylabel('Spatial resolution (km)');
xlabel('Grounding line position (km)');
grid on;


