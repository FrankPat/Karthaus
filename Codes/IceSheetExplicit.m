% Simple ice sheet model on staggered grid

close all;
clear;

m=51;
L=1.5e6; % length of domain (m)
h=zeros(m,1);
time_end=20000.;
dt=2; % time step (years)
anim=200;
gridx=L/(m-1); % grid spacing (m)
x=linspace(0,L,m)'; % horizontal distance

% Constants
A=1e-16; % Ice flow parameter
rho_ice=900.; % Ice density
grav=9.81; % Gravitation constant
n=3; % Glen index
a=0.3; % Surface mass balance

% time stepping
time_lapse = round(time_end/dt)+1;
for time_count=1:time_lapse
    time=(time_count-1)*dt;
    slope=diff(h)/gridx; % slope on staggered grid
    slope(m)=0;
    hstag=(circshift(h,[-1 0])+h)/2;
    dif=2/(n+2.)*A*(hstag.^2).*(slope.^2).*(rho_ice*grav*hstag).^n;
    h=h+(dif.*slope-circshift(dif,[1 0]).*circshift(slope,[1 0]))*dt/gridx+a*dt;
    h(1)=0;
    h(m)=0;
    
    % Plot results
    if rem(time,anim)==0
        plot(x,h);hold on;
        grid on;
        pause(0.1);
        axis([0 +Inf 0 4000]);
    end
end %% end of time stepping



%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%


