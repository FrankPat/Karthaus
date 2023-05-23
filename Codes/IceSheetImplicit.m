% Simple ice sheet model on staggered grid

close all;
clear;

m=51;
L=1.5e6; % length of domain (m)
h=zeros(m,1);
time_end=20000.;
dt=10.; % time step (years)
anim=200;
gridx=L/(m-2); % grid spacing (m)

% Vectors
x=linspace(0,L,m)'; % horizontal distance
b=x*0; % bedrock elevation

up=zeros(m,1);
dn=zeros(m,1);
cen=zeros(m,1);
dif=zeros(m,1);
con=zeros(m,1);
f=zeros(m,1);
g=zeros(m,1);
hstag=zeros(m,1);

% Constants
A=1e-16; % Ice flow parameter
rho_ice=910.; % Ice density
grav=9.81; % Gravitation constant
n=3; % Glen index
a=0.3; % Surface mass balance

% time stepping
time_lapse = round(time_end/dt)+1;
tc=0;

for time_count=1:time_lapse
    time=(time_count-1)*dt;
	s=b+h;
    slope=diff(s)/gridx; % slope on staggered grid
    for j=1:m-1
        hstag(j,1)=(h(j+1)+h(j))/2.; % ice thickness mapped on staggered grid
    end
    for j=1:m-1
        dif(j)=2./(n+2.)*A*hstag(j)^2*slope(j)^2*(rho_ice*grav*hstag(j))^n;
    end
    
    % Arrange staggered grid 
    dtdx=dt/(gridx*gridx);
    for j=2:m-1
        up(j)=dif(j)*dtdx;
        dn(j)=dif(j-1)*dtdx;
        cen(j)=1.+up(j)+dn(j);
        con(j)=h(j)+a*dt+dn(j)*b(j-1)-(cen(j)-1)*b(j)+up(j)*b(j+1);
    end

    % Tridiagonal matrix solution    
    f(1)=0;
    g(1)=0;
    f(n)=0;
    g(n)=0;
    for j=2:m-1
        f(j)=up(j)/(cen(j)-dn(j)*f(j-1));
        g(j)=(con(j)+dn(j)*g(j-1))/(cen(j)-dn(j)*f(j-1));
    end
    
    h(m)=0; % ice sheet BC
    for j=m-1:-1:1
        h(j)=g(j)+f(j)*h(j+1);
    end
    
    % Plot results
    if rem(time,anim)==0
        plot(x,h+b);hold on;
        plot(x,b);
        grid on;
        axis([0 +Inf 0 4000]);
        pause(0.1);
    end
end %% end of time stepping



%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%


