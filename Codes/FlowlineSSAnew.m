% SSA ice-sheet model
% Staggered grid for u
% Upstream differences, Crank-Nicholson scheme and Runge-Kutta RK-2
%
% Frank PATTYN, Laboratoire de Glaciologie, ULB, 2018
% Version 2.0

function [gpos,h]=FlowlineSSAnew(maxx,dx,dt,A,h,time_end)

upstream=1;
omega=2.5; % Crank-Nicholson scale factor (0=explicit; 1=implicit; >1 over-implicit)
RKorder=2; % 1: Euler; 2: Runge-Kutta-2

% Vectors
L=(maxx-2)*dx;
x=linspace(-dx,L,maxx)'; % horizontal distance
gridx=dx*1e3;
b=-x;
b(1)=b(3); % symmetric ice divide (extra grid point)
u=zeros(maxx,1)+0.1;

% Constants
secperyear=365*24*3600;
A=A*secperyear; % Ice flow parameter
rho_ice=910.; % Ice density
rho_sea=1000.; % Sea water density
grav=9.81; % Gravitation constant
n=3; % Glen index
mb=0.3; % Surface mass balance
sea_level=0; % height of eustatic sea level
m=1./3.; % basal sliding exponent
c=5e5; % basal friction coefficient 1e6

time_lapse = round(time_end/dt)+1;
bstag=(b+circshift(b,[-1 0]))/2;

for time_count=1:time_lapse
    for kk=1:RKorder % RK-2 iteration on half time step
        if kk==1
            h0=h;
        end
        RKdt=dt;
        RKdt(RKorder==2 & kk==1)=dt/2;
        % floating condition for ice sheet geometry determination
        haf=b-sea_level+h*rho_ice/rho_sea; % height above floating
        hb=b;
        hb(haf<0)=sea_level-rho_ice*h(haf<0)/rho_sea;
        s=hb+h;

        % floating condition on u-grid
        hstag=(h+circshift(h,[-1 0]))/2;
        haf=bstag-sea_level+hstag*rho_ice/rho_sea; % height above floating
        grl=ones(maxx,1)*2; % initialize grl-vector
        grl(haf>0)=0; % grounded
        grl(haf<0)=2; % floating

        % Longitudinal and driving stresses
        slope=diff(s)/gridx; % slope on u-grid
        slope(maxx)=0;
        taud=-rho_ice*grav.*hstag.*slope; % driving stress on u-grid
        txx=.25*rho_ice*grav*h*(1.-rho_ice/rho_sea); % on h-grid
        u=shelfystream(u,h,grl,taud,txx,A,c,gridx,secperyear,maxx,m,n);

        % Arrange staggered grid 
        dtdx2=RKdt/(2.*gridx);
        umin1=circshift(u,[1 0]);
        umin2=circshift(u,[2 0]);
        hmin1=circshift(h,[1 0]);
        hmin2=circshift(h,[2 0]);
        hplus1=circshift(h,[-1 0]);

        if upstream==1
            B=zeros(maxx,4);
            B(:,1)=dtdx2*umin2;
            B(:,2)=dtdx2*(-3*umin1-umin2);
            B(:,3)=dtdx2*(2*u+umin1);

            B(2,1)=0; % i=2 (ice divide)
            B(2,2)=-dtdx2*umin1(2); 
            B(2,3)=dtdx2*(u(2)-umin1(2));
            B(2,4)=dtdx2*u(2);

            B(1,:)=0; % i=1
            B(1,3)=1;
            B(maxx,:)=0; % i=maxx (shelf front)
            B(maxx,2)=-1;
            B(maxx,3)=0;
        else
            B=zeros(maxx,4);
            B(:,1)=0;
            B(:,2)=-dtdx2*umin1; 
            B(:,3)=dtdx2*(u-umin1);
            B(:,4)=dtdx2*u;

            B(1,:)=0; % i=1
            B(1,3)=1;
            B(maxx,:)=0; % i=maxx (shelf front)
            B(maxx,2)=-1;
            B(maxx,3)=0;
        end
        
        C=zeros(maxx,4);
        C(1:maxx-2,1)=B(3:maxx,1)*omega;
        C(1:maxx-1,2)=B(2:maxx,2)*omega;
        C(1:maxx,3)=B(:,3)*omega+1;
        C(2:maxx,4)=B(1:maxx-1,4)*omega;
        
        D=mb*RKdt+h0-h*(1-omega).*B(:,3)-hmin2*(1-omega).*B(:,1)- ...
            hmin1*(1-omega).*B(:,2)-hplus1*(1-omega).*B(:,4);
        D(1)=h(3);
        D(maxx)=(1-omega)*h(maxx-1);

        Am=spdiags(C,[-2 -1 0 1],maxx,maxx);
        h=Am\D;
        h(1)=h(3);
    end
    
    % plotting
    
    if rem(time_count,10)==1
        plot(x(2:end),h(2:end)+hb(2:end)); hold on;
        plot(x(2:end),b(2:end),'linewidth',2);
        plot(x(2:end),hb(2:end));
        hold off;
        grid on;
        pause(0.00001);
    end
    
end %% end of time stepping

% Determination of grounding line position on u-grid
grlj=maxx;
for j=1:maxx-2
    if grl(j)<.5 && grl(j+1)>1.5
        grl(j)=1; % last grounded grid point grlj
        grlj=j;
    end
end
        
xlabel('Horizontal distance (km)');
ylabel('Elevation (m a.s.l.)');
disp(x(grlj));
gpos=x(grlj);

save toto;

end

%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%


function [u] = shelfystream(u,h,grl,taud,txx,A,c,gridx,secperyear,maxx,m,n)

eeff=zeros(maxx,1);
f=zeros(maxx,1);
g=zeros(maxx,1);

for i=1:50
    eeff(2:maxx,1)=(diff(u)/gridx).^2; % eeff on h-grid
    eeff(1)=eeff(3);
    eeff(maxx)=eeff(maxx-1);
    eeff(eeff>0)=eeff(eeff>0).^((1-n)/(2*n));
    eeff(eeff<1)=1;
    mu=0.5*h*A^(-1/n).*eeff; % mu on h-grid
    beta2=c*abs(u).^(m-1)/(secperyear^m);
    beta2(grl==2)=0;
    beta2(1)=beta2(2);
    mu1=circshift(mu,[-1 0]);
    dn=-2*(mu+mu1)/(gridx^2)+2*(mu1-mu)/(gridx^2);
    cen=-4*(mu+mu1)/(gridx^2)-beta2;
    up=-2*(mu+mu1)/(gridx^2)-2*(mu1-mu)/(gridx^2);
    con=-taud;
    
    % Tridiagonal matrix solution    
    f(maxx)=1.;
    g(maxx)=(A*txx(maxx)^n)*gridx; % boundary condition on u-grid
    for j=maxx-1:-1:2
        f(j)=dn(j)/(cen(j)-up(j)*f(j+1));
        g(j)=(con(j)+up(j)*g(j+1))/(cen(j)-up(j)*f(j+1));
    end
    u0=u;
    u(1)=-g(2)/(1.+f(2));
    for j=2:maxx-1
        u(j)=g(j)+f(j)*u(j-1);
    end
    duxs=abs(u-u0);
    if mean(duxs)<1 % Limit on convergence
    	break;
    end
end

end

