% Finite difference model of a 2D ice sheet

clear all;
close all;

imax=31;
jmax=31;
g=9.81;
rho=910.;
secperyear=31556926;
Ao=1.e-16/secperyear;
a=0.3/secperyear;
L=1500e3/2.;
Z=(5*a*L^4/(2*Ao*(rho*g)^3))^(1/8);

nsteps=500;
dt=50*secperyear*a/Z;
nodes=imax*jmax;
sn=zeros(imax,jmax);
s=zeros(nodes,1);
node=zeros(imax,jmax);
row=zeros(4325,1);
col=zeros(4325,1);
value=zeros(4325,1);
delta=2/(imax-1);
R=zeros(nodes,1);
d=zeros(imax,jmax);

counter=0;
for j=1:jmax
    for i=1:imax
        counter=counter+1;
        node(i,j)=counter;
    end
end

for n=1:nsteps
    for j=1:jmax-1
        for i=1:imax-1
            d(i,j)=(1/4*(sn(i,j)+sn(i+1,j)+sn(i+1,j+1)+sn(i,j+1)))^5 ...
                /(4*delta^2)*( ...
                (sn(i+1,j)-sn(i,j)+sn(i+1,j+1)-sn(i,j+1))^2 ...
                +(sn(i,j+1)-sn(i,j)+sn(i+1,j+1)-sn(i+1,j))^2);
        end
    end
    % construct the finite-difference stiffness matrix
    count=0;
    for j=2:jmax-1
        for i=2:imax-1
            count=count+1;
            row(count)=node(i,j);
            col(count)=node(i,j);
            value(count)=1/dt+1/delta^2*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1));
            
            count=count+1;
            row(count)=node(i,j);
            col(count)=node(i,j+1);
            value(count)=-1/(2*delta^2)*(d(i,j)+d(i-1,j));
            
            count=count+1;
            row(count)=node(i,j);
            col(count)=node(i,j-1);
            value(count)=-1/(2*delta^2)*(d(i,j-1)+d(i-1,j-1));
            
            count=count+1;
            row(count)=node(i,j);
            col(count)=node(i+1,j);
            value(count)=-1/(2*delta^2)*(d(i,j)+d(i,j-1));
            
            count=count+1;
            row(count)=node(i,j);
            col(count)=node(i-1,j);
            value(count)=-1/(2*delta^2)*(d(i-1,j)+d(i-1,j-1));
            
            R(node(i,j))=1+sn(i,j)/dt;
        end
    end
    j=1;
    for i=2:imax-1
        count=count+1;
        row(count)=node(i,j);
        col(count)=node(i,j);
        value(count)=1;
        R(node(i,j))=0;
    end
    j=jmax;
    for i=2:imax-1
        count=count+1;
        row(count)=node(i,j);
        col(count)=node(i,j);
        value(count)=1;
        R(node(i,j))=0;
    end
    i=1;
    for j=2:jmax-1
        count=count+1;
        row(count)=node(i,j);
        col(count)=node(i,j);
        value(count)=1;
        R(node(i,j))=0;
    end
    i=imax;
    for j=2:jmax-1
        count=count+1;
        row(count)=node(i,j);
        col(count)=node(i,j);
        value(count)=1;
        R(node(i,j))=0;
    end
    count=count+1;
    row(count)=node(1,1);
    col(count)=node(1,1);
    value(count)=1;
    R(node(1,1))=0;
    
    count=count+1;
    row(count)=node(1,jmax);
    col(count)=node(1,jmax);
    value(count)=1;
    R(node(1,jmax))=0;
    
    count=count+1;
    row(count)=node(imax,jmax);
    col(count)=node(imax,jmax);
    value(count)=1;
    R(node(imax,jmax))=0;
    
    count=count+1;
    row(count)=node(imax,1);
    col(count)=node(imax,1);
    value(count)=1;
    R(node(imax,1))=0;
    
    % construct sparse matrix
    A=sparse(row,col,value);
    % Cholesky factor and solve
    s=A\R;
    sn=s(node);
    if rem(n,10)==1
        mesh(Z*sn); 
        axis([0 +Inf 0 +Inf 0 3500]);
        pause(0.1);
    end
end

figure;
spy(A);


