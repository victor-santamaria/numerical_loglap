clear

%%% Size of the domain \Omega=(-L,L)
L=1;

%%% Number of discrete points, mesh and mesh size
Nval=[50];

for N=Nval

    xi=linspace(-L,L,N+2);
    h=xi(2)-xi(1);

    Alog=LoglapRigidity(L,N); %Computing the stiffness matrix

    %%% Right-hand side of the problem and projection over finite elements
    f = @(x) 1+0*x; %%% Torsion
    F=projection(xi,f,h);

    %%% Computing the solution
    sol_log=Alog\F;

    plotsol(xi, sol_log); %Plot solution taking into account the FE basis
    
end

%%% Auxiliary functions

function [F] = projection(xi,f,h)

N=size(xi,2);
Phi = @(x) 1-abs(x);
F = zeros(N,1);

for i=1:N
    if i~=1 && i~=N
        xx = linspace(xi(i)-h,xi(i)+h,N+1);
        xx = 0.5*(xx(2:end)+xx(1:end-1));
        B1 = f(xx).*Phi((xx-xi(i))/h);
        F(i) = ((2*h)/N)*sum(B1);
    end

    if i==1
        xx = linspace(xi(i),xi(i)+h,N+1);
        xx = 0.5*(xx(2:end)+xx(1:end-1));
        B1 = f(xx).*Phi((xx-xi(i))/h);
        F(i) = ((h)/N)*sum(B1);
    end

    if i==N
        xx = linspace(xi(i)-h,xi(i),N+1);
        xx = 0.5*(xx(2:end)+xx(1:end-1));
        B1 = f(xx).*Phi((xx-xi(i))/h);
        F(i) = ((h)/N)*sum(B1);
    end

end

end

function A = LoglapRigidity(L,N)

x = linspace(-L,L,N+2);
h = x(2)-x(1);

A = zeros(N+2,N+2);
emc=-psi(1); %%Euler-Mascheroni constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filling the upper triangle of FEM matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Entry a(\phi_0,\phi_j) with j=2
A(1,3)=h*(-1/3);
A(N+2,N)=A(1,3);

%%% Entry a(\phi_0,\phi_j) with j\in\inter{3,N}
for j=4:N+1
    A(1,j)=h*(1/6*(-3*j^3*log(j)+6*j^2*log(j)+(j-2)^3*(-log(j-2))...
            +3*(j-1)^2*(j-2)*log(j-1)+(j+1)^2*(j-2)*log(j+1)-2));
    A(N+2,N+3-j)= A(1,j);
end

%%% Entry a(\phi_0,\phi_j) with j=N+1
A(1,N+2)= h*(1/6*(2*(N-3)*N^2*log(N)-N-(N-1)^3*log(N-1)-(N+1)*((N-4)*N+1)*log(N+1)-3));


%%% Case: inner upper triangle 

for i=2:N-1
    for j=i+2:N+1
        k = j-i;

        if k~=2

            PP= 1/6*h*(-(-2 + k)^3*log(-2 + k) + 4*(-1 + k)^3*log(-1 + k) - ...
                6*k^3*log(k) + 4*(1 + k)^3*log(1 + k) - (2 + k)^3*log(2 + k));
            A(i,j) = PP;
        else
            PPbis = -(1/3)*h*(-24*log(3) - 12*log(4) - 16*log(9) + 27*log(16) + log(144));
            A(i,j) = PPbis;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filling the superdiagonal of FEM matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Entry a(\phi_0,\phi_j) with j=1
A(1,2)= h*(-emc/3+5/9-2/3*log(h)-2/3*log(2)+1/3*psi(1/2))/2;
A(N+2,N+1)=A(1,2);

%%% Inner superdiagonal
for i=2:N
    A(i,i+1) = -(emc*h)/6 + (11*h)/18 - 1/3*h*log(h) ...
        + 1/12*h*(64*log(2) - 54*log(3)) + 1/3*h*log(2) + 1/6*h*psi(1/2);
end

A = A+A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filling the diagonal of FEM matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Entry a(\phi_0,\phi_0) 

A(1,1)=h*(-(1*emc)/3+8/9-2/3*log(h)+2/3*log(2)+1/3*psi(1/2));
A(N+2,N+2)=A(1,1);

for i=2:N+1
    A(i,i) = -2/9*h*(6*log(h) + 3*emc - 11 + log(64) - 3*psi(1/2));
end

end

function M = MassMatrix(x,hx)

Nx = length(x);

M=zeros(Nx,Nx);

M(1,2)=1/6;
M(Nx,Nx-1)=M(1,2);

for i=2:Nx-2
    M(i,i+1)=1/6;
end

M=M+M';

M(1,1)=2/3;
M(Nx,Nx)=M(1,1);
 
for i=2:Nx-1
    M(i,i)=2/3;
end

M = sparse(hx*M);

end

function plotsol(xi,sol_log)

    sol_log(1)=sqrt(2)*sol_log(1);

    sol_log(end)=sqrt(2)*sol_log(end);

    plot(xi, sol_log,xi,sol_log*0,'LineWidth',1);
    
end

