clear

%%% Fractional power
s=0.1;

%%% Size of the domain \Omega=(-L,L)
L=1;

c_1 = 2*(s*2^(2*s-1)*gamma(0.5*(1+2*s)))/(sqrt(pi)*gamma(1-s));

%%% Number of discrete points, mesh and mesh size
Nval=[50,100,200,400,800];
dif_norm=[];
dif_L2_norm_matlab=[];
dif_L2_norm=[];
step=[];
pendiente=[];

for N=Nval

    xi=linspace(-L,L,N+2);
    h=xi(2)-xi(1);

    As=flRigidity_w0(s,L,N);
    mass = MassMatrix(xi,h);

    %%% Right-hand side of the problem and projection over finite elements
    f = @(x) 1+0*x;
    %f = @(x) (c_1/(2*s))*((1+x).^(-2*s)+(1-x).^(-2*s));
    %f = @(x) sin(pi*x);
    F=projection(xi,f,h);

    sol_frac=As\F;

    %%% Exact solution
    exsol=@(x,s,L) 2^(-2*s)*sqrt(pi)/(gamma((1+2*s)/2)*gamma(1+s)).*(L^2-x.^2).^s;
    %exsol=@(x,s,L) 1+0.*x;

    figure(1)
    plot(xi, sol_frac,xi,exsol(xi,s,L),'x',xi,sol_frac*0,'LineWidth',1);
    %plot(xi, sol_frac,'LineWidth',1);

    err_temp=exsol(xi,s,L)'-sol_frac;
    
    e_norm = error_fl(h,s,sol_frac);

    err_L2_matlab=norm(err_temp,2);
    err_L2=sqrt(err_temp'*mass*err_temp);

    dif_norm=[dif_norm;e_norm];
    dif_L2_norm_matlab=[dif_L2_norm_matlab;err_L2_matlab];
    dif_L2_norm=[dif_L2_norm;err_L2];

    step=[step;h]


    if length(dif_norm)== 1
        pendiente=[pendiente,NaN];
    else
        pendiente=[pendiente,log(dif_norm(end)/dif_norm(end-1))/log(step(end)/step(end-1))];
    end

end

figure(2)
loglog(step,dif_norm,'LineWidth',2.5)
sl_dif=log(dif_norm(1)/dif_norm(end))/log(step(1)/step(end));
legend("Slope: "+num2str(sl_dif)); title('error H^s norm')

figure(3)
loglog(step,dif_L2_norm_matlab,'LineWidth',2.5);
sl_quad=log(dif_L2_norm_matlab(1)/dif_L2_norm_matlab(end))/log(step(1)/step(end));
legend("Slope: "+num2str(sl_quad)); title('error L2 matlab')

figure(4)
loglog(step,dif_L2_norm,'LineWidth',2.5);
sl_quad=log(dif_L2_norm(1)/dif_L2_norm(end))/log(step(1)/step(end));
legend("Slope: "+num2str(sl_quad)); title('error L^2')

[xx,B0] = plt_sol_nearboundary(xi,sol_frac,h);



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
        B1 = 2*f(xx).*Phi((xx-xi(i))/h);
        %B1 = f(xx).*Phi((xx-xi(i))/h);
        F(i) = ((h)/N)*sum(B1);
    end

    if i==N
        xx = linspace(xi(i)-h,xi(i),N+1);
        xx = 0.5*(xx(2:end)+xx(1:end-1));
        B1 = 2*f(xx).*Phi((xx-xi(i))/h);
        %B1 = f(xx).*Phi((xx-xi(i))/h);
        F(i) = ((h)/N)*sum(B1);
    end

end

end

function A = flRigidity_w0(s,L,N)


x = linspace(-L,L,N+2);
x = x(2:end-1);
h = x(2)-x(1);

A = zeros(N+2,N+2);

c = (s*2^(2*s-1)*gamma(0.5*(1+2*s)))/(sqrt(pi)*gamma(1-s));

% % %%% For taking into account the basis functions at the end sides
%%% a(\phi_0,\phi_j) 
for j=3:N+2
    k = j;

    den = 1/(4*s*(1-2*s));
        p = 2-2*s;
        q = 3-2*s;

        B1=den*(2*k^(1-2*s)-((k+1)^p-(k-1)^p)/(1-s)...
            -(2*k^q-(k+1)^q-(k-1)^q)/((1-s)*(3-2*s)));

        B3=den*(-2*k^(1-2*s)+(2*k^p-2*(k-1)^p)/(1-s)...
            +(2*(k-1)^q-k^q-(k-2)^q)/((1-s)*(3-2*s)));

    if j==N+2
        A(1,j) = (2*h^(1 - 2*s)*N^(-2*s)*(1 + N)^(-2*s)*(2 + N)^(-2*s)*(-N^(2*s)*(1 + N)^(2*s)*(2 + N)*...
            (1 + (1 + N)^2 + 4*(1 + N)*(-1 + s) - 6*s +... 
            4*s^2) + (2 + N)^(2*s)*(-N^3*(1 + N)^(2*s) + 2*N^(2*s)*(1 + N)^2*(-2 + N + 2*s))))/((-1 + s)*s*...
            (-3 + 2*s)*(-1 + 2*s));
     else
        A(1,j) = -4*h^(1-2*s)*(B1+B3);
    end
    A(N+2,N+3-j)= A(1,j);
end

%%% Parece que aqui esta el problema
%%% a(\phi_{0},\phi_{1}) & a(\phi_{N},\phi_{N+1})
A(1,2)=-((2^(1 - 2*s)*(-2 + 4^s)*h^(1 - 2*s))/((-1 + s)*s*(-3 + 2*s)));
A(N+2,N+1)=A(1,2);

%%% The usual fractional laplacian matrix

for i=2:N-1
    for j=i+2:N+1
        k = j-i;
        den = 1/(4*s*(1-2*s));
        p = 2-2*s;
        q = 3-2*s;

        B1=den*(2*k^(1-2*s)-((k+1)^p-(k-1)^p)/(1-s)...
            -(2*k^q-(k+1)^q-(k-1)^q)/((1-s)*(3-2*s)));

        B2=den*(-2*k^(1-2*s)+(2*(k+1)^p-2*k^p)/(1-s)...
            +(2*(k+1)^q-k^q-(k+2)^q)/((1-s)*(3-2*s)));

        B3=den*(-2*k^(1-2*s)+(2*k^p-2*(k-1)^p)/(1-s)...
            +(2*(k-1)^q-k^q-(k-2)^q)/((1-s)*(3-2*s)));

        B4=den*(2*k^(1-2*s)-((k+1)^p-(k-1)^p)/(1-s)...
            -(2*k^q-(k+1)^q-(k-1)^q)/((1-s)*(3-2*s)));

        A(i,j) = -2*h^(1-2*s)*(B1+B2+B3+B4);
    end
end

for i=2:N
    Q2=h^(1-2*s)*(2^p+2*s-3)/(s*(2-2*s)*(3-2*s));
    Q3=h^(1-2*s)*(13-5*2^(3-2*s)+3^(3-2*s)+s*(2^(4-2*s)-14)+4*s^2)...
        /(2*s*(1-2*s)*(1-s)*(3-2*s));
    Q4=-(h^(1-2*s))/((1-s)*(3-2*s));
    A(i,i+1) = 2*Q2+Q3+Q4;
end

A = A+A';

for i=2:N+1
    R2 = h^(1-2*s)*(4*s-6+2^(3-2*s))/(s*(1-2*s)*(1-s)*(3-2*s));
    
    R41 = h^(1-2*s)/((1-s)*(3-2*s));
    R42 = h^(1-2*s)*((2*s^2-5*s+4-2^(2-2*s))/(s*(1-2*s)*(1-s)*(3-2*s)));

    A(i,i) = 2*R2+2*(R41+R42);
end
 

%%% a(\phi_0,\phi_0) & a(\phi_{N+1},\phi_{N+1})
A(1,1) = (8*h^(1 - 2*s))/(s*(3 + 4*(-2 + s)*s));
A(N+2,N+2)=A(1,1);

A = c*A;

end

function M = MassMatrix(x,hx)

Nx = length(x);

M=zeros(Nx,Nx);

M(1,2)=1/3;
M(Nx,Nx-1)=M(1,2);

for i=2:Nx-2
    M(i,i+1)=1/6;
end

M=M+M';

M(1,1)=4/3;
M(Nx,Nx)=M(1,1);
 
for i=2:Nx-1
    M(i,i)=2/3;
end

M = sparse(hx*M);

end

function [e_norm] = error_fl(h,s,sol)

val = pi/(2^(2*s)*gamma(s+0.5)*gamma(s+1.5));
valnum = h*sum(sol);
e_norm = sqrt(val-valnum);

end

function [xx,B0]=plt_sol_nearboundary(xi,ui,h)

N=size(xi,2);

for i=1:N

    if i==1
        xx = linspace(xi(i),xi(i)+h,N+1);
        B0=ui(i)*(2*(h-1)/h-2*xx/h)+ui(i+1)*(xx/h+1/h);
        
    end
end

end



