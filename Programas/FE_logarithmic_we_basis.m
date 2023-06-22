clear

%%% Size of the domain \Omega=(-L,L)
L=1;

%%% Number of discrete points, mesh and mesh size
Nval=[50,100,200,400,800,1600];
%Nval=[100];
step=[];
dif_norm=[];
dif_L2=[];
dif_linfinity=[];
pendiente=[]; alpha=[];
slope_L2=[];

for N=Nval

    xi=linspace(-L,L,N+2);
    h=xi(2)-xi(1);

    Alog=LoglapRigidity(L,N);
    mass = MassMatrix(xi,h);

    %%% Right-hand side of the problem and projection over finite elements
    %f = @(x) 1+0*x; %%% Torsion
    %f = @(x) (-3*x.^2+1)+(log(1./(1-x.^2))+(2*log(2)+psi(1/2)+psi(1))).*(1-x.^2); %%% (1-x^2)
    f = @(x) log(1./(1^2-x.^2))+(2*log(2)+psi(1/2)+psi(1)); %%% characteristic
    %f = @(x)  (2+log(1./(1-x.^2))+(2*log(2)+psi(1/2)+psi(1))).*x; %%% x

    F=projection(xi,f,h);

    exsol=@(x) 1+0*x+0*(1-x.^2);

    sol_log=Alog\F;

    err_temp=exsol(xi)'-sol_log;

    error_norm=sqrt(E_norm(xi(1:end),err_temp,1));
    err_linf=norm(err_temp,Inf);
    %err_L2=norm(err_temp,2);
    err_L2=sqrt(err_temp'*mass*err_temp);

    dif_norm=[dif_norm;error_norm];
    dif_linfinity=[dif_linfinity;err_linf];
    dif_L2=[dif_L2;err_L2];
    
    step=[step;h];

    if length(dif_norm)== 1
        pendiente=[pendiente,NaN];
        alpha=[alpha,NaN];
        slope_L2=[slope_L2,NaN];
    else
        pendiente=[pendiente,log(dif_norm(end)/dif_norm(end-1))/log(step(end)/step(end-1))];
        alpha=[alpha,abs(log(dif_norm(end)/dif_norm(end-1)))/abs(log(abs(log(step(end))/abs(log(step(end-1))))))];
        slope_L2=[slope_L2,log(dif_L2(end)/dif_L2(end-1))/log(step(end)/step(end-1))];
    end

    figure(1)
    plot(xi, sol_log,xi,exsol(xi),'x',xi,sol_log*0,'LineWidth',1);
    %plot(xi, sol_log,'LineWidth',1);
    
end

figure(2)
loglog(step,dif_linfinity,'LineWidth',2.5);
sl_dif=log(dif_linfinity(1)/dif_linfinity(end))/log(step(1)/step(end));
legend("Slope: "+num2str(sl_dif)); title("Error L infinity norm")

figure(3)
loglog(step,dif_norm,'LineWidth',2.5);
sl_dif=log(dif_norm(1)/dif_norm(end))/log(step(1)/step(end));
legend("Slope: "+num2str(sl_dif)); title("Error norm quadrature")

figure(4)
loglog(step,dif_L2,'LineWidth',2.5);
sl_dif=log(dif_L2(1)/dif_L2(end))/log(step(1)/step(end));
legend("Slope: "+num2str(sl_dif)); title("Error L^2")


%%% Descomentar si queremos guardar la solución
write_sol(xi,sol_log,false);

write_numsol_realsol(xi,sol_log,exsol(xi),exsol,false);

[xx,B0]=plt_sol_nearboundary(xi,f(xi),h);


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
        F(i) = ((h)/N)*sum(B1);
    end

    if i==N
        xx = linspace(xi(i)-h,xi(i),N+1);
        xx = 0.5*(xx(2:end)+xx(1:end-1));
        B1 = 2*f(xx).*Phi((xx-xi(i))/h);
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
A(1,3)=h*(-2/3);
A(N+2,N)=A(1,3);

%%% Entry a(\phi_0,\phi_j) with j\in\inter{3,N}
for j=4:N+1
    A(1,j)=h*(1/3*(-3*j^3*log(j)+6*j^2*log(j)+(j-2)^3*(-log(j-2))...
            +3*(j-1)^2*(j-2)*log(j-1)+(j+1)^2*(j-2)*log(j+1)-2));
    A(N+2,N+3-j)= A(1,j);
end

%%% Entry a(\phi_0,\phi_j) with j=N+1
A(1,N+2)= h*(2/3*(2*(N-3)*N^2*log(N)-N-(N-1)^3*log(N-1)-(N+1)*((N-4)*N+1)*log(N+1)-3));


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
A(1,2)= h*(-emc/3+5/9-2/3*log(h)-2/3*log(2)+1/3*psi(1/2));
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

A(1,1)=h*(-(4*emc)/3+32/9-8/3*log(h)+8/3*log(2)+4/3*psi(1/2));
A(N+2,N+2)=A(1,1);

for i=2:N+1
    A(i,i) = -2/9*h*(6*log(h) + 3*emc - 11 + log(64) - 3*psi(1/2));
end

end


function [out] = E_norm(xi,u,r)
eps=10e-3;
%eps = 1e-3*(xi(2)-xi(1));
yi = xi;
[Ux, Uy] = meshgrid(u,u);
[X,Y] = meshgrid(xi,yi);

W = 1/2*(Ux-Uy).^2./abs(X-Y);
W(~isfinite(W)) = 0;
Z = (abs(X-Y)<=1).*W;
Z_new = zeros(length(Z),1);
for i = 1:length(Z)
    Z_new(i) = trapz(yi,Z(:,i));
end
out_1 = trapz(xi,Z_new);
killing_measure = ((r + eps - abs(xi))<=1).*(-log(r + eps -abs(xi)));
out_2 = trapz(xi,u.^2.*killing_measure');
out = out_1 + 0*out_2;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Funciones para guardar soluciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_sol(xi,solution,flag)
if flag==true
    [dt,ht,mt]=obtain_date();

    L=max(abs(xi));
    h=xi(2)-xi(1);
    N=2*L/h-1;

    name_result = strcat('num_results/','sol-log_',...
        dt,'_',ht,'h',mt,'.org');
    file_result_temp=name_result;
    outs_file_temp=fopen(file_result_temp,'w');

    fprintf(outs_file_temp,'%s %s \n','#domain: ', ...
        strcat('(-',num2str(L),',',num2str(L),')'));
    fprintf(outs_file_temp,'%s %s \n','#N: ', num2str(N));
    fprintf(outs_file_temp,'%s %4.4e \n','#h: ', 2*L/(N+1));

    fprintf(outs_file_temp,'%4.4f %4.4e \n',[xi',[2*solution(1);solution(2:end-1);2*solution(end)]].');
    ST=fclose(outs_file_temp);
end
end

function write_numsol_realsol(xi,numsol,realsol,exsol_expr,flag)
if flag==true
    [dt,ht,mt]=obtain_date();

    L=max(abs(xi));
    h=xi(2)-xi(1);
    N=2*L/h-1;

    c=func2str(exsol_expr);

    name_result = strcat('num_results/','sol-log_num-vs-exact_',...
        dt,'_',ht,'h',mt,'.org');
    file_result_temp=name_result;
    outs_file_temp=fopen(file_result_temp,'w');

    fprintf(outs_file_temp,'%s %s \n','#domain: ', ...
        strcat('(-',num2str(L),',',num2str(L),')'));
    fprintf(outs_file_temp,'%s %s \n','#Exact solution: ', c);
    fprintf(outs_file_temp,'%s %s \n','#N: ', num2str(N));
    fprintf(outs_file_temp,'%s %4.4e \n','#h: ', 2*L/(N+1));

    fprintf(outs_file_temp,'%s %s %s \n','xi','numsol','exsol');
    fprintf(outs_file_temp,'%4.4f %4.4e %4.4e \n',[xi',[2*numsol(1);numsol(2:end-1);2*numsol(end)],realsol'].');
    ST=fclose(outs_file_temp);
end
end

function [dt,ht,mt]=obtain_date()
dt = datestr(now,'dd-mm-yyyy');
ht = datestr(now,'HH');
mt = datestr(now,'MM');
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

