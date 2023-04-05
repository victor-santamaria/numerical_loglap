function  write_sol(argin)
%WRITE_SOL Guardar archivo.org con la solución
%   Esta función guarda en un archivo .org los puntos del mallado x_i y
%   y la solución y(x_i) de L_\Delta=f

if argin == true

    dt = datestr(now,'dd-mm-yyyy');
    ht = datestr(now,'HH');
    mt = datestr(now,'MM');

    disp('guardando solución...')

    name_result = strcat('num_results/','sol-log_',...
        dt,'_',ht,'h',mt,'.org');
    file_result_temp=name_result;
    outs_file_temp=fopen(file_result_temp,'w');

    %     fprintf(outs_file_temp,strcat('#',name_result,'\t - \t',donnees.nom,'\n'));
    %     fprintf(outs_file_temp,strcat('#','name: solution ode','\n'));
    %
    %     fprintf(outs_file_temp,'%4.4f %4.4e %4.4e \n',...
    %             [maillage.xi,soly_libre(1:maillage.N,end),soly_libre(maillage.N+1:end,end)].');
    
    ST=fclose(outs_file_temp);

end

end

