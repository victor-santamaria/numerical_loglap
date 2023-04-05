
%WRITE_SOL Guardar archivo.org con la solución
%   Esta función guarda en un archivo .org los puntos del mallado x_i y
%   y la solución y(x_i) de L_\Delta y = f

    dt = datestr(now,'dd-mm-yyyy');
    ht = datestr(now,'HH');
    mt = datestr(now,'MM');

    disp('guardando solución...')

    name_result = strcat('num_results/','sol-log_',...
        dt,'_',ht,'h',mt,'.org');
    file_result_temp=name_result;
    outs_file_temp=fopen(file_result_temp,'w');

    fprintf(outs_file_temp,strcat('#','name: solution logarithmic Poisson problem','\n'));
    fprintf(outs_file_temp,'%s %s \n','#domain: ', ...
        strcat('(-',num2str(L),',',num2str(L),')'));
    fprintf(outs_file_temp,'%s %s \n','#N: ', num2str(N));
    fprintf(outs_file_temp,'%s %4.4e \n','#h: ', 2*L/(N+1));
    
    fprintf(outs_file_temp,'%4.4f %4.4e \n',[xi',sol_log].');
    ST=fclose(outs_file_temp);



