diel_path = "/gscratch/chem/clairew4/dielectrics/";

l_xmin = 1;
l_xmax = 5;
l_ymin = -26;
l_ymax = -1;
l_zmin = -4;
l_zmax = 4;

s_xmin = 1;
s_xmax = 5;
s_ymin = 1;
s_ymax = 16;
s_zmin = -4;
s_zmax = 4;

A = load("temp.out");
A = sortrows(A,1);
A = sortrows(A,2);
A = sortrows(A,3);

x = A(:,1);
y = A(:,2);  
z = A(:,3); 
line = [x y z round(A(:,4)) zeros([length(x),1])]; % x y z T_rounded IA_raw
T_r = line(:,4);
IA_raw = line(:,5);

for i = 1:length(T_r)
    if (x(i) <= l_xmax) && (x(i) >= l_xmin) && ((y(i) <= l_ymax) && (y(i) >= l_ymin)) && ((z(i) <= l_zmax) && (z(i) >= l_zmin))
       IA_raw(i) = 1e5; %long rod 
       temp_longrod = line(i,4);
    elseif (x(i) <= s_xmax) && (x(i) >= s_xmin) && ((y(i) <= s_ymax) && (y(i) >= s_ymin)) && ((z(i) <= s_zmax) && (z(i) >= s_zmin)) 
       IA_raw(i) = 1e5+1; %short rod 
       temp_shortrod = line(i,4);

    else
        IA_raw(i) = T_r(i);
    end
end

IA = IA_raw - min(IA_raw) + 1;
IA_un = unique(IA);

for i = 1:length(T_r)
    if isequal(IA(i,1), 1e5 - min(IA_raw) + 1 )
    IA(i,1) = IA_un(length(IA_un)-2) + 1;
    end
    if isequal(IA(i,1), 1e5+1 - min(IA_raw) + 1 )
    IA(i,1) = IA_un(length(IA_un)-2) + 2;
    end
    
end
line = [x y z T_r IA]; % x y z T_rounded IA
% IA_new = IA;
% rows = find(line(:,1) ~= 2); 
% line(rows, :) = [];
% IA_new(rows,:) = [];
% xnew = line(:,1);
% ynew = line(:,2);
% znew = line(:,3);
% 
% IAnew = reshape(IA_new, [(max(line(:,2))-min(line(:,2))+1), (max(line(:,3))-min(line(:,3))+1)]);
% pcolor(IAnew');
% shading interp
% axis equal
% axis off
% colorbar


JA = linspace(1,length(T_r), length(T_r))';
shape = [JA x y z IA IA IA];

fileID = fopen('shape.dat_new', 'wt');
fprintf(fileID,'%s\n',' Sphere shape');
fprintf(fileID,'\t%d %s\n',length(T_r),' = number of dipoles in target');
fprintf(fileID, '%s\n', ' 1.000000 0.000000 0.000000 = A_1 vector');
fprintf(fileID, '%s\n', ' 0.000000 1.000000 0.000000 = A_2 vector');
fprintf(fileID, '%s\n', ' 1.000000 1.000000 1.000000 = (d_x,d_y,d_z)/d');
fprintf(fileID, '%s\n', ' 0.000000 0.000000 0.000000 = (x,y,z)/d');
fprintf(fileID, '%s\n', ' JA  IX  IY  IZ ICOMP(x,y,z)');
dlmwrite('shape.dat_new', shape, 'delimiter', '\t', '-append');

I = 1;
T_r_un = unique(T_r);
fid = fopen('filler_ddscat.par', 'wt');
for i = IA_un(1):IA_un(length(IA_un)-2)
    I = 0 + I;
    fprintf(fid, "%s\n", strcat("'",diel_path, 'nT_',num2str(T_r_un(i)),"K.tab' = file with refractive index ", num2str(I) ));
    I = I + 1;
end
fprintf(fid, "%s\n", strcat("'",diel_path, 'au_drude_',num2str(temp_longrod),"K.tab' = file with refractive index of long rod" ));
fprintf(fid, "%s\n", strcat("'",diel_path, 'au_drude_',num2str(temp_shortrod),"K.tab' = file with refractive index of short rod" ));



T = T_r;
for i = IA_un(1):IA_un(length(IA_un)-2)
    dT = T_r_un(i);
    n_T = 1.473 + 2.7e-4*dT;  
    w = linspace(0.5,4,101);
    n = linspace(n_T, n_T, 101);
    k = 0*w;
    Data = [1.24./w', n', k'];  
    Data = [sortrows(Data(:,1)), n', k'];
    fid = fopen(strcat(diel_path, 'nT_', num2str(dT), 'K.tab'), 'wt');
    fprintf(fid,'%s\n','nT of glycerol');
    fprintf(fid,'%s\n','1 2 3 0 0 = Specifies n or Eps');
    fprintf(fid,' %s\n','lambda  n    k');
    for j = 1:length(Data)
        fprintf(fid,' %f',Data(j,1));
        fprintf(fid,' %f',Data(j,2));
        fprintf(fid,' %f\n',Data(j,3));
    end
    fclose(fid);
end

dT = temp_longrod;
T = 298 + dT;
Eps_inf = -.0026*T+3.546;
wp = 0.0003*T + 9.0788;
gamma = (7e-5)*T+0.0337;
w = linspace(0.5,4,101);
w(1) =[];
Eps = Eps_inf - wp^2./(w.^2+1i*gamma*w);
n = sqrt((sqrt(real(Eps).^2 + imag(Eps).^2) + real(Eps))./2 );
k = sqrt((sqrt(real(Eps).^2 + imag(Eps).^2) - real(Eps))./2 );
Data = [1.24./w', n',k',real(Eps)',imag(Eps)'];
Data = sortrows(Data,1);
fid = fopen(strcat(diel_path, 'au_drude_', num2str(dT), 'K.tab'), 'wt');
fprintf(fid,' %s','drude model for au quantum beats');
fprintf(fid,'\n');
fprintf(fid,' %s','1 2 3 0 0 = Specifies n or Eps');
fprintf(fid,'\n');
fprintf(fid,' %s','lambda  n    k');
fprintf(fid,'\n');
for j = 1:length(Data)
    fprintf(fid,' %g',Data(j,1));
    fprintf(fid,' %g',Data(j,2));
    fprintf(fid,' %g',Data(j,3));
    fprintf(fid,'\n');
end
fclose(fid);

dT = temp_shortrod;
T = 298 + dT;
Eps_inf = -.0026*T+3.546;
wp = 0.0003*T + 9.0788;
gamma = (7e-5)*T+0.0337;
w = linspace(0.5,4,101);
w(1) =[];
Eps = Eps_inf - wp^2./(w.^2+1i*gamma*w);
n = sqrt((sqrt(real(Eps).^2 + imag(Eps).^2) + real(Eps))./2 );
k = sqrt((sqrt(real(Eps).^2 + imag(Eps).^2) - real(Eps))./2 );
Data = [1.24./w', n',k',real(Eps)',imag(Eps)'];
Data = sortrows(Data,1);
fid = fopen(strcat(diel_path, 'au_drude_', num2str(dT), 'K.tab'), 'wt');
fprintf(fid,' %s','drude model for au quantum beats');
fprintf(fid,'\n');
fprintf(fid,' %s','1 2 3 0 0 = Specifies n or Eps');
fprintf(fid,'\n');
fprintf(fid,' %s','lambda  n    k');
fprintf(fid,'\n');
for j = 1:length(Data)
    fprintf(fid,' %g',Data(j,1));
    fprintf(fid,' %g',Data(j,2));
    fprintf(fid,' %g',Data(j,3));
    fprintf(fid,'\n');
end
fclose(fid);
exit;
