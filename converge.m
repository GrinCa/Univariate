%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Convergence of the solution                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this script is useful to test the convergence of the method. It uses the
%.geo file and change the size of the elements, in order to change the
%number of nodes.
%%%
% Note : this algorithm use version 3.0.6 of Gmsh and version 4.1 of
% FreeFem. The algorithm is likely to fail if the version of Gmsh is higher
% than 3.0.6

clear all;

Geofile = 'Truck';

% both arrays stand for size of elastic and acoustic nodes respectively.
% For each colum of both arrays, the main script will be run. Therefore it
% is possible to put as much as values that we desire, but the size of each
% array must be the same. 
% NOTE : If the size of elastic and acoustic is different, the mesh
% elements won't be regular, which implies that the mean quadratic pressure
% is no longer a good indicator.
sizemesh_list = [1 1.2 1.5];
n_mesh = length(sizemesh_list);
% creation of the file containing the size(s) of the mesh 
system("touch sizemesh.txt");

% data to change values in the .geo file before the compilation into .msh
keyword = 'sizemesh = ';
n_key = length(keyword);

% Import existing .geo file
path = ['Geometry/',Geofile,'/'];
fid = fopen([path,Geofile,'.geo'],'rt');
dataMesh = fread(fid);
fclose(fid);
dataMesh = char(dataMesh.');

for kk=1:n_mesh
    % Find starting and end indexes of substring to replace
    startindkey = strfind(dataMesh,keyword);
    endlocate = strfind(dataMesh,';');
    endlocate_idx = find((endlocate-startindkey)>0);
    endind = endlocate(endlocate_idx(1));
    startindval = startindkey+length(keyword);
    % Test if value is to be updated, update if so
    current_val = str2double(dataMesh(startindval:endind-1));

    if current_val ~= sizemesh_list(kk)
        dataMesh = strrep(dataMesh, dataMesh(startindkey:endind), [keyword, num2str(sizemesh_list(kk)),';']);
        fid = fopen([path,Geofile, '.geo'],'wt');
        fwrite(fid,dataMesh);
        fclose (fid);
    end
    %save sizemesh to acces in the main_multi.m file.
    file_sizemesh = fopen('sizemesh.txt','a');
    fprintf(file_sizemesh,[num2str(sizemesh_list(1,kk)),'\n']);
    fclose(file_sizemesh);
    
    % compile .geo file
    disp('Compile .geo file...');
    command = ['/home/julien/Documents/Softwares/gmsh-3.0.6/bin/gmsh -3 ',path,Geofile,'.geo ','-o ',path,Geofile,'.msh']; 
    % version gmsh : 3.0.6, whiwh match with FreeFem++ 4.1. Antother
    % version could not work
    system(command);
    main_multi;
end







