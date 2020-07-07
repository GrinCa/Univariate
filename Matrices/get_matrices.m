function [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param)


disp('Reading files...');

for ii=1:length(matrix_names)
    matrix_names{ii} = strcat(mesh.file,'/',matrix_names{ii});
end

%--------------------------------------------------------------------------
%Run FreeFem++ script, IF NEEDED
%--------------------------------------------------------------------------

if (flag.rerun) % EDP updated and/or mesh updated
   t_0 = cputime;
   disp('************************');
   disp('*Rerun FreeFem++ script*');
   disp('************************');
   edpcommand = strcat('FreeFem++'," ",mesh.file,'.edp');
   system(edpcommand);
   timing.freefem = cputime-t_0;
   disp('*********************************************************');
   output = sprintf('[Get_matrices:infos] CPUtime for building of matrices %.4f s',timing.freefem);
   disp(output);
   disp('*********************************************************');
end % if

%--------------------------------------------------------------------------
% Get matrices from the files
%--------------------------------------------------------------------------

listLHS = cell(1,length(matrix_names)); 


% Matrices of the FE problem
tic;

for ii=1:length(matrix_names)
    nrows = 0;
    ncolums = 0;
    fid = fopen(matrix_names(ii),'rt');
    for jj=1:3
        line = fgets(fid);
        if jj==3
           values = str2num(strtrim(line));
           nrows = values(1);
           ncolums = values(2);
        end
    end
    fclose(fid);
    matrix_data = importdata(matrix_names(ii)," ",3);
    matrix_data = matrix_data.data;
    listLHS{ii} = sparse([matrix_data(:,1);nrows-1]+1,[matrix_data(:,2);ncolums-1]+1,[matrix_data(:,3);0]);
end
toc;
% Nodes
FEmatrices.Nodes = load(strcat(mesh.file,'/',"Nodes.txt"));
ndof = size(FEmatrices.Nodes,1);

%--------------------------------------------------------------------------
% return
%--------------------------------------------------------------------------
 
FEmatrices = build_global(FEmatrices,listLHS,param,mesh.file);

end


function FEmatrices = build_global(FEmatrices,listLHS,param,FILENAME)

% this function call the specific function for each given problem. You have
% to create a .m file for eash problem called Problem_pattern. For instance
% if Problem1 is one problem to study, then create Problem1_pattern. Then
% you need to compute it.

FEmatrices = Truck_pattern(FEmatrices,listLHS,param,FILENAME);

end

































