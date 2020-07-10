function genfolders(mesh,param)

Filename = mesh.file;

%--------------------------------------------------------------------------
% Matrices
%--------------------------------------------------------------------------

path1 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']'];

if exist(['Matrices/',Filename]) == 0
    command = ['Matrices/',Filename];
    system(['mkdir ' command]);
end

if exist(['Matrices/',Filename,'/',path1]) == 0
    command = ['Matrices/',Filename,'/',['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']'];];
    system(['mkdir ' command]);
end

for ii=1:length(param.vecfreqrange)
    if exist(['Matrices/',Filename,'/',path1,'/[' num2str(param.vecfreqrange(ii)) '][',replace(num2str(param.freqref),' ','_') ']']) == 0
    command = ['mkdir Matrices/',Filename,'/',path1,'/[' num2str(param.vecfreqrange(ii)) '][',replace(num2str(param.freqref),' ','_') ']'];
    system(command);
    end
end


%--------------------------------------------------------------------------
% Geometry : copy .msh from Geometry/filename into the main folder
%--------------------------------------------------------------------------

try
    command = ['cp ','Geometry/',Filename,'/',Filename,'.msh',' ',Filename,'.msh'];
    system(command);
catch
    disp(['[genfolders] No .msh file in the folder : /Geometry/',Filename]);
    return;
end
%--------------------------------------------------------------------------
% DataMap
%--------------------------------------------------------------------------

if exist(['DataMap/',Filename]) == 0
    command = ['DataMap/',Filename];
    system(['mkdir ' command]);
end

end