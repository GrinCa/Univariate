function convertGEO2VTK(FEmatrices,mesh,sizemesh,SOL,PARTITION,param,range)

%##########################################################################
%Please read the following lines for further informations
%##########################################################################
%this function aims to convert "SOL" array calculated with FE or WCAWE
%method into a .vtk file which Paraview can read. To understand how this
%scripts work, it is important to know what is the structure of .vtk file.
%You may find good info in the "vtk_fil_documentation.pdf" in the
%Documentation folder.
%We firstly need the Nodes file in order to get the coordinates of each
%node. Then we need SOL array (size = (i*ndof,nfreq), i=number of function). 
%It is possible thatwith the time, version compatibilities are no longer 
%working. The main probleme of this type of file(.vtk) is that Paraview 
%won't give you accurate information if it fails to read.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To sum up .vtk file contains(in this particular case):
% -header: gives the version of the vtk file
% -DATA_SET type : in our case type=UNSCTRUCTURED_GRID which enables us to
% choose as we want triangle elements for 2D or tetrahedral elements for 3D
% -POINTS : content of the Nodes file (coordinates)
% -CELLS : which refers to the .msh file. It contains the ID of the nodes
% of each elements
% -CELLS_TYPE : contains the ID of the element to use, according .vtk
% files. For instance triangles ID=5, and tetra ID=10.
% -POINT_DATA : can take several data, from scalar for pressure to tensor
% for stress field. In our case it is SCALAR id.


% FILENAME refers to .msh file and it enables us to create the path to
% store .vtk files
FILENAME = mesh.file;
connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
connectivity_table = connectivity_table(:,[1 2 3 4 5 8 6 7 9 10]);
text_field = get_text_field(FEmatrices.Nodes,connectivity_table,FILENAME);
convertVTK3(FEmatrices,text_field,SOL,PARTITION,FILENAME,param,range,sizemesh)

end


function convertVTK3(FEmatrices,text_field,SOL,PARTITION,FILENAME,param,range,sizemesh)
% -This function create as much .vtk file as we have frequencies.
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
% -index is given by the user at the call of convertGEO2VTK and contains
% the IDs of the frequencies (regarding freq array) that we want to write 
% in memory.
ndof = size(FEmatrices.Nodes,1);

SOLUTION = real(SOL);

for ii=range{1}
    disp(['*** Converting FREQ = [',num2str(param.freq(ii)),'] ***']);
    file_name = strcat('DataMap/',FILENAME,'/',FILENAME,'_sizemesh_',num2str(sizemesh),'_freq_',num2str(param.freq(ii)),'.vtk');
    fileID = fopen(file_name,'wt');
    fprintf(fileID,text_field);

    text_data = [];
    text_data = [text_data 'POINT_DATA ' num2str(ndof) '\n'];

    for n=1:length(PARTITION)
        if strcmp(PARTITION{n}{1},'scalar')
            text_data = convert_scalar(FEmatrices,text_data,SOLUTION,{ii},PARTITION{n}{2},PARTITION{n}{3},PARTITION{n}{4});
        elseif strcmp(PARTITION{n}{1},'vector')
            text_data = convert_vector(FEmatrices,text_data,SOLUTION,{ii},PARTITION{n}{2},PARTITION{n}{3},PARTITION{n}{4});
        elseif strcmp(PARTITION{n}{1},'data')
            text_data = convert_given_data(FEmatrices,text_data,real(PARTITION{n}{2}(:,ii)),PARTITION{n}{3});
        end
    end
    fprintf(fileID,text_data);
    fclose(fileID);
end
end


function text_field = get_text_field(Nodes,connectivity_table,FILENAME)
ndof = size(Nodes,1);
text_field = [];
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
text_field = [text_field ['# vtk DataFile Version 2.0\n',FILENAME,'\nASCII\n']];
% above is the header refering to the version of vtk. It might have
% changed...
text_field = [text_field 'DATASET UNSTRUCTURED_GRID\n'];       %refer to doc
text_field = [text_field ['POINTS ',num2str(ndof),' float\n']];%refer to doc
% wrinting coordinates of each node
for ii=1:ndof
    text_field = [text_field [num2str(Nodes(ii,1)),' ',...
                              num2str(Nodes(ii,2)),' ',...
                              num2str(Nodes(ii,3)),'\n']];
end

disp('*** Initialize conversion 3D ***');
text_field = [text_field ['CELLS ',num2str(size(connectivity_table,1)),' ',...
                                   num2str(11*size(connectivity_table,1)),'\n']];
for ii=1:size(connectivity_table,1)
    text_field = [text_field ['10 ',num2str(connectivity_table(ii,:))],'\n'];
end 
text_field = [text_field ['CELL_TYPES ',num2str(size(connectivity_table,1)),'\n']];
for ii=1:size(connectivity_table,1)
    text_field = [text_field '24\n'];
end


end

function text_data = convert_vector(FEmatrices,text_data,SOLUTION,index,local_nodes,global_nodes,component_name)

ndof = size(FEmatrices.Nodes,1);
indexfreq = index{1};

VECTOR = zeros(ndof,3);

VECTOR(global_nodes,1) = SOLUTION(local_nodes(1,:),indexfreq);
VECTOR(global_nodes,2) = SOLUTION(local_nodes(2,:),indexfreq);
VECTOR(global_nodes,3) = SOLUTION(local_nodes(3,:),indexfreq);

text_data = [text_data 'VECTORS ' component_name ' float\n'];
%text_data = [text_data 'LOOKUP_TABLE default\n'];

for kk=1:ndof
    text_data = [text_data num2str(VECTOR(kk,1)) ' '...
                           num2str(VECTOR(kk,2)) ' '...
                           num2str(VECTOR(kk,3)) '\n'];
end
end

function text_data = convert_scalar(FEmatrices,text_data,SOLUTION,index,local_nodes,global_nodes,component_name)

ndof = size(FEmatrices.Nodes,1);
indexfreq = index{1};

text_data = [text_data 'SCALARS ' component_name ' float 1\n'];
text_data = [text_data 'LOOKUP_TABLE default\n'];

SCALAR = zeros(ndof,1);
SCALAR(global_nodes) = SOLUTION(local_nodes,indexfreq);

for kk=1:ndof
    text_data = [text_data num2str(SCALAR(kk)) '\n'];
end


end

function text_data = convert_given_data(FEmatrices,text_data,DATA,component_name)

ndof = size(FEmatrices.Nodes,1);

text_data = [text_data 'SCALARS ' component_name ' float 1\n'];
text_data = [text_data 'LOOKUP_TABLE default\n'];

for kk=1:ndof
    text_data = [text_data num2str(DATA(kk)) '\n'];
end

end


