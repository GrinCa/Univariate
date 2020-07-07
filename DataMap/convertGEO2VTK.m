function convertGEO2VTK(mesh,Nodes,Pout,freq,index)

%##########################################################################
%Please read the following lines for further informations
%##########################################################################
%this function aims to convert "Pout" array calculated with FE or WCAWE
%method into a .vtk file which Paraview can read. To understand how this
%scripts work, it is important to know what is the structure of .vtk file.
%You may find good info in the "vtk_fil_documentation.pdf" in the
%Documentation folder.
%We firstly need the Nodes file in order to get the coordinates of each
%node. Then we need Pout array (size = (ndof,nfreq)). It is possible that
%with the time, version compatibilities are no longer working. The main
%probleme of this type of file(.vtk) is that Paraview won't give you
%accurate information if it fails to read.
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

% lines_elements: temporary cell (size=(1,nelem)) which contain the lines
% from #ELEMENTS to #ENDELEMENT in the .msh file
% elements: the same than lines_elements, same size, but the lines are now 
% arrays 
% size_line: it counts for each line (therefore each elements) the 
% numbers of integers values on this line. 
%--------------------------------------------------------------------------
lines_elements = readFile(FILENAME);
elements = cell(1,length(lines_elements));
size_line = zeros(1,length(lines_elements));
%--------------------------------------------------------------------------

%convert the string line into array of those consitutives numbers. 
for ii=1:length(lines_elements)
    tmp = str2num(strtrim(lines_elements{ii}));
    elements{ii}=tmp;
    size_line(ii)=length(elements{ii});
end

%we switch into 2D or 3D regarding the max number of integer that a line
%can contains. For .msh file in the section #ELEMENTS, a line that contains
%the nodes of a triangle is about 8 integers, while for a tetra it is about
%9 integers. Therefore, if max(size_line)==8, it refers to 2D, otherwise,
%if ...==9, it is a 3D file.
switch max(size_line)
    case 8
        text_field = get_text_field(Nodes,elements,size_line,FILENAME,'2D');
        convertVTK2(text_field,Nodes,Pout,FILENAME,freq,index);
    case 9
        text_field = get_text_field(Nodes,elements,size_line,FILENAME,'3D');
        convertVTK3(text_field,Nodes,Pout,FILENAME,freq,index);
    otherwise
        disp('[convertVTK] the file seems to be neither 2D nor 3D, stoppping the convertion');
        return;
end
end

%--------------------------------------------------------------------------
% function to be potentially optimize
%--------------------------------------------------------------------------

function dataElementStr = readFile(FILENAME)

%--------------------------------------------------------------------------
fileID = fopen(strcat(FILENAME,'.msh'),'rt');
%--------------------------------------------------------------------------

while true
    dataFile=fgets(fileID);
    if strcmp(strtrim(dataFile),'$Elements')
        nelem=str2num(strtrim(fgets(fileID)));
        dataElementStr = cell(1,nelem);
        counter=1;
        while strcmp(strtrim(dataFile),'$EndElements')~=1
            dataFile=fgets(fileID);
            dataElementStr{counter}=strtrim(dataFile);
            counter=counter+1;
        end
        break;
    end
end
fclose(fileID);
end

function convertVTK2(text_field,Nodes,Pout,FILENAME,freq,index)
ndof = size(Nodes,1);
% -This function create as much .vtk file as we have frequencies.
% -text field contains all the text data that every .vtk file needs,
% and which is the same no matter the frequency considerate. It was
% therefore useful to get it once and for all, print it in each
% file, instead of "recalculating" it every loop iteration...
% -index is given by the user at the call of convertGEO2VTK and contains
% the IDs of the frequencies (regarding freq array) that we want to write 
% in memory.
disp('***Converting DATA to 2D VTK file***');
for ii=index
    file_name = strcat('DataMap/',FILENAME,'/',FILENAME,'_freq_',...
                        num2str(freq(ii)),'Hz_.vtk');
    if exist(file_name,'file') ~= 2
        fileID = fopen(file_name,'wt');
        %write of text_field into the .vtk file
        fprintf(fileID,text_field);
        % text_data is the data that will differentiate each .vtk file, it
        % will contain the values of the pressure @freq(ii)
        text_data = [];
        for jj=1:ndof
            text_data = [text_data [num2str(real(Pout(jj,ii))) '\n']];
        end
        fprintf(fileID,text_data);
        fclose(fileID);
    end
end
end

function convertVTK3(text_field,Nodes,Pout,FILENAME,freq,index)
% this function is almost the same than the 2D one.
ndof = size(Nodes,1);
disp('***Converting DATA to 3D VTK file***');
for ii=index
    file_name = strcat('DataMap/',FILENAME,'/',FILENAME,'_freq_',num2str(freq(ii)),'Hz_.vtk');
    if exist(file_name,'file') ~= 2
        fileID = fopen(file_name,'wt');
        fprintf(fileID,text_field);
        text_data = [];
        for jj=1:ndof
            text_data = [text_data [num2str(real(Pout(jj,ii))) '\n']];
        end
        fprintf(fileID,text_data);
        fclose(fileID);
    end
end
end


function text_field = get_text_field(Nodes,elements,size_line,FILENAME,type)
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
% here we differentiate 2D to 3D cases, because of triangle and
% tetrahedron. For more details on the format of the data, please refer to
% vtk_file_documentation.
if strcmp(type,'2D')
    disp('***Initialize conversion 2D***');
    IDtriangle = find(size_line==8); 
    % find triangles elements.Indeed .msh file can contains lines, for
    % boundaries conditions, that we don't want to keep...
    text_field = [text_field [ 'CELLS ',num2str(length(IDtriangle)),' ',...
                                        num2str(4*length(IDtriangle)),'\n']];
    for ii=IDtriangle
        % Note that .vtk files are -shifted, elements  for .msh file is  for
        % .vtk file
        text_field = [text_field ['3 ',num2str(elements{ii}(6)-1),' ',...
                                       num2str(elements{ii}(7)-1),' ',...
                                       num2str(elements{ii}(8)-1),'\n']];
    end
    text_field = [text_field ['CELL_TYPES ',num2str(length(IDtriangle)),'\n']];
    % 5 refers to triangles elements according to vtk files.
    for ii=1:length(IDtriangle)
        text_field = [text_field '5\n'];
    end
elseif strcmp(type,'3D')
    disp('***Initialize conversion 3D***');
    IDtriangle = find(size_line==9);
    text_field = [text_field ['CELLS ',num2str(length(IDtriangle)),' ',...
                                       num2str(5*length(IDtriangle)),'\n']];
    for ii=IDtriangle
        text_field = [text_field ['4 ',num2str(elements{ii}(6)-1),' ',...
                                       num2str(elements{ii}(7)-1),' ',...
                                       num2str(elements{ii}(8)-1),' ',...
                                       num2str(elements{ii}(9)-1),'\n']];
    end
    text_field = [text_field ['CELL_TYPES ',num2str(length(IDtriangle)),'\n']];
    for ii=1:length(IDtriangle)
        text_field = [text_field '10\n'];
    end
end

text_field = [text_field ['POINT_DATA ' num2str(ndof) '\n']];
text_field = [text_field 'SCALARS cell_scalars float 1\n'];
text_field = [text_field 'LOOKUP_TABLE default\n'];
%after this line, convertVTK2 or 3 will then fill the data required 
end


