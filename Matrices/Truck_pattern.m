function FEmatrices = Truck_pattern(FEmatrices,listLHS,param,FILENAME)

ndof = size(FEmatrices.Nodes,1);

Hr = listLHS{1};
Hi = listLHS{2};
Qr  = listLHS{3};
Qi  = listLHS{4};

%definition of the regions
acoustic_region = 1;
PML_region = 2;

%definition of the labels
surface_label = 3;

tab_region = get_regions([acoustic_region,PML_region],ndof,FILENAME);
FEmatrices.acoustic_nodes = find(tab_region(:,1));
FEmatrices.PML_nodes = find(tab_region(:,2));

tab_labels = get_labels([surface_label],FILENAME);
FEmatrices.surf_nodes = find(tab_labels(:,1));

H = Hr+1i*Hi;
Q = Qr+1i*Qi;
FEmatrices.LHS = {H,Q};
FEmatrices.size_system = size(H,1);

end



function tab_region = get_regions(region_array,ndof,FILENAME)

connectivity_table = load(['Matrices/',FILENAME,'/connectivity_table.txt']);
region_element = load(['Matrices/',FILENAME,'/regions.txt']);

tab_region = zeros(ndof,length(region_array));

for ii=region_array
    id_elements = find(region_element == region_array(ii));
    for jj=1:length(id_elements)
        tab_region(connectivity_table(id_elements(jj),:)+1,ii) = 1;
    end
end

end


function tab_label = get_labels(label_number,FILENAME)

labels = load(['Matrices/',FILENAME,'/labels.txt']);
tab_label = zeros(length(labels(2:end,1)),length(label_number));

for ii=1:length(label_number)
    jj = find(labels(1,:) == label_number(ii));
    tab_label(:,ii) = labels(2:end,jj);
end

end
