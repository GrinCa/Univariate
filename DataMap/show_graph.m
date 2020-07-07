function show_graph(arg,data,mesh,param,folder)

if strcmp(arg,'compare_results')
    sizemesh = load('sizemesh.txt');
    if length(sizemesh) > 1
        disp('There are several sizemesh');
        disp('Please choose one amoung the following :');
        for iii=1:length(sizemesh)
            disp([num2str(iii),' -> sizemesh = ',num2str(sizemesh(iii))]);
        end
    end
    id_sample =  ['_sizemesh_',num2str(sizemesh)];
    SOLFE =      struct2cell(load(['Matrices/',mesh.file,'/',folder.path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
    SOLMDWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',folder.path2,'/SOLMDWCAWE',id_sample,'.mat']));
    SOLWCAWE =   struct2cell(load(['Matrices/',mesh.file,'/',folder.path2,'/SOLWCAWE',id_sample,'.mat']));
    
    SOLFE = SOLFE{1};
    SOLMDWCAWE = SOLMDWCAWE{1};
    SOLWCAWE = SOLWCAWE{1};
    
    
    
    if ~isempty(find(size(SOLFE)==1)) || length(size(SOLFE)) == 2
        Norm_SOLFE = zeros(param.nfreq,param.ntheta);
        Norm_SOLMDWCAWE = zeros(param.nfreq,param.ntheta);
        Norm_SOLWCAWE = zeros(param.nfreq,param.ntheta);
        if length(size(SOLFE)) == 2
            for ii=1:size(SOLFE,2)
                Norm_SOLFE(ii) = norm(SOLFE(:,ii));
                Norm_SOLMDWCAWE(ii) = norm(SOLMDWCAWE(:,ii));
                Norm_SOLWCAWE(ii) = norm(SOLWCAWE(:,ii));
            end
        else
            for ii=1:size(SOLFE,3)
                Norm_SOLFE(ii) = norm(SOLFE(:,1,ii));
                Norm_SOLMDWCAWE(ii) = norm(SOLMDWCAWE(:,1,ii));
                Norm_SOLWCAWE(ii) = norm(SOLWCAWE(:,1,ii));
            end
        end
        show2d(SOLFE,SOLMDWCAWE,SOLWCAWE,param);
    else
        Norm_SOLFE = zeros(param.nfreq,param.ntheta);
        Norm_SOLMDWCAWE = zeros(param.nfreq,param.ntheta);
        Norm_SOLWCAWE = zeros(param.nfreq,param.ntheta);
        for ii=1:size(SOLFE,2)
            for jj=1:size(SOLFE,3)
                Norm_SOLFE(ii,jj) = norm(SOLFE(:,ii,jj));
                Norm_SOLMDWCAWE(ii,jj) = norm(SOLMDWCAWE(:,ii,jj));
                Norm_SOLWCAWE(ii,jj) = norm(SOLWCAWE(:,ii,jj));
            end
        end
        show3d(Norm_SOLFE,Norm_SOLMDWCAWE,Norm_SOLWCAWE,param,folder,mesh);
    end

elseif strcat(arg,'plotTL')
    if ~isempty(find(size(data{1})==1))
        show2d(data{1},data{2},data{3},param,folder,mesh)
    else
        show3d(data{1},data{2},data{3},param,folder,mesh)
    end
end
    
end


function show2d(SOLFE_VALUES,SOLMDWCAWE_VALUES,SOLWCAWE_VALUES,param,folder,mesh)
    pathPlot = ['Matrices/' mesh.file '/' folder.path2];
    X = {param.freq,param.theta};
    sizesubset = [param.nfreq,param.ntheta];
    xlabeltxt = {"freq",'theta'};
    figure
    plot(X{find(sizesubset>1)},SOLFE_VALUES);
    xlabel(xlabeltxt{find(sizesubset>1)});
    ylabel("TL");
    title("SOLFE");
    saveas(gcf,strcat(pathPlot,'/','SOLFE','.png'));
    figure
    plot(X{find(sizesubset>1)},SOLMDWCAWE_VALUES);
    xlabel(xlabeltxt{find(sizesubset>1)});
    ylabel("TL");
    title("SOLMDWCAWE");
    saveas(gcf,strcat(pathPlot,'/','SOLMDWCAWE','.png'));
    figure
    plot(X{find(sizesubset>1)},SOLWCAWE_VALUES);
    xlabel(xlabeltxt{find(sizesubset>1)});
    ylabel("TL");
    title("SOLWCAWE");
    saveas(gcf,strcat(pathPlot,'/','SOLWCAWE','.png'));
end

function show3d(FE_VALUES,MDWCAWE_VALUES,WCAWE_VALUES,param,folder,mesh)
    
    pathPlot = ['Matrices/' mesh.file '/' folder.path2];
    [X,Y] = meshgrid(param.freq,param.theta);
    figure
    surf(X,Y,FE_VALUES');
    xlabel('freq');
    ylabel('theta');
    zlabel('TL');
    title_label = 'SOLFE';
    title(title_label);
    saveas(gcf,[pathPlot '/' title_label '.png']);
    figure
    surf(X,Y,MDWCAWE_VALUES');
    xlabel('freq');
    ylabel('theta');
    zlabel('TL');
    title_label = ['SOLMDWCAWE [' num2str(param.n_sub_range_freq*param.vecfreqrange) ',' num2str(param.n_sub_range_theta*param.vecthetarange) ']'];
    title(title_label);
    saveas(gcf,strcat(pathPlot,'/',title_label,'.png'));
    figure
    surf(X,Y,WCAWE_VALUES');
    xlabel('freq');
    ylabel('theta');
    zlabel("TL");
    title_label = ['SOLWCAWE [' num2str(param.n_sub_range_freq*param.vecfreqrange) ',' num2str(param.n_sub_range_theta*param.vecthetarange) ']'];
    title(title_label);
    saveas(gcf,strcat(pathPlot,'/',title_label,'.png'));
    figure
    surf(X,Y,abs((FE_VALUES-MDWCAWE_VALUES)./FE_VALUES)');
    xlabel('freq');
    ylabel('theta');
    zlabel('Relative error');
    title_label = ['Error MDWCAWE [' num2str(param.n_sub_range_freq*param.vecfreqrange) ',' num2str(param.n_sub_range_theta*param.vecthetarange) ']'];
    title(title_label);
    saveas(gcf,strcat(pathPlot,'/',title_label,'.png'));
    figure
    surf(X,Y,abs((FE_VALUES-WCAWE_VALUES)./FE_VALUES)');
    xlabel('freq');
    ylabel('theta');
    zlabel('Relative error');
    title_label = ['Error WCAWE [' num2str(param.n_sub_range_freq*param.vecfreqrange) ',' num2str(param.n_sub_range_theta*param.vecthetarange) ']'];
    title(title_label);
    saveas(gcf,strcat(pathPlot,'/',title_label,'.png'));
end

function showTL()

end

