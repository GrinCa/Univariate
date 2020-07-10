function show_graph(arg,VALUES,title_VALUES,mesh,param)


if strcat(arg.config,'comparison')
    show2d(VALUES,title_VALUES,param,mesh,arg)
end
    
end


function show2d(VALUES,title_VALUES,param,mesh,arg)
    %pathPlot = ['Matrices/' mesh.file '/' folder.path2];
    figure
    for ii=1:length(VALUES)
        plot(param.freq,VALUES{ii},'DisplayName',title_VALUES{ii});
        hold on
    end
    xlabel(arg.xlabel);
    ylabel(arg.ylabel);
    title(arg.title);
    legend;
    %saveas(gcf,strcat(pathPlot,'/',arg.save_name,'.png'));
end
