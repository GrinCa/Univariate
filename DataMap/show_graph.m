function show_graph(arg,VALUES,title_VALUES,param)


if strcat(arg.config,'comparison')
    show2d(VALUES,title_VALUES,param,arg)
end
    
end


function show2d(VALUES,title_VALUES,param,arg)
    figure
    if arg.log
        for ii=1:length(VALUES)
            plot(param.freq,log10(VALUES{ii}),'DisplayName',title_VALUES{ii});
            hold on
        end
    else
        for ii=1:length(VALUES)
            plot(param.freq,VALUES{ii},'DisplayName',title_VALUES{ii});
            hold on
        end
    end
    
    xlabel(arg.xlabel);
    ylabel(arg.ylabel);
    title(arg.title);
    legend;
    saveas(gcf,[arg.save_path,'/',arg.save_name,'.png']);
end
