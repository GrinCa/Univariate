function show_graph(arg,param)


if strcat(arg.config,'comparison')
    show2d(arg,param)
end
    
end


function show2d(arg,param)
    figure
    if arg.log
        for ii=1:length(arg.VALUES)
            plot(param.freq,log10(arg.VALUES{ii}),'DisplayName',arg.label{ii});
            hold on
        end
    else
        for ii=1:length(arg.VALUES)
            plot(param.freq,arg.VALUES{ii},'DisplayName',arg.label{ii});
            hold on
        end
        if arg.external_plot.is_needed
            for ii=1:length(arg.external_plot.VALUES)
                plot(param.freq,arg.external_plot.VALUES{ii},'DisplayName',arg.external_plot.label{ii});
                hold on
            end
        end
    end
    
    xlabel(arg.xlabel);
    ylabel(arg.ylabel);
    title(arg.title);
    legend;
    saveas(gcf,[arg.save_path,'/',arg.save_name,'.png']);
end
