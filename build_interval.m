function param = build_interval(param)


%build sub_range with the different freqref
index_freqref = zeros(1,param.nfreqref);
for ii=1:param.nfreqref
    [~,index_freqref_ii] = min(abs(param.freq-param.freqref(ii)));
    index_freqref(ii) = index_freqref_ii;
end

caracteristic_index = [1 index_freqref param.nfreq];
param.sub_range = cell(1,param.nfreqref);
param.n_sub_range = length(param.sub_range);
if param.nfreqref==1
    param.sub_range{1} = param.freq;
else
    for ii=1:param.nfreqref
        if ii==1
            left_edge = 1;
            right_edge = floor(mean([caracteristic_index(ii+1) caracteristic_index(ii+2)]));
            param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
        elseif ii==param.nfreqref
            left_edge = floor(mean([caracteristic_index(ii) caracteristic_index(ii+1)])) + 1;
            right_edge = param.nfreq;
            param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
        else
            left_edge = floor(mean([caracteristic_index(ii) caracteristic_index(ii+1)])) + 1;
            right_edge = floor(mean([caracteristic_index(ii+1) caracteristic_index(ii+2)]));
            param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
        end
    end
end


end