function param = build_interval(param)

%build sub_range with the different freqref
index_freqref = zeros(1,param.nfreqref);
for ii=1:param.nfreqref
    [~,index_freqref(ii)] = min(abs(param.freq-param.freqref(ii)));
end

param.sub_range = cell(1,length(param.interval_construct));
param.n_sub_range = length(param.sub_range);
id_cut = zeros(1,length(param.interval_construct)-1);
for ii=1:length(param.interval_construct)-1
    if param.interval_construct{ii}(end) == param.interval_construct{ii+1}(1)
        id_cut(ii) = index_freqref(param.interval_construct{ii}(end));
    else
        id_cut(ii) = int16(( index_freqref(param.interval_construct{ii}(end)) + ...
                             index_freqref(param.interval_construct{ii+1}(1)) )/2);
    end
end

caracteristic_index = [0 id_cut param.nfreq];

if isempty(id_cut)
    param.sub_range{1} = param.freq;
else
    for ii=1:(length(caracteristic_index)-1)
        left_edge = caracteristic_index(ii)+1;
        right_edge = caracteristic_index(ii+1);
        param.sub_range{ii} = param.freq(left_edge) : param.freqincr : param.freq(right_edge);
    end
end

end

