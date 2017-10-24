function firstspike_output = analyze_latency(ds,dirname,mycell,mycellname,triggers)


firstspike_output.stim_med_latency = [];
firstspike_output.stim_med_jitter = [];


load([getpathname(ds) filesep dirname filesep 'stims.mat']);

n = numStims(saveScript);
do = getDisplayOrder(saveScript);
stimset_reps = length(do)/n;

ON_times = [];
OFF_times = [];

for i=1:numStims(saveScript),
    p = getparameters(get(saveScript,i));
    ON_times(end+1) = p.bgpause(1);
    OFF_times(end+1) = p.bgpause(2);
end;

ON_times = unique(ON_times);
OFF_times = unique(OFF_times);

for i=1:numStims(saveScript),
    p = getparameters(get(saveScript,i));
    on_index = find(p.bgpause(1)==ON_times); 
    off_index = find(p.bgpause(2)==OFF_times);
    first_spike = zeros(stimset_reps,p.repeat);
    mean_cycleduration = mean(mean(triggers{on_index,off_index}.cycleOnset(:,2:end)-triggers{on_index,off_index}.cycleOnset(:,1:end-1)));
    for j=1:stimset_reps,
        for k=1:p.repeat,
            if k < p.repeat,
                spikedata = get_data(mycell,[(triggers{on_index,off_index}.cycleOnset(j,k)) (triggers{on_index,off_index}.cycleOnset(j,k+1))]);
            else
                spikedata = get_data(mycell,[(triggers{on_index,off_index}.cycleOnset(j,k)) (triggers{on_index,off_index}.cycleOnset(j,k)+mean_cycleduration)]);
            end
            if ~isempty(spikedata),
                cycle_spike_times = spikedata - triggers{on_index,off_index}.cycleOnset(j,k);
            else
                cycle_spike_times = NaN;
            end
            first_spike(j,k) = cycle_spike_times(1);
        end
    end
    firstspike_output.stim_med_latency(on_index,off_index) = nanmedian(nanmedian(first_spike));
end

for i=1:numStims(saveScript),
    p = getparameters(get(saveScript,i));
    on_index = find(p.bgpause(1)==ON_times); 
    off_index = find(p.bgpause(2)==OFF_times);
    for j=1:stimset_reps,
        for k=1:p.repeat,
            stim_abs_dev(j,k) = abs(first_spike(j,k)-firstspike_output.stim_med_latency(on_index,off_index));
        end
    end
    firstspike_output.stim_med_jitter(on_index,off_index) = nanmedian(nanmedian(stim_abs_dev));
end


end


