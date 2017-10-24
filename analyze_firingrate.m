function FR_output = analyze_firingrate(ds,triggers,stimset_reps,mycell,mycellname)


on_windowshift = 0.01;
off_windowshift = 0.02;
pretrain_off = 0.5;
FR_output = {};
for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        FR_output{i,j} = struct('cycle_FR',[],'mean_matchedcycle_FR',[],'norm_mean_matchedcycle_FR',[],'null_FR',[]);
    end
end
%FR_output = {};
%FR_output.cycle_FR = [];
%FR_output.mean_matchedcycle_FR = [];

for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        for k=1:stimset_reps,
            null_data = get_data(mycell,[(triggers{i,j}.cycleOnset(k,1) - pretrain_off) (triggers{i,j}.cycleOnset(k,1))]);
            FR_output{i,j}.null_FR(k,1) = (length(null_data))/pretrain_off;
            for m=1:size(triggers{i,j}.cycleOnset,2),
                spikedata = get_data(mycell,[(triggers{i,j}.cycleOnset(k,m)+on_windowshift) (triggers{i,j}.cycleOnset(k,m)+triggers{i,j}.pulseDuration+off_windowshift)]);
                if ~isempty(spikedata),
                    cycle_spike_times = spikedata - triggers{i,j}.cycleOnset(k,m);
                else
                    cycle_spike_times = NaN;
                end
                window_size = (triggers{i,j}.pulseDuration+off_windowshift)-on_windowshift;
                if ~isempty(cycle_spike_times),
                    FR_output{i,j}.cycle_FR(k,m) = length(cycle_spike_times)/window_size;
                else
                    FR_output{i,j}.cycle_FR(k,m) = 0;
                end 
            end
        end
    end
end

%matched cycle mean firing rates
for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        for k=1:size(triggers{i,j}.cycleOnset,2),
        FR_output{i,j}.mean_matchedcycle_FR = nanmean(FR_output{i,j}.cycle_FR,1);
        end
    end
end

%matched cycle mean firing rates normalized to first cycle
for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        for k=1:size(triggers{i,j}.cycleOnset,2),
        FR_output{i,j}.norm_mean_matchedcycle_FR(1,k) = FR_output{i,j}.mean_matchedcycle_FR(1,k)./FR_output{i,j}.mean_matchedcycle_FR(1,1);
        end
    end
end


end
            