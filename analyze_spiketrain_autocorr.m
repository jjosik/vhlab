function [autocorr_output,on_offduration] = analyze_spiketrain_autocorr(triggers,stimset_reps,mycell,mycellname,cycle_window,save_cellplots,...
    normalize_cellcorr,analyze_band,modulation_band,varargin)

%%MODULATION_BAND - two value cell containing lower boundary and upper
%%boundary of frequency range of band modulation to be analyzed, e.g. [25
%%50] will restrict modulation index analysis to the 25-50 Hz range

autocorr_output = {};
for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        autocorr_output{i,j} = struct('cycle_corr',[],'cycle_lags',[],'stim_corr',[],'stim_lags',[],'mean_stim_corr',[]);
    end
end

if abs(nargin) > 9,
    bin_size = varargin{1}.bin_size;
    on_windowshift = 0.01;
    off_windowshift = -(0.5-varargin{1}.total_w)+on_windowshift;
else
    bin_size = 0.001;
    on_windowshift = 0.01; %needs to be set for special windowing; consider adding to inputs
    off_windowshift = 0.02; %needs to be set for special windowing; consider adding to inputs
end
on_offduration = cell(size(triggers,1),size(triggers,2),2);
for i=1:size(triggers,1),
    for j=1:size(triggers,2), 
        for k=1:stimset_reps,
            for m=1:size(triggers{i,j}.cycleOnset,2),
                if (cycle_window == 1)&&(on_windowshift<triggers{i,j}.pulseDuration)&&(abs(off_windowshift)<triggers{i,j}.pulseDuration),
                    window_size = triggers{i,j}.pulseDuration-on_windowshift+off_windowshift;
                    spikedata = get_data(mycell,[(triggers{i,j}.cycleOnset(k,m)+on_windowshift) (triggers{i,j}.cycleOnset(k,m)+window_size)]);
                elseif (cycle_window == 1)&&(on_windowshift<triggers{i,j}.pulseDuration)&&(abs(off_windowshift)>triggers{i,j}.pulseDuration),
                    window_size = triggers{i,j}.pulseDuration;
                    spikedata = get_data(mycell,[(triggers{i,j}.cycleOnset(k,m)) (triggers{i,j}.cycleOnset(k,m)+window_size)]);
                else
                    window_size = triggers{i,j}.pulseDuration+off_windowshift;
                    spikedata = get_data(mycell,[(triggers{i,j}.cycleOnset(k,m)) (triggers{i,j}.cycleOnset(k,m)+window_size)]);
                end
                bins = 0:bin_size:window_size;
                if isempty(spikedata),
                    cycle_spike_times = NaN;
                else
                    cycle_spike_times = spikedata - triggers{i,j}.cycleOnset(k,m);
                end
                ach_bin_edges = [-(window_size):bin_size:window_size];
                N_ach = histc(cycle_spike_times,ach_bin_edges);
                ach_bin_centers = (ach_bin_edges(1:end-1)+ach_bin_edges(2:end))/2;
                N_ach = N_ach(1:end-1);
                [corr,lags] = xcorr(N_ach,N_ach,round(window_size/bin_size));
                lags = lags*bin_size;
                autocorr_output{i,j}.cycle_corr{k,m} = reshape(corr,length(corr),1);
                autocorr_output{i,j}.cycle_lags{k,m} = reshape(lags,length(lags),1);
            end
        end
        on_offduration{i,j,1} = triggers{i,j}.ONcommand;
        on_offduration{i,j,2} = triggers{i,j}.OFFcommand;
    end
end

for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        autocorr_output{i,j}.mean_stim_corr = nanmean(cell2mat(reshape(autocorr_output{i,j}.cycle_corr,1,(size(autocorr_output{i,j}.cycle_corr,1)...
            *size(autocorr_output{i,j}.cycle_corr,2)))),2);
        peak_value = max(autocorr_output{i,j}.mean_stim_corr);
        autocorr_output{i,j}.norm_mean_stim_corr = (nanmean(cell2mat(reshape(autocorr_output{i,j}.cycle_corr,1,(size(autocorr_output{i,j}.cycle_corr,1)...
            *size(autocorr_output{i,j}.cycle_corr,2)))),2))./peak_value;
        autocorr_output{i,j}.stim_lags = autocorr_output{i,j}.cycle_lags{1,1};
    end
end

if analyze_band == 1&&save_cellplots==1,
    start_lagaverage = 1/modulation_band(1,2);
    stop_lagaverage = 1/modulation_band(1,1);
    lag_background_range = [-0.1 0.1];
    for i = 1:size(autocorr_output,2),
        compare_indices_start = abs(start_lagaverage - autocorr_output{4,i}.stim_lags);
        compare_indices_stop = abs(stop_lagaverage - autocorr_output{4,i}.stim_lags);
        band_start_index = find(compare_indices_start==min(compare_indices_start));
        band_stop_index = find(compare_indices_stop==min(compare_indices_stop));
        comp_background_start = abs(lag_background_range(1,1) - autocorr_output{4,i}.stim_lags);
        comp_background_stop = abs(lag_background_range(1,2) - autocorr_output{4,i}.stim_lags);
        background_start_ind = find(comp_background_start==min(comp_background_start));
        background_stop_ind = find(comp_background_stop==min(comp_background_stop));
        mean_band(1,i) = nanmean(autocorr_output{4,i}.norm_mean_stim_corr(band_start_index:band_stop_index));
        median_background(1,i) = nanmedian(autocorr_output{4,1}.norm_mean_stim_corr(background_start_ind:background_stop_ind));
        modulation_index(1,i) = (mean_band(1,i) - median_background(1,i))/(mean_band(1,i) + median_background(1,i));
    end
else
end

if save_cellplots == 1,
    for i=1:size(triggers,1),
        for j=1:size(triggers,2),
            f = figure;
            if normalize_cellcorr == 1,
                ach_h = bar(autocorr_output{i,j}.stim_lags,autocorr_output{i,j}.norm_mean_stim_corr,'k');
            else
                ach_h = bar(autocorr_output{i,j}.stim_lags,autocorr_output{i,j}.mean_stim_corr,'k');
            end
            hold on;
            xlabel('Lags');
            ylabel('Correlation');
            if abs(off_windowshift)<triggers{i,j}.pulseDuration
                xlim([-(triggers{i,j}.pulseDuration+off_windowshift) (triggers{i,j}.pulseDuration+off_windowshift)]);
            else
                xlim([-(triggers{i,j}.pulseDuration+0.01) (triggers{i,j}.pulseDuration+0.01)]);
            end
            if normalize_cellcorr == 1,
                ylabel('Correlation');
                ylim([0 1.1]);
            else
                ylabel('Conditional rate');
                ylim([0 max(autocorr_output{i,j}.mean_stim_corr)+0.1+0.1*(max(autocorr_output{i,j}.mean_stim_corr))]);
            end
            if i==4&&(analyze_band == 1),
                yt = 0.9;
                xt = max(autocorr_output{i,j}.stim_lags)-(0.9*max(autocorr_output{i,j}.stim_lags));
                text(xt,yt,[num2str(modulation_band(1,1)),'-',num2str(modulation_band(1,2)),' Hz band MI = ',num2str(modulation_index(1,j))]);                                                                                                                                                                                                   
            else
            end
            saveas(gcf,[pwd filesep mycellname '_ACH_' num2str(triggers{i,j}.ONcommand) '_' num2str(triggers{i,j}.OFFcommand) '.fig']);
            hold off;
            close(f);
        end
    end
end
        
end
        
        
        
        
        