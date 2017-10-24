function fourier_output = analyze_fourier(ds,triggers,stimset_reps,mycell,mycellname,use_windowing,cycle_window)

on_windowshift = 0.01;
off_windowshift = 0.02;
pretrain_off = 0.5;
fourier_output = {};
for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        fourier_output{i,j} = struct('cycle_stim_amp_coeff',[],'cycle_stim_power_spect',[],'cycle_phase_spectrum',[],...
            'mean_stim_amp_coeff',[],'mean_stim_power_spect',[],'mean_phase_spectrum',[],'stimtrain_f',[],'mean_power_error',[],'cycle_spike_times',[]);
    end
end

dt = 0.001;

for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        mean_cycleduration = mean(mean(triggers{i,j}.cycleOnset(:,2:end)-triggers{i,j}.cycleOnset(:,1:end-1)));
        if (cycle_window == 1)&&(on_windowshift<triggers{i,j}.pulseDuration),
            window_size = triggers{i,j}.pulseDuration-on_windowshift+off_windowshift;
            f_times = 0:dt:window_size;
        else
            window_size = triggers{i,j}.pulseDuration+off_windowshift;
            f_times = 0:dt:window_size;
        end
        for k=1:stimset_reps,
            for m=1:size(triggers{i,j}.cycleOnset,2),
                if (cycle_window == 1)&&(on_windowshift<triggers{i,j}.pulseDuration),
                    spikedata = get_data(mycell,[(triggers{i,j}.cycleOnset(k,m)+on_windowshift) (triggers{i,j}.cycleOnset(k,m)+window_size)]);
                else
                    spikedata = get_data(mycell,[(triggers{i,j}.cycleOnset(k,m)) (triggers{i,j}.cycleOnset(k,m)+window_size)]);
                end
                if isempty(spikedata),
                    cycle_spike_times = NaN;
                else
                    cycle_spike_times = spikedata - triggers{i,j}.cycleOnset(k,m);
                end
                %optional: pre-process data with windowing (Hanning or
                %Blackman)
                if use_windowing == 1,
                    w_series = cycle_spike_times.*blackman(length(triggers{i,j}(k,m)));
                else
                    w_series = cycle_spike_times;
                end
                spike_counts = spiketimes2bins(w_series,f_times);
                [fc,freqs] = fouriercoeffs(spike_counts,dt);
                fourier_output{i,j}.cycle_stim_amp_coeff{k,m} = reshape((sqrt(fc.*conj(fc))),length(fc),1);
                fourier_output{i,j}.cycle_stim_power_spect{k,m} = reshape((fc.*conj(fc)),length(fc),1);
                fourier_output{i,j}.cycle_phase_spectrum{k,m} = reshape((atan2(imag(fc),real(fc))),length(fc),1);
                fourier_output{i,j}.stimtrain_f{k,m} = freqs;
                fourier_output{i,j}.cycle_spike_times{k,m} = cycle_spike_times;
            end
        end
    end
end
 
for i=1:size(triggers,1),
    for j=1:size(triggers,2),
        fourier_output{i,j}.mean_stim_amp_coeff = 2/length(fourier_output{i,j}.cycle_stim_amp_coeff{1,1})...
            .*nanmean(cell2mat(reshape(fourier_output{i,j}.cycle_stim_amp_coeff,1,(size(fourier_output{i,j}.cycle_stim_amp_coeff,1))...
            *(size(fourier_output{i,j}.cycle_stim_amp_coeff,2)))),2);
        fourier_output{i,j}.mean_stim_power_spect = nanmean(cell2mat(reshape(fourier_output{i,j}.cycle_stim_power_spect,1,(size(fourier_output{i,j}.cycle_stim_power_spect,1))...
            *(size(fourier_output{i,j}.cycle_stim_power_spect,2)))),2);
        fourier_output{i,j}.mean_phase_spectrum = nanmean(cell2mat(reshape(fourier_output{i,j}.cycle_phase_spectrum,1,(size(fourier_output{i,j}.cycle_phase_spectrum,1))...
            *(size(fourier_output{i,j}.cycle_phase_spectrum,2)))),2);
        fourier_output{i,j}.mean_power_error = std(cell2mat(reshape(fourier_output{i,j}.cycle_stim_power_spect,1,(size(fourier_output{i,j}.cycle_stim_power_spect,1))...
            *(size(fourier_output{i,j}.cycle_stim_power_spect,2)))),1,2)./sqrt((size(fourier_output{i,j}.cycle_stim_power_spect,1))*(size(fourier_output{i,j}.cycle_stim_power_spect,2)));
    end
end

        
        
        
  
                
                
                