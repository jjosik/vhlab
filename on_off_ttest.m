function [ sub_P,m_ON_cycle_FR_inst,m_OFF_cycle_FR_inst ] = on_off_ttest( on_windowshift,off_windowshift,all_spike_cycles,on_offduration,c_pre_window )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%INPUTS
%ON_WINDOWSHIFT - positive time in seconds indicates beginning of
%observational window over which the cycle "on" FR is being assessed with
%respect to the beginning of stimulus (light) onset.
%OFF_WINDOWSHIFT - positive time in seconds indicates end of observational
%window for previous "on" cycle FR and beginning of observation window for
%"off" cycle spontaneous FR with respect to the time of stimulus off.
%**see getblinkingstim_stacked_cycles.m for variables ALL_SPIKE_CYCLES, 
% LIGHT_ON, and LIGHT_OFF


ON_cycle_FR_inst = cell((size(all_spike_cycles,2)*size(all_spike_cycles,1)),size(all_spike_cycles,3));
OFF_cycle_FR_inst = cell((size(ON_cycle_FR_inst)));
for i = 1:size(all_spike_cycles,3),     %for n stim. on/off combos
    %on_windowshift = ceil(stim_lat(i,1)-stim_jitter(i,1));  %*ALT. - declare stim_lat to
    %inputs to make shift sensitive to measured median latency for the
    %stim.
    %off_windowshift = ceil(stim_lat(i,1)+stim_jitter(i,1)); %*
    for j = 1:size(all_spike_cycles,2), %for # reps of each stim. seq.
        for k = 1:size(all_spike_cycles,1), %for # on/off cycles in each seq.
            spikes = all_spike_cycles{k,j,i};
            ON_xLimit = [on_windowshift,(on_offduration(i,1)+off_windowshift)];  %column 1 is when observation starts; column 2 when it stops
            OFF_xLimit = [(on_offduration(i,1)+off_windowshift),(on_offduration(i,1)+on_offduration(i,2))]; 
            if OFF_xLimit(1,2)-OFF_xLimit(1,1) <= off_windowshift,
                OFF_xLimit = [-c_pre_window, 0];
            end
            dt = 0.0001;  %time bin resolution for rate measurement
            t_vector_ON = ON_xLimit(1,1):dt:ON_xLimit(1,2);
            t_vector_OFF = OFF_xLimit(1,1):dt:OFF_xLimit(1,2);
            gauss_width = 200;  %width of gaussian smoothing window in units of dt
            Yn_ON = zeros(length(t_vector_ON)-1,1);
            Yn_OFF = zeros(length(t_vector_OFF)-1,1);
            Xn_ON = [];
            Xn_OFF = [];
            for m = 2:length(t_vector_ON),
                Xn_ON(end+1) = mean([t_vector_ON(m) t_vector_ON(m-1)]);
            end
            for n = 2:length(t_vector_OFF),
                Xn_OFF(end+1) = mean([t_vector_OFF(n) t_vector_OFF(n-1)]);
            end
            if ~isempty(spikes),
                spikes = reshape(spikes,1,length(spikes));
                %if on_offduration(i,1)<=on_windowshift, 
                %    spikes_ON = (spikes(find(spikes>0)));
                %else
                    spikes_ON = (spikes(find(spikes>t_vector_ON(1)&spikes<t_vector_ON(end))));
                %end
                if on_offduration(i,2)<=off_windowshift,
                    %spikes_OFF = (spikes(find(spikes<OFF_xLimit(1,2))));
                    spikes_OFF = (all_spike_cycles{1,j,i}(find(all_spike_cycles{1,j,i}>OFF_xLimit(1,1)&all_spike_cycles{1,j,i}<OFF_xLimit(1,2))));
                else
                    spikes_OFF = (spikes(find(spikes>t_vector_OFF(1)&spikes<t_vector_OFF(end))));
                end
                if ~isempty(spikes_ON),        
                    for q1 = 1:length(spikes_ON),
                        range_ON = length(find(t_vector_ON<spikes_ON(1,q1)));
                        if range_ON <= 1,
                            Yn_ON(length(find(t_vector_ON<spikes_ON(1,q1))),1) = 1;
                        else
                            Yn_ON(length(find(t_vector_ON<spikes_ON(1,q1)))-1,1) = 1;
                        end
                    end
                else
                    Yn_ON = zeros(length(t_vector_ON)-1,1);
                end
                if ~isempty(spikes_OFF),
                    for q2 = 1:length(spikes_OFF),
                        range_OFF = length(find(t_vector_OFF<spikes_OFF(1,q2)));
                        if range_OFF <= 1,
                            Yn_OFF(length(find(t_vector_OFF<spikes_OFF(1,q2))),1) = 1;
                        else
                            Yn_OFF(length(find(t_vector_OFF<spikes_OFF(1,q2)))-1,1) = 1;
                        end
                    end
                else
                    Yn_OFF = zeros(length(t_vector_OFF)-1,1);
                end
                g_kernel = gausswin(gauss_width,4);
                norm_g_kernel = g_kernel./sum(g_kernel);
                ON_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i} = (1/dt).*conv(Yn_ON,norm_g_kernel,'same');
                OFF_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i} = (1/dt).*conv(Yn_OFF,norm_g_kernel,'same');
            else
                ON_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i} = zeros(length(t_vector_ON)-1,1);
                OFF_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i} = zeros(length(t_vector_OFF)-1,1);
            end
            % 2 options for analysis: m prefix indicates mean inst. FR, p
            % denotes peak inst. FR
            m_ON_cycle_FR_inst((j-1)*(size(all_spike_cycles,1))+k,i) = nanmean(ON_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i}(:),1);  
            m_OFF_cycle_FR_inst((j-1)*size(all_spike_cycles,1)+k,i) = nanmean(OFF_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i}(:),1);
            %p_ON_cycle_FR_inst((j-1)*(size(all_spike_cycles,1))+k,i) = max(ON_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i}(:),[],1);
            %p_OFF_cycle_FR_inst((j-1)*(size(all_spike_cycles,1))+k,i) = max(OFF_cycle_FR_inst{((j-1)*(size(all_spike_cycles,1)))+k,i}(:),[],1);
        end
    end
    X = m_ON_cycle_FR_inst(:,i);    %alt. p_ON_cycle_FR_inst, see line 93 above (also change in outputs)
    Y = m_OFF_cycle_FR_inst(:,i);   %alt. p_OFF_cycle_FR_inst, see line 94 above (also change in outputs)
    [H,p] = ttest(X,Y);
    sub_P(i,1) = p;
    
end
    
            
            

end

