function [ stimulus_triggers, all_spike_cycles, all_cycles_response, cycle_FR_inst, FR_cycle_lineavg, cycle_meanFR, stim_index_all,on_offduration,on_offreps ] = ...
    getblinkingstim_stacked_cycles( mycell, mycellname, c_pre_window, c_post_window, plot_it )
%GETBLINKINGSTIM_STACKED_CYCLES - further unpacks celldata under MYCELL
%for a collection of cells (listed under MYCELLNAME).  Extracts spiketimes for all cycles of a
%given on/off stimulus combination and stacks all 60 cycles in raster plot.  Analyzes
%firing rates, generates a PSTH matched to the rasterplot, and
%calculates and plots a gaussian smoothed instantaneous firing rate plot averaged over
%the cycle data.
% ** INPUTS **
%Plot and analysis window equal to length of each cycle plus a leading time, C_PRE_WINDOW,
%and a following time, C_POST_WINDOW.  If desired analysis window is only
%the cycle itself, both of these should be set to '0'.  Set PLOT_IT = 0 if plots should be
%saved and not displayed, '1' otherwise.
% ** OUTPUTS **
%       STIMULUS_TRIGGERS - stimulus trigger times
%       ALL_SPIKE_CYCLES - spiketimes for entire dataset, organized into
%       multi-dimensional cell arrays.  DIM 1: for cycle number 1-10; DIM
%       2: for unique stim rep. 1-6; DIM 3: for stim on/off combination
%       1-16
%       ALL_CYCLES_RESPONSE - for all cycles, returns '1' if any spikes are
%       present, '0' if empty
%       CYCLE_FR_INST - 
%       FR_CYCLE_LINEAVG - 
%       CYCLE_MEANFR - 
%       STIM_INDEX_ALL - Each row shows all locations of stim (on/off
%       combo) # n in order of full stimulus sequence.
%       ON_OFFDURATION - column 1 is time light is 'on', column 2 is time
%       light is 'off' for all 16 (rows) stim. combos
%       ON_OFFREPS - number of cycles in single stim. train 


%cd(dirname);
load('stims.mat');
n = numStims(saveScript);
do = getDisplayOrder(saveScript);
stimset_reps = length(do)/n;
[mti2_,starttime_] = tpcorrectmti(MTI2,[pwd filesep 'stimtimes.txt'],1);

for i = 1:n,
    stim_index_all(i,:) = find(do==i);  %each row shows all locations of stim (on/off combo) # n
end

for j = 1:n,
    current_stim = get(saveScript,j);
    parameters = getparameters(current_stim);
    on_offduration(j,:) = [parameters.bgpause(1,1) parameters.bgpause(1,2)];
    on_offreps(j,1) = parameters.repeat;
end

for k = 1:length(mti2_),
    light_on(k,:) = mti2_{1,k}.frameTimes;
    index_stimclass = mti2_{1,k}.stimid;
    light_off(k,:) = light_on(k,:)+on_offduration(index_stimclass,1);
end
stimulus_triggers = cat(3,light_on,light_off);  

spikedata = get_data(mycell,[(light_on(1,1)-c_pre_window) (light_off(end,end)+c_post_window)]);

all_spike_cycles = cell(on_offreps(1,1),stimset_reps,n);
for p = 1:n,
    for pp = 1:stimset_reps,
        loc = stim_index_all(p,pp);
            for b = 2:on_offreps(1,1)+1,
                for c = 1:length(spikedata),
                    if b < on_offreps(1,1)+1,
                            if spikedata(c)>=(light_on(loc,b-1)-c_pre_window) && spikedata(c)<=(light_on(loc,b)+c_post_window), %c and r windows: cycle and rep.
                                all_spike_cycles{b-1,pp,p}(end+1) = spikedata(c)-light_on(loc,b-1);
                            end
                    else if b == on_offreps(1,1)+1,
                            if spikedata(c)>=(light_on(loc,b-1)-c_pre_window) && spikedata(c)<=(light_off(loc,b-1)+on_offduration(p,2)+c_post_window),
                                all_spike_cycles{b-1,pp,p}(end+1) = spikedata(c)-light_on(loc,b-1);
                            end
                    end
                    end
                end
            end
    end
end
plot_it = 0;
%%
%if plot_it == 1,    %**this line added 12/3/14; see line 224 -- further
%updated 12/30/14 - see line 156, 120, etc.  --REMOVE THESE ADDITIONS AT LATER DATE--
if ~(exist('blink_raw_figures','dir')),
    mkdir('blink_raw_figures');
end
cd 'blink_raw_figures';
        
if exist('all_spike_cycles','var')==1,
    bins = [0.005;0.005;0.01;0.02];  %bin size designated for each on duration
    bins_adjusted = [];
    for ii = 1:length(bins),
        bins_adjusted = [bins_adjusted;repmat(bins(ii,1),length(unique(on_offduration(:,2))),1)];  %fill out bin size list for all on/off combinations
    end
    cycle_FR_inst = cell(n,(stimset_reps*on_offreps(1,1)));
    FR_grand_avg = cell(n,1);
    FR_cycle_lineavg = cell(n,on_offreps(1,1));
    Xn_array = cell(n,1);
    all_cycles_response = zeros(n,stimset_reps,on_offreps(1,1));
    for u = 1:n,    %generate figures for each on/off stim. combination
        %if plot_it == 1,        %update 12/30/14*****
        f = figure;
        %else                    %update 12/30/14*****
        %end
        % bins = 0.02;  %adjusted above to allow variable bin size
        % (date of change: 8/20/14) for stims. .02,.04,.1,.5
        xLimit = (on_offduration(u,1)+on_offduration(u,2)); %cycle limits in which to fit PSTH bins
        newXLimit = (floor(xLimit*(1/bins_adjusted(u,1))))/(1/bins_adjusted(u,1)); % adjust to ensure integer # of bins at chosen size fits interval 
        edges = [0:bins_adjusted(u,1):newXLimit];
        psth = zeros(length(edges),1); %initialize new PSTH at zero counts
        dt = 0.0001;    %time bin resolution for rate measurement
        t_vector = 0:dt:newXLimit;  %vector of rate time bin edges
        gauss_width = 200;  % width of gaussian smoothing window in units of dt
        %if plot_it == 1,        %update 12/30/14*****
        ax_upper = subplot(2,1,1);  %open plot for stacked raster
        %else                    %update 12/30/14*****
        %end
        count = 0;  % initialize count of raster stacks
        stim_cycle_response = 0;
        for v = 1:stimset_reps,
            locA = stim_index_all(u,v);
            for z = 1:on_offreps(1,1),  %generate stacked raster and psth
                count = count + 1;
                spikes = all_spike_cycles{z,v,u}(:);
                psth = psth + reshape(histc(spikes,edges),length(psth),1);
                tick_ht = 1.0;
                %if plot_it == 1,        %update 12/30/14***** 
                hold on;
                for ii = 1:length(spikes),
                    raster = line([spikes(ii) spikes(ii)],[(tick_ht*count)-tick_ht,(tick_ht*count)]);
                    set(raster,'Color',[0 0 0],'LineWidth',1.5);
                end
                %else                    %update 12/30/14*****
                %end
                %compute inst. FR for this train
                Yn = zeros(length(t_vector)-1,1);
                Xn = [];
                for qq = 2:length(t_vector),
                    Xn(end+1) = mean([t_vector(qq) t_vector(qq-1)]);
                end
                if ~isempty(spikes),
                    all_cycles_response(u,v,z) = 1;
                    stim_cycle_response = stim_cycle_response + 1; %monitors count of cycles w/ response across stim. category
                    for q = 1:length(spikes),
                        spikes = reshape(spikes,1,length(spikes));
                        range = length(find(t_vector<spikes(1,q))); %* Added 9/30/14 to accommodate scenario where length(find(t_vector<spikes))=1
                        if range <= 1,                              %*
                            Yn(length(find(t_vector<spikes(1,q))),1) = 1;   %*
                        else
                            Yn((length(find(t_vector<spikes(1,q)))-1),1) = 1;
                        end
                    end
                    g_kernel = gausswin(gauss_width,4);  %changed due to
                    %update;  works with 2010a with working signal
                    %processing toolbox
                    %n = (-(gauss_width-1)/2):((gauss_width-1)/2);        %added 11/13
                    %alpha = 4;                                          %added 11/13
                    %g_kernel = exp((-1/2) * (alpha * n/(gauss_width/2)) .^ 2)';   %added 11/13
                    norm_g_kernel = g_kernel./sum(g_kernel);
                    cycle_FR_inst{u,count} = (1/dt).*conv(Yn,norm_g_kernel,'same');  %cycle FR vector
                else 
                    cycle_FR_inst{u,count} = zeros(length(t_vector)-1,1);
                end
            end
        end
        cycle_response_rate(u,1) = (round(stim_cycle_response/count*100))/100; %percentage of cycles w/ evoked response from stim. category
        %if plot_it == 1,        %update, 12/30/14*****
        Xn_array{u,1} = Xn;
        %create light-on blocks and format plots
        on_trig_line = line([0 0],[-0.1 (stimset_reps*on_offreps(1,1)*tick_ht)]);
        set(on_trig_line,'LineStyle','--');
        a_x = [0;light_off(locA,1)-light_on(locA,1)]; 
        a_y = [(count*tick_ht+(0.5*tick_ht));(count*tick_ht+(0.5*tick_ht))];
        area_h = area(a_x,a_y,-(0.5*tick_ht),'LineStyle','none');
        set(area_h,'FaceColor','b');
        alpha(.1);
        box off;
        ylim([-0.1 (tick_ht*(stimset_reps*on_offreps(1,1)+1))]);
        set(gca,'YTick',[]);
        ylabel('Cycle');
        xlabel('Time (seconds)');
        if c_pre_window > 0 || c_post_window > 0,
            xlim([-(c_pre_window) (((on_offduration(u,1)+on_offduration(u,2))*on_offreps(1,1)))+(c_post_window)]);
        else
            xlim([(0-(0.1*newXLimit)) newXLimit+(0.1*newXLimit)]);
        end
        title(['Stim. on/off: ',num2str(on_offduration(u,1)),'/',num2str(on_offduration(u,2))]);
        hold off;
        ax_lower = subplot(2,1,2);  %open subplot and make PSTH
        hold on;
        bar((edges+(bins_adjusted(u,1)/2)),psth);
        linkaxes([ax_upper ax_lower],'x');  %can try matchaxes (vhtools) instead
        xlabel('Time (seconds)');
        ylabel('Count');
        stat_str = {['Cycle rate: ', num2str(cycle_response_rate(u,1))],'Med. jitter: ','Med. latency: '};
        annotation('textbox',[0.7 0.34 0.2 0.11],'String',stat_str);
        hold off;
        if plot_it == 0,
            saveas(gcf,[pwd filesep mycellname '_stack60_' num2str(on_offduration(u,1)) '_' num2str(on_offduration(u,2)) '.fig']);
            close(f);
        else
        end
        %else                   %update, 12/30/14*****
        %end
    end
    for d = 1:n,
        fn = figure;
        for e = 1:length(cycle_FR_inst),
            if e == 1,
                FR_vector = zeros(length(cycle_FR_inst{d,e}(:)),1);
            else
                FR_vector = FR_vector + cycle_FR_inst{d,e}(:);
            end
        end
        for ee = 1:on_offreps(1,1),
            for g = 1:stimset_reps,
                if g == 1,
                    FR_cyc_vector = zeros(length(cycle_FR_inst{d,ee}(:)),1);
                end
                FR_cyc_vector = FR_cyc_vector + cycle_FR_inst{d,((g-1)*on_offreps(1,1)+ee)}(:);
            end
            FR_cycle_lineavg{d,ee} = FR_cyc_vector./stimset_reps;
            cycle_meanFR(d,ee) = mean2(FR_cycle_lineavg{d,ee}(:));
        end
        FR_grand_avg{d,1} = FR_vector./(length(cycle_FR_inst));
        plot(Xn_array{d,1}(:),FR_grand_avg{d,1}(:),'LineWidth',1.5);
        hold on;
        xlabel('Time (sec.)');
        ylabel('Spikes/sec.)');
        title(['Cell-averaged cycle instantaneous firing rate for on/off: ',num2str(on_offduration(d,1)),'/',num2str(on_offduration(d,2))]);
        hold off;
        if plot_it == 0,
            saveas(gcf,[pwd filesep mycellname '_stack60_instFR_' num2str(on_offduration(d,1)) '_' num2str(on_offduration(d,2)) '.fig']);
            close(fn);
        else 
        end
    end
end
%else   %**these lines added for line 80 'if' statement 12/3/14 to bypass plot generation when not needed****WARNING****also bypasses computations essential to inst. rate
%end    %**
        
end

