function [f1] = grandAverage_autocorr(plot_it,celltypes,cycle_window,save_cellplots,normalize,rerun_celllist,score_grand_ACH,varargin);
% GRANDAVERAGE_AUTOCORR - Collectively averages cycle-averaged spike train
% autocorrelation for all "blinking stim" experiment cells specified under date folders as selected in the
% prompt box.  Home directory should be set to directory containing all
% experiment dates ('/Users/vhlab/Desktop/test' in the case of Colonel).  
%%%
%%%INPUTS
%%%PLOT_IT - Set to '1' if generation and saving of grand average
%%%autocorrelograms is desired, '0' otherwise.
%%%CELLTYPES - Allows cells analyzed to be restricted to categories
%%%accessed under analyzechr2blinkingresponses.m;  typically 'g' or 'mu'
%%%indicates 'good/single unit' or 'multiunit', respectively.
%%%CYCLE_WINDOW - Set to '1' if a subset of the analysis cycles is desired;
%%%e.g., if you want to analyze the last 400 ms of a 500 ms cycle window,
%%%cycle_window must be equal to '1'.  Note, the timing of the
%%%corresponding windowing must be set in the analyze_spiketrain_autocorr.m
%%%function.
%%%SAVE_CELLPLOTS - Set to '1' if saving of cycle-averaged
%%%autocorrelograms by individual cell is desired.
%%%NORMALIZE - Set to '1' if grand average correlogram should be normalized
%%%(i.e. central peak amplitude is 1).
%%%RERUN_CELLLIST - When set to '1', will recompile the array of
%%%'grand_cell' names to rerun grand average code without reanalyzing
%%%individual cells.  Must be set to '0' to analyze cells for the first
%%%time or to reanalyze/change the current value of the associate.
%*********************************************************************
%%%OUTPUTS
%%%F1 - autocorrelogram figure for the current stim. combination
%***********************************************************************
%***********************************************************************

if nargin > 7,
    modulation_band = varargin{1};
    f_min = modulation_band(1,1);
    f_max = modulation_band(1,2);
else
    %default to low-gamma
    f_min = 25;
    f_max = 50;
end

addpath('/Desktop/test');
pathFolder = uigetdir();

d = dir(pathFolder);
isub = [d(:).isdir];
test_subFold = {d(isub).name}';
test_subFold(ismember(test_subFold,{'.','..'})) = [];

[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',test_subFold);

date_directories = test_subFold(s);
grand_cells = {};
current_count = 0;

if score_grand_ACH == 1,
    %currently only determining O-score for 500 ms-ON stims. (hence index
    %from 1 to 4)
    grand_osc_score = struct();
    for j=1:4,
        grand_osc_score{1,j} = struct('score',[],'total_w',[],'sigma_fast',[],'band_peak',[]);
    end
    %for ACH
    bin_size = 0.001;
    cutoff_f = 1/bin_size;
    criterion_range = [log2((3*cutoff_f)/f_min),log2(cutoff_f/4)];
    w = 2^(floor(max(criterion_range))+1);
    total_w = 2*w;

    ach_info = struct('ach_w',[],'f_min',[],'f_max',[],'bin_size',[]);
    ach_info.ach_w = [-(w), w-1];
    ach_info.total_w = total_w;
    ach_info.f_min = f_min;
    ach_info.f_max = f_max;
    ach_info.bin_size = bin_size;
    use_FFTwindowing = 1;   %consider setting these with varargin argument
    zero_pad = 0;           %consider setting these with varargin argument
else
    ach_info = {};  %need to check that this works when not running osc. score
end

for i=1:length(date_directories),
    tic,
    if rerun_celllist == 1,
        rep_segment = 0; pre_window = 0; post_window = 0;
    else
        if score_grand_ACH == 1,
            rep_segment = 6; pre_window = 0; post_window = 0;
        else
            rep_segment = 5; pre_window = 0; post_window = 0;
        end
    end
    [~,mycell_locations,on_offduration,~]=analyzechr2blinkingresponses(date_directories{i},rep_segment,celltypes,plot_it,pre_window,post_window,ach_info);
    for j=1:size(mycell_locations,1),
        grand_cells{end+1,1} = date_directories{i};
        current_count = current_count + 1;
    end
    column_start = (current_count - size(mycell_locations,1))+1;
    for j2=1:size(mycell_locations,1),
        grand_cells{column_start+(j2-1),2} = mycell_locations{j2,1};
        grand_cells{column_start+(j2-1),3} = mycell_locations{j2,2};
    end
    toc,
end

for k=1:size(grand_cells,1),
    cd(grand_cells{k,1});
    cellexp_dir = pwd;
    ds_g = dirstruct(cellexp_dir);
    database_g = getexperimentfile(ds_g);
    mydata_g = load(database_g,'-mat');
    current_mycell = char(grand_cells{k,3});
    celldata_g = load2celllist(getexperimentfile(ds_g),strtrim(current_mycell),'-mat');
    if score_grand_ACH == 1,
        [A,I_A] = findassociate(celldata_g{1},'blink_ACH_output_Oscore','','');
        [B,I_B] = findassociate(celldata_g{1},'blink_command_Oscore','','');
    else
        [A,I_A] = findassociate(celldata_g{1},'blink_ACH_output','','');
        [B,I_B] = findassociate(celldata_g{1},'blink_command','','');
    end
    on_offduration = B.data;
    for i=1:size(A.data,1),
        for j=1:size(A.data,2),
            coll_stim_lags{i,j,k} = A.data{i,j}.stim_lags;
            coll_stim_corr{i,j,k} = A.data{i,j}.mean_stim_corr;
        end
    end
    %sub_P = alt_on_off_ttest(G_cycleFR);
    %keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

stim_corr_array = cell(size(coll_stim_corr,1),size(coll_stim_corr,2));
gmean_stim_corr = cell(size(coll_stim_corr,1),size(coll_stim_corr,2));
for i=1:size(coll_stim_corr,1),
    for j=1:size(coll_stim_corr,2),
        for k=1:size(grand_cells,1),
            stim_corr_match{1,k} = coll_stim_corr{i,j,k};
        end
        stim_corr_array{i,j} = cell2mat(stim_corr_match);
        gmean_stim_corr{i,j} = nanmean(stim_corr_array{i,j},2);
        peak_value = max(gmean_stim_corr{i,j});
        norm_gmean_stim_corr{i,j} = gmean_stim_corr{i,j}./peak_value;
    end
end

if score_grand_ACH == 1,
    %smoothing parameters
fast_param = [2,(134/(1.5*f_max))];
sigma_fast = (min(fast_param))*(cutoff_f/1000);
sigma_slow = 2*(134/(1.5*f_min))*(cutoff_f/1000);
   
    for n = 1:size(norm_gmean_stim_corr,2),
        g_kernel_tf = -(round(3*sigma_fast)):round(3*sigma_fast);
        g_kernel_ts = -(round(3*sigma_slow)):round(3*sigma_slow);
        fast_g_kernel = (1/(sqrt(2*pi)*sigma_fast))*(exp(-((g_kernel_tf*g_kernel_tf)/(2*sigma_fast*sigma_fast))));
        slow_g_kernel = (1/(sqrt(2*pi)*sigma_slow))*(exp(-((g_kernel_ts*g_kernel_ts)/(2*sigma_slow*sigma_slow))));
        if zero_pad == 1,
            fast_pad = zeros(length(g_kernel_tf),1);
            fast_padded_stim_corr = [fast_pad;reshape(norm_gmean_stim_corr{4,n}(:),length(norm_gmean_stim_corr{4,n}(:)),1);fast_pad];
            slow_pad = zeros(length(g_kernel_ts),1);
            slow_padded_stim_corr = [slow_pad;reshape(norm_gmean_stim_corr{4,n}(:),length(norm_gmean_stim_corr{4,n}(:)),1);slow_pad];
        else
            %use mirroring
            f_lead_mirror = flipud(reshape(norm_gmean_stim_corr{4,n}(1:length(g_kernel_tf)),length(g_kernel_tf),1));
            f_tail_mirror = flipud(reshape(norm_gmean_stim_corr{4,n}((end-length(g_kernel_tf)):end),length(g_kernel_tf),1));
            fast_padded_stim_corr = [f_lead_mirror;reshape(norm_gmean_stim_corr{4,n}(:),length(norm_gmean_stim_corr{4,n}(:)),1);f_tail_mirror];
            s_lead_mirror = flipud(reshape(norm_gmean_stim_corr{4,n}(1:length(g_kernel_ts)),length(g_kernel_ts),1));
            s_tail_mirror = flipud(reshape(norm_gmean_stim_corr{4,n}(1:length(g_kernel_ts)),length(g_kernel_ts),1));
            slow_padded_stim_corr = [s_lead_mirror;reshape(norm_gmean_stim_corr{4,n}(:),length(norm_gmean_stim_corr{4,n}(:)),1);s_tail_mirror];
        end
        f_smoothed_corr = conv(fast_padded_stim_corr,fast_g_kernel,'same');
        fast_smoothed_corr = f_smoothed_corr(length(g_kernel_tf)+1:end-(length(g_kernel_tf)));
        s_smoothed_corr = conv(slow_padded_stim_corr,slow_g_kernel,'same');
        slow_smoothed_corr = s_smoothed_corr(length(g_kernel_ts)+1:end-(length(g_kernel_ts)));
        %compute ACH slope at all bins
        scaling = total_w/(max(s_smoothed_corr)-min(s_smoothed_corr));
        half_slowACH = slow_smoothed_corr(1:w+1);
        slope = (half_slowACH(2:end)-half_slowACH(1:end-1))*scaling;
        left_cutoff_index = find(slope < (tan((pi*10)/180)),1,'last');
        right_cutoff_index = length(s_smoothed_corr)-left_cutoff_index;
        nopeak_ACH = fast_smoothed_corr(:);
        nopeak_ACH(left_cutoff_index:right_cutoff_index) = fast_smoothed_corr(left_cutoff_index);
        %compute FFT
        if use_FFTwindowing == 1,
            nopeak_ACH = nopeak_ACH.*blackman(length(nopeak_ACH));
        else
        end
        [fc,~] = fouriercoeffs(nopeak_ACH,bin_size);
        ACH_power_spect = reshape((fc.*conj(fc)),length(fc),1);
        ACH_power_spect = ACH_power_spect(1:(1/bin_size)/2);
        band_peak_magnitude = max(ACH_power_spect(f_min:f_max));
        mean_spect_magnitude = mean(ACH_power_spect);
        grand_osc_score.band_peak(1,n) = band_peak_magnitude;
        grand_osc_score.score(1,n) = band_peak_magnitude/mean_spect_magnitude;
    end

else
end
    
    
if plot_it == 1,
    for i=1:size(gmean_stim_corr,1),
        for j=1:size(gmean_stim_corr,2),
            f1 = figure;
            if normalize == 1,
                ach_g = bar(coll_stim_lags{i,j,1},norm_gmean_stim_corr{i,j},'k');
            else
                ach_g = bar(coll_stim_lags{i,j,1},gmean_stim_corr{i,j},'k');
            end
            hold on;
            xlabel('Lags');
            if normalize == 1,
                ylabel('Correlation');
                ylim([0 max(norm_gmean_stim_corr{i,j})+0.1*(max(norm_gmean_stim_corr{i,j}))]);
            else
                ylabel('Conditional rate');
                ylim([0 max(gmean_stim_corr{i,j})+0.1*(max(gmean_stim_corr{i,j}))]);
            end
            if i==4,
                yt = 0.9;
                ytt = 0.7;
                xt = max(coll_stim_lags{i,j,1})-(0.9*max(coll_stim_lags{i,j,1}));
                text(xt,yt,[num2str(f_min),'-',num2str(f_max),' Hz band Osc. score = ',num2str(grand_osc_score.score(1,j))]);
                text(xt,ytt,['Band peak = ', num2str(grand_osc_score.band_peak(1,j))]);
            else
            end
            saveas(gcf,['grand_average_' num2str(on_offduration{i,j,1}) '_' num2str(on_offduration{i,j,2}) '.fig']);
            hold off;
            close(f1);
        end
    end
end

    
        


    