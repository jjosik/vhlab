function [osc_score_output] = score_oscillation(bin_size,f_min,f_max,zero_pad,use_FFTwindowing)

osc_score_output = {};
for j=1:size(triggers,2),
    osc_score_output{1,j} = struct('score',[],'total_w',[],'sigma_fast',[]);
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

%for ACH
cutoff_f = 1/bin_size;
criterion_range = [log2((3*cutoff_f)/f_min),log2(cutoff_f/4)];
w = 2^(floor(max(criterion_range))+1);
total_w = 2*w;

ach_info = struct('ach_w',[],'f_min',[],'f_max',[],'bin_size',[]);
ach_info.ach_w = [-(w), w-1];
ach_info.f_min = f_min;
ach_info.f_max = f_max;
ach_info.bin_size = bin_size;
for i=1:length(date_directories),
    tic,
    if rerun_celllist == 1,
        rep_segment = 0; pre_window = 0; post_window = 0;
    else
        rep_segment = 6; pre_window = 0; post_window = 0;
    end
    [mycell,mycell_locations,on_offduration,all_spike_cycles]=analyzechr2blinkingresponses(date_directories{i},rep_segment,plot_it,pre_window,post_window,ach_info);
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
    [A,I_A] = findassociate(celldata_g{1},'blink_ACH_output','','');
    [B,I_B] = findassociate(celldata_g{1},'blink_command','','');
    on_offduration = B.data;
    for j=1:size(A.data,2),
        coll_stim_lags{4,j,k} = A.data{4,j}.stim_lags;
        coll_stim_corr{4,j,k} = A.data{4,j}.mean_stim_corr;
    end
    %sub_P = alt_on_off_ttest(G_cycleFR);
    %keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

%smoothing parameters
fast_param = [2,(134/(1.5*f_max))];
sigma_fast = (min(fast_param))*(cutoff_f/1000);
sigma_slow = 2*(134/(1.5*f_min))*(cutoff_f/1000); 

for m = 1:size(grand_cells,1),
    for n = 1:size(coll_stim_corr,2),
        g_kernel_tf = -(round(3*sigma_fast)):round(3*sigma_fast);
        g_kernel_ts = -(round(3*sigma_slow)):round(3*sigma_slow);
        fast_g_kernel = (1/(sqrt(2*pi)*sigma_fast))*(exp(-((g_kernel_tf*g_kernel_tf)/(2*sigma_fast*sigma_fast))));
        slow_g_kernel = (1/(sqrt(2*pi)*sigma_slow))*(exp(-((g_kernel_ts*g_kernel_ts)/(2*sigma_slow*sigma_slow))));
        if zero_pad == 1,
            fast_pad = zeros(length(g_kernel_tf),1);
            fast_padded_stim_corr = [fast_pad;reshape(coll_stim_corr{4,n,m}(:),length(coll_stim_corr{4,n,m}(:)),1);fast_pad];
            slow_pad = zeros(length(g_kernel_ts),1);
            slow_padded_stim_corr = [slow_pad;reshape(coll_stim_corr{4,n,m}(:),length(coll_stim_corr{4,n,m}(:)),1);slow_pad];
        else
            %use mirroring
            f_lead_mirror = flipud(reshape(coll_stim_corr{4,n,m}(1:length(g_kernel_tf)),length(g_kernel_tf),1));
            f_tail_mirror = flipud(reshape(coll_stim_corr{4,n,m}((end-length(g_kernel_tf)):end),length(g_kernel_tf),1));
            fast_padded_stim_corr = [f_lead_mirror;reshape(coll_stim_corr{4,n,m}(:),length(coll_stim_corr{4,n,m}(:)),1);f_tail_mirror];
            s_lead_mirror = flipud(reshape(coll_stim_corr{4,n,m}(1:length(g_kernel_ts)),length(g_kernel_ts),1));
            s_tail_mirror = flipud(reshape(coll_stim_corr{4,n,m}(1:length(g_kernel_ts)),length(g_kernel_ts),1));
            slow_padded_stim_corr = [s_lead_mirror;reshape(coll_stim_corr{4,n,m}(:),length(coll_stim_corr{4,n,m}(:)),1);s_tail_mirror];
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
            nopeak_ACH = nopeak_ACH.*blackman(length(nopeak_ACH);
        else
        end
        [fc,freqs] = fouriercoeffs(nopeak_ACH,bin_size);
        ACH_power_spect = reshape((fc.*conj(fc)),length(fc),1);
        ACH_power_spect = ACH_power_spect(1:(1/bin_size)/2);
        band_peak_magnitude = max(ACH_power_spect(f_min:f_max));
        mean_spect_magnitude = mean(ACH_power_spect);
        osc_score_output.score(n,m) = band_peak_magnitude/mean_spect_magnitude;
    end
end
        
        
        
        
        
        
        
