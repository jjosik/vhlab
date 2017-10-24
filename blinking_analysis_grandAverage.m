function [ f,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 ] = blinking_analysis_grandAverage( gen_plot_it,celltypes,gen_pre_window,gen_post_window )
%UNTITLED Summary of this function goes here
%   INPUTS
%           CELLTYPES - set to contain a string indicating which cell
%           categories to analyze according to the unit_quality.txt file;
%           options include 'g' and 'e' for "good" and "excellent" single units
%           and 'mu' for multi unit cells
%**IMPORTANT: BEFORE RUNNING (1/12/2015)**
%NEED TO EXPLICITLY DECLARE/CHANGE: 1.) celltypes ('mu' or 'g', here and in analyzechr2..), 2.) c_pre_window,
%3.) c_post_window (both should be 0.5 for on_off_ttest.m), 4.) plot_it = 0,
%5.) save_it = 1, 6.) do_fourier = 1, 7.) do_alt_fourier = 0, 8.) use_windowing = 0,
%9.) error_on = 0 or 1 (optional for fourier plots), 10.) fit_it = 0 or 1 (optional for 
%cycle success plots), can confirm rep_segment = 1, pre_window = 0, and
%post_window = 0, though not necessary if checked in code (set in this function, line 42).

%global plot_it;
%global pre_window;
%global post_window;
%pre_window = gen_pre_window;
%post_window = gen_post_window;
%global gen_celltypes
plot_it = 0;
c_pre_window = 0.5;
c_post_window = 0.5;
%for t test (*note: might consider linking these values to diode test somehow)
on_windowshift = 0.02;
off_windowshift = 0.02;

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
%cellload_errors = {};
current_count = 0;

for i=1:length(date_directories),
    rep_segment = 1; pre_window = 0; post_window = 0;
    %try
        [mycell,mycell_locations,on_offduration,all_spike_cycles]=analyzechr2blinkingresponses(date_directories{i},rep_segment,plot_it,pre_window,post_window);
    %catch exception
    %    cellload_errors = [cellload_errors;{date_directories{i}}];
    %    continue
    %end
    for j=1:size(mycell_locations,1),
        grand_cells{end+1,1} = date_directories{i};
        current_count = current_count + 1;
    end
    column_start = (current_count - size(mycell_locations,1))+1;
    for j2=1:size(mycell_locations,1),
        grand_cells{column_start+(j2-1),2} = mycell_locations{j2,1};
        grand_cells{column_start+(j2-1),3} = mycell_locations{j2,2};
    end
end
%lines in this block were temporarily edited out 11/13/14 to run analysis
%on E and F only -- reinstate with proper conditionals - REVERSED 11/26
G_spectral_peaks(size(grand_cells,1),1) = struct('primary',[],'secondary',[],'tertiary',[],'quarternary',[]);
G_spectral_amp(size(grand_cells,1),1) = struct('amplitudes',[]);
G_spectral_amp_freq(size(grand_cells,1),1) = struct('freqs',[]);
G_spectral_power(size(grand_cells,1),1) = struct('power',[]);
G_spectral_power_freq(size(grand_cells,1),1) = struct('p_freqs',[]);
G_meanband(size(grand_cells,1),1) = struct('gamma',[],'trans_gamma',[],'sub_gamma',[]);
G_cycle_success = cell(size(grand_cells,1),1);
for k=1:size(grand_cells,1),
    cd(grand_cells{k,1});
    cellexp_dir = pwd;
    ds_g = dirstruct(cellexp_dir);
    database_g = getexperimentfile(ds_g);
    mydata_g = load(database_g,'-mat');
    current_mycell = char(grand_cells{k,3});
    celldata_g = load2celllist(getexperimentfile(ds_g),strtrim(current_mycell),'-mat');
    %[A,I_A] = findassociate(celldata_g{1},'blink_cycle_response','','');    %blocked during histogram plots
    %[B,I_B] = findassociate(celldata_g{1},'blink_med_latency','','');       %blocked during histogram plots
    %[C,I_C] = findassociate(celldata_g{1},'blink_cycle_meanFR','','');      %blocked during histogram plots
    %[D,I_D] = findassociate(celldata_g{1},'blink_med_jitter','','');        %blocked during histogram plots
    [E,I_E] = findassociate(celldata_g{1},'stim_med_latency','','');   %blocked during grand average plot
    [F,I_F] = findassociate(celldata_g{1},'stim_med_jitter','','');     %blocked during grand average plto
    [H,I_H] = findassociate(celldata_g{1},'test_spectral_amplitudes','','');
    [J,I_J] = findassociate(celldata_g{1},'test_spectral_power','','');
    [M,I_M] = findassociate(celldata_g{1},'test_spectral_phase','','');
    [N,I_N] = findassociate(celldata_g{1},'test_spectral_peaks','','');
    [P,I_P] = findassociate(celldata_g{1},'test_gamma_strength','','');
    [Q,I_Q] = findassociate(celldata_g{1},'cycle_on_success_rate','','');
    %G_matched_cycle_response(:,:,k) = A.data;  %blocked during hist plots                             
    %G_cycle_med_lat(:,:,k) = B.data;        %blocked during hist plots
    %G_cycle_meanFR(:,:,k) = C.data;         %blocked during hist plots
    %G_delta_MJ(:,:,k) = D.data;             %blocked during hist plots
    G_stim_lat(:,k) = E.data;       %blocked during grand average plot
    G_stim_jitter(:,k) = F.data;    %blocked during grand average plot
    G_spectral_peaks(k).primary = N.data.primary;
    G_spectral_peaks(k).secondary = N.data.secondary;
    G_spectral_peaks(k).tertiary = N.data.tertiary;
    G_spectral_peaks(k).quarternary = N.data.quarternary;
    G_spectral_amp(k).amplitudes = H.data.amplitudes;
    G_spectral_amp_freq(k).freqs = H.data.frequencies;
    G_spectral_power(k).power = J.data.log_power;
    G_spectral_power_freq(k).p_freqs = J.data.frequencies;
    G_meanband(k).gamma = P.data.g_power;
    G_meanband(k).trans_gamma = P.data.trans_g_power;
    G_meanband(k).sub_gamma = P.data.sub_g_power;
    G_cycle_success{k,1} = Q.data;
    [sub_P,m_ON_cycle_FR_inst,m_OFF_cycle_FR_inst] = on_off_ttest(on_windowshift,off_windowshift,all_spike_cycles,on_offduration,c_pre_window);
    keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end
keep_index = find(keep_grandcell==1);
new_G_spectral_pri = cell(length(keep_index),1);
new_G_spectral_sec = cell(length(keep_index),1);
new_G_spectral_ter = cell(length(keep_index),1);
new_G_spectral_quart = cell(length(keep_index),1);
new_G_spectral_amp = cell(length(keep_index),1);
new_G_spectral_amp_freq = cell(length(keep_index),1);
new_G_spectral_power = cell(length(keep_index),1);
new_G_spectral_power_freq = cell(length(keep_index),1);
new_G_meanband_gamma = cell(length(keep_index),1);
new_G_meanband_trans_gamma = cell(length(keep_index),1);
new_G_meanband_sub_gamma = cell(length(keep_index),1);
for k2 = 1:length(keep_index),
    new_G_spectral_pri{k2,1} = G_spectral_peaks(keep_index(k2,1),1).primary;
    new_G_spectral_sec{k2,1} = G_spectral_peaks(keep_index(k2,1),1).secondary;
    new_G_spectral_ter{k2,1} = G_spectral_peaks(keep_index(k2,1),1).tertiary;
    new_G_spectral_quart{k2,1} = G_spectral_peaks(keep_index(k2,1),1).quarternary;
    new_G_spectral_amp{k2,1} = G_spectral_amp(keep_index(k2,1),1).amplitudes;
    new_G_spectral_amp_freq{k2,1} = G_spectral_amp_freq(keep_index(k2,1),1).freqs;
    new_G_spectral_power{k2,1} = G_spectral_power(keep_index(k2,1),1).power;
    new_G_spectral_power_freq{k2,1} = G_spectral_power_freq(keep_index(k2,1),1).p_freqs;
    new_G_meanband_gamma{k2,1} = G_meanband(keep_index(k2,1),1).gamma;
    new_G_meanband_trans_gamma{k2,1} = G_meanband(keep_index(k2,1),1).trans_gamma;
    new_G_meanband_sub_gamma{k2,1} = G_meanband(keep_index(k2,1),1).sub_gamma;
end

%%Collected histograms of cell med. latencies and jitters
%entire section from lines 145 to 185 blocked during grand average plot
Dist_med_lat_low = G_stim_lat(1,(~isnan((G_stim_lat(1,:)))));
Dist_med_jitter_low = G_stim_jitter(1,(~isnan((G_stim_jitter(1,:)))));
Dist_med_lat_high = G_stim_lat(16,(~isnan((G_stim_lat(16,:)))));
Dist_med_jitter_high = G_stim_jitter(16,(~isnan((G_stim_lat(16,:)))));

bin_edges = [0:0.001:0.08];  %in end used 0.001 bins for jitter and 0.002 for lat.;
N1 = histc(Dist_med_lat_low,bin_edges);
N2 = histc(Dist_med_jitter_low,bin_edges);
N3 = histc(Dist_med_lat_high,bin_edges);
N4 = histc(Dist_med_jitter_high,bin_edges);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
N1 = N1(1:end-1);
N2 = N2(1:end-1);
N3 = N3(1:end-1);
N4 = N4(1:end-1);

f1 = figure;
bar(bin_centers,N1,'b');
hold on;
bar(bin_centers,N3,'r');
alpha(0.5);
xlabel('Seconds');
ylabel('Count');
legend('20ms/20ms','500ms/500ms');
title('Median latency distribution');
ylim([0 (max(max([N1;N3]))+3)]);
xlim([0 0.08]);
hold off;

f2 = figure;
bar(bin_centers,N2,'b');
hold on;
bar(bin_centers,N4,'r');
alpha(0.5);
xlabel('Seconds');
ylabel('Count');
legend('20ms/20ms','500ms/500ms');
title('Median jitter distribution');
ylim([0 (max(max([N2;N4]))+3)]);
xlim([0 0.04]);
hold off;

%%Graph of cycle "on" success rate averaged over all cells
G_success_mat = cell2mat(G_cycle_success);
G_mean_success = sum(G_success_mat,1)./size(G_success_mat,1);
G_success_error = std(G_success_mat)./sqrt(size(G_success_mat,1));
X_bar = [1:length(unique(on_offduration(:,1)))];
f11 = figure;
b = bar(X_bar,G_mean_success);
hold on;
eb = errorbar(G_mean_success,G_success_error,'k.');
set(gca,'XTick',[1,2,3,4]);
set(gca,'XTickLabel',{'0.02','0.04','0.1','0.5'});
xlabel('On duration (seconds)');
ylabel('Success rate');
ylim([0 1]);
N_string = {['N = ',num2str(length(G_success_mat))]};
AN = annotation('textbox',[0.2 0.75 0.1 0.1],'String',N_string,'EdgeColor','None');
set(AN,'FontSize',15);
hold off;

%%Grand Average 4x4 plot (stim. cycle med. latency, FR, jitter, success
%%rate)
% following section blocked during hist plots
%GAvg_matched_cycle_response = nanmean(G_matched_cycle_response,3);
%e1 = nanstd(G_matched_cycle_response,1,3)/sqrt(size(G_matched_cycle_response,3));
%GAvg_cycle_med_lat = nanmean(G_cycle_med_lat,3);
%e2 = nanstd(G_cycle_med_lat,1,3)/sqrt(size(G_cycle_med_lat,3));
%GAvg_cycle_meanFR = nanmean(G_cycle_meanFR,3);
%e3 = nanstd(G_cycle_meanFR,1,3)/sqrt(size(G_cycle_meanFR,3));
%GAvg_delta_MJ = nanmean(G_delta_MJ,3);
%e4 = nanstd(G_delta_MJ,1,3)/sqrt(size(G_delta_MJ,3));
%
%f = figure;
%AN_col1 = annotation('textbox',[0.15 0.95 0.08 0.04],'String',collabels{1},'EdgeColor','None');
%AN_col2 = annotation('textbox',[0.37 0.95 0.08 0.04],'String',collabels{2},'EdgeColor','None');
%AN_col3 = annotation('textbox',[0.58 0.95 0.08 0.04],'String',collabels{3},'EdgeColor','None');
%N_col4 = annotation('textbox',[0.8 0.95 0.08 0.04],'String',collabels{4},'EdgeColor','None');
%set(AN_col1,'FontSize',20,'FontWeight','bold');
%set(AN_col2,'FontSize',20,'FontWeight','bold');
%set(AN_col3,'FontSize',20,'FontWeight','bold');
%set(AN_col4,'FontSize',20,'FontWeight','bold');
%on_label = unique(on_offduration(:,1));
%
%plots grand-averaged delta median sequence-matched latency
%for i1 = 1:length(on_label),
%    subplot(4,4,i1);
%    prop_name_M = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set_M = {'k','b','y','m'};
%    prop_values_M = repmat(color_set_M',length(unique(on_offduration(:,1))),1);
%    cycle_num1 = 1:size(GAvg_cycle_med_lat,2);
%    for j1 = 1:length(unique(on_offduration(:,2))),
%        m = (i1-1)*length(unique(on_offduration(:,1)))+j1;
%        hh(m) = errorbar(cycle_num1,GAvg_cycle_med_lat(m,:),e2(m,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(hh(m),prop_name_M(m),prop_values_M(m));
%        hold on;
%    end
%    hLeg(i1) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
%    if i1 < length(on_label),
%        set(hLeg(i1),'Visible','Off');
%    else
%        set(hLeg(i1),'Location','NortheastOutside');
%    end
%    xlim([0 10]);
%    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
%    ylabel('Delta median latency','FontSize',12,'FontWeight','bold');
%    hold off;
%end
%plots grand-averaged mean sequenced-matched cycle FR
%for i2 = 1:length(unique(on_offduration(:,1))),
%    subplot(4,4,4+i2);
%    prop_name = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set = {'k','b','y','m'};
%    prop_values = repmat(color_set',length(unique(on_offduration(:,1))),1);
%    cycle_num2 = 1:size(GAvg_cycle_meanFR,2);
%    for j2 = 1:length(unique(on_offduration(:,2))),
%        n = (i2-1)*length(unique(on_offduration(:,1)))+j2;
%        h(n) = errorbar(cycle_num2,GAvg_cycle_meanFR(n,:),e3(n,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(h(n),prop_name(n),prop_values(n));
%    hold on;    
%    end
%    hLeg1(i2) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
%    if i2 < length(on_label),
%        set(hLeg1(i2),'Visible','Off');
%    else
%        set(hLeg1(i2),'Location','NortheastOutside');
%    end
%    xlim([0 10]);
%    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
%    ylabel('Normalized mean FR','FontSize',12,'FontWeight','bold');
%    hold off;
%end

%plots sequence-matched median absolute deviation (i.e. median jitter)
%for i3 = 1:length(on_label),
%    subplot(4,4,8+i3);
%    prop_name_MJ = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set_MJ = {'k','b','y','m'};
%    prop_values_MJ = repmat(color_set_MJ',length(unique(on_offduration(:,1))),1);
%    cycle_num3 = 1:size(GAvg_delta_MJ,2);
%    for j3 = 1:length(unique(on_offduration(:,2))),
%        p = (i3-1)*length(unique(on_offduration(:,1)))+j3;
%        hhh(p) = errorbar(cycle_num3,GAvg_delta_MJ(p,:),e4(p,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(hhh(p),prop_name_MJ(p),prop_values_MJ(p));
%        hold on;
%    end
%    hLeg2(i3) = legend('0.02 off','0.04 off','0.1 Off','0.5 Off');
%    if i3 < length(on_label),
%        set(hLeg2(i3),'Visible','Off');
%    else
%        set(hLeg2(i3),'Location','NortheastOutside');
%    end
%    xlim([0 10]);
%    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
%    ylabel('Delta median jitter','FontSize',12,'FontWeight','bold');
%    hold off;
%end
%
%plots grand-averaged sequence-matched cycle response rate
%for i4 = 1:length(on_label),
%    subplot(4,4,12+i4);
%    prop_name_R = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set_R = {'k','b','y','m'};
%    prop_values_R = repmat(color_set_R',length(unique(on_offduration(:,1))),1);
%    cycle_num4 = 1:size(GAvg_matched_cycle_response,2);
%    for j4 = 1:length(unique(on_offduration(:,2))),
%        q = (i4-1)*length(unique(on_offduration(:,1)))+j4;
%        hhhh(q) = errorbar(cycle_num4,GAvg_matched_cycle_response(q,:),e1(q,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(hhhh(q),prop_name_R(q),prop_values_R(q));
%        hold on;
%    end
%    hLeg3(i4) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
%    if i4 < length(on_label),
%        set(hLeg3(i4),'Visible','Off');
%    else
%        set(hLeg3(i4),'Location','NortheastOutside');
%    end
%    xlim([0 10]);
%    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
%    ylabel('Cycle Response rate','FontSize',12,'FontWeight','bold');
%    hold off;
%end
%if plot_it == 0,
%    saveas(gcf,[test_subFold(s) filesep celltypes '_GAvg_seq_analysis' '.fig']);
%    close(f);
%else
%end
%%
%%Collected histogram of cell Fourier responses (peak frequencies)
accum_peaks_primary = [];
accum_peaks_secondary = [];
accum_peaks_tertiary = [];
accum_peaks_quarternary = [];
selected_stim = 16;  %this analysis to be based on 0.5 on, 0.5 off stim.
for k3 = 1:length(keep_index),
    accum_peaks_primary(end+1) = G_spectral_peaks(k3).primary{selected_stim,1}(:);
    accum_peaks_secondary(end+1) = G_spectral_peaks(k3).secondary{selected_stim,1}(:);
    accum_peaks_tertiary(end+1) = G_spectral_peaks(k3).tertiary{selected_stim,1}(:);
    accum_peaks_quarternary(end+1) = G_spectral_peaks(k3).quarternary{selected_stim,1}(:);
end
f_bin_edges = [0:5:80];
Nf1 = histc(accum_peaks_primary,f_bin_edges);
Nf2 = histc(accum_peaks_secondary,f_bin_edges);
Nf3 = histc(accum_peaks_tertiary,f_bin_edges);
Nf4 = histc(accum_peaks_quarternary,f_bin_edges);
f_bin_centers = (f_bin_edges(1:end-1)+f_bin_edges(2:end))/2;
Nf1 = Nf1(1:end-1);
Nf2 = Nf2(1:end-1);
Nf3 = Nf3(1:end-1);
Nf4 = Nf4(1:end-1);

f3 = figure;
bar(f_bin_centers,Nf1,'r');
hold on;
xlabel('Frequency (Hz)');
ylabel('Number of sites');
legend('primary frequency','2nd peak','3rd peak','4th peak');
xlim([0 80]);
ylim([0 (max(max([Nf1;Nf2;Nf3;Nf4]))+3)]);
hold off;

f4 = figure;
bar(f_bin_centers,Nf2,'m');
hold on;
xlabel('Frequency (Hz)');
ylabel('Number of sites');
legend('secondary frequency','2nd peak','3rd peak','4th peak');
xlim([0 80]);
ylim([0 (max(max([Nf1;Nf2;Nf3;Nf4]))+3)]);
hold off;

f5 = figure;
bar(f_bin_centers,Nf3,'b');
hold on;
xlabel('Frequency (Hz)');
ylabel('Number of sites');
legend('tertiary frequency','2nd peak','3rd peak','4th peak');
xlim([0 80]);
ylim([0 (max(max([Nf1;Nf2;Nf3;Nf4]))+3)]);
hold off;

f6 = figure;
bar(f_bin_centers,Nf4,'c');
hold on;
xlabel('Frequency (Hz)');
ylabel('Number of sites');
legend('quarternary frequency','2nd peak','3rd peak','4th peak');
xlim([0 80]);
ylim([0 (max(max([Nf1;Nf2;Nf3;Nf4]))+3)]);
hold off;

%%Alt. scatterplot of relative gamma power strength

for k5 = 1:length(keep_index),
    cell_gamma(k5,1) = new_G_meanband_gamma{k5,1}(selected_stim,1);
    cell_trans_gamma(k5,1) = new_G_meanband_trans_gamma{k5,1}(selected_stim,1);
    cell_sub_gamma(k5,1) = new_G_meanband_sub_gamma{k5,1}(selected_stim,1);
end
f9 = figure;
h_s = scatter(cell_gamma,cell_trans_gamma,'k');
hold on;
xlabel('Site mean gamma power');
ylabel('Site mean power >70 Hz');
xy_M = max(cat(1,cell_gamma,cell_trans_gamma))+(0.1*(max(cat(1,cell_gamma,cell_trans_gamma))));
xlim([0 xy_M]);
ylim([0 xy_M]);
unity_line = line([0 xy_M],[0 xy_M],'LineStyle','--');
hold off;

f10 = figure;
h_ss = scatter(cell_gamma,cell_sub_gamma,'k');
hold on;
xlabel('Site mean gamma power');
ylabel('Site mean power 10-25 Hz');
xya_M = max(cat(1,cell_gamma,cell_sub_gamma))+(0.1*(max(cat(1,cell_gamma,cell_sub_gamma))));
xlim([0 xya_M]);
ylim([0 xya_M]);
unity_line_ = line([0 xya_M],[0 xya_M],'LineStyle','--');
hold off;


%%Grand average Fourier plot
%grand average amplitude 
extract_amps = cell(length(keep_index),1);
amp_extract_matched = cell(length(keep_index),1);
amp_f_matched = cell(length(keep_index),1);
power_extract_matched = cell(length(keep_index),1);
power_f_matched = cell(length(keep_index),1);
for k4 = 1:length(keep_index),
    amp_extract_matched{k4,1} = new_G_spectral_amp{k4,1}{selected_stim,1}(:)';
    amp_f_matched{k4,1} = new_G_spectral_amp_freq{k4,1}(:,:,selected_stim);
    power_extract_matched{k4,1} = new_G_spectral_power{k4,1}{selected_stim,1}(:)';
    power_f_matched{k4,1} = new_G_spectral_power_freq{k4,1}(:,:,selected_stim);
end

amp_array = cell2mat(amp_extract_matched);
amp_f_array = amp_f_matched{1,1}{1,1}(:);
power_array = cell2mat(power_extract_matched);
power_f_array = power_f_matched{1,1}{1,1}(:);
M_G_spectral_amp = sum(amp_array,1)./size(amp_array,1);
M_G_spectral_power = sum(power_array,1)./size(power_array,1);
G_error_amp = std(amp_array,1)./sqrt(size(amp_array,1));
G_error_power = std(power_array,1)./sqrt(size(power_array,1));

f7 = figure;
h_gf = plot(amp_f_array(:,1),M_G_spectral_amp(1,:),'k');
hold on;
set(h_gf,'LineWidth',1.0);
h_gff = plot(amp_f_array(:),(M_G_spectral_amp(:)+G_error_amp(:)),'b');
set(h_gff,'LineWidth',0.5);
h_gfff = plot(amp_f_array(:),(M_G_spectral_amp(:)-G_error_amp(:)),'b');
set(h_gfff,'LineWidth',0.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
xlim([0 500]);
ylim([0 2*max(M_G_spectral_amp)]);
hold off;

f8 = figure;
h_gfp = plot(power_f_array(:),M_G_spectral_power(:),'k');
hold on;
set(h_gfp,'LineWidth',1.0);
h_gffp = plot(power_f_array(:),(M_G_spectral_power(:)+G_error_power(:)),'b');
set(h_gffp,'LineWidth',0.5);
h_gfffp = plot(power_f_array(:),(M_G_spectral_power(:)-G_error_power(:)),'b');
set(h_gfffp,'LineWidth',0.5);
xlabel('Frequency (Hz)');
ylabel('Log power');
xlim([0 500]);
ylim([0 2*max(M_G_spectral_power)]);
hold on;
    


end

