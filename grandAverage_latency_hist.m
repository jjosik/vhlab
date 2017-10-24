function [f1,f2,f3] = grandAverage_latency_hist(plot_it,celltypes,rerun_celllist);

plot_it = 0;

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

for i=1:length(date_directories),
    tic,
    if rerun_celllist == 1,
        rep_segment = 0; pre_window = 0; post_window = 0;
    else
        rep_segment = 3; pre_window = 0; post_window = 0;
    end
    [mycell,mycell_locations,on_offduration,all_spike_cycles]=analyzechr2blinkingresponses(date_directories{i},rep_segment,plot_it,pre_window,post_window);
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
    [A,I_A] = findassociate(celldata_g{1},'new_blink_firstspike','','');
    [B,I_B] = findassociate(celldata_g{1},'new_blink_firingrate','','');
    G_stim_lat{1,1,k} = A.data.stim_med_latency;
    G_stim_jitter{1,1,k} = A.data.stim_med_jitter;
    G_cycleFR(:,:,k) = B.data;
    sub_P = alt_on_off_ttest(G_cycleFR);
    keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

keep_index = find(keep_grandcell==1);
for k2 = 1:length(keep_index),
    ttest_G_stim_lat(:,:,k2) = G_stim_lat(:,:,keep_index(k2,1));
    ttest_G_stim_jitter(:,:,k2) = G_stim_jitter(:,:,keep_index(k2,1));
end
   

%new_G_stim_lat = nanmedian(cell2mat(G_stim_lat),3);
%for i=1:size(G_stim_lat,3),
%    temp_G_stim_jitter(:,:,i) = abs(minus(new_G_stim_lat(:,:),cell2mat(G_stim_lat(:,:,i))));
%end
%new_G_stim_jitter = nanmedian(temp_G_stim_jitter,3);
new_G_stim_lat = cell2mat(ttest_G_stim_lat);
new_G_stim_jitter = cell2mat(ttest_G_stim_jitter);

Dist_med_lat_low = squeeze(new_G_stim_lat(2,2,:));
Dist_med_jitter_low = squeeze(new_G_stim_jitter(2,2,:));
Dist_med_lat_high = squeeze(new_G_stim_lat(4,4,:));
Dist_med_jitter_high = squeeze(new_G_stim_jitter(4,4,:));

bin_edges = [0:0.001:0.08];  
N1 = histc(Dist_med_lat_low,bin_edges);
N2 = histc(Dist_med_jitter_low,bin_edges);
N3 = histc(Dist_med_lat_high,bin_edges);
N4 = histc(Dist_med_jitter_high,bin_edges);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
N1 = N1(1:end-1);
N2 = N2(1:end-1);
N3 = N3(1:end-1);
N4 = N4(1:end-1);

%"actual" n for these plots (due to presence of NaNs and very long latency
%cells (~100-200 ms))
n_1 = sum(N1);
n_2 = sum(N2);
n_3 = sum(N3);
n_4 = sum(N4);

f1 = figure;
bar(bin_centers,(N1./n_1),'b');
hold on;
bar(bin_centers,(N3./n_3),'r');
alpha(0.5);
xlabel('Seconds');
ylabel('Fraction of cells');
legend(['20ms/20ms, N=',num2str(n_1)],['500ms/500ms, N=',num2str(n_3)]);
title('Median latency distribution');
ylim([0 1]);
xlim([0 0.08]);
hold off;

f2 = figure;
bar(bin_centers,(N2./n_2),'b');
hold on;
bar(bin_centers,(N4./n_4),'r');
alpha(0.5);
xlabel('Seconds');
ylabel('Fraction of cells');
legend(['20ms/20ms, N=',num2str(n_2)],['500ms/500ms, N=',num2str(n_4)]);
title('Median jitter distribution');
ylim([0 1]);
xlim([0 0.08]);
hold off;

%generate all-cell latency/jitter relationship scatterplot
f3 = figure;
scatter(new_G_stim_jitter(4,4,:),new_G_stim_lat(4,4,:));
hold on;
xlabel('Pulse 1st spike jitter (ms)');
ylabel('Pulse 1st spike latency (ms)');
title('Cell-by-cell latency/jitter relationship');
ylim([0 0.08]);
xlim([0 0.08]);
X_m = nanmedian(new_G_stim_jitter(4,4,:));
Y_m = nanmedian(new_G_stim_lat(4,4,:));
h2 = line([X_m X_m],[0 0.08]);
h3 = line([0 0.08],[Y_m Y_m]);
hold off;
        


end