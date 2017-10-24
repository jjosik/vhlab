function [f1,f2] = grandAverage_fourier(plot_it,celltypes,cycle_window,selectcells_byjitter,rerun_celllist);

plot_it = 0;
cycle_window = 1;

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
        rep_segment = 4; pre_window = 0; post_window = 0;
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
    [A,I_A] = findassociate(celldata_g{1},'blink_fourier_output','','');
    [Z,I_Z] = findassociate(celldata_g{1},'new_blink_firstspike','','');
    G_stim_jitter{1,1,k} = Z.data.stim_med_jitter;
    if cycle_window == 1,
        [B,I_B] = findassociate(celldata_g{1},'blink_fourier_output_newwindow','','');
        for ii=1:size(B.data,1),
            for jj=1:size(B.data,2),
                G_power_spect{ii,jj,k} = B.data{ii,jj}.mean_stim_power_spect;
                %G_amp_spect{ii,jj,k} = B.data{i,j}.mean_stim_amp_coeff;
                %G_phase_spect{ii,jj,k} = B.data{i,j}.mean_phase_spectrum;
                G_stimtrain_f{ii,jj} = cell2mat(B.data{ii,jj}.stimtrain_f(1,1))';
            end
        end
    else
        for i=1:size(A.data,1),
            for j=1:size(A.data,2),
                G_power_spect{i,j,k} = A.data{i,j}.mean_stim_power_spect;
                %G_amp_spect{i,j,k} = A.data{i,j}.mean_stim_amp_coeff;
                %G_phase_spect{i,j,k} = A.data{i,j}.mean_phase_spectrum;
                G_stimtrain_f{i,j} = cell2mat(A.data{i,j}.stimtrain_f(1,1))';
            end
        end
    end
    %sub_P = alt_on_off_ttest(G_cycleFR);
    %keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

power_array = cell(size(G_power_spect,1),size(G_power_spect,2));
mean_power_array = cell(size(G_power_spect,1),size(G_power_spect,2));
power_array_error = cell(size(G_power_spect,1),size(G_power_spect,2));
if selectcells_byjitter == 1,
    include_byjitter = cellfun(@(G_stim_jitter) find(G_stim_jitter<=0.01),G_stim_jitter,'UniformOutput',false);
    keep_indices = find(cell2mat(cellfun(@(include_byjitter) size(include_byjitter,1),include_byjitter,'UniformOutput',false)));
    for i=1:size(G_power_spect,1),
        for j=1:size(G_power_spect,2),
            for k=1:size(keep_indices,1),
                power_extract_match{1,k} = G_power_spect{i,j,keep_indices(k)};
            end
            power_array{i,j} = cell2mat(power_extract_match);
            mean_power_array{i,j} = nanmean(power_array{i,j},2);
            power_array_error{i,j} = std(power_array{i,j},1,2)./size(power_array{i,j},2);
        end
    end
else
    for i=1:size(G_power_spect,1),
        for j=1:size(G_power_spect,2),
            for k=1:size(grand_cells,1),
                power_extract_match{1,k} = G_power_spect{i,j,k};
            end
            power_array{i,j} = cell2mat(power_extract_match);
            mean_power_array{i,j} = nanmean(power_array{i,j},2);
            power_array_error{i,j} = std(power_array{i,j},1,2)./size(power_array{i,j},2);
        end
    end
end

if plot_it == 1,
    f1 = figure;
    h_gfp = plot(G_stimtrain_f{4,4}(:),mean_power_array{4,4}(:),'k');
    hold on;
    set(h_gfp,'LineWidth',1.0);
    h_gffp = plot(G_stimtrain_f{4,4}(:),(mean_power_array{4,4}(:)+power_array_error{4,4}(:)),'b');
    set(h_gffp,'LineWidth',0.5);
    h_gfffp = plot(G_stimtrain_f{4,4}(:),(mean_power_array{4,4}(:)-power_array_error{4,4}(:)),'b');
    set(h_gfffp,'LineWidth',0.5);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    xlim([0 500]);
    ylim([0 1.5*max(mean_power_array{4,4}(:))]);
    hold off;
    
else
end

end
