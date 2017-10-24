function [f] = plot_coherence(choose_stim)

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
    for i=1:size(A.data,1),
        for j=1:size(A.data,2),
            G_power_spect{i,j,k} = A.data{i,j}.mean_stim_power_spect;
            G_stimtrain_f{i,j} = cell2mat(A.data{i,j}.stimtrain_f(1,1))';
        end
    end
    cd ..
end

coherence = {};
for k1=1:size(grand_cells,1),
    for k2 = k1:size(grand_cells,1)-1,
        cd(grand_cells{k1,1});
        cellexp_dir1 = pwd;
        ds_g1 = dirstruct(cellexp_dir1);
        database_g1 = getexperimentfile(ds_g1);
        mydata_g1 = load(database_g1,'-mat');
        mycell1 = char(grand_cells{k1,3});
        celldata_g1 = load2celllist(getexperimentfile(ds_g1),strtrim(mycell1),'-mat');
        [A,I_A] = findassociate(celldata_g1{1},'blink_fourier_output','','');
        coll_cell1_data = A.data{choose_stim(1,1),choose_stim(1,2)}.cycle_spike_times;
        cd ../grand_cells{k2,1};
        cellexp_dir2 = pwd;
        ds_g2 = dirstruct(cellexp_dir2);
        database_g2 = getexperimentfile(ds_g2);
        mydata_g2 = load(database_g2,'-mat');
        mycell2 = char(grand_cells{k2+1,3});
        celldata_g2 = load2celllist(getexperimentfile(ds_g2),strtrim(mycell2),'-mat');
        [B,I_B] = findassociate(celldata_g2{1},'blink_fourier_output','','');
        coll_cell2_data = B.data{choose_stim(1,1),choose_stim(1,2)}.cycle_spike_times;
        cross_spect_output = compute_crossspect(coll_cell1_data,coll_cell2_data,f_times);
        coherence{k1,k2+1} = cross_spect_output;
    end
end
        