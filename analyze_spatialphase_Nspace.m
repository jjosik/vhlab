function [ f ] = analyze_spatialphase_Nspace( test_frequencies, choose_stim, check_prior,fit_spatialphase, position_data, power_setting )

plot_it = 0;
c_pre_window = 0.5;
c_post_window = 0.5;
on_windowshift = 0.02;
off_windowshift = 0.02;

addpath('/Users/vhlab/Desktop/test');
pathFolder = uigetdir();

d = dir(pathFolder);
isub = [d(:).isdir];
test_subFold = {d(isub).name}';
test_subFold(ismember(test_subFold,{'.','..'})) = [];

[s,v] = listdlg('PromptString','Select folders:',...
    'SelectionMode','multiple',...
    'ListString',test_subFold);

date_directories = test_subFold(s);
collected_cells = {};
current_count = 0;

for i = 1:length(date_directories),
    rep_segment = 1; plot_it = 0; pre_window = 0; post_window = 0;
    if check_prior == 1,        %*********CHECK_PRIOR BLOCK STILL UNDER CONSTRUCTION - need to incorporate all_spike_cycles into associates and search
        mypath = ['/Users/vhlab/Desktop/test/'];
        check_dir = date_directories{i};
        ds = dirstruct([mypath check_dir]);
        database = getexperimentfile(ds);
        mydata = load(database,'-mat');
        [celldata,cellnames]=load2celllist(getexperimentfile(ds),'cell*','-mat');
        N_array = cell(length(celldata),1);
        for i2 = 1:length(celldata),
            [N,I_N] = findassociate(celldata{i2},'test_spectral_phase','','');
            [Z,I_Z] = findassociate(celldata{i2},'all_spike_cycles','',''); %need to save this associate during analysis
            N_array = N;
            Z_array = I_Z;
        end
        if ~isempty(N_array)&&~isempty(Z_array),
            remain_ = cellnames;
            for i5=1:4,
                [cellname_unit,remain_] = strtok(remain_,'_');
            end
            cellname_unit_values = str2num(cell2mat(cellname_unit));
            indices_selected = find(cellname_unit_values > 400 & cellname_unit_values < 500);
            cellnames_selected = cellnames(indices_selected);
            celldata_selected = celldata(indices_selected);
            cd(char(check_dir));
            fileID = fopen('testdirinfo.txt');
            experiment_list = textscan(fileID,'%s%s');
            blank_values = strfind(experiment_list{2},'Blink');
            blank_index = find(not(cellfun('isempty',blank_values)));
            files = experiment_list{1}(blank_index(:),1);
            fclose(fileID);
            fileID = fopen('unitquality.txt');
            unit_quality = textscan(fileID,'%s%s%s%s');
            for j2=1:4,
                unit_quality{j2}(1:2) = [];  %drop the array headers
            end
            quality_list = find((strcmp(unit_quality{4},'mu'))&(strcmp(unit_quality{4},'g')));
            quality_units = reshape(quality_list,length(quality_list),1);
            if isempty(quality_units),
                quality_units = reshape((1:length(unit_quality{1})),length(unit_quality{1}),1);
            else
            end
            fclose(fileID);
            for i3 = 1:size(files),
                for i4 = 1:length(quality_units),
                    mycell_locations{i4,2} = cellnames_selected{quality_units(i4,1),1};
                    mycell_locations{i4,1} = files{i3};
                end
            end
            cd ..
        else
        end
    else
        [mycell,mycell_locations,on_offduration,all_spike_cycles,spatial_cell_map]=analyzechr2blinkingresponses(date_directories{i},...
            rep_segment,plot_it,pre_window,post_window);
    end
    for j=1:size(mycell_locations,1),
        collected_cells{end+1,1} = date_directories{i};
        current_count = current_count + 1;
    end
    column_start = (current_count - size(mycell_locations,1))+1;
    for j2=1:size(mycell_locations,1),
        collected_cells{column_start+(j2-1),2} = mycell_locations{j2,1};
        collected_cells{column_start+(j2-1),3} = mycell_locations{j2,2};
    end
end

G_spectral_phase(size(collected_cells,1),1) = struct('phase',[],'freqs',[]);
for k=1:size(collected_cells,1),
    cd(collected_cells{k,1});
    cellexp_dir = pwd;
    ds_c = dirstruct(cellexp_dir);
    database_c = getexperimentfile(ds_c);
    mydata_c = load(database_c,'-mat');
    current_mycell = char(collected_cells{k,3});
    celldata_c = load2celllist(getexperimentfile(ds_c),strtrim(current_mycell),'-mat');
    [A,I_A] = findassociate(celldata_c{1},'test_spectral_phase','','');
    %[Z,I_Z] = findassociate(celldata_c{1},'all_spike_cycles','','');
    G_spectral_phase(k).phase = A.data.phase;
    G_spectral_phase(k).freqs = A.data.frequencies;
    [sub_P,m_ON_cycle_FR_inst,m_OFF_cycle_FR_inst] = on_off_ttest(on_windowshift,off_windowshift,all_spike_cycles,on_offduration,c_pre_window);
    keep_cells(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

keep_indices = find(keep_cells==1);
new_G_spectral_phase = cell(size(keep_indices,1),1);
new_G_freq = cell(size(keep_indices,1),1);
reselected_cellnames = cell(size(keep_indices,1),1);
for kk= 1:length(keep_indices),
    new_G_spectral_phase{kk,1} = G_spectral_phase(keep_indices(kk,1),1).phase;
    new_G_freq{kk,1} = G_spectral_phase(keep_indices(kk,1),1).freqs;
    reselected_cellnames{kk,1} = collected_cells(keep_indices(kk,1),3);
end

%neural surface distance info
posstim_n = length(position_data.fitcoordinates)/2;
posdata_arraystart = (posstim_n-(posstim_n-1))+(power_setting)*posstim_n;
power_subset = position_data.fitcoordinates(1,posdata_arraystart:end);
posdata_cellorder = cell(length(power_subset),1);
posdata_cellcoord = cell(length(power_subset),1);
for m = 1:length(power_subset),
    posdata_cellorder{m,1} = power_subset(1,m).cellname;
    posdata_cellcoord{m,1} = power_subset(1,m).coordinates;
end
posstim_XY = cell2mat(reshape(posdata_cellcoord,1,length(posdata_cellcoord)));
phase_indices = cell(length(posdata_cellorder),1);
for m1 = 1:length(posdata_cellorder),
    phase_indices = find(strcmp(reselected_cellnames,posdata_cellorder{m1,1}));
end
        
%make 3-D plot w/ 2 spatial dimensions from position stim. data against
%phase

