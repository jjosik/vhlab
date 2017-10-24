function [ f ] = alt_analyze_spatialphase_Nspace( test_frequencies, choose_stim, fit_spatialphase, position_data, power_setting )
%*****IMPORTANT*****CHOOSE ONLY ONE DATE DIRECTORY IN THE INITIAL DIALOG
%BOX - attempts to analyze more than one date will result in erroneous
%analysis

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
    mypath = ['/Users/vhlab/Desktop/test/'];
    check_dir = date_directories{i};
    ds = dirstruct([mypath check_dir]);
    database = getexperimentfile(ds);
    mydata = load(database,'-mat');
    [celldata,cellnames]=load2celllist(getexperimentfile(ds),'cell*','-mat');
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
    ds_g = dirstruct(cellexp_dir);
    database_g = getexperimentfile(ds_g);
    mydata_g = load(database_g,'-mat');
    current_mycell = char(collected_cells{k,3});
    celldata_g = load2celllist(getexperimentfile(ds_g),strtrim(current_mycell),'-mat');
    [A,I_A] = findassociate(celldata_g{1},'test_spectral_phase','','');
    [B,I_B] = findassociate(celldata_g{1},'new_blink_firingrate','','');
    G_spectral_phase(k).phase = A.data.phase;
    G_spectral_phase(k).freqs = A.data.frequencies;
    G_cycleFR(:,:,k) = B.data;
    sub_P = alt_on_off_ttest(G_cycleFR);
    keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

keep_indices = find(keep_grandcell==1);
new_G_spectral_phase = cell(size(keep_indices,1),1);
new_G_freq = cell(size(keep_indices,1),1);
reselected_cellnames = cell(size(keep_indices,1),1);
for kk= 1:length(keep_indices),
    new_G_spectral_phase{kk,1} = G_spectral_phase(keep_indices(kk,1),1).phase;
    new_G_freq{kk,1} = G_spectral_phase(keep_indices(kk,1),1).freqs;
    reselected_cellnames{kk,1} = cell2mat(collected_cells(keep_indices(kk,1),3));
end
for k=1:length(keep_indices),
    current_set = new_G_freq{k,1}(1,1,choose_stim);
    frequency_set_ = cell2mat(current_set);
    frequency_set = reshape(frequency_set_,size(frequency_set_,2),1);
    for k2=1:length(test_frequencies),
        frequency_diffs{k2,1} = abs(frequency_set-test_frequencies(k2,1));
        current_diffs = cell2mat(frequency_diffs(k2,1));
        new_spectral_set = new_G_spectral_phase{k,1}(choose_stim,1);
        new_spectral_array_ = cell2mat(new_spectral_set);
        new_spectral_array = reshape(new_spectral_array_,size(new_spectral_array_,2),1);
        test_phase(k,k2) = new_spectral_array(find(current_diffs==min(current_diffs)));
    end
end

%neural surface distance info
posstim_n = length(position_data.fitcoordinates)/2;
posdata_arraystart = (posstim_n-(posstim_n-1))+((power_setting-1)*posstim_n);
posdata_arrayend = posdata_arraystart+(posstim_n-1);
power_subset = position_data.fitcoordinates(1,posdata_arraystart:posdata_arrayend);
posdata_cellorder = cell(length(power_subset),1);
posdata_cellcoord = cell(length(power_subset),1);
for m = 1:length(power_subset),
    posdata_cellorder{m,1} = cell2mat(power_subset(1,m).cellname);
    posdata_cellcoord{m,1} = power_subset(1,m).coordinates;
end
posstim_XY = cell2mat(reshape(posdata_cellcoord,1,length(posdata_cellcoord)));
[TF,LOC] = ismember(posdata_cellorder,reselected_cellnames); 

%make 3-D plot w/ 2 spatial dimensions from position stim. data against
%phase

for k3 = 1:length(test_frequencies),
    f = figure;
    pl = scatter3(posstim_XY(1,:),posstim_XY(2,:),test_phase(LOC,k3));
    hold on;
    xlabel('neural space, x (microns)');
    ylabel('neural space, y (microns)');
    zlabel('Phase (radians)');
    zlim([-pi pi]);
    title(['Phase response by cell location at ',num2str(test_frequencies(k3,1)),' Hz, ',mat2str(cell2mat(date_directories(1)))]);
    hold off;
    if fit_spatialphase == 1,
    else
    end
end
    
  
end
        