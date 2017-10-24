function [ mycell,mycell_locations,on_offduration,all_spike_cycles,spatial_cell_map ] = analyzechr2blinkingresponses( data_dir, rep_segment,celltypes, plot_it, pre_window, post_window, varargin )
%ANALYZECHR2BLINKINGRESPONSES - Accesses appropriate databases from
%requested date folder and selects appropriate cells for analysis.  Directs
%celldata for each to appropriate analysis function depending on whether window of
%interest is stim. rep. (all 10 short cycles), or cycle rep.  Finally generates associates
%derived from that analysis.  
% ** INPUTS **
%       DATA_DIR - date folder to be analyzed
%       REP_SEGMENT - Set to '1', if window interest is stim. rep., '2' if 
%       window of interest is cycle rep.
%       CELLTYPES - Set to 'g'(and/or 'e') to select cells identified as
%       "good single units" in unit_quality.txt files; set to 'mu' to
%       analyze cells labeled "multiunit"  (DOESN'T CURRENTLY WORK)
%       PLOT_IT - set to '1' if plots are to be diplayed, '0' otherwise.
%       PRE_WINDOW = select in seconds how much time before stimulus
%       trigger should be included in plot displays
%       POST_WINDOW = select in seconds how much time AFTER stimulus ends
%       should be included in plot displays

% ***IMPORTANT -- This version needs to be in directory Desktop/test on
% Colonel to run.***

home_dir = pwd;
mypath = ['/Users/vhlab/Desktop/test/'];
ds = dirstruct([mypath data_dir]);
database = getexperimentfile(ds);
mydata = load(database,'-mat');
[celldata,cellnames]=load2celllist(getexperimentfile(ds),'cell*','-mat');
%extract relevant info from cell filenames
remain = cellnames;
for i=1:4,
    [cellname_unit,remain] = strtok(remain,'_');
end
cellname_unit_values = str2num(cell2mat(cellname_unit));
indices_selected = find(cellname_unit_values > 400 & cellname_unit_values < 500);
cellnames_selected = cellnames(indices_selected);   %***consider adding sort; there are scenarios where list order gets scrambled
celldata_selected = celldata(indices_selected);

cd(data_dir);
%access testdirinfo files to identify 'Blink' directories
%try
fileID = fopen('testdirinfo.txt');
%catch err
%    cd([pwd filesep data_dir]);
%    fileID = fopen('testdirinfo.txt');
%catch 
%    dname = uigetdir('/Desktop/test');
%    cd(dname);
%    fileID = fopen('testdirinfo.txt');
%end
experiment_list = textscan(fileID,'%s%s');
blank_values = strfind(experiment_list{2},'Blink');
blank_index = find(not(cellfun('isempty',blank_values)));
files = experiment_list{1}(blank_index(:),1);
fclose(fileID);
%access unitquality files to further narrow cell selection
fileID = fopen('unitquality.txt');
unit_quality = textscan(fileID,'%s%s%s%s');
for j=1:4,
    unit_quality{j}(1:2) = [];  %drop the array headers: columns are 1. unit #, 2. unit letter, 3. good files, and 4. unit quality 
end
%if exists(celltypes),
%    quality_units = find(strcmp(unit_quality{4},celltypes));
%else
quality_list = find((strcmp(unit_quality{4},'mu'))&(strcmp(unit_quality{4},'g')));
quality_units = reshape(quality_list,length(quality_list),1);  
%end
if isempty(quality_units),
    quality_units = reshape((1:length(unit_quality{1})),length(unit_quality{1}),1);
else
end
fclose(fileID);
do_spatialphase = 0;            %*******************change
if do_spatialphase == 1,        %this 'if' block added 2/20/15 to accommodate analyze_spatialphase.m function
    cd(char(files));
    fileID = fopen('vhlvanaloginput.vlh');
    ai_info = textscan(fileID,'%s%s');
    ai_field = cell2str(ai_info{2}(1));
    [~,remainder] = strtok(ai_field,'i');
    ai_range = strtrim(regexprep(remainder,'[[i][''][}]]',''));
    ai_range_vec = eval(ai_range)';
    fclose(fileID);
    fileID = fopen('vhlv_channelgrouping.txt');
    grouping_list = textscan(fileID,'%s%s%s');
    for p=1:3,
        grouping_list{p}(1) = [];
    end
    nolabel_grouping_list{1} = grouping_list{2};
    nolabel_grouping_list{2} = grouping_list{3};
    fclose(fileID);
    %*********per Rig 1 mapping*****************
    left_depth_indicator = [32;34;36;38;33;47;42;41;43;40;44;39;45;37;46;35];  %order of ai channels on left side of NeuroNexus probe from bottom to top
    right_depth_indicator = [63;61;59;57;62;58;48;50;49;52;51;54;53;56;55;60]; %order of ai channels on right side of NeuroNexus probe from bottom to top
    %*******************************************
    cellnum_ai_map = cell(length(nolabel_grouping_list{1}),2);
    hold_for_conversion = [];
    for q = 1:length(nolabel_grouping_list{1}),
        cellnum_ai_map{q,1} = nolabel_grouping_list{1}(q);
        hold_for_conversion = strtrim(regexprep(nolabel_grouping_list{2}(q),'[[[][]]]',''));
        cellnum_ai_map{q,2} = ai_range_vec(cellfun(@str2num,hold_for_conversion));
    end
    fifo_cellorder_left = cell(length(cellnum_ai_map),1);
    fifo_cellorder_right = cell(length(cellnum_ai_map),1);
    spatial_cell_map(1,1) = struct('directory',[],'left_channellist',[],'right_channellist',[]);
    ai_num = cell2mat(cellnum_ai_map(:,2));
    for r = 1:length(left_depth_indicator),
        fifo_cellorder_left{r,1} = cellnum_ai_map(find(ai_num==left_depth_indicator(r)),1);  %again, cell #s top-to-bottom in array column represent channel locations bottom-to-top
        fifo_cellorder_right{r,1} = cellnum_ai_map(find(ai_num==right_depth_indicator(r)),1);
    end
    spatial_cell_map(1).directory = files;
    spatial_cell_map(1).left_channellist = fifo_cellorder_left;
    spatial_cell_map(1).right_channellist = fifo_cellorder_right;
    cd ..
else 
    spatial_cell_map = [];
end
    

if rep_segment == 1,  %if window of interest is by stim. rep.********REP_SEGMENT == 1 CODE BLOCK NEEDED IF RUNNING HIST PLOT GRAND ANALYSIS
    for j = 1:size(files),
        cd(files{j});
        pre_window = 0; post_window = 0; plot_it = 0;
        dirname = files(j);
        for k = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(k,1)};
            mycellname = cellnames_selected{quality_units(k,1),1};
            [stimulus_triggers,all_spike_cycles,all_cycles_response, cycle_FR_inst,FR_cycle_lineavg,cycle_meanFR,stim_index_all,on_offduration,on_offreps] =...
                getblinkingstim_stacked_cycles(mycell,...
                mycellname,pre_window,post_window,plot_it);
            cd ..
            %fullfilename = [home_dir filesep data_dir filesep dirname filesep mycellname];
            %fprintf('%s \n',['Analyzing file ', fullfilename]);
            [matched_cycle_response_rate,delta_median_cycle_lat,norm_cycle_meanFR,delta_matched_cycle_MJ,stim_median_firstspike,...
                stim_median_jitter,norm_matched_meanFR,matched_cycle_meanFR,C_prime,F_prime,mycell] = ...
                blinking_cycletrain_info(all_spike_cycles,all_cycles_response,cycle_FR_inst,FR_cycle_lineavg,cycle_meanFR,stimulus_triggers,...
                stim_index_all,on_offduration,on_offreps,mycell,mycellname,plot_it,ds);
            %after analysis function has run
            typeA = 'stim_med_latency'; typeB = 'stim_meanFR'; typeC = 'stim_med_jitter'; typeD = 'fourier_coefficients'; typeE = 'fourier_freqs';
            owner = 'jo/analyzechr2blinkingresponses.m';
            description = 'see blinking_cycletrain_info.m Section 3. for data analysis methods';
            mycell = associate(mycell,typeA,owner,stim_median_firstspike,description);
            mycell = associate(mycell,typeC,owner,stim_median_jitter,description);
            mycell = associate(mycell,typeD,owner,C_prime,description);
            mycell = associate(mycell,typeE,owner,F_prime,description);
            mycell_locations{k,2} = cellnames_selected{quality_units(k,1),1};
            mycell_locations{k,1} = files{j};
            saveexpvar(ds,mycell,mycellname,0);
        end
        cd ../../
    end
elseif rep_segment == 2, %if window of interest is by matched cycle
    for jj = 1:size(files),
        cd(files{jj});
        pre_window = 0; post_window = 0; plot_it = 0;
        %dirname = files(jj);
        for kk = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(kk)};
            mycellname = cellnames_selected{quality_units(kk),1};
            [stimulus_triggers,all_spike_cycles,all_cycles_response, cycle_FR_inst,FR_cycle_lineavg,cycle_meanFR,stim_index_all,on_offduration,on_offreps] =...
                getblinkingstim_stacked_cycles(mycell,...
                mycellname,pre_window,post_window,plot_it);
            cd ..
            %fullfilename = [home_dir filesep data_dir filesep dirname filesep mycellname];
            %fprintf('%s \n',['Analyzing file ', fullfilename]);
            [matched_cycle_response_rate,delta_median_cycle_lat,norm_cycle_meanFR,delta_matched_cycle_MJ,stim_median_firstspike,...
                stim_median_jitter,norm_matched_meanFR,matched_cycle_meanFR,C_prime,F_prime,mycell] = ...
                blinking_cycletrain_info(all_spike_cycles,all_cycles_response,cycle_FR_inst,FR_cycle_lineavg,cycle_meanFR,stimulus_triggers,...
                stim_index_all,on_offduration,on_offreps,mycellname,plot_it);
            type1 = 'blink_cycle_response'; type2 = 'blink_med_latency'; type3 = 'blink_cycle_meanFR'; type4 = 'blink_med_jitter'; type5 = 'cycle_meanFR_nonnorm';
            owner = 'jo/analyzechr2blinkingresponses.m';
            description = 'see blinking_cycletrain_info.m Section 2. for data analysis methods';
            mycell = associate(mycell,type1,owner,matched_cycle_response_rate,description);
            mycell = associate(mycell,type2,owner,delta_median_cycle_lat,description);
            mycell = associate(mycell,type3,owner,norm_matched_meanFR,description); %switched from norm_cycle_meanFR to norm_matched_meanFR 4/9/15 (first uses inst. FR, 2nd simple cycle ave.)
            mycell = associate(mycell,type4,owner,delta_matched_cycle_MJ,description);
            mycell = associate(mycell,type5,owner,matched_cycle_meanFR,description);
            mycell_locations{kk,2} = cellnames_selected{quality_units(kk),1};  %column two gives full cellname
            mycell_locations{kk,1} = files{jj};  %column one gives t folder location
            saveexpvar(ds,mycell,mycellname,0);
        end
        cd ../../
    end
elseif rep_segment == 3,
    for jj = 1:size(files),
        cd(files{jj});
        dirname = files{jj};
        for kk = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(kk)};
            mycellname = cellnames_selected{quality_units(kk),1};
            [triggers,~,stimset_reps] = getblinkstimtriggers_(ds,mycell,mycellname,dirname);
            firstspike_output = analyze_latency(ds,dirname,mycell,mycellname,triggers);
            FR_output = analyze_firingrate(ds,triggers,stimset_reps,mycell,mycellname);
            type1 = 'new_blink_firstspike'; type2 = 'new_blink_firingrate';
            owner = 'jo/analyzechr2blinkingresponses.m';
            description = 'see analyze_latency.m and analyze_firingrate.m for analysis methods';
            mycell = associate(mycell,type1,owner,firstspike_output,description);
            mycell = associate(mycell,type2,owner,FR_output,description);
            mycell_locations{kk,2} = cellnames_selected{quality_units(kk),1};  %column two gives full cellname
            mycell_locations{kk,1} = files{jj};
            saveexpvar(ds,mycell,mycellname,0);
        end
        cd ../../
    end
    on_offduration = [];
    all_spike_cycles = [];
    spatial_cell_map = [];
elseif rep_segment == 4,
    for jj = 1:size(files),
        cd(files{jj});
        dirname = files{jj};
        use_windowing = 0;
        cycle_window = 1;
        for kk = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(kk)};
            mycellname = cellnames_selected{quality_units(kk),1};
            [triggers,~,stimset_reps] = getblinkstimtriggers_(ds,mycell,mycellname,dirname);
            fourier_output = analyze_fourier(ds,triggers,stimset_reps,mycell,mycellname,use_windowing,cycle_window);
            type1 = 'blink_fourier_output'; type2 = 'blink_fourier_output_newwindow';
            owner = 'jo/analyzechr2blinkingresponses.m';
            description = 'see analyze_fourier.m for analysis methods';
            if cycle_window == 1,
                mycell = associate(mycell,type2,owner,fourier_output,description);
            else
                mycell = associate(mycell,type1,owner,fourier_output,description);
            end
            mycell_locations{kk,2} = cellnames_selected{quality_units(kk),1};  %column two gives full cellname
            mycell_locations{kk,1} = files{jj};
            saveexpvar(ds,mycell,mycellname,0);
        end
        cd ../../
    end
    on_offduration = [];
    all_spike_cycles = [];
    spatial_cell_map = [];
elseif rep_segment == 5,
    for jj = 1:size(files),
        cd(files{jj});
        dirname = files{jj};
        cycle_window = 1;
        save_cellplots = 1;
        normalize_cellcorr = 1;
        analyze_band = 1;
        modulation_band = [25 50];
        for kk = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(kk)};
            mycellname = cellnames_selected{quality_units(kk),1};
            [triggers,~,stimset_reps] = getblinkstimtriggers_(ds,mycell,mycellname,dirname);
            [autocorr_output,f,on_offduration] = analyze_spiketrain_autocorr(triggers,stimset_reps,mycell,mycellname,cycle_window,save_cellplots,...
                normalize_cellcorr,analyze_band,modulation_band);
            type1 = 'blink_ACH_output'; type2 = 'blink_command';
            owner = 'jo/analyzechr2blinkingresponses.m';
            description = 'see analyze_spiketrain_autocorr.m for analysis methods';
            mycell = associate(mycell,type1,owner,autocorr_output,description);
            mycell = associate(mycell,type2,owner,on_offduration,description);
            mycell_locations{kk,2} = cellnames_selected{quality_units(kk),1};  %column two gives full cellname
            mycell_locations{kk,1} = files{jj};
            saveexpvar(ds,mycell,mycellname,0);
        end
        cd ../../
    end
    on_offduration = [];
    all_spike_cycles = [];
    spatial_cell_map = [];
elseif rep_segment == 6,
    for jj = 1:size(files),
        cd(files{jj});
        dirname = files{jj};
        cycle_window = 1;
        save_cellplots = 0;
        normalize_cellcorr = 0;
        %if abs(nargin) > 6,
            ach_w = varargin{1};
            modulation_band = [ach_w.f_min(1,1) ach_w.f_max(1,1)];
            analyze_band = 1;
        %else
        %    ach_w = {};
        %    analyze_band = 0;
        %end
        for kk = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(kk)};
            mycellname = cellnames_selected{quality_units(kk),1};
            [triggers,~,stimset_reps] = getblinkstimtriggers_(ds,mycell,mycellname,dirname);
            [autocorr_output,on_offduration] = analyze_spiketrain_autocorr(triggers,stimset_reps,mycell,mycellname,cycle_window,save_cellplots,...
                normalize_cellcorr,analyze_band,modulation_band,ach_w);
            type = 'blink_ACH_output_Oscore'; type_ = 'blink_command_Oscore';
            owner = 'jo/analyzechr2blinkingresponses.m';
            description = 'see analyze_spiketrain_autocorr.m for analysis methods';
            mycell = associate(mycell,type,owner,autocorr_output,description);
            mycell = associate(mycell,type_,owner,on_offduration,description);
            mycell_locations{kk,2} = cellnames_selected{quality_units(kk),1};  %column two gives full cellname
            mycell_locations{kk,1} = files{jj};
            saveexpvar(ds,mycell,mycellname,0);
            end
        cd ../../
    end
    on_offduration = [];
    all_spike_cycles = [];
    spatial_cell_map = [];
elseif rep_segment == 0,
    for jj = 1:size(files),
        cd(files{jj});
        dirname = files{jj};
        for kk = 1:length(quality_units),
            mycell = celldata_selected{1,quality_units(kk)};
            mycellname = cellnames_selected{quality_units(kk),1};
            mycell_locations{kk,2} = cellnames_selected{quality_units(kk),1};  %column two gives full cellname
            mycell_locations{kk,1} = files{jj};
        end
        cd ../../
    end
    on_offduration = [];
    all_spike_cycles = [];
    spatial_cell_map = [];
    
            
end
    
    


end

