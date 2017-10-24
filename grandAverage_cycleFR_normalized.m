function [f1,f2] = grandAverage_cycleFR_normalized(plot_it,celltypes,rerun_celllist);

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

collected_mycell_array = [];
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
    [B,I_B] = findassociate(celldata_g{1},'new_blink_firingrate','','');
    G_cycleFR(:,:,k) = B.data;
    sub_P = alt_on_off_ttest(G_cycleFR);
    keep_grandcell(k,1) = isempty(find(sub_P >= 0.05));
    cd ..
end

keep_index = find(keep_grandcell==1);
for k2 = 1:length(keep_index),
    ttest_G_cycleFR(:,:,k2) = G_cycleFR(:,:,keep_index(k2,1));
end
    
new_G_cycleFR_norm = cell(size(ttest_G_cycleFR,1),size(ttest_G_cycleFR,2));
Err_G_cycleFR_norm = cell(size(ttest_G_cycleFR,1),size(ttest_G_cycleFR,2));
for i=1:size(ttest_G_cycleFR,1),
    for j=1:size(ttest_G_cycleFR,2),
        for k=1:size(ttest_G_cycleFR,3),
            temp_G_cycleFR_norm{k,:} = ttest_G_cycleFR{i,j,k}.norm_mean_matchedcycle_FR; 
        end
        new_G_cycleFR_norm{i,j} = nanmean(cell2mat(temp_G_cycleFR_norm),1);
        Err_G_cycleFR_norm{i,j} = nanstd(cell2mat(temp_G_cycleFR_norm),1,1)./sqrt(size(temp_G_cycleFR_norm,1));
    end
end

f1 = figure;

collabels = {'0.02 on';'0.04 on';'0.1 on';'0.5 on'};
AN_col1 = annotation('textbox',[0.15 0.95 0.08 0.04],'String',collabels{1},'EdgeColor','None');
AN_col2 = annotation('textbox',[0.37 0.95 0.08 0.04],'String',collabels{2},'EdgeColor','None');
AN_col3 = annotation('textbox',[0.58 0.95 0.08 0.04],'String',collabels{3},'EdgeColor','None');
AN_col4 = annotation('textbox',[0.8 0.95 0.08 0.04],'String',collabels{4},'EdgeColor','None');
set(AN_col1,'FontSize',20,'FontWeight','bold');
set(AN_col2,'FontSize',20,'FontWeight','bold');
set(AN_col3,'FontSize',20,'FontWeight','bold');
set(AN_col4,'FontSize',20,'FontWeight','bold');

for i = 1:size(G_cycleFR,1),
    subplot(1,4,i);
    prop_name = repmat({'Color'},size(G_cycleFR,1),1);
    prop_values = {'k','b','g','r'};
    cycle_num = 1:length(new_G_cycleFR_norm{1,1});
    for j = 1:size(G_cycleFR,2),
        h(j) = errorbar(cycle_num,new_G_cycleFR_norm{i,j},Err_G_cycleFR_norm{i,j},'-k',...
            'LineWidth',1.1,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',5);
        set(h(j),prop_name(j),prop_values(j));
    hold on;    
    end
    hLeg1(i) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
    if i < size(G_cycleFR,1),
        set(hLeg1(i),'Visible','Off');
    else
        set(hLeg1(i),'Location','NortheastOutside');
    end
    xlim([0 10]);
    ylim([0.5 1.2]);
    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
    ylabel('Normalized mean FR','FontSize',12,'FontWeight','bold');
    hold off;
end

firstcycle_meanFR = cell(size(G_cycleFR,1),size(G_cycleFR,2));
for i=1:size(G_cycleFR,1),
    for j=1:size(G_cycleFR,2),
        for k=1:size(G_cycleFR,3),
            firstcycle_meanFR{i,j}(end+1,1) = G_cycleFR{i,j,k}.mean_matchedcycle_FR(1,1); 
        end
        G_firstcycle_meanFR(i,j) = nanmean(firstcycle_meanFR{i,j}(:,1)); 
    end
end
G_firstcycle_meanFR_vsON = nanmean(G_firstcycle_meanFR,2);
G_firstcycle_FRerror_vsON = std(G_firstcycle_meanFR,1,2)./sqrt(size(G_firstcycle_meanFR,2));

f2 = figure;

X_set = [1:length(G_firstcycle_meanFR_vsON)];
h1 = plot(X_set,G_firstcycle_meanFR_vsON,'k');
set(h1,'LineWidth',1.1);
hold on;
eb = errorbar(G_firstcycle_meanFR_vsON,G_firstcycle_FRerror_vsON,'k.');
set(gca,'XTick',[1,2,3,4]);
set(gca,'XTickLabel',{'0.034','0.050','0.1','0.5'});
xlabel('On duration (seconds)');
ylabel('Mean Firing Rate');
xlim([0 5]);
ylim([0 (max(G_firstcycle_meanFR_vsON)+10)]);
N_string = {['N = ',num2str(length(keep_index))]};
AN = annotation('textbox',[0.2 0.75 0.1 0.1],'String',N_string,'EdgeColor','None');
set(AN,'FontSize',15);
t = title('Mean First Cycle Firing Rate vs. ON duration');
set(t,'FontSize',15);
hold off;



end


        
    