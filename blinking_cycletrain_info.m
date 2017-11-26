function [ matched_cycle_response_rate,delta_median_cycle_lat,norm_cycle_meanFR,delta_matched_cycle_MJ,stim_median_firstspike,stim_median_jitter,C_prime,F_prime,mycell ] = ...
    blinking_cycletrain_info( all_spike_cycles,all_cycles_response,cycle_FR_inst,FR_cycle_lineavg,cycle_meanFR,stimulus_triggers,...
    stim_index_all,on_offduration,on_offreps,mycell,mycellname,plot_it,ds )
%BLINKING_CYCLETRAIN_INFO - Computes a series of spiketrain stats at 3
%levels:  
%       1.) raw observations (e.g. spike count, time to first spike) relevant to all
%       cycles.
%       2.) across sequence-matched cycles, i.e. mean counts/rates/jitter,etc.
%       tallied separately over on/off Cycle 1, Cycle 2, Cycle 3 etc.  
%       Generates a summary plot of key stats over all ON/OFF combinations.
%       3.) across stim.-matched cycles, i.e. cycle order is disregarded and
%       stats are tallied over all cycles of the same on/off timings.
%       Also generates plots of Fourier amplitude and phase spectra.
%
%       In all cases, the observational time window is the length of the
%       on-off cycle itself. 

fit_it = 1;
error_on = 0;
save_it = 1;
use_windowing = 0;
do_fourier = 1;
do_alt_fourier = 0;
plot_it_f = 1;  %NOTE: use of plot_it_f is opposite of plot_it usage in other functions, i.e. necessary to make plot, even if saving and not displaying
do_autocor = 1;

% 1.) for all cycles
pos_spikes = [];
for i = 1:size(all_spike_cycles,3),  %for n stim. on/off combos.
    for j = 1:size(all_spike_cycles,2), %for # reps. each stim sequence
        for k = 1:size(all_spike_cycles,1), % for # of on/off cycles in single sequence
            cycle_counts(k,j,i) = length(find(all_spike_cycles{k,j,i}(:)));
            if ~isempty(all_spike_cycles{k,j,i}),
                pos_spikes = all_spike_cycles{k,j,i}(find(all_spike_cycles{k,j,i}>0));
                if ~isempty(pos_spikes),
                    cycle_firstspike(k,j,i) = pos_spikes(1);
                else
                    cycle_firstspike(k,j,i) = NaN;
                end
            end
        end
    end
end

%2.)for sequence-matched cycles (and stim. matched)
for ii = 1:size(cycle_firstspike,3),
    for jj = 1:size(cycle_firstspike,1),
        matched_cycle_jitter(jj,ii) = nanstd(cycle_firstspike(jj,:,ii)); %std of first spike lat. in sequence-matched cycles
    end
end

matched_cycle_response_rate = squeeze(mean(all_cycles_response,2));
mean_cycle_counts = squeeze(sum(cycle_counts,2)./size(all_spike_cycles,2));    %mean count of sequence-matched cycles
%mean_cycle_latency = squeeze(sum(cycle_firstspike,2)./size(all_spike_cycles,2)); %mean time to first spike of sequence-match cycles
mean_cycle_latency = squeeze(nanmean(cycle_firstspike,2));
median_cycle_latency = squeeze(nanmedian(cycle_firstspike,2));
for a = 1:size(cycle_firstspike,3),
    for b = 1:size(cycle_firstspike,1),
       matched_cycle_MAD(b,a) = nanmedian(abs(cycle_firstspike(b,:,a)-median_cycle_latency(b,a)));
    end
end
% ALT. METHOD (not recommmended)
%for b = 1:size(mean_cycle_counts,2),
%    for a = 1:size(mean_cycle_counts,1),
%        cycle_firingrate(a,b) = mean_cycle_counts(a,b)./(on_offduration(b,1)+on_offduration(b,2));
%    end
%end

%plot cycle FRs as function of time (normalized wrt first cycle FR) (16
%plots x 10 data points each)
%count = 0;
%f1 = figure;

% ********** info from instantaneous firing rate *************
%for ii = 1:size(FR_cycle_lineavg,1),
%    for kk = 1:size(FR_cycle_lineavg,2),
%        peak_cycle_FR(ii,kk) = max(FR_cycle_lineavg{ii,kk}(:));
%    end
%    f = figure;
%    cycle_num = 1:size(FR_cycle_lineavg,2);
%    plot(cycle_num,peak_cycle_FR(ii,:),'-ko',...
%        'LineWidth',1.5,...
%        'MarkerEdgeColor','k',...
%        'MarkerFaceColor',[1 1 1],...
%        'MarkerSize',5);
%    hold on;
%    xlabel('Cycle number');
%    ylabel('Peak firing rate (spikes/sec.)');
%    title(['Mean peak cycle response; Stim on/off: ',num2str(on_offduration(ii,1)),'/',num2str(on_offduration(ii,2))]);
%    hold off;
%    %saveas(gcf,[pwd filesep mycellname '_FRbyCycleSeq_' num2str(on_offduration(ii,1)) '/' num2str(on_offduration(ii,2)) '.fig']);
%    %close(f);
%end

if ~(exist('blink_analysis_data','dir')),
    mkdir('blink_analysis_data');
end
cd('blink_analysis_data');

f1 = figure;
%t_string = [data_dir];
%AN = annotation('textbox',[0.4 0.95 0.2 0.04],'String',t_string);
collabels = {'0.02 on';'0.04 on';'0.1 on';'0.5 on'};
AN_col1 = annotation('textbox',[0.15 0.95 0.08 0.04],'String',collabels{1},'EdgeColor','None');
AN_col2 = annotation('textbox',[0.37 0.95 0.08 0.04],'String',collabels{2},'EdgeColor','None');
AN_col3 = annotation('textbox',[0.58 0.95 0.08 0.04],'String',collabels{3},'EdgeColor','None');
AN_col4 = annotation('textbox',[0.8 0.95 0.08 0.04],'String',collabels{4},'EdgeColor','None');
set(AN_col1,'FontSize',20,'FontWeight','bold');
set(AN_col2,'FontSize',20,'FontWeight','bold');
set(AN_col3,'FontSize',20,'FontWeight','bold');
set(AN_col4,'FontSize',20,'FontWeight','bold');
on_label = unique(on_offduration(:,1));
%plots mean sequence-matched cycle latency
%for i1 = 1:length(on_label),
%    %f2 = figure;
%    subplot(4,4,i1);
%    prop_name_L = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set_L = {'k','b','y','m'};
%    prop_values_L = repmat(color_set_L',length(unique(on_offduration(:,1))),1);
%    cycle_num1 = 1:size(mean_cycle_latency,1);
%    for j1 = 1:length(unique(on_offduration(:,2))),
%        m = (i1-1)*length(unique(on_offduration(:,1)))+j1;
%        for k1 = 1:size(mean_cycle_latency,1),
%            norm_mean_cycle_lat(m,k1) = mean_cycle_latency(k1,m)/mean_cycle_latency(1,m);
%        end
%        hh(m) = plot(cycle_num1,norm_mean_cycle_lat(m,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(hh(m),prop_name_L(m),prop_values_L(m));
%        hold on;
%    end
%    hLeg(i1) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
%    if i1 < length(on_label),
%        set(hLeg(i1),'Visible','Off');
%    else
%        set(hLeg(i1),'Location','NortheastOutside');
%    end
%    %ylim([min(norm_mean_cycle_lat)-(0.1*min(norm_mean_cycle_lat))  max(norm_mean_cycle_lat)+(0.1*max(norm_mean_cycle_lat))]);
%    xlabel('Cycle number');
%    ylabel('Relative mean latency');
%    %title(['Normalized mean latency across Sequence-matched cycles; Stim. on: ',num2str(on_label(i1))]);
%    hold off;
%end

%plots delta median sequence-matched latency
for i1 = 1:length(on_label),
    %f2 = figure;
    subplot(4,4,i1);
    prop_name_M = repmat({'Color'},length(on_offduration(:,1)),1);
    color_set_M = {'k','b','y','m'};
    prop_values_M = repmat(color_set_M',length(unique(on_offduration(:,1))),1);
    cycle_num1 = 1:size(median_cycle_latency,1);
    for j1 = 1:length(unique(on_offduration(:,2))),
        m = (i1-1)*length(unique(on_offduration(:,1)))+j1;
        for k1 = 1:size(mean_cycle_latency,1),
            delta_median_cycle_lat(m,k1) = median_cycle_latency(k1,m) - median_cycle_latency(1,m);
        end
        hh(m) = plot(cycle_num1,delta_median_cycle_lat(m,:),'-ko',...
            'LineWidth',1.5,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',5);
        set(hh(m),prop_name_M(m),prop_values_M(m));
        hold on;
    end
    hLeg(i1) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
    if i1 < length(on_label),
        set(hLeg(i1),'Visible','Off');
    else
        set(hLeg(i1),'Location','NortheastOutside');
    end
    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
    ylabel('Delta median latency','FontSize',12,'FontWeight','bold');
    hold off;
end

%plots mean sequenced-matched cycle FR
for i2 = 1:length(unique(on_offduration(:,1))),
    %f1 = figure;
    subplot(4,4,4+i2);
    prop_name = repmat({'Color'},length(on_offduration(:,1)),1);
    color_set = {'k','b','y','m'};
    prop_values = repmat(color_set',length(unique(on_offduration(:,1))),1);
    cycle_num2 = 1:size(cycle_meanFR,2);
    for j2 = 1:length(unique(on_offduration(:,2))),
        n = (i2-1)*length(unique(on_offduration(:,1)))+j2;
        for k2 = 1:size(cycle_meanFR,2),
            norm_cycle_meanFR(n,k2) = ...
                cycle_meanFR(n,k2)/cycle_meanFR(n,1);  % see line 165 in getblinkingstim_stacked_cycles.m for cycle_meanFR definition
        end
        h(n) = plot(cycle_num2,norm_cycle_meanFR(n,:),'-ko',...
            'LineWidth',1.5,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',5);
        set(h(n),prop_name(n),prop_values(n));
    hold on;    
    end
    hLeg1(i2) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
    if i2 < length(on_label),
        set(hLeg1(i2),'Visible','Off');
    else
        set(hLeg1(i2),'Location','NortheastOutside');
    end
    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
    ylabel('Normalized mean FR','FontSize',12,'FontWeight','bold');
    %title(['Mean Firing rate across Sequence-matched cycles, Normalized to first cycle response; Stim. on: ',num2str(on_label(i2))]);
    hold off;
end

%plots sequence-matched cycle jitter
%for i3 = 1:length(on_label),
%    %f3 = figure;
%    subplot(4,4,8+i3);
%    prop_name_J = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set_J = {'k','b','y','m'};
%    prop_values_J = repmat(color_set_J',length(unique(on_offduration(:,1))),1);
%    cycle_num3 = 1:size(matched_cycle_jitter,1);
%    for j3 = 1:length(unique(on_offduration(:,2))),
%        p = (i3-1)*length(unique(on_offduration(:,1)))+j3;
%        for k3 = 1:size(matched_cycle_jitter,1),
%            norm_matched_cycle_jitter(p,k3) = matched_cycle_jitter(k3,p)/matched_cycle_jitter(1,p);
%        end
%        hhh(p) = plot(cycle_num3,norm_matched_cycle_jitter(p,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(hhh(p),prop_name_J(p),prop_values_J(p));
%        hold on;
%    end
%    hLeg2(i3) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
%    if i3 < length(on_label),
%        set(hLeg2(i3),'Visible','Off');
%    else
%        set(hLeg2(i3),'Location','NortheastOutside');
%    end
%    xlabel('Cycle number','FontSize',12,'FontWeight',bold);
%    ylabel('Normalized Jitter','FontSize',12,'FontWeight',bold);
%    %title(['Sequence-matched cycle jitter, normalized to first cycle; Stim.on: ',num2str(on_label(i3))]);
%    hold off;
%end

%plots sequence-matched median absolute deviation (i.e. median jitter)
for i3 = 1:length(on_label),
    %f3 = figure;
    subplot(4,4,8+i3);
    prop_name_MJ = repmat({'Color'},length(on_offduration(:,1)),1);
    color_set_MJ = {'k','b','y','m'};
    prop_values_MJ = repmat(color_set_MJ',length(unique(on_offduration(:,1))),1);
    cycle_num3 = 1:size(matched_cycle_jitter,1);
    for j3 = 1:length(unique(on_offduration(:,2))),
        p = (i3-1)*length(unique(on_offduration(:,1)))+j3;
        for k3 = 1:size(matched_cycle_MAD,1),
            delta_matched_cycle_MJ(p,k3) = matched_cycle_MAD(k3,p) - matched_cycle_MAD(1,p);
        end
        hhh(p) = plot(cycle_num3,delta_matched_cycle_MJ(p,:),'-ko',...
            'LineWidth',1.5,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',5);
        set(hhh(p),prop_name_MJ(p),prop_values_MJ(p));
        hold on;
    end
    hLeg2(i3) = legend('0.02 off','0.04 off','0.1 Off','0.5 Off');
    if i3 < length(on_label),
        set(hLeg2(i3),'Visible','Off');
    else
        set(hLeg2(i3),'Location','NortheastOutside');
    end
    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
    ylabel('Delta median jitter','FontSize',12,'FontWeight','bold');
    hold off;
end

%plots sequence-matched cycle response rate
for i4 = 1:length(on_label),
    %f4 = figure;
    subplot(4,4,12+i4);
    prop_name_R = repmat({'Color'},length(on_offduration(:,1)),1);
    color_set_R = {'k','b','y','m'};
    prop_values_R = repmat(color_set_R',length(unique(on_offduration(:,1))),1);
    cycle_num4 = 1:size(matched_cycle_response_rate,2);
    for j4 = 1:length(unique(on_offduration(:,2))),
        q = (i4-1)*length(unique(on_offduration(:,1)))+j4;
        hhhh(q) = plot(cycle_num4,matched_cycle_response_rate(q,:),'-ko',...
            'LineWidth',1.5,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 1 1],...
            'MarkerSize',5);
        set(hhhh(q),prop_name_R(q),prop_values_R(q));
        hold on;
    end
    hLeg3(i4) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
    if i4 < length(on_label),
        set(hLeg3(i4),'Visible','Off');
    else
        set(hLeg3(i4),'Location','NortheastOutside');
    end
    xlabel('Cycle number','FontSize',12,'FontWeight','bold');
    ylabel('Cycle Response rate','FontSize',12,'FontWeight','bold');
    %title(['Sequence-matched cycle response rate; Stim. on: ',num2str(on_label(i4))]);
    hold off;
end
if plot_it == 0,
    saveas(gcf,[pwd filesep mycellname 'seq_analysis' '.fig']);
    close(f1);
else
end

%extra subplot showing median sequence-matched cycle latency
%f5 = figure;
%collabels = {'0.02 on';'0.04 on';'0.1 on';'0.5 on'};
%AN_col4 = annotation('textbox',[0.15 0.95 0.08 0.04],'String',collabels{1},'EdgeColor','None');
%AN_col5 = annotation('textbox',[0.37 0.95 0.08 0.04],'String',collabels{2},'EdgeColor','None');
%AN_col6 = annotation('textbox',[0.58 0.95 0.08 0.04],'String',collabels{3},'EdgeColor','None');
%AN_col7 = annotation('textbox',[0.8 0.95 0.08 0.04],'String',collabels{4},'EdgeColor','None');
%for i5 = 1:length(on_label),
%    subplot(1,4,i5);
%    prop_name_A = repmat({'Color'},length(on_offduration(:,1)),1);
%    color_set_A = {'k','b','y','m'};
%    prop_values_A = repmat(color_set_A',length(unique(on_offduration(:,1))),1);
%    cycle_num5 = 1:size(median_cycle_latency,1);
%    for j5 = 1:length(unique(on_offduration(:,2))),
%        mm = (i5-1)*length(unique(on_offduration(:,1)))+j5;
%        for k5 = 1:size(median_cycle_latency,1),
%            norm_median_cycle_lat(mm,k5) = median_cycle_latency(k5,mm)/median_cycle_latency(1,mm);
%        end
%        h5(mm) = plot(cycle_num5,norm_median_cycle_lat(mm,:),'-ko',...
%            'LineWidth',1.5,...
%            'MarkerEdgeColor','k',...
%            'MarkerFaceColor',[1 1 1],...
%            'MarkerSize',5);
%        set(h5(mm),prop_name_A(mm),prop_values_A(mm));
%        hold on;
%    end
%    hLeg(i5) = legend('0.02 off','0.04 off','0.1 off','0.5 off');
%    if i5 < length(on_label),
%        set(hLeg(i5),'Visible','Off');
%    else
%        set(hLeg(i5),'Location','NortheastOutside');
%    end
%    xlabel('Cycle number');
%    ylabel('Relative median latency');
%    %title(['Normalized median latency across Sequence-matched cycles; Stim. on: ',num2str(on_label(i1))]);
%    hold off;
%end
%%
%3.)for only stim.-matched cycles (insensitive to location in sequence)

for iii = 1:size(cycle_firstspike,3)
    all_cycle_jitter(iii,1) = std2(cycle_firstspike(:,:,iii));
    mean_firstspike(iii,1) = nanmean(nanmean(cycle_firstspike(:,:,iii)));
    stim_median_firstspike(iii,1) = nanmedian(nanmedian(cycle_firstspike(:,:,iii)));
    for a1 = 1:size(cycle_firstspike,2),
        for b1 = 1:size(cycle_firstspike,1),
            median_dev(b1,a1,iii) = abs(cycle_firstspike(b1,a1,iii)-(stim_median_firstspike(iii,1)));
        end
    end
    stim_median_jitter(iii,1) = nanmedian(nanmedian(median_dev(:,:,iii)));
end
%firstspike_skew = abs(mean_firstspike-stim_median_firstspike);
mean_stim_counts = squeeze(sum(sum(cycle_counts)));
for bb = 1:length(mean_stim_counts),
    stim_firingrate(bb,1) = mean_stim_counts(bb)./(on_offduration(bb,1)+on_offduration(bb,2));
end

%fig1 = figure;
%colormap('hot');
%[X1,Y1] = meshgrid(0.5:1:(0.5+length(unique(on_offduration(:,1)))),0.5:1:(0.5+length(unique(on_offduration(:,2)))));
%A = reshape(mean_firstspike,[length(unique(on_offduration(:,1))) length(unique(on_offduration(:,2)))]);
%latency_plot = imagesc(A);
%direction = [1 0 0];
%rotate(latency_plot,direction,180);
%set(gca,'XTick',[1,2,3,4]);
%set(gca,'YTick',[1,2,3,4]);
%set(gca,'XTickLabel',{'0.02','0.04','0.1','0.5'});
%set(gca,'YTickLabel',{'0.5','0.1','0.04','0.02'});
%xlabel('On duration (secs.)');
%ylabel('Off duration (secs.)');
%colorbar;
%title('Mean latency across all cycles by stim. on/off combos');
%if plot_it == 0,
%    saveas(gcf,[pwd filesep mycellname '_MeanLat_AC' '.fig']);
%    close(fig1);
%else
%end

%fig2 = figure;
%colormap('hot');
%[X1,Y1] = meshgrid(0.5:1:(0.5+length(unique(on_offduration(:,1)))),0.5:1:(0.5+length(unique(on_offduration(:,2)))));
%C = reshape(all_cycle_jitter,[length(unique(on_offduration(:,1))) length(unique(on_offduration(:,2)))]);
%jitter_plot = imagesc(C);
%direction = [1 0 0];
%rotate(jitter_plot,direction,180);
%set(gca,'XTick',[1,2,3,4]);
%set(gca,'YTick',[1,2,3,4]);
%set(gca,'XTickLabel',{'0.02','0.04','0.1','0.5'});
%set(gca,'YTickLabel',{'0.5','0.1','0.04','0.02'});
%xlabel('On duration (secs.)');
%ylabel('Off duration (secs.)');
%colorbar;
%title('First spike jitter across cycles by stim. on/off combos');
%if plot_it == 0,
%    saveas(gcf,[pwd filesep mycellname '_Jitter_AC' '.fig']);
%    close(fig2);
%else
%end

%fig3 = figure;
%colormap('hot');
%[X1,Y1] = meshgrid(0.5:1:(0.5+length(unique(on_offduration(:,1)))),0.5:1:(0.5+length(unique(on_offduration(:,2)))));
%B = reshape(stim_median_firstspike,[length(unique(on_offduration(:,1))) length(unique(on_offduration(:,2)))]);
%median_latplot = imagesc(B);
%direction = [1 0 0];
%rotate(median_latplot,direction,180);
%set(gca,'XTick',[1,2,3,4]);
%set(gca,'YTick',[1,2,3,4]);
%set(gca,'XTickLabel',{'0.02','0.04','0.1','0.5'});
%set(gca,'YTickLabel',{'0.5','0.1','0.04','0.02'});
%xlabel('On duration (secs.)');
%ylabel('Off duration (secs.)');
%colorbar;
%title('First spike median latency across cycle by stim. on/off combos');
%if plot_it == 0,
%    saveas(gcf,[pwd filesep mycellname '_medianLat_AC' '.fig']);
%    close(fig3);
%else
%end

%fig4 = figure;
%colormap('hot');
%[X1,Y1] = meshgrid(0.5:1:(0.5+length(unique(on_offduration(:,1)))),0.5:1:(0.5+length(unique(on_offduration(:,2)))));
%D = reshape(stim_median_jitter,[length(unique(on_offduration(:,1))) length(unique(on_offduration(:,2)))]);
%med_jitterplot = imagesc(D);
%direction = [1 0 0];
%rotate(med_jitterplot,direction,180);
%set(gca,'XTick',[1,2,3,4]);
%set(gca,'YTick',[1,2,3,4]);
%set(gca,'XTickLabel',{'0.02','0.04','0.1','0.5'});
%set(gca,'YTickLabel',{'0.5','0.1','0.04','0.02'});
%xlabel('On duration (secs.)');
%ylabel('Off duration (secs.)');
%colorbar;
%title('First spike median jitter across cycles by stim. on/off combos');
%if plot_it == 0,
%    saveas(gcf,[pwd filesep mycellname '_medianJitter_AC' '.fig']);
%    close(fig4);
%else
%end

%fig5 = figure;
%colormap('hot');
%[X1,Y1] = meshgrid(0.5:1:(0.5+length(unique(on_offduration(:,1)))),0.5:1:(0.5+length(unique(on_offduration(:,2)))));
%E = reshape(firstspike_skew,[length(unique(on_offduration(:,1))) length(unique(on_offduration(:,2)))]);
%skew_plot = imagesc(E);
%direction = [1 0 0];
%rotate(skew_plot,direction,180);
%set(gca,'XTick',[1,2,3,4]);
%set(gca,'YTick',[1,2,3,4]);
%set(gca,'XTickLabel',{'0.02','0.04','0.1','0.5'});
%set(gca,'YTickLabel',{'0.5','0.1','0.04','0.02'});
%xlabel('On duration (secs.)');
%ylabel('Off duration (secs.)');
%colorbar;
%title('Skew of first spike distribution across all cycles by stim. on/off combo');
%if plot_it == 0,
%    saveas(gcf,[pwd filesep mycellname '_skew_AC' '.fig']);
%    close(fig5);
%else 
%end
%T5 = cell2table(E,...
%    'VariableNames',{[unique(on_offduration(:,2))]},'RowNames',{[unique(on_offduration(:,1))]});
%writetable(T5,'FirstSpike_SkewTable.dat','WriteRowNames',true);

%calculate and plot success per cycle  (**to be used in paper**)
stim_success_count = sum(sum(all_cycles_response,3),2);
stim_success_rate = sum(sum(all_cycles_response,3),2)./(size(all_cycles_response,3)*size(all_cycles_response,2));
cycle_on_success_rate = sum(reshape(stim_success_count',length(unique(on_offduration(:,1))),length(unique(on_offduration(:,2)))),1)./...
    (size(all_cycles_response,3)*size(all_cycles_response,2)*length(unique(on_offduration(:,1))));
cycle_on_error = std(reshape(stim_success_rate',length(unique(on_offduration(:,1))),length(unique(on_offduration(:,2)))),0,1)...
    ./sqrt(length(unique(on_offduration(:,1))));
figA = figure;
if fit_it == 1,     %********************REMEMBER TO SET FIT_IT********************
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0],...
        'Upper',[Inf,Inf],...
        'Startpoint',[1,1]);
    ft = fittype('1-a.*(exp(-(b.*x)))','options',s);    %saturating exponential model **see alt. fit below
    [c,gof] = fit(on_offduration(:,1),stim_success_rate(:,1),ft);
    g = errorbar((unique(on_offduration(:,1))),cycle_on_success_rate,cycle_on_error,'.');
    hold on;
    f_curve = plot(c);
    set(f_curve,'Color','red');
    legend('hide');  %**check if this works
else
    g = errorbar((unique(on_offduration(:,1))),cycle_on_success_rate,cycle_on_error,'ko-',...
        'LineWidth',1.5);
end
hold on;
g2 = line([0 0.55],[1 1]);
set(g2,'LineStyle','--');
xlim([0 0.55]);
ylim([0 1.1]);
xlabel('On Time (in secs.)');
ylabel('Success rate');
an = annotation('textbox',[0.3,0.25,0.2,0.1],'String',['Adjusted R^2: ',num2str(gof.adjrsquare)]); 
set(an,'FontWeight','bold');
hold off;
if save_it == 1,
    saveas(gcf,[pwd filesep mycellname '_all_cycle_success' '.fig']);
    close(figA);
else
end

%Alt. fit (Naka-Rushton)
%ft = fittype('(x^a)./((x^a)+(b^a))','options',s);
%*NOTE - for 6-27-14, cell 6 both models provided comparable gofs: adj.
%R-squares of .6421 and .6422 respectively
        
%Section 3.2
%************************************************************************%
%************************************************************************%
%Fourier analysis - by stim.-matched cycle  (observation window is one
%cycle; dt time bins averaged across stacked cycles)
collected_coeff = cell(1,length(on_offduration));
if do_fourier == 1
    dt = 0.001;
    for j2 = 1:length(on_offduration),
        f_times = 0:dt:(on_offduration(j2,1)+on_offduration(j2,2));
        for j3 = 1:size(all_spike_cycles,2),
            for j4 = 1:size(all_spike_cycles,1),
                %optional: pre-process data with Hanning window
                if use_windowing == 1,
                    w_series = all_spike_cycles{j4,j3,j2}(:).*hanning(length(all_spike_cycles{j4,j3,j2}));
                else
                    w_series = all_spike_cycles{j4,j3,j2}(:);
                end
                spike_counts = spiketimes2bins(w_series,f_times);
                [fc,freqs] = fouriercoeffs(spike_counts,dt);
                stim_amp_coeff{j4,j3,j2} = reshape(sqrt(fc.*conj(fc)),length(fc),1);
                stim_power_coeff{j4,j3,j2} = reshape((fc.*conj(fc)),length(fc),1);
                phi{j4,j3,j2} = reshape((atan2(imag(fc),real(fc))),length(fc),1);
                stimtrain_f{j4,j3,j2} = freqs;
            end
        end
    end
    new_coeff = squeeze(reshape(stim_amp_coeff,(size(stim_amp_coeff,1)*size(stim_amp_coeff,2)),1,length(on_offduration)));
    C_prime = cell(length(on_offduration),1);
    e = cell(length(on_offduration),1);
    %derive mean amplitude spectrum
    for k1 = 1:length(on_offduration),
        f_times_re = 0:dt:(on_offduration(k1,1)+on_offduration(k1,2));
        for k2 = 1:length(f_times_re),
            for k3 = 1:size(new_coeff,1),
                C(k3,k2) = new_coeff{k3,k1}(k2);
            end
        end
        C_prime{k1,1} = 2.*(abs(sum(C,1)./size(C,1)));
        if error_on == 1,   %block contains calculations for errorbars
            e{k1,1} = std(C,1)./sqrt(size(C,1));
            %alpha = 0.05;
            %nu = 2*on_offreps(1,1);
            %error_low = nu/chi2inv((alpha/2),nu);
            %error_high = nu/chi2inv(1-(alpha/2),nu);
        else
        end
        C = [];
    end
    %derive mean power spectrum
    new_power_c = squeeze(reshape(stim_power_coeff,(size(stim_power_coeff,1)*size(stim_power_coeff,2)),1,length(on_offduration)));
    C_p_prime = cell(length(on_offduration),1);
    tmean_A = zeros(size(new_power_c,1),length(on_offduration));
    tmean_B = zeros(size(new_power_c,1),length(on_offduration));
    for k4 = 1:length(on_offduration),
        f_times_re2 = 0:dt:(on_offduration(k4,1)+on_offduration(k4,2));
        for k5 = 1:length(f_times_re2),
            for k6 = 1:size(new_power_c,1),
                C_p(k6,k5) = new_power_c{k6,k4}(k5);
            end
        end
        C_p_prime{k4,1} = sum(C_p,1)./size(C_p,1);
        if error_on == 1,   %block contains calculations for chi-square-based errorbars
            alpha = 0.05;
            nu = 2*on_offreps(1,1);
            error_low = nu/chi2inv((alpha/2),nu);
            error_high = nu/chi2inv(1-(alpha/2),nu);
        else
        end
        %Trial-averaged band power, 15-50 Hz band comparison to >50 Hz band
        current_f_set = cell2mat(stimtrain_f(1,1,k4));
        min_band_A = 15;
        max_band_A = 50;
        min_band_B = max_band_A;
        max_band_B = max(current_f_set);
        diff_minA = (abs(min_band_A-current_f_set));
        diff_maxA = (abs(max_band_A-current_f_set));
        diff_minB = (abs(min_band_B-current_f_set));
        diff_maxB = (abs(max_band_B-current_f_set));
        nearest_minA_match = find(min(diff_minA)==diff_minA);
        nearest_maxA_match = find(min(diff_maxA)==diff_maxA);
        nearest_minB_match = find(min(diff_minB)==diff_minB);
        nearest_maxB_match = find(min(diff_maxB)==diff_maxB);
        band_indices_A = [nearest_minA_match:nearest_maxA_match];
        band_indices_B = [nearest_minB_match:nearest_maxB_match];
        for r = 1:size(C_p,1),
            tmean_A(r,k4) = (sum(C_p(r,band_indices_A(:))))/length(band_indices_A);
            tmean_B(r,k4) = (sum(C_p(r,band_indices_B(:))))/length(band_indices_B);
        end
                    
        C_p = [];
    end
    %generate trial-averaged band comparison plots
    for r1 = 1:length(on_offduration),
        alpha_t = 0.05;
        [H,p] = ttest(tmean_A(:,r1),tmean_B(:,r1));
        fg = figure;
        hc = scatter(tmean_A(:,r1),tmean_B(:,r1),'k');
        hold on;
        xlabel(['Site trial-averaged power (' num2str(min_band_A) '-' num2str(max_band_A) ' Hz)']);
        ylabel(['Site trial-averaged power (>' num2str(min_band_B) ' Hz)']);
        title(['Matched trial-averaged power,' num2str(min_band_A) '-' num2str(max_band_A) ' Hz vs. >' num2str(max_band_A) ' Hz']);
        xy_max = max(cat(1,tmean_A(:,r1),tmean_B(:,r1)))+(0.1*(max(cat(1,tmean_A(:,r1),tmean_B(:,r1)))));
        xlim([0 xy_max]);
        ylim([0 xy_max]);
        unity_line = line([0 xy_max],[0 xy_max],'LineStyle','--');
        if p > alpha_t,
            N_string = {['p = ' num2str(p)]};
        else
            N_string = {['**p = ' num2str(p)]};
        end
        AN = annotation('textbox',[0.2 0.75 0.1 0.1],'String',N_string,'EdgeColor','None');
        set(AN,'FontSize',15);
        if save_it == 1,
            saveas(gcf,[pwd filesep mycellname 'trial_av_pcomp_',num2str(on_offduration(r1,1)),'_',num2str(on_offduration(r1,2)), '.fig']);
            close(fg);
        else
        end
    end
        
        
    
end
F_prime = [];
if plot_it_f == 1, %lines blocked when used with multicell_AverageAmpSpect.m
    f2 = figure;  %make amplitude spectrum figure
    for j5 = 1:length(on_offduration),
        F = stimtrain_f{1,1,j5}(:);
        F_prime{j5,1} = F;
        subplot(length(unique(on_offduration(:,1))),length(unique(on_offduration(:,2))),j5);
        if error_on == 1,
            h_f = plot(stimtrain_f{1,1,j5}(:),C_prime{j5,1}(:));
            set(h_f,'Color',[0 0 0]);       %IMPORTANT:  these errorbars still need to be fixed (not entirely happy with this method)
            %alt. method to errorbar function:
            hold on;
            for j_e = 1:length(e{j5,1}),
                h_e = line([(stimtrain_f{1,1,j5}(1,j_e)) stimtrain_f{1,1,j5}(1,j_e)],...
                    [(C_prime{j5,1}(1,j_e)-e{j5,1}(1,j_e)),(C_prime{j5,1}(1,j_e)+e{j5,1}(1,j_e))]); 
                set(h_e,'LineWidth',3);
                h_e.Color(4) = 0.3;
                w = warning('query','last');
                id = w.identifier;
                warning('off',id);
            end 
            xlim([0 (1/dt*1/2)]);
            ylim([0 2*max(C_prime{j5,1}(:))]);
            xlabel('freq (Hz)','FontSize',12,'FontWeight','bold');
            ylabel('Amplitude (V)','FontSize',12,'FontWeight','bold');
            warning('on',id);
        else
        h_f = plot(stimtrain_f{1,1,j5}(:),C_prime{j5,1}(:));
        xlim([0 (1/dt*1/2)]);
        if max(C_prime{j5,1})>0,
            ylim([0 2*max(C_prime{j5,1})]);
        else 
            ylim([0 (2*max(C_prime{j5,1})+1)]);
        end
        xlabel('freq (Hz)','FontSize',12,'FontWeight','bold');
        ylabel('Amplitude (V)','FontSize',12,'FontWeight','bold');
        end
    end
    if save_it == 1,
        saveas(gcf,[pwd filesep mycellname 'amplitude_spectrum' '.fig']);
        close(f2);
    else
    end
else
end
FF_prime = [];
if plot_it_f == 1,
    f3 = figure;    %make power spectrum figure
    for j6 = 1:length(on_offduration),
        FF = stimtrain_f{1,1,j6}(:);
        FF_prime{j6,1} = FF;
        subplot(length(unique(on_offduration(:,1))),length(unique(on_offduration(:,2))),j6);
        if error_on == 1,
            h_ff = semilogy(stimtrain_f{1,1,j6}(:),C_p_prime{j6,1}(:),...
                [stimtrain_f{1,1,j6}(:) stimtrain_f{1,1,j6}(:)],[error_low.*C_p_prime{j6,1}(:) error_high.*C_p_prime{j6,1}(:)]);
            property_array = {'LineWidth','Color'};
            value_array = {0.5,'black';0.5,'blue';0.5,'blue'};
            set(h_ff,property_array,value_array);
            xlim([0 (1/dt*1/2)]);
            ylim([0 10]);
            xlabel('freq (Hz)','FontSize',12,'FontWeight','bold');
            ylabel('Log power (dB)','FontSize',12,'FontWeight','bold');
        else
        h_ff = plot(stimtrain_f{1,1,j6}(:),C_p_prime{j6,1}(:));
        xlim([0 (1/dt*1/2)]);
        if max(C_p_prime{j6,1})>0,
            ylim([0 2*max(C_p_prime{j6,1})]);
        else 
            ylim([0 (2*max(C_p_prime{j6,1})+1)]);
        end
        xlabel('freq (Hz)','FontSize',12,'FontWeight','bold');
        ylabel('Power (V^2)','FontSize',12,'FontWeight','bold');
        end
    end
    if save_it == 1,
        saveas(gcf,[pwd filesep mycellname 'power_spectrum' '.fig']);
        close(f3);
    else
    end
else
end
%phase spectrum mean
new_phi = squeeze(reshape(phi,(size(phi,1))*(size(phi,2)),1,length(on_offduration)));
ph_prime = cell(length(on_offduration),1);
for k7 = 1:length(on_offduration),
    f_times_re3 = 0:dt:(on_offduration(k7,1)+on_offduration(k7,2));
    for k8 = 1:length(f_times_re3),
        for k9 = 1:size(new_phi,1),
            ph(k9,k8) = new_phi{k9,k7}(k8);
        end
    end
    ph_prime{k7,1} = sum(ph,1)./size(ph,1);
    ph = [];
end
F3_prime = [];
if plot_it_f == 1,
    f4 = figure;    %make phase spectrum figure
    for j7 = 1:length(on_offduration),
        F3 = stimtrain_f{1,1,j7}(:);
        F3_prime{j7,1} = F3;
        subplot(length(unique(on_offduration(:,1))),length(unique(on_offduration(:,2))),j7);
        h_f3 = plot(stimtrain_f{1,1,j7}(:),ph_prime{j7,1}(:));
        xlim([0 (1/dt*1/2)]);
        ylim([-pi pi]);
        xlabel('freq (Hz)','FontSize',12,'FontWeight','bold');
        ylabel('Phase (rads)','FontSize',12,'FontWeight','bold');
    end
    if save_it == 1,
        saveas(gcf,[pwd filesep mycellname 'phase_spectrum' '.fig']);
        close(f4);
    else
    end
else
end

%find power spectral peaks**********************
theta = 0.18; %average power spectrum threshold
for j8 = 1:length(on_offduration),
    [peaks{j8,1},locs{j8,1}] = findpeaks(C_p_prime{j8,1}(:),'MinPeakHeight',(theta/(2*pi*mean(cycle_meanFR(j8,:)))),'MinPeakDistance',10);  %changed to include thresholding 2/5/15
end
primary_peaks = cell(length(on_offduration),1);
secondary_peaks = cell(length(on_offduration),1);
tertiary_peaks = cell(length(on_offduration),1);
quarternary_peaks = cell(length(on_offduration),1);
for j9 = 1:length(on_offduration),
    if length(locs{j9,1})< 1,
        primary_peaks{j9,1} = NaN;
    else
        primary_peaks{j9,1} = stimtrain_f{1,1,j9}(locs{j9,1}(1,1));
    end
    if length(locs{j9,1})< 2,
        secondary_peaks{j9,1} = NaN;
    else
        secondary_peaks{j9,1} = stimtrain_f{1,1,j9}(locs{j9,1}(1,2));
    end
    if length(locs{j9,1})< 3,
        tertiary_peaks{j9,1} = NaN;
    else 
        tertiary_peaks{j9,1} = stimtrain_f{1,1,j9}(locs{j9,1}(1,3));
    end
    if length(locs{j9,1})< 4,
        quarternary_peaks{j9,1} = NaN;
    else 
        quarternary_peaks{j9,1} = stimtrain_f{1,1,j9}(locs{j9,1}(1,4));
    end
end 
%alt. gamma assessment method (2/6/2015) - mean power over gamma range and
%extra-gamma range; depict in scatterplot under grand ave. code to show 
%relative gamma versus "non-gamma" power
for j9 = 1:length(on_offduration),
    gamma_index_range = find(stimtrain_f{1,1,j9}>=25&stimtrain_f{1,1,j9}<=70);
    trans_gamma_index_range = find(stimtrain_f{1,1,j9}>70&stimtrain_f{1,1,j9}<100);
    sub_gamma_index_range = find(stimtrain_f{1,1,j9}>=10&stimtrain_f{1,1,j9}<25);
    gamma_power(j9,1) = mean(C_p_prime{j9,1}(gamma_index_range(1):gamma_index_range(end)));
    trans_gamma_power(j9,1) = mean(C_p_prime{j9,1}(trans_gamma_index_range(1):trans_gamma_index_range(end)));
    sub_gamma_power(j9,1) = mean(C_p_prime{j9,1}(sub_gamma_index_range(1):sub_gamma_index_range(end)));
end
    
    
%save Fourier analysis data to disk
amplitude_spectrum(1,1) = struct('cellname',[],'amplitudes',[],'frequencies',[]);
amplitude_spectrum(1).amplitudes = C_prime; 
amplitude_spectrum(1).frequencies = stimtrain_f; 
amplitude_spectrum(1).cellname = mycellname;
power_spectrum(1,1) = struct('cellname',[],'log_power',[],'frequencies',[]);
power_spectrum(1).log_power = C_p_prime; 
power_spectrum(1).frequencies = stimtrain_f; 
power_spectrum(1).cellname = mycellname;
phase_spectrum(1,1) = struct('cellname',[],'phase',[],'frequencies',[]);
phase_spectrum(1).phase = ph_prime; 
phase_spectrum(1).frequencies = stimtrain_f; 
phase_spectrum(1).cellname = mycellname;
peak_frequencies(1,1) = struct('cellname',[],'primary',[],'secondary',[],'tertiary',[],'quarternary',[]);
peak_frequencies(1).primary = primary_peaks;
peak_frequencies(1).secondary = secondary_peaks;
peak_frequencies(1).tertiary = tertiary_peaks;
peak_frequencies(1).quarternary = quarternary_peaks;
peak_frequencies(1).cellname = mycellname;
mean_bandpower(1,1) = struct('cellname',[],'g_power',[],'trans_g_power',[],'sub_g_power',[]);
mean_bandpower(1).cellname = mycellname;
mean_bandpower(1).g_power = gamma_power;
mean_bandpower(1).trans_g_power = trans_gamma_power;
mean_bandpower(1).sub_g_power = sub_gamma_power;
type_a = 'test_spectral_amplitudes'; type_p = 'test_spectral_power'; type_ph = 'test_spectral_phase'; type_k = 'test_spectral_peaks'; type_g = 'test_gamma_strength';
type_sr = 'cycle_on_success_rate';
owner = 'jo/blinking_cycletrain_info.m';
description = 'see blinking_cycletrain_info.m Section 3.2 (~line 513) for data analysis methods';
mycell = associate(mycell,type_a,owner,amplitude_spectrum,description);
mycell = associate(mycell,type_p,owner,power_spectrum,description);
mycell = associate(mycell,type_ph,owner,phase_spectrum,description);
mycell = associate(mycell,type_k,owner,peak_frequencies,description);
mycell = associate(mycell,type_g,owner,mean_bandpower,description);
mycell = associate(mycell,type_sr,owner,cycle_on_success_rate,description);
saveexpvar(ds,mycell,mycellname,0);




%ALT. Fourier analysis - by stim.-matched cycle sequence (observation
%window is length of full 10-cycle trial; dt time bins averaged across
%trials) ***********NOTE: THIS SECTION IS INCOMPLETE; MUST COMMENT OUT IF
%RUNNING FUNCTION OR SELECT DO_ALT_FOURIER = 0 ****************************************************
if do_alt_fourier == 1,
    dt_ = 0.0025;
    %first need to reorganize all_spike_cycles into 10-cycle sequence
    current_train = [];
    for n = 1:size(all_spike_cycles,3),
        for n1 = 1:size(all_spike_cycles,2),
            for n2 = 1:size(all_spike_cycles,1),
                current_train = [current_train;all_spike_cycles{n2,n1,n}(:)];
            end
                all_seq_train{n1,n} = current_train;
                current_train = [];
        end
    end
    for j6 = 1:length(on_offduration),
        f_times_ = 0:dt_:((on_offduration(j6,1)+on_offduration(j6,2))*size(all_spike_cycles,1));
        for j7 = 1:size(all_seq_train,2),
            %note: can consider windowing function, but may not be
            %necessary in long data (10-cycle) form, i.e. this is an alt. to windowed FFT. 
            spike_counts_alt = spiketimes2bins(all_seq_train{j7,j6}(:),f_times_);
            [fc_alt,freqs_alt] = fouriercoeffs(spike_counts_alt,dt_);
            long_stimamp_coeff{j7,j6} = reshape(fc_alt,length(fc_alt),1);
            long_stimpower_coeff{j7,j6} = reshape((fc_alt.*conj(fc_alt)),length(fc_alt),1);
            long_phi{j7,j6} = reshape((atan2(imag(fc_alt),real(fc_alt))),length(fc_alt),1);
            long_stimtrain_f{j7,j6} = freqs_alt;
        end
    end
    alt_new_coeff = squeeze(reshape(long_stimamp_coeff,(size(long_stimamp_coeff,1)*size(long_stimamp_coeff,2)),1,length(on_offduration)));
    C_prime = cell(length(on_offduration),1);
    %derive mean amplitude spectrum
     for kk1 = 1:length(on_offduration),
        f_times_re = 0:dt:(on_offduration(kk1,1)+on_offduration(kk1,2));
        for kk2 = 1:length(f_times_re),
            for kk3 = 1:size(alt_new_coeff,1),
                C(kk3,kk2) = alt_new_coeff{kk3,kk1}(kk2);
            end
        end
        C_prime{kk1,1} = abs(sum(C,1)./size(C,1));
        C = [];
     end
    %derive mean power spectrum
    alt_new_power_c = squeeze(reshape(long_stimpower_coeff,(size(long_stimpower_coeff,1)*size(long_stimpower_coeff,2)),1,length(on_offduration)));
    C_p_prime = cell(length(on_offduration),1);
    for kk4 = 1:length(on_offduration),
        f_times_re2 = 0:dt:(on_offduration(kk4,1)+on_offduration(kk4,2));
        for kk5 = 1:length(f_times_re2),
            for kk6 = 1:size(alt_new_power_c,1),
                C_p(kk6,kk5) = alt_new_power_c{kk6,kk4}(kk5);
            end
        end
        C_p_prime{kk4,1} = sum(C_p,1)./size(C_p,1);
        C_p = [];
    end
    %mean phase spectrum
else
    
end

%************************************************************************%
%************************************************************************%
%Autocorrelogram
if do_autocor == 1,
    ints_per_stim = cell((size(all_spike_cycles,1)*size(all_spike_cycles,2)),length(on_offduration));
    ints_per_train = [];
    cat_ints_per_stim = cell(1,length(on_offduration));
    ach_dt = 0.002;     %*********time step may need adjustment, consider trying 1 ms (alt.: consider variable time step for different stim. sets)
    for jj1 = 1:length(on_offduration),
        for jj2 = 1:size(all_spike_cycles,2),
            for jj3 = 1:size(all_spike_cycles,1),
                s_series = all_spike_cycles{jj3,jj2,jj1}(:);
                for a = 1:length(s_series),
                    tr_series = repmat(s_series(a),length(s_series),1);
                    ints_per_train = tr_series-s_series(1:end);
                end
                ints_per_train = reshape(ints_per_train,length(ints_per_train),1);
                ints_per_stim{((jj2-1)*(on_offreps(1,1))+jj3),jj1} = ints_per_train(:);
            end
        end
        cat_ints_per_stim{1,jj1} = vertcat(ints_per_stim{:,jj1});
    end
    
    %plot autocorrelation histogram for each stim set
    for jj4 = 1:length(on_offduration),
        half_shift_range = on_offduration(jj4,1)+on_offduration(jj4,2);
        ach_bin_edges = [-(half_shift_range):ach_dt:half_shift_range];
        N_ach = histc(cat_ints_per_stim{1,jj4},ach_bin_edges);
        ach_bin_centers = (ach_bin_edges(1:end-1)+ach_bin_edges(2:end))/2;
        N_ach = N_ach(1:end-1);
        f5 = figure;
        ach_h = bar(ach_bin_centers,N_ach,'k');
        hold on;
        inv_bin_centers = (-1).*(ach_bin_centers);
        ach_i = bar(inv_bin_centers,N_ach,'k');
        xlabel('Shift (secs.)');
        ylabel('Counts');
        xlim([-(half_shift_range) half_shift_range]);
        ylim([0 (max(N_ach)+1)]);
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep mycellname 'ACH_' num2str(on_offduration(jj4,1)) '_' num2str(on_offduration(jj4,2)) '.fig']);
            close(f5);
        else
        end
    end
else
end
        
        

%Compute oscillation score (Muresan et al. 2008, J. Neurophys.)
%Step 1. autocorrelation histogram (see above)
%Step 2. Smoothing the autocorrelogram - use fast Gaussian kernel (sigma
%window 1 or 2 ms for 0.5-1 kHz correlogram)



cd ..

end




