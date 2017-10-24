function [triggers,do,stimset_reps] = getblinkstimtriggers_(ds, mycell, mycellname, dirname);

triggers = {};

load([getpathname(ds) filesep dirname filesep 'stims.mat']);


n = numStims(saveScript);
do = getDisplayOrder(saveScript);
stimset_reps = length(do)/n;
[mti2_,starttime_] = tpcorrectmti(MTI2,[getpathname(ds) filesep dirname filesep 'stimtimes.txt'],1);

ON_times = [];
OFF_times = [];

for i=1:numStims(saveScript),
    p = getparameters(get(saveScript,i));
    ON_times(end+1) = p.bgpause(1);
    OFF_times(end+1) = p.bgpause(2);
end;

ON_times = unique(ON_times);
OFF_times = unique(OFF_times);

for i=1:length(ON_times), 
    for j=1:length(OFF_times),
        triggers{i,j} = struct('trainOnset',[],'ONcommand',[],'OFFcommand',[],'pulseDuration',[],'cycleOnset',[]);
    end;
end;

frame_rate = 60;
projector_frame = 1000/frame_rate;
for i=1:numStims(saveScript),
    p = getparameters(get(saveScript,i));
    on_index = find(p.bgpause(1)==ON_times); 
    off_index = find(p.bgpause(2)==OFF_times);
    triggers{on_index,off_index}.ONcommand = ON_times(on_index);
    triggers{on_index,off_index}.OFFcommand = OFF_times(off_index);
    triggers{on_index,off_index}.pulseDuration = (ceil(1000*ON_times(on_index)/projector_frame)*projector_frame)/1000;
    trainOnsetIndexes = find(do==i);
    
    for j=1:length(trainOnsetIndexes),
        triggers{on_index,off_index}.trainOnset(j) = mti2_{trainOnsetIndexes(j)}.frameTimes(1);
        for k=1:p.repeat,
            triggers{on_index,off_index}.cycleOnset(j,:) = mti2_{trainOnsetIndexes(j)}.frameTimes(1:end);
        end;
    end;
    
end;

end



