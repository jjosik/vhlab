function  sub_P = alt_on_off_ttest( G_cycleFR )

on_rate = cell(size(G_cycleFR,1),size(G_cycleFR,2));
off_rate = cell(size(G_cycleFR,1),size(G_cycleFR,2));
for i = 1:size(G_cycleFR,3),
    for j = 1:size(G_cycleFR,1),
        for k = 1:size(G_cycleFR,2),
            on_rate{j,k}(:,1) = nanmean(G_cycleFR{j,k,i}.cycle_FR,2);
            off_rate{j,k}(:,1) = G_cycleFR{j,k,i}.null_FR;
        end
    end
    array_on_rate = cell2mat(on_rate);
    on_dist = reshape(array_on_rate,(size(array_on_rate,1)*size(array_on_rate,2)),1);
    array_off_rate = cell2mat(off_rate);
    off_dist = reshape(array_off_rate,(size(array_off_rate,1)*size(array_off_rate,2)),1);
    [H,p] = ttest(on_dist,off_dist);
    sub_P(i,1) = p;
end


    
end
   
            
            

