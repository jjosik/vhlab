function cross_spectrum_output = compute_crossspect(mycell1,mycell2,f_times)

if size(mycell1,1)==size(mycell2,1)&&size(mycell1,2)==size(mycell2,2),
    check_match = 1;
else
    check_match = 0;
    mycell2 = reshape(mycell2,size(mycell1,1),size(mycell1,2));
end

cross_spect_output = {};
cross_spect_output{i,j} = struct('cross_spect',[],'cohere_num',[],'coherence',[]);

for i = 1:size(mycell1,1),
    for j = 1:size(mycell1,2),
        spike_counts1 = spiketimes2bins(mycell1{i,j},f_times);
        spike_counts2 = spiketimes2bins(mycell2{i,j},f_times);
        [fc1,freqs1] = fouriercoeffs(spike_counts1,dt);
        [fc2,freqs2] = fouriercoeffs(spike_counts2,dt);
        cross_spectrum_output.cross_spect{i,j} = reshape(fc1.*conj(fc2)/length(fc1),length(fc1),1);
        power_spect1{i,j} = reshape(fc1.*conj(fc1)/length(fc1),length(fc1),1);
        power_spect2{i,j} = reshape(fc2.*conj(fc2)/length(fc2),length(fc2),1);
    end
end

mean_cross_spect = nanmean(cell2mat(reshape(cross_spect,1,size(cross_spect,1)*size(cross_spect,2))),2);
mean_power_spect1 = nanmean(cell2mat(reshape(power_spect1,1,size(power_spect1,1)*size(power_spect1))),2);
mean_power_spect2 = nanmean(cell2mat(reshape(power_spect2,1,size(power_spect2,1)*size(power_spect2))),2);
cross_spectrum_output.cohere_num = abs(mean_cross_spect).^2;
cross_spectrum_output.coherence = (cross_spectrum_ouput.cohere_num)./(mean_power_spect1.*mean_power_spect2);
        

