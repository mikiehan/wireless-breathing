%% parameters
% Fs : sampling rate in Hz (ex 41100)
% fmin : min frequency in KHz (ex 11)
% B : Bandwidth in Hz (ex 10000)
% on : on interval in millisec (ex 30)
% off: off interval in millisec (ex 0)
% start_distance_range in m (ex 0.5) 
% filename : the recorded wav filename (without .wav) (ex 'rajat_after_run_1'
% inputWavDir : the directory that has input wav files (ex '/Users/profhan/Downloads/Wireless/Data/FMCW Data/')
% outputResultDir : the directory to save the output results (ex '/Users/profhan/Downloads/Wireless_Networking_UT-master/fft_results')
% dist_range in m (ex 1) max range of twice of distance
function rx_FMCW_param(Fs, fmin, B, on, off, start_distance_range, dist_range, filename, inputWavDir, outputResultDir)
    tic 
    close all
    vs=340;

    mic = 1;

    %most important new data %1.sin18k  - 120 ms time %2.batmapper - 30 ms time - %3.fmcw_11 - 60 ms time 
    %bandwidth - 13k and min 10k - chirp generate
    %fmin = 10; %3; %11;fmax = 22; %21;B = 12000; %10000
    %importnat
    %batmapper_1, handone - 10 k, 8 B, 60 ms - older - fmcw_11_1 - 3k, 10B , 60 ms - 
    %bat - 17 27 
    %% Variables for experiments
    fmax = fmin + B/1000;
    K=(on/1000)*Fs; % num samples of on interval
    step_size=1;

    %movement, normal breathing, angle, music 
    %% Estimate FFT for each repetition of chirp with step_size 
    disp(strcat(inputWavDir,'/',filename,'.wav'))
    [Rx,fs] = audioread(strcat(inputWavDir,'/',filename,'.wav'));

    %Applying FFT and removing a few initial data
    y = fftFilter(Rx(:,mic)',Fs,fmin*1000,fmax*1000,500);
    y = y(2*Fs:length(y));

    %mic2y = Rx(:,2)';
    mic1y = Rx(:,1)';

    % TX signal
    TX = genChirp(Fs,on/1000,off/1000,fmin*1000,B,length(Rx)).';

    %synchronizing
    yR0 = genChirp(Fs,on/1000,0,fmin*1000,B,1);

    %synchronizing to remove the intial noise
    maxIndex=syncFMCWSymbol(y,0,yR0,length(yR0));
    disp(maxIndex);

    y = y(maxIndex:length(y));

    % windowing RX signal
    w = 0.5-0.5*cos(2*pi/K*(0:length(y)-1));
    y = y.*w;

    disp('y size')
    disp(size(y))
    len_chirp = Fs*(on+off)/1000;
    Ne = floor(length(y)/(len_chirp));
    results = [];

    disp(Ne-1)
    for i = 1:step_size:Ne-1 %Ne is number of chirps % for ith chirp
        y0 = y(1+i*len_chirp:1+(i+1)*len_chirp); 
        x0 = TX(1+i*len_chirp:1+(i+1)*len_chirp).';
        prod = y0.*x0; 
        [f,fft1] = plotFFT(prod,Fs);

        %apply filter
        filter= min(fft1)*(ones(1,size(f,2)));filter(:,1000:1400) = 1; filter(:,446:458) = 1;    
        results = [results; fft1];
    end
    %% estimate distance peaks in ffts
    dist = vs*f*on/B;
    avg_series = mean(results);
    figure; % Figure 1 plots avg amplitude over all chirps  
    plot(dist/2, avg_series)
    xlim([start_distance_range/2 dist_range/2])
    xlabel('Distance (m)')
    new_results1 = results ;
    [s,d] = size(results);

    % Write this to file.
    fft_results_dir = strcat(outputResultDir); %, '/fft_results');
    %% Distance estimation to a file  
    %% each row is one chirp
    filepath = strcat(fft_results_dir,'/','results_',filename,'.txt');
    disp(size(results))
    fid = fopen(filepath,'wt');
    for ii = 1:size(results,1)
        fprintf(fid,'%g\t',results(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
%     %%Write x-axis values of distance estimation plot - not so important
%     filepath = strcat(fft_results_dir,'/','dist_arr_',filename,'.txt');
%     disp(size(dist))
%     fid = fopen(filepath,'wt');
%     for ii = 1:size(dist,2)
%         fprintf(fid,'%g\t',dist(ii));
%         fprintf(fid,'\n');
%     end
%     fclose(fid)
% 
%     %%plot the distance estimation profile for jth time-step 
%     dist_peak = [];
%     figure1  = figure;
%     for j = 100:100 % (Mikie) is this 100 just randomly chosen????!!
%         index_count = 1;cur_dist = 0.0;
%         avg_series =  results(j,:); % get 100th row from results
%         %[rows, cols] = size(avg_series);
%         %[r, c] = size(dist);
%         while cur_dist < start_distance_range
%             cur_dist = dist(1, index_count);
%             index_count = index_count + 1;
%         end
% 
%         last_index_count = index_count + 1;
%         cur_dist = 0.0;
% 
%         while cur_dist < dist_range
%             cur_dist = dist(1, last_index_count);
%             last_index_count = last_index_count + 1;
%         end
%         dist_new = dist(1,index_count:last_index_count);
% 
%         %figure;
%         hold on
%         norm_series = avg_series(1,index_count:last_index_count)./max(avg_series(1,index_count:last_index_count));
%         %norm_series(norm_series_old  < 0.5) = 0;
%         %plot(dist_new,avg_series(1,index_count:last_index_count));
%         plot(dist_new./2, norm_series);
% 
%         hold on;
%         [pks, locs] = findpeaks(avg_series(1,index_count:last_index_count),'MinPeakProminence',0.7e-4,'Annotate','extents');
% 
%         plot(dist(1,index_count:last_index_count)./2.0,norm_series,dist_new(1,locs)./2.0,pks,'o');
%         dist_val = min(dist_new(1,locs)./2.0);
% 
%         title(strcat(filename), 'interpreter', 'latex');
% 
%         xlabel('Distance (m)')
%         ylabel('Z - Magnitude spectrum of X*Y')
%     end
% 
%     disp('time elapsed')
%     timeElapsed = toc
% 
%     saveas(figure1,strcat(outputResultDir, '/figures', filename, '.jpg'))
%     disp('The index point and last point for dista array')
%     index_count, last_index_count
end
