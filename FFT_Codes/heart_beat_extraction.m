
filepathPrefix = strcat('/Users/profhan/Downloads/Wireless_Networking_UT-master/');
filepathRX = strcat(filepathPrefix, 'Passive_Recording');
filename = strcat('heart_on_skin_external_mic_8000fs');
%filename = strcat('heart_on_skin_external_mic_8000fs_2');
%filename = strcat('heart_on_skin_2min');
%filename = strcat('left_wrist_8000fs_1');
%filename = strcat('left_wrist_no_breath');
%filename = strcat('heart_1cm_external_mic_8000fs');

[Rx,Fs] = audioread(strcat(filepathRX,'/',filename,'.wav'));
Rx = Rx(2*Fs:length(Rx)); % remove a few initial points

% Step 1-1. Apply 5th Order Butterworth Low Pass Filter 
fc = 25; % cut off frequency in Hz
Wn = fc/(Fs/2); % Fs/2 is the nyquist frequency
[b,a] = butter(5,Wn,'low');

Rx = filter(b,a,Rx);
Fs, Wn

figure;
spectrogram(Rx, Fs, [], [], Fs, 'yaxis');
ylim([0 0.1]);
%xlim([0.2 1]);

% Step 1-2. Downsampling by n
n = 40; % downsampling by 40 (8000/40 = 200 Hz)
Rx = downsample(Rx, n);
fs = Fs/n; % decreased sampling rate
Rx = Rx(1:floor(size(Rx)/fs) * fs);

% spectrogram after downsampling
figure;
spectrogram(Rx, fs, [], [], fs, 'yaxis');

figure;
plot(0:size(Rx)-1, Rx);
xlim([0 fs*3]); % plot just 10 seconds 

% Step 1-3. K-means method to remove artifacts
block = 5; % signal block of 5 seconds
num_blocks = floor(size(Rx)/fs/block);
Rx = Rx(1:num_blocks * fs * block);
Rx_trimmed = [];
for i = 1:num_blocks % for each block
    y = Rx((i-1) * fs * block + 1 : i * fs * block);
    % divide into 5 equal parts
    max_amp_x = zeros(5,1);
    std_y = zeros(5,1);
    for j = 1:5
        y_j = y((j-1) * fs + 1: j * fs); % one signal part
        max_amp_x(j, 1) = max(y_j); % max_amplitude as x coord
        std_y(j, 1) = std(y_j); % standard deviation as y coord        
    end
    [idx,centroids] = kmeans([max_amp_x std_y], 2);
    x_d = abs(centroids(1,1) - centroids(2,1))/min(centroids(1,1), centroids(2,1)); % delta x
    C = 1;
    if (centroids(1,2) > centroids(2,2)) %
        C = 2;
    end
    Sn = zeros(5,1);
    for j = 1:5
        y_j = y((j-1) * fs + 1: j * fs);
        if(x_d >= 0.5)
            if(idx(j) ~= C) % skip this signal part
                Sn(j) = 1; 
            else % do not skip 
                Rx_trimmed = [Rx_trimmed; y_j];
            end 
        else % do not skip
            Rx_trimmed = [Rx_trimmed; y_j];
        end 
    end    
end

Rx_trimmed = Rx_trimmed';
total_samples = length(Rx_trimmed);
duration = total_samples/fs;

figure;
plot(Rx_trimmed);
findpeaks(Rx_trimmed,'Annotate','extents');
%xlim([0 fs * 5]);

total_samples = length(Rx_trimmed);
duration = total_samples/fs;
max_num_peaks = duration/0.30; % min beat-to-beat interval is 300 milisec 
[pks, locs]= findpeaks(Rx_trimmed);
[sorted_pks, sorted_idx] = sort(pks, 'descend');
sorted_pks_s = sorted_pks(1:max_num_peaks);
sorted_idx_s = sorted_idx(1:max_num_peaks);

figure;
hold on;
plot([1:total_samples], Rx_trimmed);
plot(locs(sorted_idx_s), sorted_pks_s, 'o');
xlim([0 fs*10]); % 10 sec

figure;
hold on;
plot([1:total_samples]/fs, Rx_trimmed);
plot(locs(sorted_idx_s)/fs, sorted_pks_s, 'o');
xlim([0 5]);
%plot(sorted_pks,'o'); 

% spectrogram after K-means
figure;
spectrogram(Rx_trimmed, fs, [], [], fs, 'yaxis');

% Step 2. S1 sound extraction from clean acoustic pulse signal
% Step 2-1. Short time Fourier Transform with blackman window of 32 samples
% with 50 percent overlab 
% (with 200Hz this provides 160 millisec time resolution)
% (paper suggests 150 millisec time resolution)
M = 32; 
win = blackman(M, 'periodic'); % blackman window
L = 200; % fft length
[ss, ff, tt] = stft(Rx_trimmed, fs, 'Window', win,'OverlapLength',(M/2), 'FFTLength' , L);
% make single sided amplitude
ff = ff(L/2:end); 
ss = abs(ss/L);
ss = ss(L/2:end,:);
ss(2:end-1,:) = 2 * ss(2:end-1, :); 
% translate into dB
psd_db = mag2db(ss); 
[max_mean, freq_row] = max(mean(psd_db,2)); % calculate avg psd across time and select freq with max
freq_selected = ff(freq_row); 

max_psd = max(max(psd_db));

% spectrogram with blackman window 
figure;
spectrogram(Rx_trimmed, win, M/2, L, fs, 'yaxis', 'onesided');
xlim([0 5]);

% Step 2-2. Extract grids with P ≥ Pmax − Pt, where Pt∈[5,10]dB such that m∈[3,14]
% where m is the number of groupings in 4 sec duration (in actual paper it
% is 5 seconds)
block = 4; 
time_grid = (M/2)/fs;
n_blocks = floor(tt(end)/block);
m_min = ceil(40/60*block); % 3
m_max = ceil(200/60*block); % 14

groups_all = [];
for db_thresh = 5:12
    b = 1;
    valid = true;
    groups_thresh = [];
    while (valid && b <= n_blocks)
        groups_block = [];
        % the interval is from [(b-1)* block , b * block] seconds
        t_start = (b-1) * block/time_grid + 1;
        t_end = b * block/time_grid;
       
        time_grids_idx = find(psd_db(freq_row,t_start:t_end) >= max_psd - db_thresh); 
        time_grids_idx(end + 1) = -1; % adds new endpoint to very end so code picks up end of last group of consecutive values
        diff_idx = find(diff(time_grids_idx)~=1);
        if(~isempty(diff_idx))
            i_start = 1;
            for i = diff_idx(1:end)
                if (i_start ~= i) 
                    groups_block = [groups_block ; db_thresh (time_grids_idx(i_start) + t_start - 1) (time_grids_idx(i) + t_start - 1)];
                end
                i_start = i + 1;
            end 
        end
        n_groups = size(groups_block, 1);
        if(n_groups >= m_min && n_groups <= m_max)
            groups_thresh = [groups_thresh; groups_block];
        else
            disp(groups_block);
            groups_thresh = []; % clear out groups for this threshold
            valid = false;
        end
        b = b + 1; % increment b 
    end
    if(valid)
        groups_all = [groups_all ; groups_thresh];
    end
end 
% % % 
% % % 
% % % 
% % % figure;
% % % plot(f, psd);
% % % % Extract grids with P ≥ Pmax − Pt, where Pt∈[5,10]dB such that m∈[4,17].
% % % % How to get dB? 
% % % 
% % % % We know HR ranges from 40 to 200 bpm which corresponds to 
% % % % 1500 milisec - 300 milisec beat-to-beat interval 
% % % min_b2b_interval = 300; 
