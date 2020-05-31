
filepathPrefix = strcat('/Users/profhan/Downloads/Wireless_Networking_UT-master/');
filepathRX = strcat(filepathPrefix, 'Passive_Recording');
%filename = strcat('heart_on_skin_external_mic_8000fs');
%filename = strcat('heart_on_skin_2min');
filename = strcat('left_wrist_8000fs_1');
%filename = strcat('left_wrist_no_breath');


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
figure;
plot(Rx_trimmed);
xlim([0 fs * 5]);

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
ss = ss(1:floor(L/2)+1,:);
ss(2:end-1,:) = 2 * ss(2:end-1, :); 
% translate into dB
psd_db = mag2db(ss); 
max_psd = max(max(psd_db));

% spectrogram with blackman window 
figure;
spectrogram(Rx_trimmed, win, M/2, L, fs, 'yaxis', 'onesided');
xlim([0 5]);

% Step 2-2. Extract grids with P ≥ Pmax − Pt, where Pt∈[5,10]dB such that m∈[4,17]
% where m is the number of groupings in 5 sec duration
% for db_thresh = 5:10 
%     % how to choose freq? (freq with max psd? only a few freq range?) 
%     for f = 1:length(ff) % just do all freq??
%         t_idx = find(psd_db(f,1:fs*5) >= (max_psd(f,:) - db_thresh)); % contains index's with higher psd value
%         % extract t_start_m t_end_m pairs 
%         for m = 4:17
%         end 
%     end 
% end 
% % 
% % 
% % 
% % figure;
% % plot(f, psd);
% % % Extract grids with P ≥ Pmax − Pt, where Pt∈[5,10]dB such that m∈[4,17].
% % % How to get dB? 
% % 
% % % We know HR ranges from 40 to 200 bpm which corresponds to 
% % % 1500 milisec - 300 milisec beat-to-beat interval 
% % min_b2b_interval = 300; 
