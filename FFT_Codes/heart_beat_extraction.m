
filepathPrefix = strcat('/Users/profhan/Downloads/Wireless_Networking_UT-master/');
filepathRX = strcat(filepathPrefix, 'Passive_Recording');
%filename = strcat('heart_on_skin_2min');
filename = strcat('left_wrist_no_breath');


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
n = 100; % downsampling by 100 (411 Hz)
Rx = downsample(Rx, n);
fs = Fs/n; % decreased sampling rate
Rx = Rx(1:floor(size(Rx)/fs) * fs);

figure;
plot(0:size(Rx)-1, Rx);
xlim([0 fs*3]); % plot just 10 seconds 

% Step 1-2.5 (?) Apply Blackman window (sample size of 256??). 
% Should it be here before K-means?
% Paper says "The power spectrum of the acoustic pulse signal with a downsampled frequency of 210 Hz 
% is calculated in the algorithm using a Blackman window of 32 samples (approximately 150 milliseconds) 
% with an overlap of 50% between successive frames."
M = length(Rx); 
w = .42-.5*cos(2*pi*(0:M-1)/(M-1))+.08*cos(4*pi*(0:M-1)/(M-1));
Rx = Rx.*w;

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
                Rx_trimmed = [Rx_trimmed y_j];
            end 
        else % do not skip
            Rx_trimmed = [Rx_trimmed y_j];
        end 
    end    
end

figure;
plot(Rx_trimmed);
xlim([0 fs * 3]);

% Step 2. S1 sound extraction from clean acoustic pulse signal
% Step 2-1. Short time Fourier Transform
[f, psd] = plotSTFT(Rx_trimmed, fs); % psd single-sided amplitue
p_max = max(psd); % 2.5e-24 (seems too small)
figure;
plot(f, psd);

