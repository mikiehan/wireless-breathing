

filepathPrefix = strcat('/Users/profhan/Downloads/Wireless_Networking_UT-master/');
filepathRX = strcat(filepathPrefix, 'Passive_Recording');
filename = strcat('heart_on_skin_2min');

[Rx,fs] = audioread(strcat(filepathRX,'/',filename,'.wav'));
 Rx = Rx(2*fs:length(Rx)); % remove a few initial points

% Step 1-1. Apply 5th Order Butterworth Low Pass Filter 
fc = 25; % cut off frequency in Hz
Wn = fc/(fs/2); % Fs/2 is the nyquist frequency
[b,a] = butter(5,Wn,'low');

Rx = filter(b,a,Rx);
fs, Wn

figure;
spectrogram(Rx, 2048, [], [], fs, 'yaxis');
ylim([0 0.1]);
%xlim([0.2 1]);

% Step 1-2. Downsampling by n
n = 100; % downsampling by 100
Rx = downsample(Rx, n);
fs = fs/n;
figure;
plot(0:size(Rx)-1, Rx);
xlim([0 441*10]); % plot just 10 seconds 

% Step 1-3. K-means method to remove artifacts
