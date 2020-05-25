[data, Fs] = audioread("/Users/ethanzh/UT/Research/Ethan-Mikie-Recordings/breathing/20cm.wav");

% get the length of the audio file
N = length(data)
audio_length = N/Fs

data = data(1:Fs*audio_length);

% PARAMETERS
% spectrogram(song, windowSize, windowOverlap, freqRange, fs, 'yaxis');

% Larger Window Size value increases frequency resolution
% Smaller Window Size value increases time resolution
% Specify a Frequency Range to be calculated for using the Goertzel function
% Specify which axis to put frequency

%% USAGE
spectrogram(data, 256, [], [], Fs, 'yaxis');
