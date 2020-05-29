date = "2020-05-26"
upper_bound = 21
path = "/Users/ethanzh/UT/Research/Ethan-Mikie-Recordings/breathing/"

all_wavs = strcat(path, date, "/*.wav")

wavfiles = dir(all_wavs)
for file = wavfiles'
    disp(file.name)
    
    [data, Fs] = audioread(strcat(path, date, "/", file.name));

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
    spectrogram(data, 2048, [], [], Fs, 'yaxis');
    title(file.name);
    
    % if this is a chirp, set upper bound = 21, else set to 4
    is_chirp = contains(file.name, "chirp"); 
    if is_chirp == 0
        ylim([0 4]); % limit y axis to breathing range (0kHz-4kHz)
    else
        ylim([0 21]); % limit y axis to max frequency
    end
    
    disp(file.name)
    png_name = file.name(1:end-4)

    saved_image_name = strcat(path, date, "/", png_name, ".png");
    disp(saved_image_name);

    saveas(gcf, saved_image_name);
end