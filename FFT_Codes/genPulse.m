
fnx = @(x,freq, snr, att) awgn(sin(2 * pi * freq * x).* exp(att * abs(x) + normrnd(0, 0.01, size(x))), snr, 'measured');

total_duration = 60; % seconds
fs = 8000; % sampling rate 
t = 0:1/fs:total_duration;

heart_rate = 59; % 59 bpm 
heartbeat_interval = 60/heart_rate;

d2 = [];
% make the heartbeat interval varies (offset) 
% also make the gain varies for each pulse (gain)
for d = 0:heartbeat_interval:total_duration
    offset = floor(normrnd(d, 0.1)/(1/fs)) * (1/fs);
    gain = normrnd(1, 0.2); 
    d2 = [d2; offset  gain];
end

% the last three params are 
fn = 20; % 1) freq of tone (e.g. 20 Hz)
noise = 10; % 2) snr of gausian white noise to be added (e.g 10) 
atn = -30; % 3) attuation factor for the sin wave (e.g. -30) 
y = pulstran(t, d2, fnx, fn , noise , atn);


figure;
hold on;
plot(t,y);
plot(d2(:, 1), zeros(size(d2)), 'x');
xlim([0 10]);
xlabel('Time (s)')
ylabel('Waveform')

%sound(y, fs);
filename = "synth-" + "hr" + heart_rate + "-snr" + noise + "-atn" + atn;
audiowrite(strcat(filename, '.wav'), y, fs); 

fid = fopen(strcat(filename, '.txt'),'wt');
for ii = 1:size(d2,1)
    fprintf(fid,'%g\t',d2(ii,1));
    fprintf(fid,'\n');
end
fclose(fid)

