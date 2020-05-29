% S : signal
% Fs : sampling rate

function [f,P1] = plotSTFT(S, Fs, flag)
    if nargin < 3  
        flag = false; 
    end
    L = length(S); % num samples
    f = Fs*(0:(L/2))/L;
    %f = f/1000; 

    Y = stft(S);
   
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
   
    P1(2:end-1) = 2*P1(2:end-1);
    
    if (flag)
        figure; plot(f,P1)
        title('Single-Sided Amplitude Spectrum of S(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
    end
end