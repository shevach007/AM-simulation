%creating a simulation of Am modulation - creating waveform from bits,
%modulating the waveform, sending several signals through the same medium, 
%preforming de-modulation and passing the signal through LPF and recover the bits  
clear; close all; clc;

T = 20;
samples = 10000;
t = linspace(0, T, samples);
%data bits that we wish to modulate and then recover, you can add any
%10-bit vector to the matrix and the program will modulate, plot, and
%recover the waveform (remember to add the carrier frequency at line 45)
signals_bits = [0, 1, 0, 1, 1, 0, 0, 1, 1, 1; %bit vector 1
                1, 1, 0, 0, 1, 0, 1, 0, 1, 1; %bit vector 2
                0, 0, 1, 1, 0, 0, 1, 1, 0, 0; %bit vector 3
                1, 0, 1, 1, 0, 1, 1, 1, 0, 0]; %bit vector 4

num_of_signals = size(signals_bits, 1);
length_of_signal = size(signals_bits, 2);
time_per_symbol = samples/length(signals_bits(1,:));

%-------map bits to time vector---------
x_signals = zeros(num_of_signals, samples);
for j=1:num_of_signals
    for i=1:length(signals_bits(1,:))
        start_idx = round((i-1)*time_per_symbol) + 1;
        end_idx = round(i*time_per_symbol);
        x_signals(j, start_idx:end_idx) = signals_bits(j, i);
    end
end
figure;
title("signals in time domain")
for i = 1:num_of_signals
    subplot(num_of_signals, 1, i);
    plot(t, x_signals(i,:));
    title(['Signal ', num2str(i)]);
    ylim([-0.2, 1.2]); % Set y-axis limits for consistency
    for j=1:length_of_signal
        txt = signals_bits(i, j);
        txt = text(j*2 - 1, 0, num2str(txt));
        txt.Color = [1, 0, 0];
    end
end
sgtitle('unmodulated signals in time domain')
%--------create carriers in different frequencies---------
freq_vec = [5, 20, 40, 60];% add frequencies here is you add signal
carriers = zeros(num_of_signals, samples);
for j=1:num_of_signals
    carriers(j,:) = cos(2*pi*freq_vec(j)*t);
end

%-------plot carriers-----------
figure;
for i = 1:num_of_signals
    subplot(num_of_signals, 1, i);
    plot(t, carriers(i,:));
    title(['carrier ', num2str(i)]);
    ylim([-1.2, 1.2]); 
end
sgtitle('carrier signals')

%--------create modulated signals--------
modulated_signals = zeros(num_of_signals, samples);
for i = 1:num_of_signals
    modulated_signals(i,:) = x_signals(i,:) .* carriers(i,:);
end
figure;
for i = 1:num_of_signals
    subplot(num_of_signals, 1, i);
    plot(t, modulated_signals(i,:));
    title(['modulated signal ', num2str(i)]);
    ylim([-1.2, 1.2]); 
    for j=1:length_of_signal
        txt = signals_bits(i, j);
        txt = text(j*2 - 1,0,num2str(txt));
        txt.Color = [1, 0, 0];
    end
end
sgtitle('modulated signals in time domain')

%-------sum all the signals (simulating transmition throgh a common
%medium)-------

sum_all_trans_signals = zeros(1, samples);
for i = 1:num_of_signals
    sum_all_trans_signals(1,:) = modulated_signals(i, :) + sum_all_trans_signals(1,:);
end

figure;
plot(t, sum_all_trans_signals);
title("all transmited signals through common medium");

%-----demodulating signals, if we multiply the modulated signal with the 
% same carrier frequency we can rocover the origital bits, the data will
% be at 0 Hz------
demodulated_signals = zeros(3, samples);
for i = 1:num_of_signals
    demodulated_signals(i, :) = sum_all_trans_signals(1,:).* carriers(i, :);
end
figure;
for i = 1:num_of_signals
    subplot(num_of_signals, 1, i);
    plot(t, demodulated_signals(i,:));
    title(['de-modulated signal ', num2str(i)]);
    ylim([-1.2, 1.2]); 
end
%showing the signals in frequency domain
sum_all_trans_signals_freq = fft(sum_all_trans_signals);
fs = 1 / (t(2) - t(1)); 
f = linspace(0, fs/2, floor(samples/2) + 1);

modulated_signals_freq = zeros(3, samples);
for i = 1:num_of_signals
    modulated_signals_freq(i,:) = fft(modulated_signals(i,:));
end

% plot(f, abs(sum_all_trans_signals_freq(1:floor(samples/2) + 1)))
% title("transmited signals in frequency domain");

figure;
hold on
for i = 1:num_of_signals
    hz_var = freq_vec(i);
    plot(f, abs(modulated_signals_freq(i,1:floor(samples/2) + 1)), 'DisplayName' ,[num2str(hz_var), ' Hz carrier signal']);
    title('all signals in frequency domain');
    xlabel("frequency (Hz)");
    ylabel("Amplitude")
    %ylim([-1.2, 1.2]); 
end
hold off
legend
%constructing the LPF to recover the actual data bits we want at 0 Hz
cutoff_freq = 1; %cutoff frequency in Hz
filter_order = 100; %order of the filter
[b, a] = fir1(filter_order, cutoff_freq / (0.5 * fs), 'low');

filtered_signals = zeros(num_of_signals, samples);
for i = 1:num_of_signals
    filtered_signals(i,:) = filter(b, a, demodulated_signals(i,:));
end

figure;
for i = 1:num_of_signals
    subplot(num_of_signals, 1, i);
    plot(t, filtered_signals(i,:));
    title(['filtered signals ', num2str(i)]);
    ylim([-1.2, 1.2]); 
    for j=1:length_of_signal
        txt = signals_bits(i, j);
        txt = text(j*2 - 1,0,num2str(txt));
        txt.Color = [1, 0, 0];
    end
end
sgtitle('all the filtered signals')

figure;
for i = 1:num_of_signals
    subplot(num_of_signals, 2, 2*i-1);
    plot(t, x_signals(i,:));
    title(['Signal ', num2str(i)]);
    ylim([-0.2, 1.2]); 
    for j=1:length_of_signal
        txt = signals_bits(i, j);
        txt = text(j*2 - 1,0,num2str(txt));
        txt.Color = [1, 0, 0];
    end
    subplot(num_of_signals, 2, 2*i);
    plot(t, filtered_signals(i,:));
    title(['Recovered Signal ', num2str(i)]);
    ylim([-0.2, 1.2]); 
    for j=1:length_of_signal
        txt = signals_bits(i, j);
        txt = text(j*2 - 1,0,num2str(txt));
        txt.Color = [1, 0, 0];
    end
end
sgtitle('Comparing the original waveform and the recovered waveform')

