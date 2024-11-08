%%
% This code:
% 
% Visualizes raw EEG data.
% a) Applies stop ban filters at specific frequencies.
% b) Performs band-pass filtering.
% c) Computes FFT and power distribution over 30s epochs.
% d) Detects sleep stages by analyzing frequency band power in each epoch 
% by comparing them with defined thresholds for each sleep stage.
% e) Normalizes data using a manual Box-Cox transformation.
%%

close all;
clear all;

%% a) Preparation of removing EEG artefacts at f=0.99 Hz and f=2 Hz - Visualization of Unmodified EEG

data = load('data_eeg.mat');

%Extract sampling frequency (fs) and EEG signal
fs = data.dati.fs;   %Sampling frequency
eeg = data.dati.eeg; %EEG signal

%Calculate the total duration of the signal
N = length(eeg);     %Number of samples in the signal
t = (0:N-1) / fs;    %Time vector

%Raw EEG
figure;
plot(t, eeg, 'b', 'LineWidth', 1.2); %Blue color and thicker line for better visibility
title('Raw EEG Signal');
xlabel('Time [s]');
ylabel('Amplitude');
grid on; 

%The figure shows that there is a lot of noise in the signal, making it unusable as-it is

%% Filtering at 0.99 Hz Frequency (IIR)
load("coefficients_SOS_G.mat");%loading coefficients of IIR1

eeg_iir = sosfilt(SOS, eeg);

figure;

%Plot raw EEG and filtered EEG in subplots with a unified x-axis
sb(1) = subplot(2, 1, 1);
plot(t, eeg, 'b', 'LineWidth', 1); %Original EEG in blue
title('Raw EEG Signal');
ylabel('Amplitude');
grid on;

sb(2) = subplot(2, 1, 2);
plot(t, eeg_iir, 'r', 'LineWidth', 1); %Filtered EEG in red
title('Filtered EEG Signal (0.99 Hz)');
xlabel('Time [s]');
ylabel('Amplitude');
grid on;

% Link the x-axes of the subplots for synchronized zoom and pan
linkaxes(sb, 'x');

%% Filtering at 2 Hz Frequency (IIR)
load("coefficients_SOS2_G2.mat");
eeg_iir2 = sosfilt(SOS2, eeg_iir);
%% 

% Spectrum of the signal without artefacts
N = length(eeg_iir2);
frequencies_filt_iir = (0:N-1) * (fs / N); % Frequency vector for the filtered signal
IIR_filt = fft(eeg_iir2);
P2_filt = abs(IIR_filt / N); % Double-sided spectrum
P1_filt = P2_filt(1:N/2+1); % Single-sided spectrum
P1_filt(2:end-1) = 2 * P1_filt(2:end-1);
f_filt = frequencies_filt_iir(1:N/2+1);

% Plot the spectrum of the filtered signal
figure;
plot(f_filt, P1_filt, 'g', 'LineWidth', 1.2); % Green color and thicker line for visibility
xlabel('Frequency (Hz)');
ylabel('|P1(filt)(f)|');
title('Spectrum of EEG Signal without artifacts');
grid on; % Add grid for better readability

%% b) Band-pass filter between 0.5 Hz and 30 Hz
low_cutoff = 0.5;      %Lower cutoff frequency (Hz)
high_cutoff = 30;    %Upper cutoff frequency (Hz)
[b, a] = butter(4, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass'); %4th-order filter
filtered_eeg = filtfilt(b, a, eeg_iir2); %Apply the filter

figure; plot(t,filtered_eeg); title('2x Filtered EEG'); xlabel('time [s]'); ylabel('Amplitude');

% %% Compressed Spectral Array (CSA)
% it takes a long time to run so we just comment it 
%%%uncomment until "end of comment" to visualize csa
%%%(estimating time 2min)
% window_length = fs * 5;  %5-second window
% shift = window_length * 0.5; % 50% overlap
% nfft = 2048;% FFT size
% noverlap = fs * 0.5;% Overlap for Welch
% 
% i = 1;  %Starting index
% counter = 1;
% N = length(filtered_eeg);
% 
% while i + window_length < N
%     eeg_segment = detrend(filtered_eeg(i:i + window_length)); % Detrend each window
%     [PSD, f] = pwelch(eeg_segment, hamming(fs), noverlap, nfft, fs);
%     PSD = PSD(f < 30); % Limit to frequencies below 50 Hz
%     CSA(counter, :) = PSD;
%     time(counter) = t(i + window_length / 2) / 60; % Time in minutes
%     i = i + shift;
%     counter = counter + 1;
% end
% 
% % CSA visualization adjustment
% [TimeGrid, FreqGrid] = meshgrid(time, f(f < 30));
% figure;
% waterfall(TimeGrid, FreqGrid, CSA');
% xlabel('Time [minutes]');
% ylabel('Frequency [Hz]');
% zlabel('Power [\muV^{2}/Hz]');
% title('Compressed Spectral Array (CSA)');


%EEG_normalized = filtered_eeg;

%% c) Epoch and Frequency Band Analysis
epoch_length = 30 * fs;  %30-second epoch in samples (3000 samples per epoch)
num_epochs = floor(length(filtered_eeg) / epoch_length);  %Number of full 30-second epochs

% d) Frequency band definitions (in Hz)
bands = {
    'Delta', 1, 4;
    'Theta', 4, 7;
    'Alpha', 7, 10;
    'Beta', 10, 30;
};

%Preallocate cell arrays to store FFT results and power percentages per band
fft_results = cell(1, num_epochs);       %To store FFT result for each epoch
band_percentages = cell(1, num_epochs);  %To store power percentage for each band
sleep_stages = cell(1, num_epochs);      %To store detected sleep stage per epoch

% Preallocate arrays to store power percentages for Alpha, Beta, Theta, and Delta bands
alpha_power = zeros(1, num_epochs);
beta_power = zeros(1, num_epochs);
theta_power = zeros(1, num_epochs);
delta_power = zeros(1, num_epochs);


%% Adjust Thresholds for Better Detection of Sleep Stages

alpha_threshold = 30;  % Alpha percentage threshold for Wakefulness
theta_threshold = 2;  % Theta percentage threshold for N1
delta_threshold_N3 = 63;  % Minimum delta power for N3
beta_threshold_REM = 3.3;  % Minimum beta power for REM
for i = 1:num_epochs
    % Extract the current epoch
    epoch_data = filtered_eeg((i-1)*epoch_length + 1 : i*epoch_length);
    
    % Compute FFT
    fft_result = fft(epoch_data);
    N = length(epoch_data);
    freqs = (0:N-1) * (fs/N);  % Frequency axis for FFT
    half_freqs = freqs(1:floor(N/2));  % Only positive frequencies
    half_fft = abs(fft_result(1:floor(N/2)));  % FFT magnitude

    % Total power for normalization
    total_power = sum(half_fft.^2);

    % Calculate the percentage power for each band
    band_percent = zeros(size(bands, 1), 1);  % To store percentage for each band
    for b = 1:size(bands, 1)
        band_idx = find(half_freqs >= bands{b, 2} & half_freqs <= bands{b, 3});
        band_power = sum(half_fft(band_idx).^2);
        band_percent(b) = (band_power / total_power) * 100;
    end

    % Store FFT and band percentage
    fft_results{i} = fft_result;
    band_percentages{i} = band_percent;

    %% Detect Sleep Stage
    alpha_percentage = band_percent(3);
    theta_percentage = band_percent(2);
    delta_percentage = band_percent(1);
    beta_percentage = band_percent(4);

    % Store power percentage for each band
    alpha_power(i) = alpha_percentage;
    beta_power(i) = beta_percentage;
    theta_power(i) = theta_percentage;
    delta_power(i) = delta_percentage;
    
    % Define sleep stages
    if alpha_percentage > alpha_threshold
        sleep_stage = 'Wakefulness';
    elseif theta_percentage > theta_threshold && alpha_percentage < alpha_threshold && delta_percentage < 5
        sleep_stage = 'N1';
    elseif theta_percentage > theta_threshold && delta_percentage < delta_threshold_N3 && alpha_percentage < alpha_threshold
        sleep_stage = 'N2';
    elseif delta_percentage > delta_threshold_N3 && alpha_percentage < alpha_threshold && beta_percentage < beta_threshold_REM
    sleep_stage = 'N3';
    elseif    beta_percentage > beta_threshold_REM &&theta_percentage > theta_threshold 
         sleep_stage = 'REM';
    else
        sleep_stage = 'unkown';
    end
    
    % Store detected sleep stage
    sleep_stages{i} = sleep_stage;
    %% plot percentage distribution in a chosen epoch
    if i==240 || i==700 || i==100
        figure;
        hold on;

        % Plot FFT magnitude for each frequency band with specific colors
        %colors = ['0.3, 0.3, 0.6', '0.1, 0.5, 0.8', 'c', '0.9, 0.1, 0.7', 'b', 'y', 'k'];  % Colors for each band
        colors = [
            0.3, 0.3, 0.6;  % Color for Delta
            0.1, 0.5, 0.8;  % Coulor for Theta
            0, 1, 1;        % Cyan for Alpha
            0.9, 0.1, 0.7;  % Magenta for Beta
        ];

        for b = 1:size(bands, 1)
            band_idx = find(half_freqs >= bands{b, 2} & half_freqs <= bands{b, 3});
            plot(half_freqs(band_idx), half_fft(band_idx), 'Color', colors(b, :), 'LineWidth', 2); 
        end


        % Add legend and labels
        legend(bands{:, 1});
        title(['Epoch ' num2str(i) ' - Sleep Stage: ' sleep_stages{i}]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');

        % Add percentages of each frequency band in the title
        band_labels = strcat(bands(:, 1), ': ', num2str(band_percent, '%.2f%%'));
        annotation('textbox', [0.15, 0.7, 0.1, 0.1], 'String', band_labels, 'FitBoxToText', 'on');

        hold off;
    end

end

%% Visualization of different EEG waves
% Create the visualization figure for EEG wave bands
figure;

% Loop through each band to display its filtered signal
for i = 1:size(bands, 1)
    band_name = bands{i, 1};
    
    % Designing a bandpass filter for each band
    [b, a] = butter(2, [bands{i, 2}, bands{i, 3}] / (fs / 2), 'bandpass');
    signal = filtfilt(b, a, filtered_eeg); % Application of the filter
    filtered_signals.(band_name) = signal; % Filtered signal storage
    
    t = (0:length(signal)-1) / fs;
    
    subplot(size(bands, 1), 1, i);
    plot(t, signal, 'LineWidth', 1.5, 'Color', colors(i, :)); 
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(['EEG Signal - Wave ', band_name]);
    grid on;
end

sgtitle('Visualization of the Different EEG Waves');

%% Plot the power in Alpha, Beta, Theta, and Delta bands across all epochs
figure;
hold on;

% Plot each frequency band's power across epochs
plot(alpha_power, '-*', 'DisplayName', 'Alpha Power (%)', 'Color', 'c');
plot(beta_power, '-*', 'DisplayName', 'Beta Power (%)', 'Color', '0.9, 0.1, 0.7');
plot(theta_power, '-*', 'DisplayName', 'Theta Power (%)', 'Color', '0.1, 0.5, 0.8');
plot(delta_power, '-*', 'DisplayName', 'Delta Power (%)', 'Color', '0.3, 0.3, 0.6');

% Add legend and labels
legend;
xlabel('Epoch');
ylabel('Power Percentage (%)');
title('Power in Alpha, Beta, Theta, and Delta Frequency Bands Across Epochs');
hold off;

%% Calculate and Print Sleep Stage Durations

% Define the unique sleep stages
unique_stages = unique(sleep_stages);

% Initialize a structure to store durations for each stage
sleep_durations = struct();

% Calculate the duration of each sleep stage
for s = 1:length(unique_stages)
    stage = unique_stages{s};
    
    % Count the number of epochs corresponding to the current sleep stage
    num_epochs_in_stage = sum(strcmp(sleep_stages, stage));
    
    % Total duration in seconds (each epoch is 30 seconds)
    total_duration_sec = num_epochs_in_stage * 30;
    
    % Convert to minutes and hours
    total_duration_min = total_duration_sec / 60;
    total_duration_hr = total_duration_min / 60;
    
    % Store the durations in the structure
    sleep_durations.(stage) = struct('Seconds', total_duration_sec, 'Minutes', total_duration_min, 'Hours', total_duration_hr);
    
    % Print the results
    disp(['Sleep Stage: ', stage]);
    disp(['  Duration: ', num2str(total_duration_sec), ' seconds']);
    disp(['  Duration: ', num2str(total_duration_min), ' minutes']);
    disp(['  Duration: ', num2str(total_duration_hr), ' hours']);
end

%% Hypnogram

function hypnogram(sleep_stages, window_size)
    unique_stages = {'Wakefulness', 'REM', 'NREM'};
    stage_values = zeros(1, length(sleep_stages));
    
    % Color for each stage
    colors = [
        0.9, 0.1, 0.7;  % Pink for Wakefulness
        1.0, 0.6, 0.1; % Orange for REM
        0.2, 0.5, 1.0; % Blue for NREM
    ];

    % Conversion of stages sleep in group (N1, N2, N3, N4 regroup in NREM)
    for i = 1:length(sleep_stages)
        if ismember(sleep_stages{i}, {'N1', 'N2', 'N3', 'N4'})
            stage_values(i) = find(strcmp(unique_stages, 'NREM'));
        elseif strcmp(sleep_stages{i}, 'REM')
            stage_values(i) = find(strcmp(unique_stages, 'REM'));
        else
            stage_values(i) = find(strcmp(unique_stages, sleep_stages{i}));
        end
    end

    % Application of the sliding majority window to smooth the hypnogram
    smoothed_values = stage_values;
    window_by_three = floor(window_size / 2.6);
    for i = 1:length(stage_values)
        % Windows terminals 
        start_idx = max(1, i - window_by_three);
        end_idx = min(length(stage_values), i + window_by_three);
        
        % Extraction of the values
        window_values = stage_values(start_idx:end_idx);
        
        % Calculation of the majority stage in the window
        majority_stage = mode(window_values);
        
        % Replace values corresponding to minority stages by the majority stage
        smoothed_values(i) = majority_stage;
    end

    % Tracing
    figure;
    plot(1:length(smoothed_values), smoothed_values, '-', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5]);
    
    hold on;
    for i = 1:length(unique_stages)
        plot(find(smoothed_values == i), smoothed_values(smoothed_values == i), 'o', 'LineWidth', 2, ...
             'Color', colors(i, :));
    end
    hold off;
    
    % Reversing the order of the y-axis so that wakefulness appears at the top
    ylim([0.5, length(unique_stages) + 0.5]);
    set(gca, 'YDir', 'reverse');  % Reversing of the y-axis
    
    yticks(1:length(unique_stages));
    yticklabels(unique_stages);
    
    xlabel('Epochs','FontSize',14);
    ylabel('Sleep Stage','FontSize',14);
    title('Hypnogram','FontSize',18);
    grid on;
end

hypnogram(sleep_stages, 20)

%% e) Box-Cox Transformation
eeg_shifted = filtered_eeg - min(filtered_eeg) + 1; % Make data positive
lambda = 0.3; %Lambda value for Box-Cox transformation
%Add a small constant before the transformation
small_constant = 0.01; %Adjust based on data range
eeg_adjusted = eeg_shifted + small_constant;

%Box-Cox Transformation
if lambda == 0
    transformed_eeg = log(eeg_adjusted);
else
    transformed_eeg = (eeg_adjusted.^lambda - 1) / lambda;
end

disp(['Lambda Parameter for Box-Cox: ', num2str(lambda)]);

%Visualization of Distribution after Box-Cox Transformation
figure;
histogram(transformed_eeg, 50); %Histogram with 50 bins
hold on;
mu_transformed = mean(transformed_eeg);
sigma_transformed = std(transformed_eeg);
x_values_transformed = linspace(min(transformed_eeg), max(transformed_eeg), 100);
y_values_transformed = (1 / (sigma_transformed * sqrt(2 * pi))) * exp(-0.5 * ((x_values_transformed - mu_transformed) / sigma_transformed).^2);
y_values_transformed = y_values_transformed * numel(transformed_eeg) * (max(transformed_eeg) - min(transformed_eeg)) / 50; % Normalize for histogram
plot(x_values_transformed, y_values_transformed, 'r', 'LineWidth', 2); % Red line for normal distribution curve
xlabel('Transformed EEG Signal Value');
ylabel('Frequency');
title('Histogram and Normal Distribution after Box-Cox Transformation');
hold off;

%Normality Test, Skewness and Kurtosis
skewness_transformed = skewness(transformed_eeg);  % Skewness calculation after transformation
skewness_original = skewness(filtered_eeg);        % Skewness before transformation
kurtosis_transformed = kurtosis(transformed_eeg) - 3;  % Kurtosis (flattening) after transformation
kurtosis_original = kurtosis(filtered_eeg) - 3;        % Kurtosis before transformation

%Display skewness and kurtosis results
disp(['Skewness before Box-Cox Transformation: ', num2str(skewness_original)]);
disp(['Skewness after Box-Cox Transformation: ', num2str(skewness_transformed)]);
disp(['Kurtosis before Box-Cox Transformation: ', num2str(kurtosis_original)]);
disp(['Kurtosis after Box-Cox Transformation: ', num2str(kurtosis_transformed)]);
