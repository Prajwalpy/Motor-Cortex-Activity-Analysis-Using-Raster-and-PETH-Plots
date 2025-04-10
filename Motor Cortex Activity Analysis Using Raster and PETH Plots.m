% Clear Workspace and Close Figures
close all; clc; clear;

%% Loading Data Files
% Load neural data 
NeuralData = load('NeuralData.mat');

% Load behavioral data
BehaviorData = load('BehaviorData.mat');

% Load All data
data = load('All.mat');
AllData = data.BehavTaskMarkers;

%% Define Common Parameters and Neuron Selections
selected_channels = {'Chan40', 'Chan30', 'Chan54'};  
selected_units = [1, 1, 1];                          
nbins = 50;                                         
bins = linspace(-1, 1, nbins+1);                     
spacing = bins(2) - bins(1);                                                



%% Question 1: Using BME_526_HW2_BehaviorData.mat

% Defining movement types and corresponding event times in order
movement_types1 = {'Index', 'Index & Middle', 'Middle', 'Thumb', 'Thumb & Index', 'Thumb, Index & Middle', 'Thumb & Middle'};
event_times = {BehaviorData.i_times, BehaviorData.im_times, BehaviorData.m_times, BehaviorData.t_times, BehaviorData.ti_times, BehaviorData.tim_times, BehaviorData.tm_times};

% Looping over each neuron
for neuron_idx = 1:length(selected_channels)
    figure;
    set(gcf, 'Position', [100, 100, 1400, 700]);
    sgtitle(['Neuron Response - ', selected_channels{neuron_idx}, ' (Unit ', num2str(selected_units(neuron_idx)), ') [BME-526-HW2-BehaviorData.mat]']);
    
    % Extract spike times for the current neuron
    spike_times = NeuralData.Channels.(selected_channels{neuron_idx});
    neuron_spike_times = spike_times(spike_times(:,2) == selected_units(neuron_idx), 3);
    
    % Loop through each movement type
    for mov_idx = 1:length(movement_types1)
        current_event_times = event_times{mov_idx};
        numTrials = length(current_event_times);
        
        %% Raster Plot
        subplot(2, length(movement_types1), mov_idx);
        hold on;
        AP_perBin = zeros(numTrials, nbins);
        all_spikes = [];
        for trial = 1:numTrials
            cueTime = current_event_times(trial);
            trial_spikes = neuron_spike_times(neuron_spike_times >= (cueTime - 1) & neuron_spike_times <= (cueTime + 1)) - cueTime;
            AP_perBin(trial, :) = histcounts(trial_spikes, bins);
            all_spikes = [all_spikes; trial_spikes(:)];
            plot(trial_spikes, trial * ones(size(trial_spikes)), 'k.', 'MarkerSize', 5);
        end
        title(movement_types1{mov_idx});
        xlabel('Time (s)');
        ylabel('Trial');
        xlim([-1,1]);
        ylim([0, numTrials+1]);
        
        %% PSTH Plot
        thresholds = [60, 35, 30];
        subplot(2, length(movement_types1), mov_idx + length(movement_types1));
        hold on;
        std_perBin = std(AP_perBin);
        counts = histcounts(all_spikes, bins) / (numTrials * spacing); 
        smoothed_counts = smoothdata(counts, 'gaussian', 5);
        

        x_vals = bins(1:end-1);
        fill_x = [x_vals, fliplr(x_vals)];
        fill_y = [smoothed_counts - 2*std_perBin, fliplr(smoothed_counts + 2*std_perBin)];
        fill(fill_x, fill_y, [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(x_vals, smoothed_counts, 'b', 'LineWidth', 2);
        yline(thresholds(neuron_idx), '--r', 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Firing Rate (Hz)');
        ylim([0, max(smoothed_counts + 2*std_perBin) + 5]);
    end
end

%% 1a. Response Fields:
% The field of the data in this assignment refers to be the motor response field, which means it's about measuring how neurons fire in response to specific finger movements, 
% like moving your index finger or combining thumb and index movements.

%1b. Individual Neuron Responses:
%Based on the raster plots and PSTHs we can get insights into the response fields of the three chosen neurons:
% Channel 40: This neuron shows the strongest response to moving the index and middle fingers together, with a peak firing rate of about 90 Hz, which is notably higher than other conditions.
% Channel 30: This neuron seems to respond most to moving the thumb, index, and middle fingers together, with a peak firing rate of about 50 Hz.
% Channel 54: This neuron appears to respond most to moving the index finger alone, with a peak firing rate of about 40 Hz.


%% 2a. Decoding Finger Movements:
%The selection of channels 40, 30, and 54 proves to be strategically advantageous for decoding finger movements due to their distinct response fields:
%    - Chan54 (Unit 1): This neuron helps identify when the index finger is moved individually, with a peak response of 35 Hz and lower activity for other movements. 
%    - Chan30 (Unit 1): This neuron effectively detects when the index finger is moved alone or in combination with the thumb, showing strong activity around 45 Hz. 
%    - Chan40 (Unit 1): This neuron reliably detects when the index and middle fingers are moved together, with a peak response of 90 Hz, which is significantly higher than for other actions.
    
% 2b. Action Potential Firing Rate Thresholds:
%    - Chan54:  30 Hz threshold is good for detecting index finger movement alone, as its peak is 35 Hz, above other movements (20-25 Hz).
%    - Chan30:  Threshold of about 35 Hz works to catch index-related movements, as its peaks for these are higher (40-45 Hz) than others (around 20-35 Hz).
%    - Chan40:  60 Hz threshold is ideal for spotting when you move index and middle fingers together, as its peak is 90 Hz, well above others (up to 45 Hz).
% Each neuron’s threshold differs because their baselines and peak responses vary.



%% ------------------------------- EXTRA CREDITS---------------------------------------------------

%% Question 1.Electrode Mapping

% Load Electrode Mapping Data 
load('ElectrodeMap.mat');

% Extract channel numbers from the selected channel names 
channel_numbers = str2double(extractAfter(selected_channels, 'Chan'));

% Define a grid layout for electrode positions 
[X_grid, Y_grid] = meshgrid(1:10, 1:10);

%% 1(a): Plot Mapping for Cerebus Channels
figure;
hold on;
title('Electrode Mapping on Monkey Cortex (Cerebus Channels)');
xlabel('Electrode X Position');
ylabel('Electrode Y Position');
grid on;
axis equal;
for i = 1:length(channel_numbers)
    idx = find(MapStruct.CerebusChannel == channel_numbers(i));
    if ~isempty(idx)
        electrode_index = MapStruct.CerebusElectrode(idx);
        elec_x = X_grid(electrode_index);
        elec_y = Y_grid(electrode_index);
        scatter(elec_x, elec_y, 150, 'r', 'filled');
        text(elec_x + 0.2, elec_y, ['Chan' num2str(channel_numbers(i))], 'FontSize', 12);
    end
end
hold off;

%% 1(b): Plot Mapping for TDT Channels
figure;
hold on;
title('Electrode Mapping on Monkey Cortex (TDT Channels)');
xlabel('Electrode X Position');
ylabel('Electrode Y Position');
grid on;
axis equal;
for i = 1:length(channel_numbers)
    idx = find(MapStruct.TDTChannel == channel_numbers(i));
    if ~isempty(idx)
        % Here we use the CerebusElectrode index for layout (adjust if necessary)
        electrode_index = MapStruct.CerebusElectrode(idx);
        elec_x = X_grid(electrode_index);
        elec_y = Y_grid(electrode_index);
        scatter(elec_x, elec_y, 150, 'b', 'filled');
        text(elec_x + 0.2, elec_y, ['Chan' num2str(channel_numbers(i))], 'FontSize', 12);
    end
end
hold off;



%% Question 2: Using BME_526_HW2_All.mat

% Converting fields from AllData structure to numeric arrays
finger_data = [AllData.finger];
time_data = [AllData.TimeStampSec];

% Defining movement types and corresponding finger values
movement_types2 = {'Thumb', 'Index', 'Middle', 'Thumb & Index', 'Thumb, Index & Middle', 'Thumb & Middle'};
finger_values = [1, 2, 3, 4, 5, 6];

% Extract event times based on finger values
event_times2 = cell(1, numel(movement_types2));
for i = 1:numel(movement_types2)
    finger_rows = (finger_data == finger_values(i));
    event_times2{i} = time_data(finger_rows);
end

% Loop over each selected neuron
for neuron_idx = 1:length(selected_channels)
    figure;
    set(gcf, 'Position', [100, 100, 1400, 700]);
    sgtitle(['Neuron Response - ', selected_channels{neuron_idx}, ' (Unit ', num2str(selected_units(neuron_idx)), ') [BME-526-HW2-All.mat]']);
    
    % Extract spike times for the current neuron
    spike_times = NeuralData.Channels.(selected_channels{neuron_idx});
    neuron_spike_times = spike_times(spike_times(:,2) == selected_units(neuron_idx), 3);
    
    % Loop through each movement type
    for mov_idx = 1:length(movement_types2)
        current_event_times = event_times2{mov_idx};
        numTrials = numel(current_event_times);
        
        %% Raster Plot
        subplot(2, length(movement_types2), mov_idx);
        hold on;
        AP_perBin = zeros(numTrials, nbins);
        all_spikes = [];
        for trial = 1:numTrials
            cueTime = current_event_times(trial);
            % Select spikes in the window [-1,1] and align them
            trial_spikes = neuron_spike_times(neuron_spike_times >= (cueTime - 1) & neuron_spike_times <= (cueTime + 1)) - cueTime;
            AP_perBin(trial, :) = histcounts(trial_spikes, bins);
            all_spikes = [all_spikes; trial_spikes(:)];
            plot(trial_spikes, trial * ones(size(trial_spikes)), 'k.', 'MarkerSize', 5);
        end
        title(movement_types2{mov_idx});
        xlabel('Time (s)');
        ylabel('Trial');
        xlim([-1,1]);
        ylim([0, numTrials+1]);
        
        %% PSTH Plot
        thresholds = [35, 30, 20];
        subplot(2, length(movement_types2), mov_idx + length(movement_types2));
        hold on;
        std_perBin = std(AP_perBin);
        counts = histcounts(all_spikes, bins) / (numTrials * spacing); 
        smoothed_counts = smoothdata(counts, 'gaussian', 5);
        
        x_vals = bins(1:end-1);
        fill_x = [x_vals, fliplr(x_vals)];
        fill_y = [smoothed_counts - 2*std_perBin, fliplr(smoothed_counts + 2*std_perBin)];
        fill(fill_x, fill_y, [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(x_vals, smoothed_counts, 'b', 'LineWidth', 2);
        yline(thresholds(neuron_idx), '--r', 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Firing Rate (Hz)');
        ylim([0, max(smoothed_counts + 2*std_perBin) + 5]);
    end
end


%% Individual Neuron Responses:
%The selection of channels 40, 30, and 54 proves to be strategically advantageous for decoding finger movements due to their distinct response fields:
% Chan54: This neuron seems to respond most strongly when moving the thumb or middle finger individually, suggesting its role is in controlling these specific single finger movements.
% Chan30: This neuron shows the strongest response when moving the thumb, index, and middle fingers together, indicating it’s involved in coordinating multiple fingers at once.
% Chan40: This neuron responds most when the thumb and index fingers are moving together, pointing to its role in controlling this specific combination of fingers.

% Decoding Finger Movements: 
% Chan54: This neuron is particularly good at detecting individual thumb or middle finger movements, with strong firing rates of 20–30 Hz, making it useful for decoding these specific movements. 
% Chan30: This neuron shows its highest activity, around 30–35 Hz, when the thumb, index, and middle fingers are moved together, so it’s well-suited for decoding this combined movement. 
% Chan40: This neuron also peaks at 35–40 Hz for that same combined movement, and it responds strongly to individual thumb movement (30–35 Hz), which could help decode both.


%  Action Potential Firing Rate Thresholds:
%  Chan54: A threshold of 20 Hz seems good for detecting individual thumb or middle finger movements, as its firing rate for these is above this level, while combined movements are usually lower.
%  Chan30: Setting the threshold at 30 Hz works well for detecting the combined thumb, index, and middle finger movement, as this is where it shows the strongest response.
%  Chan40: A 35 Hz threshold is suitable for detecting the same combined movement, though it might also catch some individual thumb movements, which can reach up to 35 Hz.
%%


