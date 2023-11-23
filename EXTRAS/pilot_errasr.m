% Load the data - normal spikes and then add alil noisy pulse 
% Also pass it through high pass filter and also make it 2d
% 
nwb = nwbRead('/Users/rishabhsinghal/lab_data/NWB_files/Timon_1711101_stim.nwb');

% get the artifacts from 16 channels in the nwb file

% step1 - get the indices first - two columns corrosponding to start and
% stop indices of artifacts

ts = nwb.acquisition.get('Behavior_Timestamp');
ts = ts.toTable();

stim = ts.('stimTrial');
stim_ts = ts.('stimOnTime');

stim_ts = stim_ts(stim==1);
stim_ts = stim_ts.*40000;

pre_window = 10000;
post_window = 29999;

% intialise a matrics trial x time x channels
%%

spike_raw = nwb.acquisition.get('FEF_recording').data(:,:);
% plot(spike_raw)

%% grab few trials
stim_ts = stim_ts(1:2);
data_n = NaN(length(stim_ts),45000,16);

for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    data_n(i,:,:) = spike_raw(:,takeout:takeafter)';
end
%% Reshape the matirx

% Permute the dimensions
data_n = permute(data_n, [3, 1, 2]);

% Reshape the matrix
data_n = reshape(data_n, 16, []);

yOffset = 6;

% Create the plot iteratively
figure;
hold on;
for channel = 1:10
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ data_n(channel,:));
end

xlabel('time');
ylabel('channels');

title('Iterative Plot with Separation');
hold off;

%% Highpass filter


% Design a high-pass filter
Fs = 40000; % Sample frequency (adjust based on your data)
Fcutoff = 200; % Cutoff frequency (adjust as needed)

d = designfilt('highpassiir', 'FilterOrder', 4, ...
               'PassbandFrequency', Fcutoff, ...
               'PassbandRipple', 0.2, ...
               'SampleRate', Fs);

% Apply the high-pass filter to each row
filteredMatrix = zeros(size(data_n));
for i = 1:size(data_n, 1)
    filteredMatrix(i, :) = filtfilt(d, data_n(i, :));
end

Fs = 40000; % Sample frequency (adjust based on your data)
Fcutoff = 200; % Cutoff frequency (adjust as needed)

d = designfilt('highpassiir', 'FilterOrder', 4, ...
               'PassbandFrequency', Fcutoff, ...
               'PassbandRipple', 0.2, ...
               'SampleRate', Fs);

% Apply the high-pass filter to each row
dataCleaned0 = zeros(size(dataCleaned));
for i = 1:size(dataCleaned, 1)
    dataCleaned0(i, :) = filtfilt(d, dataCleaned(i, :));
end
%%

% Create a figure with the desired size
figure('Units', 'inches', 'Position', [1, 1, 8 ,6]);
x = [1:length(dataCleaned0(1,10700:11000,14))]/40000;
plot(x,filteredMatrix(1,10700:11000,14),'color','black','LineWidth',2);
hold on; 
% plot(x,dataCleaned0(1,10700:11000,14),'color','r','LineWidth',0.8);
% Create a line plot

yticks([])
xticks([])
% Set the font size for the axes
set(gca, 'FontSize', 18);

% Remove the top and right axes lines
box off; % Turns off the box surrounding the plot
ax = gca; % Get current axes
ax.YAxisLocation = 'left';
ax.XAxisLocation = 'bottom';
ax.Box = 'off'; % Turn off the box surrounding the plot

% Set the remaining axes to have a larger font size
ax.FontSize = 18;
% lgd = legend('raw', 'clean');
% set(lgd, 'Box', 'off', 'FontSize', 18);
xlabel('Time');
ylabel('Voltage')
% Optionally, if you want to explicitly remove the ticks from the top and right axes
% (not usually necessary with 'box off' command):
% ax.XAxis.TickDir = 'out'; % Ensures ticks are only outside the bottom axis
% ax.YAxis.TickDir = 'out'; % Ensures ticks are only outside the left axis

% Adjust the figure's PaperSize property to make sure the saved figure has the correct dimensions
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6, 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0, 0, 6, 5]);

% Now, if you want to save the figure, you can use the print function
print('pulse','-dsvg'); % Saves the figure at 300 dpi as a PNG file
%%
% Example data


% Create a figure with the desired size
figure('Units', 'inches', 'Position', [1, 1, 8 ,6]);
x = [1:length(dataCleaned0(1,:,14))]/40000;
plot(x,filteredMatrix(1,:,14),'color','black','LineWidth',0.8);
hold on; 
plot(x,dataCleaned0(1,:,14),'color','r','LineWidth',0.8);
% Create a line plot


% Set the font size for the axes
set(gca, 'FontSize', 18);

% Remove the top and right axes lines
box off; % Turns off the box surrounding the plot
ax = gca; % Get current axes
ax.YAxisLocation = 'left';
ax.XAxisLocation = 'bottom';
ax.Box = 'off'; % Turn off the box surrounding the plot

% Set the remaining axes to have a larger font size
ax.FontSize = 18;
lgd = legend('raw', 'clean');
set(lgd, 'Box', 'off', 'FontSize', 18);
xlabel('Time');
ylabel('Voltage')
% Optionally, if you want to explicitly remove the ticks from the top and right axes
% (not usually necessary with 'box off' command):
% ax.XAxis.TickDir = 'out'; % Ensures ticks are only outside the bottom axis
% ax.YAxis.TickDir = 'out'; % Ensures ticks are only outside the left axis

% Adjust the figure's PaperSize property to make sure the saved figure has the correct dimensions
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6, 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0, 0, 6, 5]);

% Now, if you want to save the figure, you can use the print function
print('raw_clean','-dsvg'); % Saves the figure at 300 dpi as a PNG file


%%
% Display the filtered matrix (optional)


% ADD artificial stimulation

% Assuming flattenedMatrix is of size 16x400

% Define timepoints for the pulse pattern
t1 = 501;  % Example start timepoint
t2 = 10500; % Example end timepoint

% Duration of a single pulse
pulseDuration = (t2 - t1 + 1) /20;

% Generate a modulated sine wave for a single pulse
t = linspace(0, 2*pi, pulseDuration);
singlePulse = sin(t)*0.15;

% Repeat the pulse for the desired number of pulses (3 in this example)
pulsePattern = repmat(singlePulse, 1, 20);
X = filteredMatrix;
X_stim = filteredMatrix;
% Add the pulse pattern to each row
X_stim(:, t1:t2) = filteredMatrix(:, t1:t2)+repmat(pulsePattern, size(filteredMatrix, 1), 1);

% second pulse 
t1 = 40501;  % Example start timepoint
t2 = 50500; % Example end timepoint

% Duration of a single pulse
pulseDuration = (t2 - t1 + 1) /20;

% Generate a modulated sine wave for a single pulse
t = linspace(0, 5*pi, pulseDuration);
singlePulse = cos(t)*0.15;

% Repeat the pulse for the desired number of pulses (3 in this example)
pulsePattern = repmat(singlePulse, 1, 20);

% Add the pulse pattern to each row
X_stim(:, t1:t2) = filteredMatrix(:, t1:t2)+repmat(pulsePattern, size(filteredMatrix, 1), 1);


yOffset = 6;

% Create the plot iteratively
figure;
hold on;
for channel = 1:10
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ X_stim(channel,:));
end

xlabel('time');
ylabel('channels');

title('Iterative Plot with Separation');
hold off;
figure;
hold on;
for channel = 1:10
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ X(channel,:));
end

xlabel('time');
ylabel('channels');

title('Iterative Plot with Separation');
hold off;

%% ERRASR PART - Clean it and look at the PCA

[V,explained] = pca_errasr(X',X_stim',16,0.1);
%%
% Assuming you have X_centered and V from the PCA

% Compute scores
scores = X_stim' * V;

% Compute the variance explained by each PC
variance_explained = var(scores);

% Sort PCs by explained variance in descending order
[~, order] = sort(variance_explained, 'descend');
V_sorted = V(:, order);

% Plot sorted explained variance
figure;
bar(variance_explained(order) / sum(variance_explained) * 100);
title('Explained Variance by Each PC');
xlabel('Principal Component');
ylabel('Percentage of Total Variance Explained');


%%
% X_stim = X_stim - mean(X_stim, 1);
X_stim = X_stim';
X = X';
%% PLot again
% coeff = pca_manifold(X_stim);
%coeff = PCA_stochastic(X_stim,10);
[coeff,~,~,~,explained1] = pca(X_stim);
%%
Fs = 40000; % Sample frequency (adjust based on your data)
Fcutoff = 200; % Cutoff frequency (adjust as needed)

d = designfilt('highpassiir', 'FilterOrder', 4, ...
               'PassbandFrequency', Fcutoff, ...
               'PassbandRipple', 0.2, ...
               'SampleRate', Fs);

% Apply the high-pass filter to each row
X = load('X.mat');
X = X.data;
X = squeeze(X(1,:,:));
data = squeeze(X);

for i = 1:size(data, 1)
    X(i, :) = filtfilt(d, data(i, :));
end

%%
Xstim = load('Xstim.mat');
Xstim = Xstim.pcaMatOverTrials;
X_stim = Xstim(:,:,1);

%%

[coeff1,explained] = Becket_pca(X_stim,X,0.05);
explained = explained*100;
%%
figure(1);
[coeff,~,~,~,explained1] = pca(X_stim);
projected_data = X_stim * coeff;

plot(projected_data(:,8));
title('Normal PCA')
%%
figure(2);
projected_data = X_stim * coeff1;
plot(projected_data(:,8));
title('ERRASR PCA')

%%
dataCleaned = pcaCleaned;
for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    spike_raw(:,takeout:takeafter) = squeeze(dataCleaned(:,:,i));
end