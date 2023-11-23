% the script to use the nwb file to clean the artifacts

% read the file


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
stim_ts = stim_ts(1:4);
data_n = NaN(length(stim_ts),40000,16);
% data_n = NaN(10,40000,16);

for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    data_n(i,:,:) = spike_raw(:,takeout:takeafter)';
end
%% Grab spontaneous data as well
data = NaN(40000,16,length(stim_ts));
for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)- 50000);
    takeafter = floor(stim_ts(i) - 10000-1);
    data(:,:,i) = spike_raw(:,takeout:takeafter)';
end

% % Plotting parameters`
yOffset = 6;
% data = data_n(1:20,:,:);
% Create the plot iteratively
figure;
hold on;
for channel = 1:16
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ data(:,channel,1));
end
% 
% for row = 1:10
%     yValues = row*0.1 + yOffset;  % Add offset to y-values
%     plot(yValues+ data_n(1,:,row));
% end
% 
% 
% % Adjust plot appearance
% xlabel('time');
% ylabel('channels');
% 
% title('Iterative Plot with Separation');
% hold off;

%%

opts = ERAASR.Parameters();
opts.Fs = 40000; % samples per second
Fms = opts.Fs / 1000; % multiply to convert ms to samples

opts.thresholdHPCornerHz = 250;
opts.thresholdChannel = 8;
opts.thresholdValue = 3;

opts.alignChannel = 1;
opts.alignUpsampleBy = 12;
opts.alignWindowPre = Fms * 0.5;
opts.alignWindowDuration = Fms * 15;

% 60 ms stim, align using 20 ms pre start to 110 post
opts.extractWindowPre = Fms * 30;
opts.extractWindowDuration = Fms * 560;
opts.cleanStartSamplesPreThreshold = Fms * 1;
        
opts.cleanHPCornerHz = 10; % light high pass filtering at the start of cleaning
opts.cleanHPOrder = 4; % high pass filter order 
opts.cleanUpsampleBy = 1; % upsample by this ratio during cleaning
opts.samplesPerPulse = Fms * 5; % 3 ms pulses - changes to 5
opts.nPulses = 100;

opts.nPC_channels = 12;
opts.nPC_trials = 3;
opts.nPC_pulses = 6;

opts.omit_bandwidth_channels = 3;
opts.omit_bandwidth_trials = 1;
opts.omit_bandwidth_pulses = 1;

opts.alignPulsesOverTrain = true; % do a secondary alignment within each train, in case you think there is pulse to pulse jitter. Works best with upsampling
opts.pcaOnlyOmitted = true; % if true, build PCs only from non-omitted channels/trials/pulses. if false, build PCs from all but set coefficients to zero post hoc
opts.lamda = 0.1;
opts.cleanOverChannelsIndividualTrials = true;
opts.cleanOverPulsesIndividualChannels = false;
opts.cleanOverTrialsIndividualChannels = false;

opts.cleanPostStim = false; % clean the post stim window using a single PCR over channels

opts.showFigures = false; % useful for debugging and seeing well how the cleaning works
opts.plotTrials = 1; % which trials to plot in figures, can be vector
opts.plotPulses = 1; % which pulses to plot in figures, can be vector
opts.figurePath = [pwd '/timon_171101_postbug'];% folder to save the figures
% if ~exist(opts.figurePath, 'dir')
%     mkdir(opts.figurePath)
% else
%     disp('The folder does exist.');
% end

opts.saveFigures = false; % whether to save the figures
opts.saveFigureCommand = @(filepath) print('-dpng', '-r300', [filepath '.png']); % specify a custom command to save the figure

%% Do alignment and cleaning procedure

[dataCleaned, extract] = ERAASR.cleanTrials(data_n,data, opts);
%
% Plot the cleaned traces on a single trial

% data = data_n(1:20,:,:);
%% Create the plot iteratively
figure(1);
hold on;
yOffset = 5;
for channel = 1:16
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ dataCleaned(1,:,channel));
end
%%
figure(1);
hold on;
yOffset = 5;
for channel = 1:16
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ data_n(1,:,channel));
end

%%

figure();
plot(squeeze(dataCleaned(1, :, 1)),'-');
box off;

%% PUt the data back into the field
for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    spike_raw(:,takeout:takeafter) = squeeze(dataCleaned(i,:,:))' ;
end

%% write the same spike raw to nwb as a processed module but same manner
electrical_series = types.core.ElectricalSeries(...
    'description', 'Post_ERRASR_data_no_HP',...
    'starting_time', 0, ... % seconds
    'starting_time_rate', 40000, ... % Hz
    'data', spike_raw, ...
    'electrodes', nwb.acquisition.get('FEF_recording').electrodes.loadAll, ...
    'data_unit', 'volts'); 

ecephys_module = types.core.ProcessingModule(...
    'description', 'extracellular electrophysiology post errasr');
% 
ecephys_module.nwbdatainterface.set('ERRASR', electrical_series);
% 
% nwb.processing.set('ecephys', ecephys_module);
nwb.processing.set('ERRASR', ecephys_module);

%% save the file


nwbExport(nwb, '/Users/rishabhsinghal/lab_data/NWB_files/Timon_1711101_stim_new.nwb')