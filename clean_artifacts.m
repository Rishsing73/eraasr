% the script to use the nwb file to clean the artifacts

% read the file
filename = '/Users/rishabhsinghal/lab_data/kirk171026.nwb';
nwb = nwbRead(filename); %CHNAGE

ts = nwb.acquisition.get('Behavior_Timestamp');
ts = ts.toTable();

stim = ts.('stimTrial');
stim_ts = ts.('stimOnTime');

stim_ts = stim_ts(stim==1); %all trials are stimulation trials
stim_ts = stim_ts.*40000;

pre_window = 10000;
post_window = 29999;

% intialise a matrics trial x time x channels
spike_raw = nwb.acquisition.get('FEF_recording').data(:,:);


%% grab few trials - select all trials - CHNAGE
% n_trials = 4; %uncommment if you want to run short data only
% stim_ts = stim_ts(1:n_trials);
data_n = NaN(length(stim_ts),40000,16);
for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    data_n(i,:,:) = spike_raw(:,takeout:takeafter)';
end


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

opts.omit_bandwidth_channels = 2; % omit these channels PCs
opts.omit_bandwidth_trials = 1;
opts.omit_bandwidth_pulses = 1;

opts.alignPulsesOverTrain = true; % do a secondary alignment within each train, in case you think there is pulse to pulse jitter. Works best with upsampling
opts.pcaOnlyOmitted = false; % if true, build PCs only from non-omitted channels/trials/pulses. if false, build PCs from all but set coefficients to zero post hoc
opts.lamda = 0.1;
opts.cleanOverChannelsIndividualTrials = false;
opts.cleanOverPulsesIndividualChannels = false;
opts.cleanOverTrialsIndividualChannels = false;

opts.cleanPostStim = false; % clean the post stim window using a single PCR over channels

opts.showFigures = false; % useful for debugging and seeing well how the cleaning works - CHANGE
opts.plotTrials = 1; % which trials to plot in figures, can be vector
opts.plotPulses = 1; % which pulses to plot in figures, can be vector
opts.figurePath = [pwd '/timon_171101_postbug'];% folder to save the figures
% if ~exist(opts.figurePath, 'dir')
%     mkdir(opts.figurePath)
% else
%     disp('The folder does exist.');
% end

opts.saveFigures = false; % whether to save the figures - CHANGE
opts.saveFigureCommand = @(filepath) print('-dpng', '-r300', [filepath '.png']); % specify a custom command to save the figure

%% Do alignment and cleaning procedure

[dataCleaned, extract] = ERAASR.cleanTrials(data_n, opts); % we have removed 0.5 ms of data from each pulse - 50 ms of data from 500ms

%% Create the plot iteratively - plot for cleaned data vs uncleaned
figure(2);  
hold on;
trial = 2;
yOffset = 10;
for channel = 1:16
    yValues = channel* yOffset;  % Add offset to y-values
    plot(yValues+ dataCleaned(trial,:,channel),'r');
    plot(yValues+ data_n(trial,:,channel),'b');
end
ylabel('voltage')
xlabel('time')



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
nwbExport(nwb, filename) %CHANGE