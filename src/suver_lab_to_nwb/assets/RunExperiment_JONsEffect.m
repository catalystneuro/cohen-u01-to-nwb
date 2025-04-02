%%
% run experiment with 3 cameras, and 1 AM systems 2400 headstage
%   - dorsal camera captures both antenna, while lateral cameras capture
%   dorsal/ventral movements of singular antenna
%   - headstage simultaneously records flucations in membrane potential of
%   patched cell
% 
%% last edited by Olivia on 5/27/2024
% ****** 'filteredI' was changed to 'filteredVm' on 8/30/24 (only on
%           successfully patched APN2 cells)
%       (labeled wrong in the begining because I didnt kow how to read the
%       amplifier correctly, lol #growth)
%    double cheked by Marie
%
%% Typical use:
% RunExperiment_JONsEffect('glued', '2-4dpe', '24C06', 1)
%%          

function [] = RunExperiment_JONsEffect(condition, age, genotype, expNumber)

AM_SYSTEMS = 0;
FPS = 100; % frame rate
numTrials = 25; 
Fs = 10000; % sample rate
trialLength = 10; %in seconds
duration = trialLength; 

if nargin < 3
   display('Please enter the following inputs:')
   display('1: condition (string) - is antennae glued or free?')
   display('2: age (string) - input age range of this fly')  
   display('3: genotype (string) - input genotype for this experiment')   
   display('4: exptNumber (int) - number associated with this set of trials (will append to previous if file exists with this exptNumber already)')
   return
end

formatOut = 'yyyy_mm_dd';
dateStr = datestr(date, formatOut);
dirStr = ['E:\Data\' dateStr];
if ~isfolder(dirStr)
    mkdir(dirStr)
end

%% Open data structure and count trials

% check whether a saved data file exists with today's date
saveStr = ['C:\Users\suver lab\Vanderbilt\SuverLab - E.phys Rig\Olivia\Data\' dateStr, '\' dateStr,'_E',num2str(expNumber)]
directory = dir([saveStr, '.mat']);
if isempty(directory)
    % if no saved data exists, then this is the first trial
    nn = 0;
else
    load(saveStr, 'data'); %load current data file
    nn = length(data);     %most recent piece of data (will append to this)
    display(['Appending to existing data structure for this expt (num trials saved= ' num2str(nn) ')']);
end

%% Reset aquisition engines
daqreset;

aio =  daq('ni');
aio.Rate = Fs;

Vm = addinput(aio, "Dev1", "ai0","Voltage"); % Signal from 10xVm port on AM Systems amp
Vm.TerminalConfig = 'SingleEnded';
I = addinput(aio, "Dev1", "ai1","Voltage"); % Signal from I port on AM Systems amp
I.TerminalConfig = 'SingleEnded';
filteredVm = addinput(aio, "Dev1", "ai2","Voltage"); % Signal from filtered Im port on AM Systems amp
filteredVm.TerminalConfig = 'SingleEnded';
puffer = addinput(aio, "Dev1", "ai3","Voltage"); % Signal from Puffer
puffer.TerminalConfig = 'SingleEnded';
tachometer = addinput(aio,"Dev1","ai4","Voltage"); % Tachometer Signal
tachometer.TerminalConfig = 'SingleEnded';

%record_cams = addinput(aio, "Dev1","ai4","Voltage"); % initiating camera output
trig_cams = addoutput(aio, "Dev1","ao1","Voltage"); % initiating camera trigger
cameraTrigger = 4; % volts (determined by product documentation on Basler's website)

%% Trigger pulse creation
n_samp = round(duration * Fs); % Total number of samples to be written to DAQ
trigger = zeros(1, n_samp); % Initialize an array of all zeros with length n_samp
interval_size = round(Fs/FPS); % Calculate the length of one square wave interval in samples (samples per frame)
samples_high = round(interval_size / 5); % Duration of the high phase of the square wave

start_r = samples_high;

% Generate the square wave signal
while start_r <= n_samp
    trigger(start_r:start_r + samples_high) = cameraTrigger; % Assign 4 to the interval corresponding to high phase
    start_r = start_r + interval_size; % Move to the next interval
end

%% Iterate through entire stimulus set
for trialNum = 1:numTrials
    nn = nn + 1
    %% Save information about this piece of data
    data(nn).date = dateStr;                                % date of experiment (embedded in filename and directory)
    data(nn).expNumber = expNumber;                         % experiment number (can be multiple for one cell/fly)
    data(nn).trial = nn;                                    % trial number within this experiment (can run multiple sets of trials within one experiment)
    data(nn).condition = condition;                         % glued antenna or not glued antenna
    data(nn).age = age;                                     % age range of experimental flies
    data(nn).genotype = genotype;                           % genotype of experimental flies
    data(nn).samplerate = Fs;                               % samplerate of stimulus and data acquired
    data(nn).fps = FPS;
    data(nn).scaleCurrent = 100;%200;                       % axon = 100 % scaling factor for picoamps (with MultiClamp 200B!)
    data(nn).scaleVoltage = 100;%10;                        % scaling factor for mV
   
    display(['Trial num= ' num2str(trialNum)])

    %% Camera Initialization in NiDAQ
    imaqreset;
    imaqmex('feature','-limitPhysicalMemoryUsage',false);

    data(nn).fps = FPS;
    data(nn).nframes = data(nn).fps*trialLength;
    
    %% set up the first camera (fly's left antenna - lateral cam)
    vid1 = videoinput('gentl', 1, 'Mono8'); % fly's right lateral camera
    vid1.ROIPosition = [0 0 640 480]; % offsetX, offsetY, width, height
    src1 = getselectedsource(vid1);
    src1.ExposureAuto = 'Off';
    src1.ExposureTime = 2000; %ring light setting: 9930; % lower values = higher frame rate
    src1.GainAuto = 'Once';
    %src1.Gamma = 0.8;
    %src1.BalanceRatio = 1.5; % couldnt figure this variable out at the time, and didnt seem to effet my video quality
    
    %configure line 4 (green) from Basler - output voltage pulse to nidaq
    %(reads trigger pulse)
    src1.LineSelector = "Line4";
    src1.LineMode = "Output";
    src1.LineInverter = "False";
    src1.LineSource = "FrameTriggerWait";

    %configure line 3 (brown) from Basler - voltage pulse from nidaq to camera
    src1.LineSelector = "Line3";
    src1.LineMode = "Input";
    src1.TriggerSelector = 'FrameStart';
    src1.TriggerSource = "Line3";
    src1.TriggerMode = "on";
    src1.TriggerActivation = "RisingEdge";
    src1.TriggerDelay = 0;
    src1.AcquisitionFrameRateEnable = "False"; % when set to false, frame rate is determined by roi, exposure, gain, and gamma values
    src1.AcquisitionFrameRate = FPS;
    src1.DeviceLinkThroughputLimitMode = 'Off';

    % save video
    vidSaveStr1 = [saveStr, '_Video_lateral_flyRight_' num2str(nn)]
    vidfile1 = vidSaveStr1;
    logfile1 = VideoWriter(vidfile1, 'Motion JPEG AVI');
    vid1.LoggingMode = 'Disk';

    logfile1.FrameRate = FPS; %
    numFrames = logfile1.FrameRate*trialLength; %
    vid1.FramesPerTrigger = numFrames; %

    vid1.DiskLogger = logfile1;
    triggerconfig(vid1,'hardware', 'DeviceSpecific', 'DeviceSpecific')
    
    %% set up the second camera (both antennae & wing envolope - ventral cam)
    vid2 = videoinput('gentl', 2, 'Mono8'); % ventral camera
    vid2.ROIPosition = [0 0 640 480]; % offsetX, offsetY, width, height
    src2 = getselectedsource(vid2);
    src2.ExposureAuto = 'Off';
    src2.ExposureTime = 300; % ring light setting: 3000; % lower values = higher frame rate
    src2.Gain = 0;
    %src2.BalanceRatio = 1.5; % couldnt figure this variable out at the time, and didnt seem to effet my video quality
    
    %configure line 4 (green) from Basler - output voltage pulse to nidaq
    %(reads trigger pulse)
    src2.LineSelector = "Line4";
    src2.LineMode = "Output";
    src2.LineInverter = "False";
    src2.LineSource = "FrameTriggerWait";

    %configure line 3 (brown) from Basler - voltage pulse from nidaq to camera
    src2.LineSelector = "Line3";
    src2.LineMode = "Input";
    src2.TriggerSelector = 'FrameStart';
    src2.TriggerSource = "Line3";
    src2.TriggerMode = "on";
    src2.TriggerActivation = "RisingEdge";
    src2.TriggerDelay = 0;
    src2.AcquisitionFrameRateEnable = "False"; % when set to false, frame rate is determined by roi, exposure, gain, and gamma values
    src2.AcquisitionFrameRate = FPS;
    src2.DeviceLinkThroughputLimitMode = 'Off';

    % save video
    vidSaveStr2 = [saveStr, '_Video_lateral_ventral_' num2str(nn)]
    vidfile2 = vidSaveStr2;
    logfile2 = VideoWriter(vidfile2, 'Motion JPEG AVI');
    vid2.LoggingMode = 'Disk';

    logfile2.FrameRate = FPS; %
    numFrames = logfile2.FrameRate*trialLength; %
    vid2.FramesPerTrigger = numFrames; %

    vid2.DiskLogger = logfile2;
    triggerconfig(vid2,'hardware', 'DeviceSpecific', 'DeviceSpecific')

    %% set up the third camera (fly's right antenna - lateral cam)
    vid3 = videoinput('gentl', 3, 'Mono8'); % fly's right lateral camera
    vid3.ROIPosition = [0 0 640 480]; % offsetX, offsetY, width, height
    src3 = getselectedsource(vid3);
    src3.ExposureAuto = 'Off';
    src3.ExposureTime =  2000; %ring light setting:9930; % lower values = higher frame rate
    src3.GainAuto = 'Once';
    %src3.Gamma = 0.8;
    %src1.BalanceRatio = 1.5; % couldnt figure this variable out at the time, and didnt seem to effet my video quality
    
    %configure line 4 (green) from Basler - output voltage pulse to nidaq
    %(reads trigger pulse)
    src3.LineSelector = "Line4";
    src3.LineMode = "Output";
    src3.LineInverter = "False";
    src3.LineSource = "FrameTriggerWait";

    %configure line 3 (brown) from Basler - voltage pulse from nidaq to camera
    src3.LineSelector = "Line3";
    src3.LineMode = "Input";
    src3.TriggerSelector = 'FrameStart';
    src3.TriggerSource = "Line3";
    src3.TriggerMode = "on";
    src3.TriggerActivation = "RisingEdge";
    src3.TriggerDelay = 0;
    src3.AcquisitionFrameRateEnable = "False"; % when set to false, frame rate is determined by roi, exposure, gain, and gamma values
    src3.AcquisitionFrameRate = FPS;
    src3.DeviceLinkThroughputLimitMode = 'Off';

    % save video
    vidSaveStr3 = [saveStr, '_Video_lateral_flyLeft' num2str(nn)]
    vidfile3 = vidSaveStr3;
    logfile3 = VideoWriter(vidfile3, 'Motion JPEG AVI');
    vid3.LoggingMode = 'Disk';

    logfile3.FrameRate = FPS; %
    numFrames = logfile3.FrameRate*trialLength; %
    vid3.FramesPerTrigger = numFrames; %

    vid3.DiskLogger = logfile3;
    triggerconfig(vid3,'hardware', 'DeviceSpecific', 'DeviceSpecific')
    
    % start video streaming & data aquisition
    preview(vid1)
    preview(vid2)
    preview(vid3)
    nn;
    start(vid1)
    start(vid2)
    start(vid3)
 
%% record data in matrix
    trialdata = readwrite(aio, trigger',"OutputFormat","Matrix"); % 100000*4 matrix of data points
    % each column of the matrix is representative of the inputs: 
    data(nn).Vm = trialdata(:,1); 
    data(nn).I = trialdata(:,2);
    data(nn).filteredVm = trialdata(:,3);
    data(nn).puffer = trialdata(:,4);
    data(nn).tachometer = trialdata(:,5);

    write(aio, 0) %one zero for one output signal (cameras) - resetting camera trigger to 0
    save(saveStr, 'data')

end

disp('taaaadaaaa! data collection complete :)')

end
