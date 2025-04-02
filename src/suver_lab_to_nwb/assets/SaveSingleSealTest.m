%%
% Run and save single seal test trial! Will save 10 s with a 10mV peak-peak
% pulse at 100 Hz :)
%
% Last edited by Marie on 05/11/2023 - signals look correct w/ AM today!
% Configured for use with AM Systems 2400, gain =10, EXT /10, seal test = 10 mV
% 10GOhm/*100MOhm* (probe gain LOW) (100mV/nA)
%%

function [] = SaveSingleSealTest(expNumber)

if nargin < 1
    disp('Please enter the following inputs:')
    disp('1: exptNumber (int) - number associated with this set of trials (will append to previous if file exists with this exptNumber already)')
    return
end
%%
formatOut = 'yyyy_mm_dd';
dateStr = datestr(date, formatOut);

parentDir = 'C:\Users\suver lab\Documents\Olivia\Data\';
dirName = ['C:\Users\suver lab\Documents\Olivia\Data\' dateStr '\'];
fileStr = [dateStr '_SealTest_' num2str(expNumber)];
saveStr = [dirName fileStr];
if exist(dirName) ~= 7 %make the folder if it's not there already for this data set
   mkdir(parentDir, dateStr);
end

%% Make a square wave for seal test
EXT_DIV = 10; % /10 for EXT input on amplifier
SAMPLERATE = 10000;
NUM_SEC_PULSE = 1;
pulseAmp = .01 % in V (10mV)
pulseFreq = 100; % in Hz
dutyCycle = 0.5; % equal high/low square wave

on = 1:SAMPLERATE/pulseFreq;
oneCycle = [zeros(1,(SAMPLERATE/pulseFreq)/2) ones(1,length(on)).*pulseAmp zeros(1,(SAMPLERATE/pulseFreq)/2)];
onePulse = EXT_DIV.*repmat(oneCycle,1,NUM_SEC_PULSE*SAMPLERATE/pulseFreq);


%% Run seal test for one "trial"
nn = 1;
%% Save information about this piece of data
data(nn).date = dateStr;          % date of experiment (embedded in filename and directory_
data(nn).expNumber = expNumber;   % experiment number (can be multiple for one cell/fly)
data(nn).trial = nn;              % trial number within this experiment (can run multiple sets of trials within one experiment)
data(nn).numSecOut = 2;           % number of seconds to save and analyze
data(nn).scaleCurrent = 100;      % scaling factor for nA (prev. picoamps?) w/ gain = 10 (100mV/nA) in voltage clamp w/ AM 2400
data(nn).scaleVoltage = 100;      % scaling factor for mV  (meas. in voltage clamp, gain = 10 & x10 output) w/ AM 2400 
data(nn).SAMPLERATE = SAMPLERATE;

%% Reset aquisition engines
daqreset;
% %% Set up and send digital output via NIDAQ board
% sa = daq.createSession('ni'); %establish a connection with the NIDAQ
% sa.Rate = SAMPLERATE; %samplerate (Hz)
% AOch = sa.addAnalogOutputChannel('Dev1',0,'Voltage'); %set up one channels for output (1st for valve control, 2nd for MFC)
% sa.queueOutputData([onePulse']); %load this trial's seal test voltage square wave signal
% 
% %% Configure input channels
% AIch = sa.addAnalogInputChannel('Dev1',0:1,'Voltage'); %set up three analog inputs (voltage, current, analogOutput command)
% AIch(1).TerminalConfig = 'SingleEnded'; %this is necessary to get absolute (rather than relative) value at this input!!
% AIch(2).TerminalConfig = 'SingleEnded';

%% Set up analog input and output via NIDAQ board
sa = daq('ni'); %establish a connection with the NIDAQ
sa.Rate = SAMPLERATE; %samplerate (Hz)
 % 
AOch = addoutput(sa,'Dev1',0,'Voltage');
analogOutSignal1 = zeros(size(onePulse(1,:)')); %load zeros to make sure light is off at start of experiment
preload(sa, analogOutSignal1); %load zero signal
 % 
AIch = addinput(sa,'Dev1',0:2,'Voltage'); %set up two analog inputs (voltage, current)
AIch(1).TerminalConfig = 'SingleEnded';
AIch(2).TerminalConfig = 'SingleEnded';
%AIch(3).TerminalConfig = 'SingleEnded';


%% Run trial using analog signal
%dataIn = sa.startForeground; %send the data to NIDAQ analog out
dataIn = readwrite(sa,onePulse',"OutputFormat","Matrix");

%% Collect data
voltage = dataIn(:,1)*data(nn).scaleVoltage;
current = dataIn(:,2)*data(nn).scaleCurrent;

data(nn).Vm = voltage; 
data(nn).I = current;  

%% Save data, avoiding re-writing over an existing file if the user is providing a repeat input!
origExptNumber = expNumber; incrFileNum = 0;
currExptNum = expNumber;
while exist([saveStr '.mat'], 'file')
    %display(['file already exists: ' saveStr])
    currExptNum = currExptNum+1;
    fileStr = [dateStr '_SealTest_' num2str(currExptNum)];
    saveStr = [dirName fileStr];
    incrFileNum = 1;
end
if incrFileNum == 1 %let the user know that a new expt number has been assigned to this seal test
    display(['**Seal test expt number ' num2str(origExptNumber) ' already exists. Saving as expt ' num2str(currExptNum) '!**'])
end
save(saveStr, 'data');

figure; hold on
plot(voltage, 'k')
plot(current, 'b')


%% Plot the seal test (raw data)
%PlotSingleTrial(data(nn), nn, 'SealTest'); 
%% Show average stats for this seal
dateStr
currExptNum
ComputeCellStats(dateStr, currExptNum);

