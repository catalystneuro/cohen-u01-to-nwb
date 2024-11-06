function rec = windtrial
%AMS - edited from "windtrial_phantom" in "Antennae_Opto_0623"
%daq setup is for opto trials - may need to modify later
%to use with new photron camera

%returns struct (rec) as fly2 file

%need to reset daq, don't need to clear all b/c you're in a function
daqreset

%% DEFINE VARIABLES HERE FOR EACH TRIAL
%fly info for save file   
flynum = 16;
flytrial = 1;    
flycross = 'TESTING'; 
%directory name UPDATE IF NEEDED
datafolder = "C:\Users\Foxlab\Dropbox\Fox_lab\Amy\Antennae_windtrials";

%stimulus info (in seconds) 

stimtime = .5; %length of stimulus 
stimint = 1; %interval between stimuli/after stimulus
stimwait = 0.25; %time before first stim
stimrep = 1; %repetitions 

%leaving opto stuff in FOR NOW because of how daq is set up!
%set stimvolt to 0 so it doesn't go off
stimvolt = 0; 

%% END VARIABLE EDITS


%% MESSAGE BOX 
%this is just copy/pasted from the help page for the inputdlg function
%not as good as a gui like xdrift b/c doesn't save overwritten variables
prompt = {'Cross', 'Fly Number','Trial Number', 'Stimulus Length', 'Inter-Stimulus Interval',...
    'Lag Time', 'Repetitions', 'Opto Volt'};
dlgtitle = 'Trial Data';
fieldsize = [1 30; 1 30; 1 30; 1 30; 1 30; 1 30; 1 30; 1 30];
%default values
definput = {flycross, num2str(flynum), num2str(flytrial), num2str(stimtime), num2str(stimint),...
    num2str(stimwait), num2str(stimrep), num2str(stimvolt)};
%creates a cell array
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
%if there's something empty then cancel trial
if isempty(answer)
    return
end

%rewrite data
flycross = answer{1};
flynum = str2double(answer{2});
flytrial = str2double(answer{3});
stimtime = str2double(answer{4});
stimint = str2double(answer{5});
stimwait = str2double(answer{6});
stimrep = str2double(answer{7});
stimvolt = str2double(answer{8});


%define total stim length
stimlength = stimwait + stimrep*(stimtime+stimint); 

%% Make DAQ objects & channels  
%INITIATE DAQ
d = daq('ni'); 
d.Rate = 10000; %10k rate

%ADD DAQ CHANNELS
fastectrig = addoutput(d, "Dev2","ao0","Voltage");
optotrig = addoutput(d,"Dev2","ao1","Voltage"); %DISCONNECTED
phantomtrig = addoutput(d,"Dev2","ao2","Voltage");
windtrig = addoutput(d, "Dev2","ao3","Voltage"); %WIND

CamSync = addinput(d, "Dev2", "ai0","Voltage"); 
CamTrigger = addinput(d, "Dev2", "ai1","Voltage"); 
OptoTrigger = addinput(d, "Dev2", "ai2","Voltage");
LWingbeatAmp = addinput(d, "Dev2", "ai3","Voltage");
RWingbeatAmp = addinput(d, "Dev2", "ai4","Voltage");
WingbeatFreq = addinput(d, "Dev2", "ai5","Voltage");
LHutchen = addinput(d, "Dev2", "ai6","Voltage");
RHutchen = addinput(d, "Dev2", "ai7","Voltage");
%there's no AI 8-15 on the daq?
PTrigger = addinput(d, "Dev2", "ai16","Voltage");
WTrigger = addinput(d, "Dev2", "ai17","Voltage");

%START SUB-FUNCTION (this would return "scandata" and not rec!)

%CREATE TRIGGER/INPUT VECTORS
%camera trigger is the stimlength multiplied by the rate, then UP about
%.25sec after? could use the stimtime for the trigger but if you want the
%stimulus to be longer then the file will be too big. 
%all input vectors need to be single columns, and all the same size!
%camera trigger - goes up a the end of the trial
%multiply by 3.3 because the fastec cameras want a voltage of 3.3
trigtime = 0.25; %time after trigger/when trig is 1! using 250 ms so 
camtrig = [zeros(stimlength*d.Rate,1); ones(trigtime*d.Rate,1)]*3.3;
%Phantom camera is TTL! needs 5V to trigger, won't trigger at 3.3
phantomtrig = [zeros(stimlength*d.Rate,1); ones(trigtime*d.Rate,1)]*5;
totallength = stimlength + trigtime;
%KEEPING THIS IN FOR NOW, OPTO LIGHT NOT IN PLACE
%opto trig, might have variable sizes and lengths so add with a for loop
otrig = [zeros(stimwait*d.Rate,1)]; %add rate time
for i=1:stimrep %COULD EDIT
    otrig = [otrig; ones(stimtime*d.Rate,1); zeros(stimint*d.Rate,1)];
end
%need to add to the end so it's the same length as the trigger
otrig = [otrig; zeros(0.25*d.Rate,1)];
%multiply by voltage if we want to vary intensity (change LED driver to
%"mod")
otrig = otrig*stimvolt;

%WIND TRIGGER
wtrig = [zeros(stimwait*d.Rate,1)]; %add rate time
for i=1:stimrep %COULD EDIT
    wtrig = [wtrig; ones(stimtime*d.Rate,1); zeros(stimint*d.Rate,1)];
end
wtrig = [wtrig; zeros(0.25*d.Rate,1)]; 
windvolt = 9; %solenoid meant to run at 7V, doesn't react until 9V
%solenoid power source should be 12V btw
wtrig = wtrig*windvolt; 


%indicator LED

%% GO!
%create output matrix, needs to be MxN matrix, where M is num data scans
%and N is num output channels
outdata = [camtrig otrig phantomtrig wtrig];
%readwrite function: inScanData = readwrite(d,outScanData); 
scandata = readwrite(d,outdata,"OutputFormat","Matrix");

%% CREATE STRUCT AND SAVE FILE
%create struct, needs to follow xdrift naming!
rec = struct;
%daq data
rec.daq.fs = d.Rate;
rec.daq.data = scandata;
rec.daq.tstamps = linspace(0,totallength,totallength*d.Rate)';
rec.daq.channelnames = ["CamSync" "CamTrigger" "OptoTrigger" "LWingBeatAmp" "RWingBeatAmp" "WingBeatFreq" "LHutchen" "RHutchen" "PTrigger" "WindTrigger"];

%need to create base file name here for save
%need to convert flycross from char vector to string
flycross = convertCharsToStrings(flycross);
ffn = flycross + '_' + string(datetime('now'),'yyMMdd_HHmmss') + '_f' + string(flynum) + '_r' + string(flytrial);

%trial info, can also be figured out from daq data
rec.trial.info = ffn; %trial info (filename)
rec.trial.duration = stimlength;
rec.trial.stimlength = stimtime;
rec.trial.interval = stimint;
rec.trial.delay = stimwait;
rec.trial.reps = stimrep;
rec.trial.flynum = flynum;
rec.trial.flytrial = flytrial;
rec.trial.stimvolt = stimvolt;
rec.trial.windvolt = windvolt;

%SAVE
%create trial folder, this will also create the cross folder if it doesn't
%exist! just make sure the datafolder is correct!
%directory file name
dfn = datafolder + filesep + flycross + filesep + ffn;
mkdir(dfn);
%save file name (fly2 file)
sfn = dfn + filesep + ffn + '.fly2';
save(sfn,'rec');

%END SUB-FUNCTION HERE?

%PLOT TO CHECK WBA
wbaL = rec.daq.data(:,4);
wbaR = rec.daq.data(:,5);
wbf = rec.daq.data(:,6);
plot(wbaL - wbaL(1))
hold on
plot(wbaR - wbaR(1))
plot(wbf - wbf(1))
hold off

end
