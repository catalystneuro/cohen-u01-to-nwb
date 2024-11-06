%INPUT: single trial folder, plot (0 or 1) 
%analyzes all data within a file and saves struct as .pcf
%can use this per trial, called as part of "folderanalysis"
%toplot and tosave should =1 if yes

function trialdata = trialanalysis_photron(trial_path, env, toplot, tosave)

%if not all the arguments, redefine
%we ain't about checking for special cases
if nargin < 3
    trial_path = uigetdir();
    toplot = 1;
    env = 40; %CHANGE ENVELOPE SIZE HERE if running single trial
    tosave = 0;
end

trialdata = struct;

%load in data
%daq
fd = dir([trial_path filesep '*.fly2']);
ddata = load(fullfile(fd.folder, filesep, fd.name),'-mat');
ddata = ddata.rec;

%haltere data
hd = dir([trial_path filesep '*_haltamp.mat']);
hdata = load(fullfile(hd.folder, filesep, hd.name));

%antenna/head data
ad = dir([trial_path filesep '*_antennatracking*.csv']);
%make a variable to know if head fixed or head free
% headfree = 0;
% if isempty(ad)
%     ad = dir([trial_path filesep '*_headtracking*.csv']);
%     headfree = 1;
% end
adata = readtable(fullfile(ad.folder, filesep, ad.name), 'NumHeaderLines',1);

%wing data
wd = dir([trial_path filesep '*_PROC.mat']);
wdata = load(fullfile(wd.folder, filesep, wd.name));

%antenna columns: 
%1 = coords (frame number), then x y likelihood for each
%2-4 = L_ant; 5-7 = R_ant; 8-10 = L_base; 11-13 = R_base
%head tracking:
%14-16 = L_OC; 17-19 = R_OC; 20-22 = L_pOR; 23-25 = R_pOR; 
%26-28 = L_aVT; 29-31 = R_aVT

%ASSIGN DATA VARIABLES
%daq
%trial info
duration = ddata.trial.duration;
stimlength = ddata.trial.stimlength;
flynum = ddata.trial.flynum;
trialnum = ddata.trial.flytrial;
volt = ddata.trial.stimvolt;
%data
stamps = ddata.daq.tstamps;
ctrig = ddata.daq.data(:,2); %fastec cam trigger
%otrig = ddata.daq.data(:,3); %opto trigger
ptrig = ddata.daq.data(:,9); %photron trigger
wtrig = ddata.daq.data(:,10); %wind trigger

%NOT CURRENTLY USING
%including in struct though so it's there if needed?
%otrig = ddata.daq.data(:,3); %opto trigger
dwbf = ddata.daq.data(:,6);
%wbf data from wingbeat amplifier needs to be multiplied by 100: 
% 1V = 100Hz
dwbf = dwbf*100; %daq wingbeat freq, added 7/17
dwbaL = ddata.daq.data(:,4); %temp
% twbaL = twbaL - mean(twbaL(1:2499)); %set to 0 for denoising
dwbaR = ddata.daq.data(:,5);
% twbaR = twbaR - mean(twbaR(1:2499));
% hutchenL = ddata.daq.data(:,7); 
% hutchenR = ddata.daq.data(:,8); 

%ALIGN CAMERAS TO DAQ - CHECK WITH MIKE AGAIN?
%find trigger times
ctrigtime = find(ctrig>3,1)/ddata.daq.fs;
ptrigtime = find(ptrig>3,1)/ddata.daq.fs;

%topcam (currently skipping sidecam!)
meta_t = fastecMetaReader(fullfile(trial_path, 'TOPCAM_000000.txt'));
ts_t = linspace(1/meta_t.fs, meta_t.numframes/meta_t.fs, meta_t.numframes);
ts_t = ts_t-(ts_t(end)-ctrigtime);
%photron/phantom camera
mp = dir([trial_path filesep 'HALTCAM*.mii']); %since the number changes every time
meta_p = photronMetaReader(fullfile(mp.folder, mp.name));
ts_p = linspace(1/meta_p.fs, meta_p.numframes/meta_p.fs, meta_p.numframes);
ts_p = ts_p-(ts_p(end)-ptrigtime);

%% ANTENNAE
%subtract height from y data because y = 0 is at top of image
Lx = adata.L_ant;
Ly = meta_t.height - adata.L_ant_1;
Rx = adata.R_ant;
Ry = meta_t.height - adata.R_ant_1;
%base data just in case
Lxb = adata.L_base;
Lyb = meta_t.height - adata.L_base_1;
Rxb = adata.R_base;
Ryb = meta_t.height - adata.R_base_1;
%find theta
Lt = atan2d((Ly-Lyb),(Lx-Lxb));
Lt = wrapTo360(Lt);
Rt = atan2d((Ry-Ryb),(Rx-Rxb));
%find distance between L and R
dant = [];
for j = 1:length(Lx)
    dant(j) = pdist2([Lx(j) Ly(j)],[Rx(j) Ry(j)]);
end  
dant = dant';

%% WING DATA
%denoise daq WBA freq data
dwbf = hampel(dwbf);
%flyalyzer data
wbL = wdata.fly.wingL.angle;
% wbLr = wdata.fly.wingL.root;
wbR = wdata.fly.wingR.angle;
%need to wrap R wing data
wbR = wrapTo180(wbR);
% wbRr = wdata.fly.wingR.root;
%create timestamps: flyalyzer tstamps aren't good, it's always every 20 
wstamps = linspace(ts_t(1),ts_t(end),length(ts_t)/20);

%zeroing removed 7/17: will be zeroed later!
%zero everything, camera upix is 1351 -> 1350/20 = 67.5 - see
%"cameratimes_phantom.m" for calculation
% wbL = wbL-mean(wbL(1:60));
% wbR = wbR-mean(wbR(1:60));

%7/17 alright we're gonna modify the freq calc bc something's wrong
%NOTE: Mike's code, this sums up image pixels to get the frequency. This
%does include the haltere because it's in the video, but the wings are much
%bigger so it shouldn't make too much difference.
%MODIFIED SUMVID
vr = VideoReader([trial_path filesep 'TOPCAM_000000.avi']);
%create array for video sums
int = zeros(1,vr.NumFrames);
%go through frames
for j = 1:vr.NumFrames
    %read image, greyscale but saves as color
    img = read(vr,j);
    %take only red channel since they're all the same
    img = img(:,:,1);
    %add up all the pixels in that image and add to time series
    int(j)=sum(sum(img));
end
%pxx is made of FFTs for each timestamp, f is freq, t is time
%zscore gets rid of DC info (low freq)
%window size is 64 w 63 overlap (or 128/127), or use 20/0 so same
%sample number as flyalyzer
%freq range is 150-300, 2000 is fps, 'reassigned' does a thing
%reassigned is kinda like contrast stretching?
[pxx,f,t]=spectrogram(zscore(int),64,63,150:300,2000,'reassigned');
%time/freq ridge is pretty much the wbf plot, has a .001 penalty to
%make sure there's not much difference between y values (so it doesn't
%jump)
%NOTE: there might be some weird harmonics if the penalty isn't small
%enough, can check by plotting t,calc_wbf 
calc_wbf = tfridge(pxx,f,.001);

%NOTE: there might be some weird harmonics if the penalty isn't small
%enough, can check by plotting t,calc_wbf over top of the spectrogram to
%check if it lines up
    % spectrogram(zscore(int),128,127,150:300,2000,'yaxis','reassigned');
    % plot(t,calc_wbf,'r')
%add NaNs back in to beginning so same length as timestamps
calc_wbf = [NaN((length(int) - length(calc_wbf)),1); calc_wbf];

%% Haltere Data
%in ampalyzer, haltere base needs to be to the right of the haltere
%base, not in the center! otherwise lots of weird artefacts
%haltere data
%get coordinates and subtract y from height of frame (y=0 at top)
p_height = meta_p.height;
%haltere and wing root
hroot = hdata.haltrpos;
hroot(2) = p_height-hroot(2);
wroot = hdata.wingrpos;
wroot(2) = p_height-wroot(2);   
hx = hdata.haltpos(:,1);
hy = p_height - hdata.haltpos(:,2);
%this bit is from Mike and some trial and error...
%find a transformation matrix that turns xy in haltere root to 0,0 and xy in
%wing root to 1,0
tform = fitgeotrans([hroot;wroot],[0,0;1,0],'nonreflectivesimilarity');
%use that transformation matrix for the tracked points
[hxn, hyn] = transformPointsForward(tform,hx,hy);
%get angles and wrap
Ht = atan2d(hyn,hxn); %PHOTRON
%get AMPLITUDE (envelope function)
%envelope of 14 looks best for phantom, about 40 for photron -DEFINED AT TOP THIS IS WHY
%WE'RE GOING THROUGH FOLDERS DAMMIT
[Htup,Htlow] = envelope(Ht,env,'peak');
%find peaks the other way... 10k fps at 200 Hz peak-to-peak should be 50,
%so half a stroke is 25, should be just above that...
[Htu_pks, Htu_locs, Htu_w] = findpeaks(Ht, 'MinPeakHeight',20, 'MinPeakDistance',30);
[Htl_pks, Htl_locs, Htl_w] = findpeaks(-Ht, 'MinPeakHeight',20, 'MinPeakDistance',30);
Htl_pks = -Htl_pks;
% %distance, normalized to distance between wing and root
% %removed 6/28 bc not a great time? kept Hd though for reference
% %alternates between top and bottom of wing stroke; every other peak has a
% %hutchen
Hd = pdist2([0 0],[hxn,hyn]);
% %dips at top of stroke....
% [Hdup,Hdlow] = envelope(Hd,env,'peak');
% [Hdu_pks, Hdu_locs, Hdu_w] = findpeaks(Hd, 'MinPeakHeight',.4, 'MinPeakDistance',30);

%get FREQUENCY
%spectrogram function basically makes an fft (power x function: pxx) for
%each time t, samples every 256 points (255 overlapping)
% 64 & 128 is too small because fps is 10k! 
%to be optimal range
%detrend gets rid of all the low freq info (dc info?)
%looking for first set of freq and no harmonics (max wbf is around 230?)
[pxx,f,~]=spectrogram(detrend(Ht),256,255,0:300,meta_p.fs);
%find time frequency ridge (kinda like the avg?)
Hfreq = tfridge(pxx,f,.001);
%since spec uses sliding window, first set of overlap is not used, add NaNs
%to make same size as timestamps
Hfreq = [NaN((length(Ht) - length(Hfreq)),1); Hfreq];

if isequal(toplot,1) %plot data
    close all
    figure('WindowState','maximized');
    t =  tiledlayout(2,3);

    nexttile
    hold on
    %yyaxis left
    plot(wstamps,wbR,'r')
    %for some reason this is dashed if you don't specify it...
    plot(wstamps,wbL,'-b') 
    % yyaxis right
    % plot(stamps,wbf)
    title('Wingbeat Amp R (red) & L (blue)')
    hold off
    
    %7/17 look at freq
    nexttile
    hold on
    yyaxis left
    ylabel('WB Amp')
    plot(stamps,dwbf)
    yyaxis right
    ylabel('Calc')
    plot(ts_t,calc_wbf)
    
    hold off
    nexttile
    hold on
    plot(ts_t,dant)
    title('Distance btwn Antennae')
    hold off

    nexttile
    hold on
    yyaxis left
    plot(ts_t,Lt,'b')
    ylabel("Ltheta")
    yyaxis right
    plot(ts_t,Rt,'r')
    ylabel("Rtheta")
    title('Antennal angle')
    hold off

    nexttile
    hold on
    plot(ts_p,Hfreq)
    title('Haltere Freq')
    hold off

    % nexttile
    % hold on
    % yyaxis left
    % plot(ts_p,(Htup-Htlow),'r')
    % ylim([150 200])
    % yyaxis right
    % plot(ts_p,(Hdup-Hdlow),'b')
    % ylim([.35 .65])
    % title('Haltere Envelope Change angle (red) & length (blue)')
    % hold off

    nexttile
    hold on
    envelope(Ht,env,'peak')
    title('Haltere angle')
    hold off

    % nexttile
    % hold on
    % plot(Ht)
    % plot(Htu_locs,Htu_pks,'o')
    % plot(Htl_locs,Htl_pks,'o')
    % title('Haltere angle peaks')
    % hold off


    % nexttile
    % hold on
    % envelope(Hd,env,'peak')
    % plot(Hdu_locs,Hdu_pks,'o')
    % plot(Hdl_locs,Hdl_pks,'o')
    % title('Haltere dist')
    % hold off

    title(t,{"Fly: " + num2str(flynum) + " Trial: " + num2str(trialnum), ...
        "Stim " + num2str(stimlength) + " Voltage: " + num2str(volt), ...
        "Envelope: " + num2str(env)});
end

 %% Add to struct
 %note: nothing is normalized to the mean yet! This is only the raw
 %data!
 %not added: hutchens
 %trial info
 trialdata.info.duration = duration;
 trialdata.info.stimlength = stimlength;
 trialdata.info.flynum = flynum;
 trialdata.info.trialnum = trialnum;
 trialdata.info.envelope = env;
 trialdata.info.volt = volt;
 %daq data
 trialdata.trial.tstamps = stamps;
 trialdata.trial.ctrig = ctrig;
 %trialdata.trial.otrig = otrig;
 trialdata.trial.wbaL = dwbaL;
 trialdata.trial.wbaR = dwbaR;
 trialdata.trial.wbf = dwbf;
 trialdata.trial.ptrig = ptrig;
 trialdata.trial.wtrig = wtrig;
 %camera data
 trialdata.ant.tstamps = ts_t';
 trialdata.ant.dist_ant = dant;
 trialdata.ant.Rt = Rt;
 trialdata.ant.Lt = Lt;
 trialdata.wing.wbL = wbL;
 trialdata.wing.wbR = wbR;
 trialdata.wing.wstamps = wstamps;
 trialdata.wing.wbf = calc_wbf;
 trialdata.halt.tstamps = ts_p';
 trialdata.halt.Ht = Ht;
 trialdata.halt.Htup = Htup;
 trialdata.halt.Htlow = Htlow;
 trialdata.halt.Hdist = Hd;
 % trialdata.halt.Hdup = Hdup;
 % trialdata.halt.Hdlow = Hdlow;
 trialdata.halt.Hfreq = Hfreq;
 trialdata.halt.Htu_pks = Htu_pks;
 trialdata.halt.Htu_locs = Htu_locs;
 trialdata.halt.Htu_w = Htu_w;
 trialdata.halt.Htl_pks = Htl_pks;
 trialdata.halt.Htl_locs = Htl_locs;
 trialdata.halt.Htl_w = Htl_w;
 % trialdata.halt.Hdu_pks = Hdu_pks;
 % trialdata.halt.Hdu_locs = Hdu_locs;
 % trialdata.halt.Hdu_w = Hdu_w;
 % trialdata.halt.Hdl_pks = Hdl_pks;
 % trialdata.halt.Hdl_locs = Hdl_locs;
 % trialdata.halt.Hdl_w = Hdl_w;

 % %head data if applicable TO DO
 % if headfree == 1
 %     trialdata.head.pitch = [];
 %     trialdata.head.yaw = [];
 %     trialdata.head.roll = [];
 % end

 if tosave == 1
     %% save
     %save analysis so it runs faster if you're playing around with the graphs
     %date & time
     %dt = string(datetime('now'),'yyMMdd_HHmmss');
     %save file name (Tsh file)
     flycross = extractBetween(trial_path,'trials\','\');
     flycross = convertCharsToStrings(flycross);
     sfn = convertCharsToStrings(trial_path) + filesep + flycross + '_f' + flynum + '_r' + trialnum + '_trial.pcf';
     save(sfn,'trialdata');
 end