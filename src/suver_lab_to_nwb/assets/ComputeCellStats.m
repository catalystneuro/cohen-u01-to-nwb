
%%
% Compute various statistics about cell access given
%  seal test signals (current + voltage).
%
% last updated on 04/22/2024 by Marie
%%

function [] = ComputeCellStats(date, expt)
if nargin < 2
   display('Please enter: (1) date (e.g. ''2016_02_10''')
   display('              (2) expt number (e.g. 2)')
   return
end
AM_SYSTEMS = 1; %if 1, will adjust a threshold (AM systems produces 10mV pulse, Axon 5mV)
TEST = 0;
dir = 'C:\Users\suver lab\Documents\Olivia\Data\';
filename = [dir '\' date '\' date '_SealTest_' num2str(expt)]
load(filename)

a = 1; %almost always have only data sample of this test
iResp = data(a).I; %recorded current trace in response to test pulses
vInj = data(a).Vm; %injected voltage step, i.e. 10mV from Vcomp signal on A-M Systems 2400


%% determine on/off locations of voltage pulse, and plot
if AM_SYSTEMS
    THRESH = 3; %this is hard-coded and works for the 10mV pulse we typically use! helps avoid wiggly start of traces we sometimes see
    SCALE = 10; %matches amp @ 100MOhm and x10 fixed output for AM 2400
else
    THRESH = 1;
    SCALE = 1000; %matches Axon
end
onInds = find(diff(vInj)>THRESH);
if TEST plot(diff(vInj)); end
offInds = find(diff(vInj)<-THRESH);
%ignore the very beginning on/off hits (wiggly start)
onInds = onInds(onInds > 50);
offInds = offInds(offInds > 50);
MIN_LEN = 70; %in samples; again, hard-coded because we sample at 10KHz and the pulse from the amp is fixed-width

%% Plot it
if TEST
    figure;
    hold on
    plot(vInj, 'k')
    plot(iResp, 'b')
    ylabel('mV'); xlabel('samples')
end

while offInds(1) < onInds(1) %first ind is off
    offInds = offInds(2:end);
end
%first pass - get rid of double-hits, which happen quite regularly with this signal (diff)
onInds(diff(onInds)<5) = nan;
offInds(diff(offInds)<5) = nan;
onInds = onInds(~isnan(onInds));
offInds = offInds(~isnan(offInds));
if TEST
    line([onInds onInds],[max(vInj) max(vInj)-10], 'Color', 'g')
    line([offInds offInds],[max(vInj) max(vInj)-10], 'Color', 'r')
end
if length(onInds) ~= length(offInds) %if there is still an uneven number, something may be off
    if length(onInds) > length(offInds)
        display('Trimming last onInd, which is extra')
        onInds = onInds(1:end-1); %last onInd is extra
    end
    if length(onInds) ~= length(offInds)
        display('Number of onInds and offInds is not equal. Investigate on/off ind detection!')
        return
    end
end

%% Plot it
if TEST
    figure;
    hold on
    plot(vInj, 'k')
    plot(iResp, 'b')
    ylabel('mV'); xlabel('samples')
    line([onInds onInds],[max(vInj) max(vInj)-10], 'Color', 'g')
    line([offInds offInds],[max(vInj) max(vInj)-10], 'Color', 'r')
end

%% grab all snippets of Im around voltage pulse
prePostSamples = 50; %grab 50 ms before and after pulse
scrsz = get(0,'ScreenSize'); %[left,bottom,width,height]
ll = scrsz(4)*0.01; bb = scrsz(4)*0.1; ww = scrsz(3)/4; hh = scrsz(4)*0.8;
figure('Color','w', 'Position',[ll,bb,ww,hh]);
ImSnip = []; VmSnip = []; jj = 1;
for ii = 1:length(onInds) %skip first and last ones, often not enough pre-/post-sample time
    if ii ~= 1 && ii ~= length(onInds)
        onInd = onInds(ii)-prePostSamples;
        offInd = offInds(ii)+prePostSamples;
        lenTraces(jj) = length(onInd:offInd);
        meanTraces(jj) = mean(iResp(onInd:offInd));
        ImSnip{jj} = iResp(onInd:offInd);
        VmSnip{jj} = vInj(onInd:offInd);
        jj = jj+1;
    end
end
%shift longer samples to left to align
for ii = 1:length(ImSnip)
    if length(ImSnip{ii}) > min(lenTraces)
        aa = circshift(ImSnip{ii},-1);        
        ImTrace(ii,:) = aa(1:end-1);
        aa = circshift(VmSnip{ii},-1);
        VmTrace(ii,:) = aa(1:end-1);
        hold on;
        plot(ImTrace(ii,:), 'c'); plot(VmTrace(ii,:), 'r');
    else
        ImTrace(ii,:) = ImSnip{ii};
        VmTrace(ii,:) = VmSnip{ii};
    end
end

%% Plot mean over individual traces
ImAvg = mean(ImTrace);
VmAvg = mean(VmTrace);
plot(ImAvg, 'b');
plot(VmAvg, 'k');
ylabel('mV or pA')

%determine steady-state portions of pulses
lenPulse = offInds(1)-onInds(1); %length of the pulse
fractionFromMid = 0.25; %middle region (fraction) of pulse to take stats from
centerPulse = lenPulse/2;
sectionBeg = prePostSamples+centerPulse-fractionFromMid*lenPulse;
sectionEnd = prePostSamples+centerPulse+fractionFromMid*lenPulse;
midSection = ImAvg(sectionBeg:sectionEnd);
ImPulse = mean(midSection);
line([sectionBeg sectionEnd],[ImPulse ImPulse], 'Color', 'm', 'LineWidth', 2)
midSection = VmAvg(sectionBeg:sectionEnd);
VmPulse = mean(midSection);
line([sectionBeg sectionEnd],[VmPulse VmPulse], 'Color', 'c', 'LineWidth', 2)

% plot baseline section
baseVmTrace = VmAvg(1:(prePostSamples-5));
baseVm = mean(baseVmTrace);
line([1 length(baseVmTrace)], [baseVm baseVm], 'Color', 'c', 'LineWidth', 2)
baseImTrace = ImAvg(1:(prePostSamples-5));
baseIm = mean(baseImTrace);
line([1 length(baseImTrace)], [baseIm baseIm], 'Color', 'm', 'LineWidth', 2)

%Compute Ra using middle portion of trace (steady state)
ImPulse = ImPulse-baseIm
VmPulse = VmPulse-baseVm

[I0 I0x] = max(ImAvg);
plot(I0x, I0, '*m')
I0 = I0-baseIm;
Iratio = I0/ImPulse;
Ga = I0/VmPulse; %in nS (pA/mV)
Ra = 1/Ga*SCALE; %convert to mOhms  
Gm = ImPulse/(VmPulse-ImPulse/Ga); %check this
Rm = 1/Gm*SCALE; %convert to mOhms %check this
Ri = VmPulse/ImPulse*SCALE ;
title({[date ' expt ' num2str(expt) ' electrode ']; ...
    ['I0 = ' num2str(I0) 'pA, ImPulse = ' num2str(ImPulse)  ' Vpulse = ' num2str(VmPulse) 'mV'];...
    ['(Ga = ' num2str(Ga) 'nS, Gm = ' num2str(Gm) ')'];['Ra = ' num2str(Ra) 'mOhm,Rinput=' num2str(Ri) 'mOhm';]}, 'Interpreter', 'none')
