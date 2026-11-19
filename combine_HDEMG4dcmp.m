%% Combine trials for HDEMG decomposition
% Treat each contraction intensity separately
% and later compare reference conditions from the two intensities.
clear; close all

%% Inputs
pname = 'C:\Users\Brent\Desktop\Research\Masters\Karla\';
tqCU = 20;    % cut-off frequency for low-pass Butterworth filter of torque
fsamp = 2000;
n_channels = 64; %nChans +2 because of 2 aux channels but not necessary
IED = 8;
P = 23; % folder number (2=1, 3=2, ..., 8=8, 9=9, ...)
dur = 36999; % duration of the file in data pts
fDecomp = {'MVC20.mat' 'MVC40.mat'};
extra = 0.5*fsamp;

%% Get folders
cd(pname);
folders = dir;
for ff = 1:length(folders)
    if contains(folders(ff).name(1),'0') ||  contains(folders(ff).name(1),'1')
        keep(ff,1) = 1;
    else
        keep(ff,1) = 0;
    end
end
folders(keep==0) = [];
clear keep
cd(pname)

%% Go to desired folder and loop through files in that folder
cd(folders(P).name)
cd('Data')
spikePath = cd;
filesAll = dir;
files(:,1) = extractfield(filesAll,'name');
MVC100 = contains(files,'MVC');
f100 = files(MVC100);
MVC20 = contains(files,'_20');
if sum(MVC20) == 0
    MVC20 = contains(files,'20_');
    MVC20(MVC100==1) = 0;
end
f20 = files(MVC20);
MVC40 = contains(files,'_40');
if sum(MVC40) == 0
    MVC40 = contains(files,'40_');
    MVC40(MVC100==1) = 0;
end
f40 = files(MVC40);
ref = contains(files,'ISO');
fRef = files(ref);
fn = [f20;f40];

%% Remove outlier trials based on 1 median absolute deviation
cd ..
load('analysed.mat','holdEMG','files');
cd('Data')
c = 0;
a = 0;
for bb = 1:6
    out = isoutlier(mean(holdEMG{1,bb},2),'median',ThresholdFactor=1);
    lose = find(out==1);
    lose = lose+c;
    bye{bb} = lose;
    if bb == 1
        byeFile = lose;
    else
        byeFile = vertcat(byeFile,lose);
    end
    c = c+size(holdEMG{1,bb},1);
    clear out lose
end

fnUse = fn;
if c < length(fn)
    cya = contains(fn,'otb+');
    fnUse(cya) = [];
    FN = fnUse;
end
fnUse(byeFile) = [];
MVC20 = contains(fnUse,'20');
f{1} = fnUse(MVC20);
MVC40 = contains(fnUse,'40');
f{2} = fnUse(MVC40);
if sum(MVC20) + sum(MVC40) ~= length(fnUse)
    error('Something went wrong with file determination.')
end

path_current = cd;
cd ..
save('files.mat',"fnUse");
cd(path_current)

% Explicitly reorder files to SHO, STR, REF
%f{3} = [f{1}(7:11); f{1}(3:6); f{1}(1:2)];
%f{1} = f{3};

for bb = 1 %1:2
    for tt = 1:length(f{bb})

        if bb == 2 && tt == 1
            clear EMG emg ref_signal
        end
        % Go to Spike2 folder
        cd(spikePath)

        % Get Spike2 data
        fLoad = char(f{bb}(tt));
        spikeData = load(fLoad,'Torque','Sync');

        % Crop Spike2 data based on sync signal
        voltage = max(spikeData.Sync.values);
        start = find(spikeData.Sync.values>=voltage*0.95,1,'first');

        % Get OTB data
        if isfolder('OTB')
            cd('OTB')
        end
        otbFile = regexprep([fLoad(1:end-4) '_EMG.mat'],'.mat','.otb+');
        [otb.SamplingFrequency, otb.data] = OpenOTBfilesConvFact(otbFile);
        otb.data = otb.data';

        % Crop OTB data based on sync signal
        otb.Sync = otb.data(:,65);
        voltageRev = min(otb.Sync);
        startOTB = find(otb.Sync<=voltageRev*10,1,'first');
        if length(otb.data) >= dur
            last = dur;
        end
        otb.Data = otb.data(startOTB:startOTB+last,:);

        %% Reference signal (torque data)
        if ~exist('last','var')
            start = 1;
            last = spikeData.Torque.length-start;
        elseif last ~= dur
            last = spikeData.Torque.length-start;
        end

        %% Crop again so just the steady state is used
        s1 = fsamp;
        s2 = s1*2;                      %2 sec (onset)
        h0 = s2+s1*4+1;                 %4 sec from onset
        % Hold phase
        h1 = h0+s1/2+s1*2;              %6.5 sec from onset
        h2 = h1+s1*6-1;                 %12.5 sec from onset
        %h1 = 1; %10001
        %h2 = dur; %29000

        %% Torque data - combine files with vertcat
        if ~exist('ref_signal','var')
            %ref_signal(:,1) = Wfilt(spikeData.Torque.values(start:start+last),tqCU,'low',fsamp);
            ref_signal(:,1) = Wfilt(spikeData.Torque.values(h1-extra:h2),tqCU,'low',fsamp); 
        elseif exist('ref_signal','var') && tt == length(f{bb})
            tq(:,1) = Wfilt(spikeData.Torque.values(h1:h2+extra),tqCU,'low',fsamp); 
            ref_signal = vertcat(ref_signal,tq);
        else
            %tq(:,1) = Wfilt(spikeData.Torque.values(start:start+last),tqCU,'low',fsamp);
            tq(:,1) = Wfilt(spikeData.Torque.values(h1:h2),tqCU,'low',fsamp);
            ref_signal = vertcat(ref_signal,tq);
        end

        %% EMG data - combine files with horzcat
        % Determine length
        %emgLength = size(otb.Data,1);
        emgLength = 12000; 
        %emgLength = dur; %19000
        % Create empty matrix
        if ~exist('EMG','var') 
            EMG = NaN(n_channels,emgLength+extra);
            check = 1;
        elseif exist('EMG','var') && tt == length(f{bb})
            emg = NaN(n_channels,emgLength+extra);
            check = 2;
        else
            emg = NaN(n_channels,emgLength);
            check = 2;
        end
        % Put EMG data in empty matrix and combine files with horzcat
        for n = 1:n_channels
            if check == 1
                %EMG(n,:) = otb.Data(:,n);
                EMG(n,:) = otb.Data(h1-extra:h2,n);
            elseif check == 2 && tt == length(f{bb})
                emg(n,:) = otb.Data(h1:h2+extra,n);
            else
                %emg(n,:) = otb.Data(:,n);
                emg(n,:) = otb.Data(h1:h2,n);
            end
        end

        if check == 2
            EMG = horzcat(EMG,emg);
        end

        h = figure;
        if ~exist('emg','var')
            emg = EMG;
        end
        plot(EMG(1,:));
        %emgPlot = abs(EMG(1,:));
        %plot(mean(emgPlot));
        waitfor(h);

        clearvars -except bb EMG extra dur f* FN P IED n_channels pname ref_signal spikePath tqCU tt
    end

    %% Store required variables
    signal.data = EMG;
    signal.fsamp = fsamp;
    signal.nChan = size(signal.data,1);
    signal.ngrid = 1;
    signal.gridname = 'GR08MM1305';
    signal.muscle = {'Tibialis Anterior'};
    signal.torque(1,:) = ref_signal;
    signal.files = f{bb};
    signal.edge = extra/fsamp;

    %% Save signal variable in "Decomp" folder
    cd(pname)
    cd(folders(P).name)
    if ~isfolder('Decomp')
        mkdir('Decomp')
    end

    cd('Decomp')
    cd('Test')

    % Save data for MUedit decomposition
    save([fDecomp{bb}(1:end-4)],'signal')

    % Save data for SCD decomposition
    %emg = signal.data;
    %save(fDecomp{bb},'emg')

    clear signal
end