%% Analysis of Karla Bosse's Master's project
clear;
close all;

%% Required functions
% Wfilt Wcu vline
% Barycentre:
% https://de.mathworks.com/matlabcentral/answers/133395-compute-centroid-of-a-matrix

%% Inputs
%P = 2; %from 2-18
pname = 'C:\Users\Brent\Desktop\Research\Masters\Karla\';
checkMVC = 0; % check MVC files
checkEMG = 0;
tqCU = 20;    % cut-off frequency for low-pass Butterworth filter of torque
angCU = 6;    % cut-off frequency for low-pass Butterworth filter of angle
c = [0 0 0 0 0 0 0 0 0]; %MVC20-05-45;20ISO-LEN-SHO;40ISO-LEN-SHO
fsamp = 2000;
IED = 8;
%x = [1];
[b,a] = butter(4,[20 400]/(fsamp/2));
%[~,D] = bandpass(x,[10 500],1/TA.interval);
%emg = filtfilt(D,x);

%% Details
%1-s ramp (2-3 s) to reach 20% or 40% MVC
%2-s rotation (5-7 s)
%ISO-LEN-SHO @ 20-40% MVC

%% Outputs
%TRef(:,1:3) = (1) max. abs. torque matching error, (2+3) var+CV of SS hold

%% Loop through participants
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

for ff = 2:18 %8=8

    cd(pname)
    cd(folders(ff).name)
    cd('Data')
    spikePath = cd;
    load('TQangFit.mat','P');

    %% Loop through files
    files = dir;
    mvcCount = 0;
    for tt = 1:length(files)
        if contains(files(tt).name,'.mat')
            keep(tt,1) = 1;
            if contains(files(tt).name,'MVC')
                mvcCount = mvcCount+1;
            end
        else
            keep(tt,1) = 0;
        end
    end
    files(keep==0) = [];
    if contains(files(end).name,'TQangFit')
        files(end) = [];
    end

    %% Alternative way
    %     files = uipickfiles;
    %     mvcCount = 0;
    %     for tt = 1:length(files)
    %         if contains(files{tt},'MVC')
    %             mvcCount = mvcCount+1;
    %         end
    %     end

    %% Load files
    mvcCount = ceil(mvcCount/3);
    mvc = NaN(mvcCount,6);
    counter = 0;

    % Use submax trials within 1 MAD only
    ct = 0;
    for xx = 1:length(files)
        if contains(files(xx).name,'MVC','IgnoreCase',true) == 1
            ct = ct+1;
            fMVC{ct,1} = files(xx).name;
        end
    end

    current_path = cd;
    cd ..
    load("files.mat");
    filesOld = files;
    clear files
    files = vertcat(fMVC,fnUse);

    cd(current_path)

    for tt = 1:length(files)

        counter = counter+1;

        if counter == 1 && isfolder('OTB')
            cd('OTB')
            otbPath = cd;
        elseif ~isempty(dir('*.otb+'))
            otbPath = cd;
        end

        %% Load data
        %fname = files(tt).name;
        fname = files{tt};
        spikeFile = [spikePath '\' fname];
        %spikeFile = files{tt};
        if isfile(spikeFile) && ~exist('otbPath','var')
            spikeData = load(spikeFile,'Torque','Angle','SQ*','Sync');
        elseif isfile(spikeFile)
            spikeData = load(spikeFile,'Torque','Angle','Sync');
        end

        % Get sync signal and crop
        voltage = max(spikeData.Sync.values);
        start = find(spikeData.Sync.values>=voltage*0.95,1,'first');

        % If OTB was used (i.e. not the talker in Spike2)
        if ~isfield(spikeData,'SQ_1')
            %otbFile = regexprep([files(tt).name(1:end-4) '_EMG.mat'],'.mat','.otb+');
            otbFile = regexprep([files{tt}(1:end-4) '_EMG.mat'],'.mat','.otb+');
            [otb.SamplingFrequency, otb.data] = OpenOTBfilesConvFact(otbFile);
            %[otb.SamplingFrequency, otb.data, ~, bits] = OpenOTBfilesConvFact(otbFile);
            %disp(bits)
            otb.data = otb.data';
            %otb.Time(:,1) = 0:1/fsamp:1/fsamp*(length(otb.data)-1/fsamp);
            otb.Sync = otb.data(:,65);
            voltageRev = min(otb.Sync);
            startOTB = find(otb.Sync<=voltageRev*10,1,'first');

            if length(otb.data) >= 38999 && ~contains(fname,'MVC')
                last = 38999;
            elseif contains(fname,'MVC')
                last = length(otb.data)-startOTB;
            end

            otb.Data = otb.data(startOTB:startOTB+last,:);

            %% Get matlab file converted using OTB (note that timestamps are
            %% different, but timestamps are not used in this script)
            %otbFile = [otbPath '\' files(tt).name(1:end-4) '_EMG.mat'];
            % If using uiputfiles:
            %[~,otbFile] = fileparts(files{tt});
            %otbFile = [otbFile '_EMG.mat'];
            %otb = load(otbFile,'Data','Time','SamplingFrequency');
            if otb.SamplingFrequency ~= fsamp
                error(['Sampling frequency was not ' num2str(fsamp) ' Hz.'])
            end
        end

        % Check sample rate
        if 1/spikeData.Torque.interval ~= fsamp || 1/spikeData.Angle.interval ~= fsamp
            error(['Sampling frequency was not ' num2str(fsamp) ' Hz.'])
        end

        %% Define analysis windows
        s1 = fsamp;
        s2 = s1*2;                      %2 sec (onset)
        h0 = s2+s1*4+1;                 %4 sec from onset (6 s in total)
        % Hold phase
        h1 = h0+s1/2+s1*2;              %6.5 sec from onset
        h2 = h1+s1*6-1;                 %12.5 sec from onset
        % Hold phase start
        h3 = h1;                        %6.5 sec from onset
        h4 = h1+s1*1-1;                 %7.5 sec from onset
        % Hold phase end
        h5 = h2-s1*1+1;                 %11.5 sec from onset
        h6 = h2;                        %12.5 sec from onset

        %% 9 conditions
        if contains(fname,'MVC20')
            aa = 1;
        elseif contains(fname,'MVC05')
            aa = 2;
        elseif contains(fname,'MVC45')
            aa = 3;
        elseif contains(fname,'ISO') && contains(fname,'20')
            aa = 4;
        elseif contains(fname,'LEN') && contains(fname,'20')
            aa = 5;
        elseif contains(fname,'SHO') && contains(fname,'20')
            aa = 6;
        elseif contains(fname,'ISO') && contains(fname,'40')
            aa = 7;
        elseif contains(fname,'LEN') && contains(fname,'40')
            aa = 8;
        elseif contains(fname,'SHO') && contains(fname,'40')
            aa = 9;
        else
            error('Filename does not match any of the conditions.')
        end

        %% Count
        c(aa) = c(aa)+1;

        %% Synchronise data
        %time(:,1) = spikeData.Sync.times(start:start+last);
        if ~exist('last','var')
            start = 1;
            last = spikeData.Torque.length-start;
        elseif last ~= 38999
            last = spikeData.Torque.length-start;
        end

        tq(:,1) = Wfilt(spikeData.Torque.values(start:start+last),tqCU,'low',fsamp);
        %tq(:,1) = interp1(spikeData.Torque.times,tqx,time,'linear','extrap');

        ang(:,1) = Wfilt(spikeData.Angle.values(start:start+last),angCU,'low',fsamp);
        %ang = interp1(spikeData.Angle.times,angx,time,'linear','extrap');

        %% Torque data
        tqPassive = polyval(P,ang); %P1 to fit lengthening correctly
        %tqPassive = mean(tq(s:e));
        tqActive = tq-tqPassive;

        %% EMG data
        n_channels = 64; %nChans +2 because of 2 aux channels but not necessary

        % Determine length
        if isfield(spikeData,'SQ_1')
            emgLength = spikeData.SQ_1.length;
        else
            emgLength = size(otb.Data,1);
        end

        % Create empty matrix
        emg = NaN(n_channels,emgLength);
        if tt == 1
            EMG = cell(1,9);
        end

        % Put EMG data in empty matrix
        for n = 1:n_channels

            if isfield(spikeData,'SQ_1')
                emg(n,:) = spikeData.( ['SQ_' num2str(n)] ).values;
            else
                emg(n,:) = otb.Data(:,n);
            end

            % Check FFT of signals
            if checkEMG
                if c(aa) == 1
                    EMG{1,aa}(n,:) = emg(n,:);
                end

                if c(aa) ~= 1
                    EMG{1,aa}(1:n_channels,1:end+emgLength) = horzcat(EMG{1,aa},emg);
                end
            end
        end

        % Single-differential of channels from 13x5 grid
        emgGrid = [emg(1:12,:); flipud(emg(13:25,:)); emg(26:38,:);...
            flipud(emg(39:51,:)); emg(52:64,:)]; % Rearrange
        if exist('x','var')
            emgGrid(x,:) = 0; % Remove invalid channels
        end
        emgBi2(1,:) = emgGrid(63,:);
        emgBi2(2,:) = emgGrid(60,:);
        emgBi = diff(emgBi2);
        emgBi = filtfilt(b,a,emgBi);
        %emgRectBi{c(aa),aa}(1,:) = rectify(emgBi,fsamp);
        emgRMSBi{c(aa),aa}(1,:) = rmsDC(emgBi, 250, 1, fsamp);
        emgDiff = diff(emgGrid);
        emgDiff([12 25 38 51],:) = [];
        emgDiff(sum(isnan(emgDiff),2)>0,:) = [];

        % Filter & rectify
        for n = 1:size(emgDiff,1)
            emgDiff(n,:) = filtfilt(b,a,emgDiff(n,:));
            %emgRect{c(aa),aa}(n,:) = rectify(emgDiff(n,:),fsamp);
            emgRMS{c(aa),aa}(n,:) = rmsDC(emgDiff(n,:),250,1,fsamp);
        end

        %% MVC files
        if aa < 4

            % Data
            [mvc(c(aa),aa),tqMaxX] = max(tqActive); % maximum torque value
            mvc(c(aa),aa) = round(mvc(c(aa),aa),1);
            mvc(c(aa),aa+3) = round(ang(tqMaxX),1); % angle at maximum torque

            % Plot
            if checkMVC
                h = figure(1); set(h,'Name',files(tt).name);
                subplot(211); hold on; plot(tq,'-'); plot(tqPassive,':'); plot(tqActive,'.');
                hline(0);
                vline(tqMaxX,'k'); %vline(emgMaxX,'r');
                %vline(s,'c'); vline(e,'c');
                subplot(212); plot(ang);
                vline(tqMaxX,'k'); %vline(emgMaxX,'r');

                waitfor(h)
            end

        else

            % Torque matching data
            tqMean = mean(tqActive(h1:h2));
            tqMatch = tqActive(h1:h2)-tqMean;
            TRef{aa-3}(c(aa),1) = max(abs(tqMatch));
            TRef{aa-3}(c(aa),2) = std(tqActive(h1:h2));
            TRef{aa-3}(c(aa),3) = std(tqActive(h1:h2))/tqMean*100;

            % Data over time
            tqCond{aa-3}(c(aa),:) = NaN(1,length(tq));
            angCond{aa-3}(c(aa),:) = NaN(1,length(ang));
            emgCond{aa-3}(c(aa),:) = NaN(1,emgLength);
            emgCondBi{aa-3}(c(aa),:) = NaN(1,emgLength);
            tqCond{aa-3}(c(aa),:) = tqActive;
            angCond{aa-3}(c(aa),:) = ang;
            emgCond{aa-3}(c(aa),:) = mean(emgRMS{c(aa),aa},'omitnan');
            emgCondBi{aa-3}(c(aa),:) = emgRMSBi{c(aa),aa};

            % Ramp data
            rampEMG{aa-3}(c(aa),:) = mean(emgRMS{c(aa),aa}(:,2*fsamp+1:3*fsamp),2);
            rampEMGbi{aa-3}(c(aa),:) = mean(emgRMSBi{c(aa),aa}(:,2*fsamp+1:3*fsamp),2);
            rampTQ(c(aa),aa-3) = mean(tqActive(2*fsamp+1:3*fsamp));

            % Rotation data
            rotEMG{aa-3}(c(aa),:) = mean(emgRMS{c(aa),aa}(:,5*fsamp+1:7*fsamp),2);
            rotEMGbi{aa-3}(c(aa),:) = mean(emgRMSBi{c(aa),aa}(:,5*fsamp+1:7*fsamp),2);
            rotTQ(c(aa),aa-3) = mean(tqActive(5*fsamp+1:7*fsamp));

            % Steady-state data
            holdEMG{aa-3}(c(aa),:) = mean(emgRMS{c(aa),aa}(:,h1:h2),2);
            holdEMGbi{aa-3}(c(aa),:) = mean(emgRMSBi{c(aa),aa}(:,h1:h2),2);
            holdEMGcv{aa-3}(c(aa),:) = std(emgRMS{c(aa),aa}(:,h1:h2),0,2)./mean(emgRMS{c(aa),aa}(:,h1:h2),2)*100;
            holdEMGcvbi{aa-3}(c(aa),:) = std(emgRMSBi{c(aa),aa}(:,h1:h2))./mean(emgRMSBi{c(aa),aa}(:,h1:h2))*100;
            holdEMGmf{aa-3}(c(aa),:) = medfreq(emgDiff(:,h1:h2)',fsamp);
            holdEMGmfbi{aa-3}(c(aa),:) = medfreq(emgBi(:,h1:h2),fsamp);
            holdSS(c(aa),aa-3) = mean(tqActive(h1:h2));
            holdSS(c(aa),aa+3) = mean(ang(h1:h2));
            holdSE(c(aa),aa-3) = mean(tqActive(h3:h4));
            holdSE(c(aa),aa+3) = mean(tqActive(h5:h6));

            % EMG RMS map using steady-state data (lanczos3)
            tmp = mean(emgRMS{c(aa),aa}(:,h1:h2),2);
            holdMap = [tmp([1 1:11],1) tmp(12:23,1) tmp(24:35,1) tmp(36:47,1) tmp(48:59,1)];
            holdMap(1,1) = 0;
            tmp2=holdMap/sum(holdMap(:)); % Normalise
            [m,n]=size(tmp2);
            [I,J]=ndgrid(1:m,1:n); % Coordinates
            holdEMGmap{aa-3}(c(aa),:) = ([dot(I(:),tmp2(:)),  dot(J(:),tmp2(:))]-1)*IED; % Centroid
            %regionprops(true(size(holdMap)), holdMap,  'WeightedCentroid'); %X & Y flipped
        end

        clearvars -except P counter otbPath spikePath ff pname checkMVC tqCU angCU c mvc P folders files TRef holdE* holdS* tqCond angCond emgCond fsamp IED EMG x checkEMG emgR* ramp* rot* b a
    end

    if exist('mvc','var')
        [T,Tx] = max(mvc(:,1));
    end

    SS = round(holdSE(:,7:12)-holdSE(:,1:6),1);

    SSx = find(abs(SS)>=1.0);
    if ~isempty(SSx)
        disp('Mean torque diff. in first and last seconds is >= 1 Nm for some trials: check SSx.')
    end

    if checkEMG
        for jj = 1:size(EMG,2)

            h = plotHDEMG(EMG{1,jj},fsamp);
            waitfor(h);

            %% Delete noisy channels
            %     prompt = "Which channels should be deleted? ";
            %     if tt == 1
            %         x = input(prompt);
            %     end
            %     EMG{:}(x,:) = NaN;
        end
    end

    % Find max EMG
    %[~,emgMaxX] = max(emgRMS{Tx,1},[],2);
    [~,emgMaxX] = max(emgRMS{Tx,1}(:,500:end),[],2);
    emgMaxX = emgMaxX+500-1;
    plotHDEMG(emgRMS{Tx,1},emgMaxX); 
    [~,emgMaxBiX] = max(emgRMSBi{Tx,1},[],2);
    Ebi(1,:) = mean(emgRMSBi{Tx,1}(emgMaxBiX-499:emgMaxBiX+499));

    for aa = 1:length(emgMaxX)
        E(1,aa) = mean(emgRMS{Tx,1}(aa,emgMaxX(aa,1)-499:emgMaxX(aa,1)+499));
    end

    holdSS(holdSS==0) = NaN;

    for bb = 1:6
        AVGrampEMG(:,bb) = mean(rampEMG{1,bb})./E*100;
        AVGrampEMGbi(:,bb) = mean(rampEMGbi{1,bb})./Ebi*100;
        AVGrampTQ(1,bb) = round(mean(rampTQ(:,bb))/T*100,2);
        AVGrotEMG(:,bb) = mean(rotEMG{1,bb})./E*100;
        AVGrotEMGbi(:,bb) = mean(rotEMGbi{1,bb})./Ebi*100;
        AVGrotTQ(1,bb) = round(mean(rotTQ(:,bb))/T*100,2);
        AVGholdEMG(:,bb) = mean(holdEMG{1,bb})./E*100;
        AVGholdEMGbi(:,bb) = mean(holdEMGbi{1,bb})./Ebi*100;
        AVGholdEMGcv(:,bb) = mean(mean(holdEMGcv{1,bb}));
        AVGholdEMGcvbi(:,bb) = mean(holdEMGcvbi{1,bb});
        AVGholdEMGmf(:,bb) = mean(mean(holdEMGmf{1,bb}));
        AVGholdEMGmfbi(:,bb) = mean(holdEMGmfbi{1,bb});
        AVGholdEMGmapRow(:,bb) = round(mean(holdEMGmap{1,bb}(:,1)),1);
        AVGholdEMGmapCol(:,bb) = round(mean(holdEMGmap{1,bb}(:,2)),1);
        AVGholdTQ(1,bb) = round(mean(holdSS(:,bb),'omitnan')/T*100,2);
    end

    cd(pname)
    cd(folders(ff).name)
    save('analysed_new.mat');

    clearvars -except pname checkMVC checkEMG tqCU angCU fsamp IED folders
    c = [0 0 0 0 0 0 0 0 0]; %MVC20-05-45;20ISO-LEN-SHO;40ISO-LEN-SHO
    [b,a] = butter(4,[20 400]/(fsamp/2));

end