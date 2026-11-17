close all
%clear

% Sample rate
fsamp = 2000;
%fsamp = 2049.122;
% Rotation duration
ssDur = round(3*fsamp);
% Rotation velocity
desiredVel = -5;
%rotDur = round(20/(abs(desiredVel)+1)*fsamp);
if desiredVel > 0
    desiredVel = desiredVel-1.5;
elseif desiredVel < 0
    desiredVel = desiredVel+1.5;
end


% PID
bb = 1;

for tt = 1
    %% Check discharge rates of each MU for each condition
    % Get trial length
    if exist('fname','var')
        load(fname{tt});
    end
    if ~exist('trialDur','var')
        trialDur = (length(signal.torque)-fsamp)/length(signal.files);
    end
    totalDur = length(signal.torque);

    % Crop trials based on condition
    if ~exist('sho','var')
        ISO = contains(signal.files,'ISO');
        iso = sum(ISO);
        LEN = contains(signal.files,'LEN');
        len = iso+sum(LEN);
        SHO = contains(signal.files,'SHO');
        sho = len+sum(SHO);
        %ISO = 2; iso = 2; LEN = 2; len = 4; SHO = 4; sho = 8;
    end
    pad = round(fsamp/2,0);
    %% Steady state only vs. dynamic phase
    %crop = [pad+1 trialDur*iso+pad; trialDur*iso+pad+1 trialDur*len+pad; trialDur*len+pad+1 trialDur*sho+pad]; % crop 7-s steady state
    % crop = [pad+1+fsamp*4 trialDur*iso+pad; trialDur*iso+pad+1+fsamp*4 trialDur*len+pad; trialDur*len+pad+1+fsamp*4 trialDur*sho+pad]; % crop 3-s steady state
    time(:,1) = 1/fsamp:1/fsamp:1/fsamp*length(signal.torque);
    vel(:,1) = Wvel(signal.angle,time);
    velSS(:,1) = find(vel<=desiredVel);
    velGaps = find(diff(velSS)>=fsamp*2);
    %velGaps = find(diff(velSS)>=rotDur);
    s = velSS([1; velGaps+1]); % start points
    e = velSS([velGaps; end]); % end points
    bye = find((e-s)<fsamp);
    s(bye) = [];
    e(bye) = [];
    %s = round((e-s)/2)+s;
    %e = round((e-s)/2)+s;
    h = figure; plot(vel); ylim([-20 20]); hold on; 
    vline(s,'r-'); vline(e,'r:');
    crop = [s e]; % crop dynamic phase
    vline(crop(:,1),'b-'); vline(crop(:,2),'g-');
    waitfor(h);

    %% Load and analyse
    % Counter
    ii = 1;

    %% Get spike train by binarizing MUPulses
    spikeTrain = binario(MUPulses,totalDur);

    % Loop through conditions
    for kk = 1:size(crop,1)

        % AND TRIALS
        if kk == 1
            trialNo = sum(ISO);
        elseif kk == 2
            trialNo = sum(LEN);
        elseif kk == 3
            trialNo = sum(SHO);
            %else
            %    error('>3 conditions in kk.')
        end

        %% Loop through each MU
        for jj = 1:length(MUPulses)

            % Find specifc MU in given condition
            MU = MUPulses{jj};

            if size(crop,1) > 3
                trialNo = 1;
            end

            for trial = 1:trialNo
                if trial == 1 % first trial
                    first = find(MU>=crop(kk,1),1,'first');
                    last = find(MU<=crop(kk,1)+ssDur-1,1,'last');
                    if jj == 1
                        firstCrop(trial,kk,jj) = crop(kk,1);
                    end
                elseif trial == trialNo % last trial
                    first = find(MU>=crop(kk,2)-ssDur+1,1,'first');
                    last = find(MU<=crop(kk,2),1,'last');
                    if jj == 1
                        firstCrop(trial,kk,jj) = crop(kk,2)-ssDur+1;
                    end
                else % all trials in between
                    first = find(MU>=(trial-1)*ssDur+crop(kk,1),1,'first');
                    last = find(MU<=(trial-1)*ssDur+crop(kk,1)+ssDur-1,1,'last');
                    if jj == 1
                        firstCrop(trial,kk,jj) = (trial-1)*ssDur+crop(kk,1);
                    end
                end
                %trial = 1; % for concatenating trials, uncomment this line and comment above lines

                if size(crop,1) > 3
                    last = find(MU<=crop(kk,2),1,'last');
                end

                first = round(first);
                last = round(last);

                if isempty(first) || isempty(last) || last <= first
                    dr(jj,trial) = NaN;
                    drnum(jj,trial) = NaN;
                    drmed(jj,trial) = NaN;
                    drcv(jj,trial) = NaN;
                    %DRbi(jj,ii) = NaN;
                else
                    DRall = 1./(diff(MU(first:last))/fsamp);
                    drnum(jj,trial) = length(MU(first:last));
                    dr(jj,trial) = mean(DRall,'omitnan');
                    drmed(jj,trial) = median(DRall,'omitnan');
                    drcv(jj,trial) = std(DRall)/mean(DRall)*100;
                    MUuse{trial,kk,jj} = MU(first:last);
                    %DRbi(jj,trial) = bimodalitycoeff(DRall);
                end

                clear DRall
            end % end of trial no. loop

            clear MU
        end % end of motor unit no. loop

        % Average values among trials for each condition (kk)
        % 1 = REF; 2 = LEN; 3 = SHO
        DRnum(:,kk) = mean(drnum,2,'omitnan'); % discharge number
        DR(:,kk) = mean(dr,2,'omitnan'); % mean discharge rate
        DRmed(:,kk) = mean(drmed,2,'omitnan'); % median discharge rate
        DRcv(:,kk) = mean(drcv,2,'omitnan'); % CV of discharge rate
        clear dr*

    end % end of condition loop

    MUunique(bb,tt) = sum(sum(isnan(DR)));
    DRnan = isnan(DR);
    for aa = 1:size(DRnan,1)
        if sum(DRnan(aa,:),2)>0
            DRrow(aa,1) = 1;
        end
    end
    if exist('DRrow','var')
        DRrem = find(DRrow==1);
        DRnum(DRrem,:) = [];
        DR(DRrem,:) = [];
        DRmed(DRrem,:) = [];
        DRcv(DRrem,:) = [];
        PNR(DRrem) = [];
        MUPulses(DRrem) = [];
        IPTs(DRrem,:) = [];
        Cost(DRrem) = [];
    end

    for jj = 1:size(DR,1)
        % DRnan is a matrix where MUs to remove are indicated by 1s
        DRnan(jj,:) = isnan(DR(jj,:));
    end

    %% Find matched MUs
    [~,MUmaxX] = max(sum(DRnan));
    DRremCol = find(DRnan(:,MUmaxX)==1);
    MU2go = DRremCol;
    MUtokeep(:,1) = 1:length(MUPulses);
    MUmatch = setdiff(MUtokeep,MU2go);
    clear DRremCol MU2go

    % Set inputs for CST
    win_len = round(0.4*fsamp,0);
    h_win = hann(win_len);
    h_winNorm = h_win/sum(h_win);
    half_win = win_len/2; % half window to remove artifacts from filtering

    % Loop through each condition separately
    for kk = 1:3

        if size(crop,1) > 3
                trialNo = 1;
        end

        % AND TRIALS
        if kk == 1
            trialNo = sum(ISO);
        elseif kk == 2
            trialNo = sum(LEN);
        elseif kk == 3
            trialNo = sum(SHO);
        %else
            %error('>3 conditions in kk.')
        end

        DRremCol = find(DRnan(:,kk)==1);
        MU2go = DRremCol;
        MU2keep = setdiff(MUtokeep, MU2go);
        % Calculate cumulative spike train by summing binarized
        % MUPulses in each condition (i.e. non-matched MUs among conditions)

        % Looop through each trial separately
        clear first last
        for trial = 1:trialNo
            if trial == 1 % first trial
                first = crop(kk,1);
                last = crop(kk,1)+ssDur-1;
            elseif trial == trialNo % last trial
                first = crop(kk,2)-ssDur+1;
                last = crop(kk,2);
            else % all trials in between
                first = (trial-1)*ssDur+crop(kk,1);
                last = (trial-1)*ssDur+crop(kk,1)+ssDur-1;
            end

            if size(crop,1) > 3
                    last = crop(kk,2);
            end

            first = round(first);
            last = round(last);

            % % Calculate the cumulative spike train (cst)
            % cst(trial,:) = sum(spikeTrain(MU2keep,first:last));
            % % Filter (and later average) the CST
            % % 'same' keeps the output vector the exact same length as the input signal
            % Fcst_trial_raw = conv(cst(trial,:), h_winNorm, 'same');
            % Fcst_trial(trial,:) = Fcst_trial_raw(half_win+1:end-half_win);

            % % Divide cst by number of active motor units and
            % % divide again by the time window, so you have a
            % % mean firing rate per active motor unit (Hz) for that condition
            % % (in spikes per second per motor unit)
            % missingMU = find(~any(spikeTrain(MU2keep,first:last), 2));
            % numMUs = length(MU2keep) - length(missingMU);
            % window_duration = (last - first + 1) / fsamp;
            % cstNormMean_trial(trial,1) = (sum(Fcst_trial(trial,:)) / numMUs) / window_duration;
            % %(sum(cst(trial,:)) / numMUs) / window_duration;
            % rawMean = mean(Fcst_trial(trial,:));
            % rawMedian = median(Fcst_trial(trial,:));
            % cstNorm_trial(trial,1) = cstNormMean_trial(trial,1) * (rawMedian / rawMean);

        end % end of trial loop

        %[AUC_delta(1,kk),AUC_alpha(1,kk),AUC_beta(1,kk),AUC_gamma(1,kk)] = intramuscular_coherence(MUuse(:,kk,MUmatch),cst,fsamp,firstCrop(:,kk));

        % % Get variability of filtered Fcst and average among trials
        % FcstSD_trial = iqr(Fcst_trial,2);
        % FcstSD(1,kk) = mean(FcstSD_trial);
        % cstNorm(1,kk) = mean(cstNorm_trial);
        % cstNormMean(1,kk) = mean(cstNormMean_trial);

        DRMU1 = find(sum(DRnan,2)==2);
        if ~isempty(DRMU1)
            if length(DRMU1)>1
                MUspecial1(:,1) = 1:length(DRMU1);
                [~,MUspecial1(:,2)] = find(DRnan(DRMU1,:)==0);
            elseif ~isempty(DRMU1)
                MUspecial1(1,1) = 1;
                [~,MUspecial1(1,2)] = find(DRnan(DRMU1,:)==0);
            end
        end
        DRMU2 = find(sum(DRnan,2)==1);
        if ~isempty(DRMU2)
            if length(DRMU2)>1
                MUspecial2(:,1) = 1:length(DRMU2)*2;
                [~,MUspecial2(:,2)] = find(DRnan(DRMU2,:)==0);
            elseif ~isempty(DRMU2)
                MUspecial2(1:2,1) = [1;2];
                [~,MUspecial2(1:2,2)] = find(DRnan(DRMU2,:)==0);
            end
        end

        clear Fcst_trial FcstSD_trial cstNorm_trial cstNormMean_trial

    end

    A(bb,:) = [size(IPTs,1) round(mean(DRcv,'omitnan'),2) round(std(DRcv),2) round(min(DRcv),2) round(max(DRcv),2)];
    % DR stats (mean, median, std, min, max)
    AA(bb,:) = [round(mean(DR,'omitnan'),2) round(mean(DRmed,'omitnan'),2) round(std(DRmed),2) round(min(DRmed),2) round(max(DRmed),2)];
    % PNR stats (mean, std, min, max)
    AAA(bb,:) = [round(mean(PNR),2) round(std(PNR),2) round(min(PNR),2) round(max(PNR),2)];

    %% AUC in frequency bands
    %bands = table(AUC_delta,AUC_alpha,AUC_beta,AUC_gamma);

    %% CST
    % cstNorm;
    % FcstSD;

    %% MU removal
    if exist('DRrem','var')
        MUgone = length(DRrem);
    else
        MUgone = 0;
    end

    %% Unique motor units
    if exist('MUspecial1','var')
        specialMU1{bb,tt} = MUspecial1;
        clear MUspecial1 DRMU1
    else
        specialMU1{bb,tt} = NaN;
    end
    if exist('MUspecial2','var')
        specialMU2{bb,tt} = MUspecial2;
        clear MUspecial2 DRMU2
    else
        specialMU2{bb,tt} = NaN;
    end

    %[length(DR) std(DR) iqr(DRmed)]

    %save([fname{tt}(1:end-4) '#COHERENCE.mat'])

    %% Clearvars
    %clear DR* first last

end % end of trial loop

%clearvars -except fsamp fname pname folder specialMU1 specialMU2  tt

% Determine unique MUs
for aa = 1:length(specialMU1)
    if ~sum(sum(isnan(specialMU1{aa,tt}))) && ~sum(sum(isnan(specialMU2{aa,tt})))
        yes{aa,1} = vertcat(specialMU1{aa,tt}(:,2),specialMU2{aa,tt}(:,2));
    elseif ~sum(sum(isnan(specialMU2{aa,tt})))
        yes{aa,1} = specialMU2{aa,tt}(:,2);
    elseif ~sum(sum(isnan(specialMU1{aa,tt})))
        yes{aa,1} = specialMU1{aa,tt}(:,2);
    else
        yes{aa,1} = 0;
    end
end

for aa = 1:length(yes)
    for bb = 1:3
        cc = find(yes{aa,1}==bb);
        MUunique(aa,bb) = length(cc);
    end
end

%[A AA AAA];