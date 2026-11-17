% Code with intramuscular coherence analysis
clear; close all

% Input
doPlot = 0;

% Sample rate
fsamp = 2000;
%fsamp = 2049;

% Files
fname = {'MVC20.mat' 'MVC40.mat'};

% Get folders
pname = 'C:\Users\Brent\Desktop\Research\Masters\Karla'; %pathname
cd(pname)
folders = dir;
folder = folders;
folder = {folders([folders.isdir]).name};
folder = folder(~ismember(folder ,{'.','..','00LJ','DFG','Old','Test','Hug','Revision'}));

if doPlot == 1

    figure(1); t = tiledlayout(6,6,'TileSpacing','tight','Padding','compact');

end

%% For doPlot
for bb = 1:length(folder)

    disp([num2str(bb), '/', num2str(length(folder))]);
    start_time = datetime('now', 'Format', 'HH:mm:ss');
    disp(char(start_time));

    cd(pname)
    cd(folder{bb})
    cd("Decomp")
    %cd('test')

    %% Check discharge rates of each MU for each condition in Karla's experiment
    for tt = 2 % 1 = MVC20; 2 = MVC40

        % Get trial length
        load(fname{tt});
        trialDur = (length(signal.torque)-fsamp)/length(signal.files);
        totalDur = length(signal.torque);

        % Crop trials based on condition
        ISO = contains(signal.files,'ISO');
        iso = sum(ISO);
        LEN = contains(signal.files,'LEN');
        len = iso+sum(LEN);
        SHO = contains(signal.files,'SHO');
        sho = len+sum(SHO);
        %% Steady state only
        crop = [1001 trialDur*iso+1000; trialDur*iso+1001 trialDur*len+1000; trialDur*len+1001 trialDur*sho+1000]; % crop steady state

        %% Load and analyse
        % Load file
        clear signal
        load([fname{tt}(1:end-4) '_dcmpUse.mat']) % saved in DEMUSE
        %load([fname{tt}(1:end-4) '_dcmp#.mat']) % to be loaded into DEMUSE

        % Counter
        ii = 1;

        %% Get spike train by binarizing MUPulses
        spikeTrain = binario(MUPulses,totalDur);

        % Loop through conditions
        for kk = 1:3

            % AND TRIALS
            if kk == 1
                trialNo = sum(ISO);
            elseif kk == 2
                trialNo = sum(LEN);
            elseif kk == 3
                trialNo = sum(SHO);
            else
                error('>3 conditions in kk.')
            end

            %% Loop through each MU
            for jj = 1:length(MUPulses)

                % Find specifc MU in given condition
                MU = MUPulses{jj};

                for trial = 1:trialNo
                    if trial == 1 % first trial
                        first = find(MU>=crop(kk,1),1,'first');
                        last = find(MU<=crop(kk,1)+trialDur-1,1,'last');
                        if jj == 1
                            firstCrop(trial,kk,jj) = crop(kk,1);
                        end
                    elseif trial == trialNo % last trial
                        first = find(MU>=crop(kk,2)-trialDur+1,1,'first');
                        last = find(MU<=crop(kk,2),1,'last');
                        if jj == 1
                            firstCrop(trial,kk,jj) = crop(kk,2)-trialDur+1;
                        end
                    else % all trials in between
                        first = find(MU>=(trial-1)*trialDur+crop(kk,1),1,'first');
                        last = find(MU<=(trial-1)*trialDur+crop(kk,1)+trialDur-1,1,'last');
                        if jj == 1
                            firstCrop(trial,kk,jj) = (trial-1)*trialDur+crop(kk,1);
                        end
                    end
                    %trial = 1; % for concatenating trials, uncomment this line and comment above lines

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

        [DRrow,DRcol] = find(DRcv>=40); %40

        if ~isempty(DRrow)
            % for aa = 1:3
            %     aaa = find(DRcol==aa);
            %     DR(DRrow(aaa),aa) = NaN;
            %     DRmed(DRrow(aaa),aa) = NaN;
            % end
            DR(DRrow,DRcol) = NaN; DRmed(DRrow,DRcol) = NaN;
            for jj = 1:length(DR)
                % DRnan is a matrix where MUs to remove are indicated by 1s
                DRnan(jj,1:3) = isnan(DR(jj,:));
            end
            DRrem = find(sum(DRnan,2)>0);
            % DRpnr is a matrix where MUs to remove are indicated by 1s
            DRrempnr(:,1) = find(round(PNR,0)<28);
            % DRcommon is a vector with shared MUs to remove based on CV
            % and PNR exclusion criteria
            DRremCommon = intersect(DRrem, DRrempnr);
            CV_rem = DRrem(~ismember(DRrem, DRrempnr));     % elements of A not present in B
            PNR_rem = DRrempnr(~ismember(DRrempnr, DRrem)); % elements of B not present in A
            DRremCommonP = length(DRremCommon)/(length(CV_rem)+length(PNR_rem))*100;
            % Remove noisy (CV + PNR based) MUs from CST calculation for
            % each condition:
            %MUtogo = [CV_rem; PNR_rem; DRremCommon];
            %MUtogo = sort(MUtogo);
            %% Find matched MUs
            [~,MUmaxX] = max(sum(DRnan));
            DRremCol = find(DRnan(:,MUmaxX)==1);
            MUtogo = [DRremCol; DRrempnr];
            MU2go = unique(MUtogo);
            MUtokeep(:,1) = 1:length(MUPulses);
            MUmatch = setdiff(MUtokeep,MU2go);
            clear DRremCol MUtogo MU2go

            % Set inputs for CST
            win_len = 0.4*fsamp;
            h_win = hann(win_len);
            h_winNorm = h_win/sum(h_win);
            half_win = win_len/2; % half window to remove artifacts from filtering

            % Loop through each condition separately
            for kk = 1:3

                % AND TRIALS
                if kk == 1
                    trialNo = sum(ISO);
                elseif kk == 2
                    trialNo = sum(LEN);
                elseif kk == 3
                    trialNo = sum(SHO);
                else
                    error('>3 conditions in kk.')
                end

                DRremCol = find(DRnan(:,kk)==1);
                MUtogo = [DRremCol; DRrempnr];
                MU2go = unique(MUtogo);
                MU2keep = setdiff(MUtokeep, MU2go);
                % Calculate cumulative spike train by summing binarized
                % MUPulses in each condition (i.e. non-matched MUs among conditions)

                % Looop through each trial separately
                clear first last
                for trial = 1:trialNo
                    if trial == 1 % first trial
                        first = crop(kk,1);
                        last = crop(kk,1)+trialDur-1;
                    elseif trial == trialNo % last trial
                        first = crop(kk,2)-trialDur+1;
                        last = crop(kk,2);
                    else % all trials in between
                        first = (trial-1)*trialDur+crop(kk,1);
                        last = (trial-1)*trialDur+crop(kk,1)+trialDur-1;
                    end

                    % Calculate the cumulative spike train (cst)
                    cst(trial,:) = sum(spikeTrain(MU2keep,first:last));
                    % Filter (and later average) the CST
                    % 'same' keeps the output vector the exact same length as the input signal
                    Fcst_trial_raw = conv(cst(trial,:), h_winNorm, 'same');
                    Fcst_trial(trial,:) = Fcst_trial_raw(half_win+1:end-half_win);

                    % Divide cst by number of active motor units and
                    % divide again by the time window, so you have a
                    % mean firing rate per active motor unit (Hz) for that condition
                    % (in spikes per second per motor unit)
                    missingMU = find(~any(spikeTrain(MU2keep,first:last), 2));
                    numMUs = length(MU2keep) - length(missingMU);
                    window_duration = (last - first + 1) / fsamp;
                    cstNormMean_trial(trial,1) = (sum(Fcst_trial(trial,:)) / numMUs) / window_duration;
                    %(sum(cst(trial,:)) / numMUs) / window_duration;
                    rawMean = mean(Fcst_trial(trial,:));
                    rawMedian = median(Fcst_trial(trial,:));
                    cstNorm_trial(trial,1) = cstNormMean_trial(trial,1) * (rawMedian / rawMean);

                end % end of trial loop

                [AUC_delta(1,kk),AUC_alpha(1,kk),AUC_beta(1,kk),AUC_gamma(1,kk)] = intramuscular_coherence(MUuse(:,kk,MUmatch),cst,fsamp,firstCrop(:,kk));

                % Get std of filtered Fcst and average among trials
                FcstSD_trial = iqr(Fcst_trial,2);
                FcstSD(1,kk) = mean(FcstSD_trial);
                cstNorm(1,kk) = mean(cstNorm_trial);
                cstNormMean(1,kk) = mean(cstNormMean_trial);

                clear Fcst_trial FcstSD_trial cstNorm_trial cstNormMean_trial

            end

            DRMU1 = find(sum(DRnan,2)==2);
            if ~isempty(DRMU1)
                bye = find(PNR(DRMU1)<28);
                DRMU1(bye) = [];
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
                bye = find(PNR(DRMU2)<28);
                DRMU2(bye) = [];
                if length(DRMU2)>1
                    MUspecial2(:,1) = 1:length(DRMU2)*2;
                    [~,MUspecial2(:,2)] = find(DRnan(DRMU2,:)==0);
                elseif ~isempty(DRMU2)
                    MUspecial2(1:2,1) = [1;2];
                    [~,MUspecial2(1:2,2)] = find(DRnan(DRMU2,:)==0);
                end
            end
            DRnum(DRrem,:) = [];
            DR(DRrem,:) = [];
            DRmed(DRrem,:) = [];
            DRcv(DRrem,:) = [];
            %DRbi(DRrem,:) = [];
            PNR(DRrem) = [];
            MUPulses(DRrem) = [];
            IPTs(DRrem,:) = [];
            Cost(DRrem) = [];
        else
            for jj = 1:length(DR)
                % DRnan is a matrix where MUs to remove are indicated by 1s
                DRnan(jj,1:3) = isnan(DR(jj,:));
            end
            DRrem = find(sum(DRnan,2)>0);
            % DRpnr is a matrix where MUs to remove are indicated by 1s
            DRrempnr(:,1) = find(round(PNR,0)<28);

            % DRcommon is a vector with shared MUs to remove based on CV
            % and PNR exclusion criteria
            DRremCommon = intersect(DRrem, DRrempnr);
            CV_rem = DRrem(~ismember(DRrem, DRrempnr));     % elements of A not present in B
            PNR_rem = DRrempnr(~ismember(DRrempnr, DRrem)); % elements of B not present in A
            DRremCommonP = length(DRremCommon)/(length(CV_rem)+length(PNR_rem))*100;
            % Remove noisy (CV + PNR based) MUs from CST calculation for
            % each condition:
            %MUtogo = [CV_rem; PNR_rem; DRremCommon];
            %MUtogo = sort(MUtogo);
            %% Find matched MUs
            [~,MUmaxX] = max(sum(DRnan));
            DRremCol = find(DRnan(:,MUmaxX)==1);
            MUtogo = [DRremCol; DRrempnr];
            MU2go = unique(MUtogo);
            MUtokeep(:,1) = 1:length(MUPulses);
            MUmatch = setdiff(MUtokeep,MU2go);
            clear DRremCol MUtogo MU2go

            % Set inputs for CST
            win_len = 0.4*fsamp;
            h_win = hann(win_len);
            h_winNorm = h_win/sum(h_win);
            half_win = win_len/2; % half window to remove artifacts from filtering

            % Loop through each condition separately
            for kk = 1:3

                % AND TRIALS
                if kk == 1
                    trialNo = sum(ISO);
                elseif kk == 2
                    trialNo = sum(LEN);
                elseif kk == 3
                    trialNo = sum(SHO);
                else
                    error('>3 conditions in kk.')
                end

                DRremCol = find(DRnan(:,kk)==1);
                MUtogo = [DRremCol; DRrempnr];
                MU2go = unique(MUtogo);
                MU2keep = setdiff(MUtokeep, MU2go);
                % Calculate cumulative spike train by summing binarized
                % MUPulses in each condition (i.e. non-matched MUs among conditions)

                % Looop through each trial separately
                clear first last
                for trial = 1:trialNo
                    if trial == 1 % first trial
                        first = crop(kk,1);
                        last = crop(kk,1)+trialDur-1;
                    elseif trial == trialNo % last trial
                        first = crop(kk,2)-trialDur+1;
                        last = crop(kk,2);
                    else % all trials in between
                        first = (trial-1)*trialDur+crop(kk,1);
                        last = (trial-1)*trialDur+crop(kk,1)+trialDur-1;
                    end

                    % Calculate the cumulative spike train (cst)
                    cst(trial,:) = sum(spikeTrain(MU2keep,first:last));
                    % Filter (and later average) the CST
                    % 'same' keeps the output vector the exact same length as the input signal
                    Fcst_trial_raw = conv(cst(trial,:), h_winNorm, 'same');
                    Fcst_trial(trial,:) = Fcst_trial_raw(half_win+1:end-half_win);

                    % Divide cst by number of active motor units and
                    % divide again by the time window, so you have a
                    % mean firing rate per active motor unit (Hz) for that condition
                    % (in spikes per second per motor unit)
                    missingMU = find(~any(spikeTrain(MU2keep,first:last), 2));
                    numMUs = length(MU2keep) - length(missingMU);
                    window_duration = (last - first + 1) / fsamp;
                    cstNormMean_trial(trial,1) = (sum(Fcst_trial(trial,:)) / numMUs) / window_duration;
                    %(sum(cst(trial,:)) / numMUs) / window_duration;
                    rawMean = mean(Fcst_trial(trial,:));
                    rawMedian = median(Fcst_trial(trial,:));
                    cstNorm_trial(trial,1) = cstNormMean_trial(trial,1) * (rawMedian / rawMean);

                end % end of trial loop

                [AUC_delta(1,kk),AUC_alpha(1,kk),AUC_beta(1,kk),AUC_gamma(1,kk)] = intramuscular_coherence(MUuse(:,kk,MUmatch),cst,fsamp,firstCrop(:,kk));

                % Get std of filtered Fcst and average among trials
                FcstSD_trial = iqr(Fcst_trial,2);
                FcstSD(1,kk) = mean(FcstSD_trial);
                cstNorm(1,kk) = mean(cstNorm_trial);
                cstNormMean(1,kk) = mean(cstNormMean_trial);

                clear Fcst_trial FcstSD_trial cstNorm_trial cstNormMean_trial

            end

            DRMU1 = find(sum(DRnan,2)==2);
            if ~isempty(DRMU1)
                bye = find(PNR(DRMU1)<28);
                DRMU1(bye) = [];
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
                bye = find(PNR(DRMU2)<28);
                DRMU2(bye) = [];
                if length(DRMU2)>1
                    MUspecial2(:,1) = 1:length(DRMU2)*2;
                    [~,MUspecial2(:,2)] = find(DRnan(DRMU2,:)==0);
                elseif ~isempty(DRMU2)
                    MUspecial2(1:2,1) = [1;2];
                    [~,MUspecial2(1:2,2)] = find(DRnan(DRMU2,:)==0);
                end
            end
            DRnum(DRrem,:) = [];
            DR(DRrem,:) = [];
            DRmed(DRrem,:) = [];
            DRcv(DRrem,:) = [];
            %DRbi(DRrem,:) = [];
            PNR(DRrem) = [];
            MUPulses(DRrem) = [];
            IPTs(DRrem,:) = [];
            Cost(DRrem) = [];
        end

        PNRrow = find(round(PNR,0)<28);
        DRnum(PNRrow,:) = [];
        DR(PNRrow,:) = [];
        DRmed(PNRrow,:) = [];
        DRcv(PNRrow,:) = [];
        PNR(PNRrow) = [];
        MUPulses(PNRrow) = [];
        IPTs(PNRrow,:) = [];
        Cost(PNRrow) = [];

        %% Summary stats
        % Number of firings and CV stats
        A(bb,:) = [size(IPTs,1) round(mean(DRcv,'omitnan'),2) round(std(DRcv),2) round(min(DRcv),2) round(max(DRcv),2)];
        % DR stats (mean, median, std, min, max)
        AA(bb,:) = [round(mean(DR,'omitnan'),2) round(mean(DRmed,'omitnan'),2) round(std(DRmed),2) round(min(DRmed),2) round(max(DRmed),2)];
        % PNR stats (mean, std, min, max)
        AAA(bb,:) = [round(mean(PNR),2) round(std(PNR),2) round(min(PNR),2) round(max(PNR),2)];

        %% AUC in frequency bands
        bands = table(AUC_delta,AUC_alpha,AUC_beta,AUC_gamma);

        %% CST
        % cstNorm;
        % FcstSD;

        %% MU removal
        MUgone = table(length(DRrem),length(DRrempnr),length(DRremCommon),DRremCommonP);

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

        %% Histogram
        if doPlot == 1
            DRdiff(:,1) = round(DR(:,2)-DR(:,1),1);
            DRdiff(:,2) = round(DR(:,3)-DR(:,1),1);

            %edges = -6.5:1:6.5;
            if min(min(DRdiff)) < -6 ||  max(max(DRdiff)) > 6
                nexttile; histogram(DRdiff(:,1),'BinWidth',1); hold on; histogram(DRdiff(:,2),'BinWidth',1); title(['P' num2str(bb),': ' num2str(length(DR)) ' MUs'],'FontWeight','normal','fontsize',12)
                box off; set(gca,'TickDir','out');
            else
                nexttile; histogram(DRdiff(:,1),'BinWidth',1); hold on; histogram(DRdiff(:,2),'BinWidth',1); title(['P' num2str(bb),': ' num2str(length(DR)) ' MUs'],'FontWeight','normal','fontsize',12)
                xlim([-6 6]); xticks([-6 -4 -2 0 2 4 6]);
                box off; set(gca,'TickDir','out');
            end

            if tt == 1
                lim(bb,1:2) = [min(min(DRdiff)) max(max(DRdiff))];
            else
                lim(bb,3:4) = [min(min(DRdiff)) max(max(DRdiff))];
            end

        end
        save([fname{tt}(1:end-4) '_coherence.mat'])

        %% Clearvars
        clear DR* first last

    end % end of trial loop

    clearvars -except doPlot fsamp fname pname folder specialMU1 specialMU2  tt
end % end of participant loop

end_time = datetime('now', 'Format', 'HH:mm:ss');
disp(char(end_time));

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

if doPlot == 1
    % Fix axes
    nexttile(1); yticks([0 3 6]); title({'\bf20% MVC\rm','P1: 9 MUs'},'fontsize',12)
    nexttile(2); yticks([0 3 6]);  title({'\bf40% MVC\rm','P1: 9 MUs'},'fontsize',12)
    nexttile(3); title({'\bf20% MVC\rm','P2: 34 MUs'},'fontsize',12)
    nexttile(4); title({'\bf40% MVC\rm','P2: 33 MUs'},'fontsize',12)
    nexttile(5); title({'\bf20% MVC\rm','P3: 35 MUs'},'fontsize',12)
    nexttile(6); title({'\bf40% MVC\rm','P3: 34 MUs'},'fontsize',12)
    nexttile(8); xlim([-6 8]); xticks([-6 -4 -2 0 2 4 6 8])
    nexttile(14); xlim([-14 6]); xticks([-14 -6 -4 -2 0 2 4 6])
    nexttile(17); xlim([-6 6]); xticks([-6 -4 -2 0 2 4 6])
    nexttile(21); xlim([-6 13]); xticks([-6 -4 -2 0 2 4 6 12])

    % Add legend
    nexttile(6); legend({'LEN'; 'SHO'}); legend('boxoff'); legend('Location','northeastoutside');

    % Add labels
    xlabel(t,'\Delta Discharge rate (Hz)','fontweight','bold','fontsize',14)
    ylabel(t,'Count (#)','fontweight','bold','fontsize',14)

    fig = gcf;  %or one particular figure whose handle you already know, or 0 to affect all figures
    set( findall(fig, '-property', 'fontsize'), 'fontsize', 10)

    cd(pname);

    DR20 = [10.3	10.3	11.6
        10.7	9	10.5
        11.2	10.3	10.7
        13.5	11.5	14.5
        14.9	14.4	13.3
        12.2	12	13.9
        12.8	10.7	12.7
        11	8.2	11.1
        14.2	12	13.2
        11.3	9.8	11
        12.1	12.1	12.7
        13.1	10.9	13
        13.1	12.4	12.7
        12.2	12.1	11.8
        15.4	12.3	15.8
        12.6	11.3	12.3
        10.8	10.3	10.6];

    DR40 = [14.3	12.5	14.4
        14.2	11.4	13.8
        12	13	13.5
        16	16.1	20.4
        17.9	16.9	18.2
        14.9	12.5	16.1
        16.6	13.8	17.3
        12.6	11.2	14.9
        16.1	15	16.5
        13.1	12.6	13.9
        15.3	14.9	15.7
        14.9	14.4	16.2
        14.1	14.6	15.7
        15.4	15	15.6
        17.8	17.4	18.8
        14	13.7	14.5
        12.7	12.2	13.6];

    DRdiff(:,1) = round(DR20(:,2)-DR20(:,1),1);
    DRdiff(:,2) = round(DR20(:,3)-DR20(:,1),1);

    nexttile(35);
    histogram(DRdiff(:,1),'BinWidth',1); hold on; histogram(DRdiff(:,2),'BinWidth',1); title(['N=17 (mean): 25 MUs'],'FontWeight','bold','fontsize',10)
    xlim([-6 6]); xticks([-6 -4 -2 0 2 4 6]);
    box off; set(gca,'TickDir','out');

    DRdiff(:,1) = round(DR40(:,2)-DR40(:,1),1);
    DRdiff(:,2) = round(DR40(:,3)-DR40(:,1),1);

    nexttile(36);
    histogram(DRdiff(:,1),'BinWidth',1); hold on; histogram(DRdiff(:,2),'BinWidth',1); title(['N=17 (mean): 21 MUs'],'FontWeight','bold','fontsize',10)
    xlim([-6 6]); xticks([-6 -4 -2 0 2 4 6]);
    box off; set(gca,'TickDir','out');

    clear

end
