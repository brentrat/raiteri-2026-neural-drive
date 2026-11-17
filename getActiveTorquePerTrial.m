%% Get steady-state active torque output of individual trials

% Clear workspace
clear; close all

%% Get paths of all folders
[~,topLevelFolder] = uigetfile('Select file in top folder of interest');
% topLevelFolder = pwd;
% Get a list of all files and folders in this folder
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory
dirFlags = [files.isdir];
% Extract only those that are directories
subFolders = files(dirFlags); % A structure with extra info
% Get only the folder names into a cell array.
subFolderNames = {subFolders(4:end-3).name};

%% Loop through each folder independently
for ff = 1:17

    pName = ['C:\Users\Brent\Desktop\Research\Masters\Karla\' subFolderNames{ff} '\'];

    % Load analysed files per participant
    load([pName 'analysed_new.mat'],'holdSS','T','TRef');

    %% Get individual trials (active torque (not normalised) in holdSS (columns 1:6)
    % outcome intensity condition pid trial

    % outcome
    % active torque
    %vals = round(holdSS(:,1:6),1); % do NOT normalise to %MVT

    % active torque steadiness (CV)
    % vals = nan(size(holdSS,1),6);
    % for cols = 1:6
    %     rows = size(TRef{cols},1);
    %     vals(1:rows,cols) = round(TRef{cols}(:,3),2);
    % end

    % active torque matching
    vals = nan(size(holdSS,1),6);
    for cols = 1:6
        rows = size(TRef{cols},1);
        vals(1:rows,cols) = round(TRef{cols}(:,1),2);
    end
    
    % outcome per trial
    outcome = vals(:);

    % intensity + condition
    trialNo = size(vals,1);
    condition = [ones(trialNo,1)*1; ones(trialNo,1)*2; ones(trialNo,1)*3];
    intensity = [ones(length(condition),1)*1; ones(length(condition),1)*2];
    condition = [condition;condition];

    % pid
    rows = size(outcome,1);
    pidKey = ['p' subFolderNames{ff}(1:2)];
    pid = repmat(pidKey,[rows,1]);

    % trial
    trialCond(:,1) = 1:trialNo;
    trial = repmat(trialCond,[6 1]);

    % remove NaNs
    bye = isnan(outcome);
    pid(bye,:) = [];
    intensity(bye) = [];
    condition(bye) = [];
    trial(bye) = [];
    outcome(bye) = [];
    tRows(:,1) = 1:length(outcome);

    % Create a table to store results for each participant
    % Do not overwrite previous data
    if exist('results','var')
        finalRow = size(results,1);
        tRows = tRows+finalRow;
    else
        finalRow = 0;
    end

    % Store in a 6-column table now
    % pid intensity condition trial outcome
    results(tRows, :) = table(pid(tRows-finalRow,:), intensity(tRows-finalRow), condition(tRows-finalRow), trial(tRows-finalRow), outcome(tRows-finalRow), ...
        'VariableNames', {'pid', 'intensity', 'condition', 'trial', 'outcome'});

    clearvars -except topLevelFolder subFolderNames pName ff results
end

% Save the results table to a .csv file for statistical (SAS/R) analysis
cd(topLevelFolder);
%writetable(results,'activeTorqueNmResultsMean_Raiteri.csv');
writetable(results,'activeTorqueCV_Raiteri.csv');

% Save the results table to a matfile
%save(fullfile(topLevelFolder, 'activeTorqueResults.mat'), 'results');