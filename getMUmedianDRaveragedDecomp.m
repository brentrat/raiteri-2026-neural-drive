%% Get output of averaged motor units from each decomposition

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

    p1 = ['C:\Users\Brent\Desktop\Research\Masters\Karla\' subFolderNames{ff} '\Decomp\'];
    p2 = ['C:\Users\Brent\Desktop\Research\Masters\Karla\' subFolderNames{ff} '\' subFolderNames{ff}(1:2) '\'];

    % Load all 4 decomposition files per participant
    for kk = 1:2
        % DEMUSE + MUedit @ 20 then 40% MVC
        if kk == 1
            load([p2 'MVC20#CHECK.mat']);
            a = DRmed;
            clear DRmed
            load([p1 'MVC20_check.mat']);
            b = DRmed;
            c = vertcat(a,b);
            outcome(:,1) = round(mean(c),1);
            intensity = 1;
        elseif kk == 2
            load([p2 'MVC40#CHECK.mat']);
            a = DRmed;
            clear DRmed
            load([p1 'MVC40_check.mat']);
            b = DRmed;
            c = vertcat(a,b);
            outcome(:,1) = round(mean(c),1);
            intensity = 2;
        else
            error('k needs to be between 1 and 4.')
        end

    %% Average individual motor unit results (median DR = DRmed)
    % Store in a 5-column table later on
    % pid decomp intensity condition outcome

    % 1) outcome stored above
    
    % 2) pid
    rows = size(outcome,1);
    pidKey = ['p' subFolderNames{ff}(1:2)];
    pid = repmat(pidKey,[rows,1]);

    % 3) decomp not needed

    % 4) intensity
    intensity = repmat(intensity,[rows,1]);

    % 5) MUnumber not needed

    % 6) condition - 3 in total
    condition = [1; 2; 3];
    
    % Create a table to store results for each participant
    % Do not overwrite previous data
    if exist('results','var')
        finalRow = size(results,1);
        tRows = (1:rows)'+finalRow;
    else
        tRows = (1:rows)';
    end

    % Store in a 5(!)-column table now
    % pid decomp intensity condition outcome
    results(tRows, :) = table(pid, intensity, condition, outcome, ...
        'VariableNames', {'pid', 'intensity', 'condition', 'outcome'});

    clearvars -except topLevelFolder subFolderNames p1 p2 kk ff results 
end
end

% Save the results table to a .csv file for statistical (SAS/R) analysis
cd(topLevelFolder);
writetable(results,'motorUnitResultsAveragedDecomp_Raiteri.csv');

% Save the results table to a matfile
save(fullfile(topLevelFolder, 'motorUnitResultsAveragedDecomp_Raiteri.mat'), 'results');