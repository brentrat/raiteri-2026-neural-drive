clear; %close all

% Load summary data
pathName = 'C:\Users\Brent\Desktop\Research\Masters\Karla';
cd(pathName)
folders = dir;
is_dir_flag = [folders.isdir];
folders = folders(is_dir_flag);
folder = folders(4:20);

for jj =  1:17

    cd(pathName)

    cd(folder(jj).name)

    load('analysed_new.mat','-regexp', '^(?!EMG)\w'); 

    %% Make a nice figure
    % 1st EMG, 2nd angle, 3rd torque
    for aa = 1:6
        C1{aa}(jj,:) = mean(emgCond{aa}/mean(E)*100,1,'omitnan');
        C2{aa}(jj,:) = mean(angCond{aa},1,'omitnan');
        C3{aa}(jj,:) = mean(tqCond{aa}/T*100,1,'omitnan');
        if jj == 7
            c1{aa} = emgCond{aa}/mean(E)*100;
            c2{aa} = angCond{aa};
            c3{aa} = tqCond{aa}/T*100;
        end
        TQerror(jj,aa) = round(mean(TRef{aa}(:,1)),1);
        TQstd(jj,aa) = round(mean(TRef{aa}(:,2)),1);
        TQcv(jj,aa) = round(mean(TRef{aa}(:,3)),2);
        TQhold(jj,aa) = AVGholdTQ(1,aa);
        TQholdAbs(jj,aa) = AVGholdTQ(1,aa)/100*T;
        EMGhold(jj,aa) = round(mean(AVGholdEMG(:,aa)),2);
        EMGholdbi(jj,aa) = round(AVGholdEMGbi(1,aa),2);
        EMGcv(jj,aa) = round(mean(AVGholdEMGcv(:,aa)),2);
        EMGcvbi(jj,aa) = round(AVGholdEMGcvbi(1,aa),2);
        EMGmf(jj,aa) = round(mean(AVGholdEMGmf(:,aa)),2);
        EMGmfbi(jj,aa) = round(AVGholdEMGmfbi(1,aa),2);
        EMGmapCol(jj,aa) = AVGholdEMGmapCol(1,aa);
        EMGmapRow(jj,aa) = AVGholdEMGmapRow(1,aa);

        % Rotation
        TQrot(jj,aa) = AVGrotTQ(1,aa);
        TQramp(jj,aa) = AVGrampTQ(1,aa);
        EMGrot(jj,aa) = round(mean(AVGrotEMG(:,aa)),2);
        EMGramp(jj,aa) = round(mean(AVGrampEMG(:,aa)),2);
    end

    trials(jj,:) = c;

    if jj == 1
        MVC = mvc;
        MVC(:,7) = ones(size(mvc,1),1)*jj;
    else
        test = mvc;
        test(:,7) = ones(size(mvc,1),1)*jj;
        MVC = vertcat(MVC,test);
    end

    clearvars -except trials C* c* TQ* EMG* pathName folder jj MVC

end

cd(pathName)

%% Prism figure
xx = 7;
for aa = 1:6
    %REF{aa}(:,2) = mean(C2{aa},1,'omitnan'); REF{aa}(:,1) = mean(C1{aa},1,'omitnan'); REF{aa}(:,3) = mean(C3{aa},1,'omitnan');
    % C stores the between-trial mean: C1 = EMG; C2 = ANG; C3 = TQ
    REF{aa}(:,2) = C2{aa}(xx,:); REF{aa}(:,1) = C1{aa}(xx,:); REF{aa}(:,3) = C3{aa}(xx,:);
end
FIGind = table(REF{1},REF{2},REF{3},REF{4},REF{5},REF{6});
TIME(:,1) = 1/2000:1/2000:1/2000*39000;
TIME = TIME-1.5;

% [f,location] = uiputfile('analysed_new_summary.mat');
% if isequal(f,0) || isequal(location,0)
%     disp('File not saved.')
% else
%     save('analysed_new_summary.mat')
% end