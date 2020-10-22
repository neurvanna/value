for sessionID = [1:7,9,10,11,13,15,17]
%for sessionID = 3

    clearvars -except sessionID;
    close all;
    
    db(1).subject = 'AL021';
    db(1).date = '2019-06-05'; %k2, k3 good 1 2
    db(2).subject = 'AL021';
    db(2).date = '2019-06-06'; %k1,k2,k3, zo good 3 4 5 6
    db(3).subject = 'AL021';
    db(3).date = '2019-06-07';%k1,k2 good
    db(4).subject = 'AL022';
    db(4).date = '2019-06-19';%k1,k2,k3,zo good
    db(5).subject = 'MW003';
    db(5).date = '2019-08-11';%k1,k2,zo good
    db(6).subject = 'MW003';
    db(6).date = '2019-08-12';%k1,k2 good
    db(7).subject = 'AL029';
    db(7).date = '2019-10-23';%k3,zo ok
    db(9).subject = 'AL026';
    db(9).date = '2019-11-01';%k1, k2
    db(10).subject = 'AL026';
    db(10).date = '2019-11-02';%k1, k2 good  
    db(11).subject = 'AL026';
    db(11).date = '2019-11-03';%k1, k2
    db(13).subject = 'AL026';
    db(13).date = '2019-11-05';%k1, k2
    db(15).subject = 'AL026';
    db(15).date = '2019-11-07';%k1, k2 ok   
    db(17).subject = 'AL026';
    db(17).date = '2019-11-13';%k ok 
    db(18).subject = 'AL016';
    db(18).date = '2019-07-11';%k1,k2
    %     db(8).subject = 'AL028';
    %     db(8).date = '2019-10-29';%k2


    %     db(12).subject = 'AL028';
    %     db(12).date = '2019-11-04';%k2

    %db(14).subject = 'AL026';
    %db(14).date = '2019-11-06';%k1, k2

%     db(16).subject = 'AL026';
%     db(16).date = '2019-11-12';%k2


    
    subject = db(sessionID).subject;
    date = db(sessionID).date;
    %% Define folders
    subjectsFolder = '\\znas.cortexlab.net\Subjects';
    % subjectsFolder = '\\basket.cortexlab.net\data\nick';
    %subjectsFolder = '\\znas.cortexlab.net\Subjects';
    alignDir = fullfile(subjectsFolder, subject, date, 'alignments');
    % driftPlotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\OpticTract\driftPlots';
    driftPlotFolder = 'C:\Users\annaL\Documents\PhD\ephys_results\driftPlots';
    probeDataPlotFolder = 'C:\Users\annaL\Documents\PhD\ephys_results\probeDataPlots';
    
    %if ~isfolder(alignDir)
    mkdir(alignDir)
    %end
    
    %% Get basic info
    
    % if ~exist(alignDir, 'dir')
    %     mkdir(alignDir);
    % end
    [tags, hasEphys] = getEphysTags(subject, date);
    % determine what exp nums exist
    [expNums, blocks, hasBlock, pars, isMpep, tl, hasTimeline] = ...
        dat.whichExpNums(subject, date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    Timeline = tl{1};
    
    useFlipper = true;
    % useFlipper = false;
    bBlocktoTL = readNPY([alignDir '\correct_block_2_to_timeline_1.npy']);
    %bProbe1toProbe2 =readNPY('Y:\AL022\2019-06-19\alignments\correct_ephys_K3_to_ephys_K1.npy');
    
    subjectsFolder = getRootDir(subject, date);
    TLdir = fullfile(subjectsFolder,'1')
    load(fullfile(subjectsFolder,'2', [date '_2_' subject '_Block.mat']))
    
    Fs = 30000;
    for t = 1:length(tags)
        
        probeTag = tags(t);
        probe_number = t;
        dataDir = char(fullfile(subjectsFolder, 'data_new_format', tags(t)));
        bTLtoMaster = readNPY(fullfile(alignDir, ...
            sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{t})));
        
        
        %% spikes
        sp = loadAllKsDir(subject, date);
        spikes.clusters = sp(t).spikeTemplates;
        spikes.times = sp(t).st;
       % s = spikes.times;
       % writeNPY(s, fullfile(dataDir, 'spikes.times.npy'))
       st = readNPY(fullfile(dataDir, 'spikes.times.npy'));
       sc = readNPY(fullfile(dataDir, 'spikes.clusters.npy'));
       ft = readNPY(fullfile(dataDir, 'trials.feedback_times.npy'));
       mt = readNPY(fullfile(dataDir, 'moves.moveOnsets.npy'));
       for i = 1:length(ft)
           mt_before_ft(i) = max(mt(find(mt<ft(i))));
       end
       c_locationAll = readtable(fullfile(dataDir, 'clusters.brainLocation.csv'));
        c_location = c_locationAll(:,4);
        Q1 = readNPY(fullfile(dataDir, 'Q1.npy'));
        Q2 = readNPY(fullfile(dataDir, 'Q2.npy'));
        h = readNPY(fullfile(dataDir, 'h.npy'));
        allClu = unique(sc);
        signif = zeros(1,length(allClu));
        st_bins = zeros(length(ft), length(allClu));
        st_bins_rotated = zeros(length(ft), length(allClu), length(ft));
        counter = 0;
        %% using correlations
        
%         realCorQ1 = zeros(1, length(allClu));
%         realCorQ2 = zeros(1, length(allClu));
%         realCorH = zeros(1, length(allClu));
%         signifQ1 = zeros(1, length(allClu));
%         signifQ2 = zeros(1, length(allClu));
%         signifH = zeros(1, length(allClu));
        rotatedCor = zeros(length(allClu), length(ft));

        for j = 1:length(allClu)
            currClu = allClu(j);
            currSt = st(find(sc==currClu));
            for i = 1:length(ft)
                st_bins(i, j) = numel(find(currSt<mt_before_ft(i) & currSt>mt_before_ft(i)-0.1)); % trial*cluster
            end
        end

            for k = 1:length(ft)
                st_bins_rotated(:,:,k) = circshift(st_bins, k); % along second axis
            end
    %         for j = 1:length(allClu)
    %             realCorQ1(j) = corr(st_bins(:,j), Q1);
    %             for k = 1:length(ft)
    %                 rotatedCorQ1(j, k) = corr(st_bins_rotated(:,j, k), Q1);
    %             end
    %             realCorQ2(j) = corr(st_bins(:,j), Q2);
    %             for k = 1:length(ft)
    %                 rotatedCorQ2(j, k) = corr(st_bins_rotated(:,j, k), Q2);
    %             end
    %             realCorH(j) = corr(st_bins(:,j), h);
    %             for k = 1:length(ft)
    %                 rotatedCorH(j, k) = corr(st_bins_rotated(:,j, k), h);
    %             end
    %             if (realCorQ1(j)>prctile(rotatedCorQ1(j,:),97.5) | realCorQ1(j)<prctile(rotatedCorQ1(j,:),2.5)) 
    %                signifQ1(j) = 1;
    %             end
    %             if ( realCorQ2(j)>prctile(rotatedCorQ2(j,:),97.5) | realCorQ2(j)<prctile(rotatedCorQ2(j,:),2.5))
    %                 signifQ2(j) = 1;
    %             end
    %             if (realCorH(j)>prctile(rotatedCorH(j,:),97.5) | realCorH(j)<prctile(rotatedCorH(j,:),2.5))
    %                 signifH(j) = 1;
    %             end  
    %         end
    %         sigQ1 = c_location(find(signifQ1==1),1);
    %         sigQ2 = c_location(find(signifQ2==1),1);
    %         sigH = c_location(find(signifH==1),1);
    %         allLoc = c_location;
    %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding', sprintf('%s_%s_%s_signifQ1.mat',subject, date, probeTag{1})), 'sigQ1')
    %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding', sprintf('%s_%s_%s_signifQ2.mat',subject, date, probeTag{1})), 'sigQ2')
    %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding', sprintf('%s_%s_%s_signifH.mat',subject, date, probeTag{1})), 'sigH')
    %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\allLocs', sprintf('%s_%s_%s_allLoc.mat',subject, date, probeTag{1})), 'allLoc')
    %         

            %% using linear regression

            fb = double([block.events.feedbackValues]);
            fb(find(fb==0)) = -1;
            fb = fb';
            rsp = double(block.events.responseValues);
            rsp = rsp';
            pastfb = double(block.events.feedbackValues);
            pastfb(find(pastfb==0)) = -1;
            pastfb = [1 pastfb(1:end-1)];
            pastfb = pastfb';
            pastrsp = double(block.events.responseValues);
            pastrsp = [1 pastrsp(1:end-1)];
            pastrsp = pastrsp';
            interaction = pastrsp.*pastfb;
            signifRegQ1 = [];
            signifRegQ2 = [];
            signifRegH = [];
            signifRegpastFb = [];
            signifRegRsp = [];
            CPD_Rsp = [];
            CPD_pastFb = [];
            CPD_Q1 = [];
            CPD_Q2 = [];
            CPD_H = [];
            CPD_Rsp_rotated = [];
            CPD_pastFb_rotated = [];
            CPD_Q1_rotated = [];
            CPD_Q2_rotated = [];
            CPD_H_rotated = [];
            signifPct = zeros(1,100);
            use_cluster = [];
            for j = 1:length(allClu)
                if sum(st_bins(:, j))>50
                    use_cluster(end+1) = j;
                end
            end
            
            for j = 1:length(use_cluster) %for each cell, compute how much adding Q1 helps in reality and in rotated datasets. If a lot, that's a significant cell

                    j/length(use_cluster)
                    len2 = floor(size(st_bins,1)/2);
                    mdl0_withoutRsp = fitglm([pastrsp(1:len2), rsp(1:len2), pastfb(1:len2), pastrsp(1:len2).*pastfb(1:len2)], st_bins(1:len2, use_cluster(j)), 'linear');
                    mdl0_Rsp = fitglm([pastrsp(1:len2), rsp(1:len2), pastfb(1:len2), pastrsp(1:len2).*pastfb(1:len2), Q1(1:len2)], st_bins(1:len2, use_cluster(j)), 'linear');
                    CPD_Rsp(j) = (mdl0_withoutRsp.Deviance - mdl0_Rsp.Deviance);

        %             mdl0_withoutQ1Q2H = fitglm([fb, rsp, pastrsp, interaction], st_bins(:, j));
        %             
        %             mdl0_Q1 = fitglm([fb, rsp, pastrsp, interaction, Q1], st_bins(:, j));
        %             CPD_Q1(j) = (mdl0_withoutQ1Q2H.Deviance - mdl0_Q1.Deviance);
        %             
        %             mdl0_Q2 = fitglm([fb, rsp, pastrsp, interaction, Q2], st_bins(:, j));
        %             CPD_Q2(j) = (mdl0_withoutQ1Q2H.Deviance - mdl0_Q2.Deviance);
        %             
        %             mdl0_H = fitglm([fb, rsp, pastrsp, interaction, h], st_bins(:, j));
        %             CPD_H(j) = (mdl0_withoutQ1Q2H.Deviance - mdl0_H.Deviance);

                    for k = 2:length(ft)/2
                        shifted = st_bins_rotated(1:len2,:,k);

                        mdl1_withoutRsp = fitglm([pastrsp(1:len2), rsp(1:len2), pastfb(1:len2), pastrsp(1:len2).*pastfb(1:len2)], shifted(:, use_cluster(j)), 'linear');

                        mdl1_Rsp = fitglm([pastrsp(1:len2), rsp(1:len2), pastfb(1:len2), pastrsp(1:len2).*pastfb(1:len2), Q1(1:len2)], shifted(:, use_cluster(j)), 'linear');
                        CPD_Rsp_rotated(j,k) = (mdl1_withoutRsp.Deviance - mdl1_Rsp.Deviance); % neurons*trial sample

        %                 mdl1_withoutQ1Q2H = fitglm([fb, rsp, pastrsp, interaction], shifted(:, j));                
        %                 
        %                 mdl1_Q1 = fitglm([fb, rsp, pastrsp, interaction, Q1], shifted(:, j));
        %                 CPD_Q1_rotated(j,k) = (mdl1_withoutQ1Q2H.Deviance - mdl1_Q1.Deviance);  
        %                 
        %                 mdl1_Q2 = fitglm([fb, rsp, pastrsp, interaction, Q2], shifted(:, j));
        %                 CPD_Q2_rotated(j,k) = (mdl1_withoutQ1Q2H.Deviance - mdl1_Q2.Deviance); 
        %                 
        %                 mdl1_H = fitglm([fb, rsp, pastrsp, interaction, h], shifted(:, j));
        %                 CPD_H_rotated(j,k) = (mdl1_withoutQ1Q2H.Deviance - mdl1_H.Deviance); 

                    end

        %             if  CPD_Q1(j)>prctile(CPD_Q1_rotated(j,:),95)  
        %                  signifRegQ1(j) = 1;
        %             end
        %             
        %             if  CPD_Q2(j)>prctile(CPD_Q2_rotated(j,:),95)
        %                  signifRegQ2(j) = 1;
        %             end
        %             
        %             if  CPD_H(j)>prctile(CPD_H_rotated(j,:),95) 
        %                  signifRegH(j) = 1;
        %             end

                      if CPD_Rsp(j)>prctile(CPD_Rsp_rotated(j,:),95)  
                           signifRegRsp(j) = 1;
                      else
                          signifRegRsp(j) = 0;
                      end
                      for t = 1:100
                          if CPD_Rsp(j)>prctile(CPD_Rsp_rotated(j,:),t)  
                           signifPct(t) = signifPct(t)+1;
                          end
                      end
                      if CPD_Rsp(j)<prctile(CPD_Rsp_rotated(j,:),1)
                          j
                      end

            end
            sigRsp = c_location(use_cluster(find(signifRegRsp==1)),1);
    %         sigQ1 = c_location(find(signifRegQ1==1),1);
    %         sigQ2 = c_location(find(signifRegQ2==1),1);
    %         sigH = c_location(find(signifRegH==1),1);
            allLoc = c_location(use_cluster,1);
            totalCellCount = length(use_cluster);
            save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\regRspQ1full', sprintf('%s_%s_%s_signifRegRsp.mat',subject, date, probeTag{1})), 'signifRegRsp')
            save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\signifPctQ1full', sprintf('%s_%s_%s_signifPct.mat',subject, date, probeTag{1})), 'signifPct')
            save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\totalCellQ1full', sprintf('%s_%s_%s_totalCellCount.mat',subject, date, probeTag{1})), 'totalCellCount')
            save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\allLocsQ1full', sprintf('%s_%s_%s_allLoc.mat',subject, date, probeTag{1})), 'allLoc')
            save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ1full', sprintf('%s_%s_%s_sigLoc.mat',subject, date, probeTag{1})), 'sigRsp')
            save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\useClusterQ1full', sprintf('%s_%s_%s_useCluster.mat',subject, date, probeTag{1})), 'use_cluster')
            %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\regQ1', sprintf('%s_%s_%s_signifRegQ1.mat',subject, date, probeTag{1})), 'sigQ1')
    %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\regQ2', sprintf('%s_%s_%s_signifRegQ2.mat',subject, date, probeTag{1})), 'sigQ2')
    %         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\regH', sprintf('%s_%s_%s_signifRegH.mat',subject, date, probeTag{1})), 'sigH')


    end
end






mat = dir('C:\Users\annaL\Documents\PhD\value_coding\sigLocs\*.mat'); 
for q = 1:length(mat) 
    signReg{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\sigLocs\' mat(q).name]);
end

mat = dir('C:\Users\annaL\Documents\PhD\value_coding\allLocs\*.mat'); 
for q =  1:length(mat) %1:length(mat) 
    allReg{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\allLocs\' mat(q).name]);
end

allSigRegions = [];
for j = 1:length(signReg)
    for k =1:height(signReg{1, j}.sigRsp)
        allSigRegions = [allSigRegions, signReg{1, j}.sigRsp{k,1}];
    end
end

allSnSRegions = [];
for j =  1:length(mat) %1:length(allReg)
    for k =1:height(allReg{1, j}.allLoc)
        allSnSRegions = [allSnSRegions, allReg{1, j}.allLoc{k,1}];
    end
end


brain_groups = {["VISa", "VISam", "VISl", "VISp", "VISpm", "VISrl"]; % visual cortex
    ["RSP"]; % RSP
    ["ACA", "AUD", "COA", "DP", "ILA", "MOp", "MOs", "OLF", "ORB", "ORBm", "PIR", "PL", "SSp", "SSs", "TT"]; % non-visual cortex except RSP
    ["CL", "LD", "LGd", "LH", "LP", "MD", "MG", "PO", "POL", "PT", "RT", "SPF", "TH", "VAL", "VPL", "VPM"]; % thalamus
    ["CA", "CA1", "CA2", "CA3", "DG", "SUB", "POST"]; % hippocampal
    ["APN", "IC", "MB", "MRN", "NB", "PAG", "RN", "ZI"]; % midbrain except SC
    ["SC"]; % SC
    ["ACB", "CP", "GPe", "LS", "LSc", "LSr", "MS", "OT", "SNr", "SI"];% basal ganglia
    ["BLA", "BMA", "EP", "EPd", "MEA"] % cortical subplate
    }

bg_counter = zeros(1, length(brain_groups));
otherRegions = {};
for i = 1:length(allSigRegions)
    flag = 0;
    for k = 1:length(brain_groups)
        for j = 1:length(brain_groups{k})
            if contains(allSigRegions{i},brain_groups{k}(j))
                bg_counter(k) = bg_counter(k)+1;
                flag = 1;
            end
        end
    end
    if flag == 0
        if ~ismember(allSigRegions{i},otherRegions)
            otherRegions{end+1} = allSigRegions{i};
        end
    end
end


bg_counter2 = zeros(1, length(brain_groups));
otherRegions2 = {};
for i = 1:length(allSnSRegions)
    flag = 0;
    for k = 1:length(brain_groups)
        for j = 1:length(brain_groups{k})
            if contains(allSnSRegions{i},brain_groups{k}(j))
                bg_counter2(k) = bg_counter2(k)+1;
                flag = 1;
            end
        end
    end
    if flag == 0
        otherRegions2{end+1} = allSnSRegions{i};
    end
end

figure;
X = categorical({'vis ctx', 'RSP', 'other ctx', 'thal', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate', 'other'});
X = reordercats(X,{'vis ctx', 'RSP', 'other ctx', 'thal', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate', 'other'});
Y = [bg_counter, numel(otherRegions)];
Y2 = [bg_counter2, numel(otherRegions2)];
bar(X,Y2)
hold on
bar(X,Y)
title('number of cells in each region')
legend('total', 'significantly coding past fb at the time of movement')

X = categorical({'vis ctx', 'RSP', 'other ctx', 'thal', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate'});
X = reordercats(X,{'vis ctx', 'RSP', 'other ctx', 'thal', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate'});
Y = bg_counter./bg_counter2;
figure; 
bar(X,Y)
hold on; plot(xlim, [0.05, 0.05], 'r')
title('percentages significant in each region')



clear signifPct
clear totalCellCount

mat = dir('C:\Users\annaL\Documents\PhD\value_coding\signifPctInter\*.mat');  
for q = 1:34%1:length(mat) 
    signifPct{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\signifPctInter\' mat(q).name]);
end
mat = dir('C:\Users\annaL\Documents\PhD\value_coding\totalCellCount\*.mat'); 
for q = 1:34%1:length(mat) 
    totalCellCount{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\totalCellCount\' mat(q).name]);
end

figure;hold on;
for i =1:34%1:length(mat)
plot(signifPct{1, i}.signifPct  /totalCellCount{1, i}.totalCellCount  )
end
hold on; plot([0,100], [1,0], 'k')
xlabel('percentile of permuted data')
ylabel('percentage of cells better than Xth percentile of randomly permuted data')
title('coding of interaction between previous feedback and action')


%%
% 1.
% -100:0 before feedback
% mdl0_withoutRsp = fitglm([pastrsp, rsp])
% mdl0_Rsp = fitglm([pastrsp, rsp, pastfb])
% C:\Users\annaL\Documents\PhD\value_coding\regRsp
% C:\Users\annaL\Documents\PhD\value_coding\signifPct 
% C:\Users\annaL\Documents\PhD\value_coding\allLocs
% C:\Users\annaL\Documents\PhD\value_coding\sigLocs
% is working reasonably well, all sessions

% 2.
% -100:0 before feedback
% mdl0_withoutRsp = fitglm([pastrsp, rsp, pastfb])
% mdl0_Rsp = fitglm([pastrsp, rsp, pastfb, pastrsp.*pastfb])
% C:\Users\annaL\Documents\PhD\value_coding\regRspInter
% C:\Users\annaL\Documents\PhD\value_coding\signifPctInter
% C:\Users\annaL\Documents\PhD\value_coding\allLocsInter
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsInter
% works for some sessions


% 3.
% -100:0 before feedback
% mdl0_withoutRsp = fitglm([pastrsp, rsp])
% mdl0_Rsp = fitglm([pastrsp, rsp, Q1])
% C:\Users\annaL\Documents\PhD\value_coding\regRspQ1
% C:\Users\annaL\Documents\PhD\value_coding\signifQ1
% C:\Users\annaL\Documents\PhD\value_coding\allLocsQ1
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ1
% works for some sessions

% 4.
% -100:0 before feedback
% mdl0_withoutRsp = fitglm([pastrsp, rsp])
% mdl0_Rsp = fitglm([pastrsp, rsp, h])
% C:\Users\annaL\Documents\PhD\value_coding\regRspH
% C:\Users\annaL\Documents\PhD\value_coding\signifH
% C:\Users\annaL\Documents\PhD\value_coding\allLocsH
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsH
% works for some sessions

% 5.
% -100:0 before movement
% mdl0_withoutRsp = fitglm([pastrsp, rsp])
% mdl0_Rsp = fitglm([pastrsp, rsp, Q1])
% C:\Users\annaL\Documents\PhD\value_coding\regRspQ1m
% C:\Users\annaL\Documents\PhD\value_coding\signifQ1m
% C:\Users\annaL\Documents\PhD\value_coding\allLocsQ1m
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ1m
% works for some sessions

% 6.
% -100:0 before movement
% mdl0_withoutRsp = fitglm([pastrsp, rsp])
% mdl0_Rsp = fitglm([pastrsp, rsp, Q2])
% C:\Users\annaL\Documents\PhD\value_coding\regRspQ2m
% C:\Users\annaL\Documents\PhD\value_coding\signifQ2m
% C:\Users\annaL\Documents\PhD\value_coding\allLocsQ2m
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ2m
% works for some sessions

% 6.
% -100:0 before movement
% mdl0_withoutRsp = fitglm([pastrsp, rsp])
% mdl0_Rsp = fitglm([pastrsp, rsp, H])
% C:\Users\annaL\Documents\PhD\value_coding\regRspHm
% C:\Users\annaL\Documents\PhD\value_coding\signifHm
% C:\Users\annaL\Documents\PhD\value_coding\allLocsHm
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsHm
% works for some sessions

%7. 
% -100:0 before movement
% mdl0_withoutRsp = fitglm([pastrsp(1:len2), rsp(1:len2), pastfb(1:len2), pastrsp(1:len2).*pastfb(1:len2)], st_bins(1:len2, use_cluster(j)), 'linear');
% mdl0_Rsp = fitglm([pastrsp(1:len2), rsp(1:len2), pastfb(1:len2), pastrsp(1:len2).*pastfb(1:len2), Q1(1:len2)], st_bins(1:len2, use_cluster(j)), 'linear');
% C:\Users\annaL\Documents\PhD\value_coding\regRspQ1full
% C:\Users\annaL\Documents\PhD\value_coding\signifQ1full
% C:\Users\annaL\Documents\PhD\value_coding\allLocsQ1full
% C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ1full              