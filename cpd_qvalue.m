% first run analyse_sides_multiple
%then modeling_workflow_likelihoods
%then cpd_qvalue


for sessionId = 1
    %for k = 9
    clearvars -except sessionId;
    close all;
    
    db(1).subject = 'AL021';
    db(1).date = '2019-06-05'; %k2, k3
    db(2).subject = 'AL021';
    db(2).date = '2019-06-06'; %k1,k2,k3, zo
    db(3).subject = 'AL021';
    db(3).date = '2019-06-07';%k1,k2
    db(4).subject = 'AL022';
    db(4).date = '2019-06-19';%k1,k2,k3,zo
    db(5).subject = 'MW003';
    db(5).date = '2019-08-11';%k1,k2,zo
    db(6).subject = 'MW003';
    db(6).date = '2019-08-12';%k1,k2
    db(7).subject = 'AL029';
    db(7).date = '2019-10-23';%k3,zo
    db(8).subject = 'AL028';
    db(8).date = '2019-10-29';%k2
    db(9).subject = 'AL026';
    db(9).date = '2019-11-01';%k1, k2
    db(10).subject = 'AL026';
    db(10).date = '2019-11-02';%k1, k2
    db(11).subject = 'AL026';
    db(11).date = '2019-11-03';%k1, k2
    db(12).subject = 'AL028';
    db(12).date = '2019-11-04';%k2
    db(13).subject = 'AL026';
    db(13).date = '2019-11-05';%k1, k2
    db(14).subject = 'AL026';
    db(14).date = '2019-11-06';%k1, k2
    db(15).subject = 'AL026';
    db(15).date = '2019-11-07';%k1, k2
    db(16).subject = 'AL026';
    db(16).date = '2019-11-12';%k2
    db(17).subject = 'AL026';
    db(17).date = '2019-11-13';%k2
    db(18).subject = 'AL016';
    db(18).date = '2019-07-11';%k1,k2
    
    
    subject = db(sessionId).subject;
    date = db(sessionId).date;
    [tags, hasEphys] = getEphysTags(subject, date);
    tagsNum = length(tags);
    %%
    for t = 1:length(tags)
    clearvars -except t sessionId tagsNum subject date hasEphys tags;
    startOfRange = 30;
    endOfRange = 60;
    load(['mousedata_' subject '.mat'])
    C = mousedata.dates_allSession;
    index = false(1, numel(C));
    if sessionId==4
        sess = 26;
    else if sessionId==18
            sess = 39;
        else
            for sess = 1:numel(C)
                index(sess) = (string(C{sess}) == date);
                sess = find(index==1);
            end
        end
    end
    
    trials_before = cumsum(mousedata.trial_number_session(1:sess-1));
    starting_trial = trials_before(end)+1;
    final_trial = trials_before(end)+mousedata.trial_number_session(sess);
    %% Define folders
    subjectsFolder = getRootDir(subject, date);
    alignDir = fullfile(subjectsFolder,  'alignments');
    codingPlotFolder = 'C:\Users\annaL\Documents\PhD\ephys_results\depth_graded';
    
    %% load the dataset
    Q1All = load(['q1_' subject '_all']);
    Q2All = load(['q2_' subject '_all']);
    habitsAll = load(['habits_' subject '_all']);
    %608 568
    Q1session = Q1All.q1(starting_trial:final_trial);
    Q2session = Q2All.q2(starting_trial:final_trial);
    habits_session = habitsAll.habits(starting_trial:final_trial);
    sp = loadAllKsDir(subject, date);
    load(fullfile(subjectsFolder,'2', [date '_2_' subject '_Block.mat']))
    
    %% Get basic info
    
    [expNums, blocks, hasBlock, pars, isMpep, tl, hasTimeline] = ...
        dat.whichExpNums(subject, date);
    TLexp = expNums(hasTimeline);
    TLexp = TLexp(end);
    useFlipper = true;
    
        probeTag = tags(t);
        probe_number = t;
        myKsDirSorting = fullfile(subjectsFolder, strcat('\ephys_', probeTag), '\sorting');
        
        bBlocktoTL = readNPY([alignDir '\correct_block_2_to_timeline_1.npy']);
        bTLtoMaster = readNPY(fullfile(alignDir, ...
            sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{t})));
        
        %%
        
        % Align to movement
        
%         Fs = numel(block.inputs.wheelTimes)/ (block.inputs.wheelTimes(end) - block.inputs.wheelTimes(1));
%         [moveOnsets, moveOffsets, moveAmps, vel,t, pos] = findWheelMoves3(block.inputs.wheelMMValues, block.inputs.wheelMMTimes, Fs, []);
%         f_neg0 = moveOnsets;
%         
%         distances = [];
%         ind = [];
%         ind2 = [];
%         % in this cycle I only find those trials where movement happened after
%         % interactive on (because AL021 had exp def without interactive delay so sometimes would start moving early)
%         % ind is the index of the moveOnset that setisfies this criteria, and ind2
%         % is the number of the corresponding trial
%         for i = 1:size(block.events.interactiveOnTimes,2)-1
%             %             if ~isempty (find(moveOnsets(:)-block.events.trialNumTimes(i)>0 & moveOnsets(:)-block.events.feedbackTimes(i)<0,1))
%             %              [ind(end+1)] = find(moveOnsets(:)-block.events.trialNumTimes(i)>0 & moveOnsets(:)-block.events.feedbackTimes(i)<0,1);
%             if ~isempty (find(moveOnsets(:)-block.events.interactiveOnTimes(i)>0 & moveOnsets(:)-block.events.feedbackTimes(i)<0,1))
%                 [ind(end+1)] = find(moveOnsets(:)-block.events.interactiveOnTimes(i)>0 & moveOnsets(:)-block.events.feedbackTimes(i)<0,1);
%                 %
%                 ind2(end+1) = i;
% 
%             end
%         end
%         
%         moves(ind2) = moveOnsets(ind);
%         moves(setdiff(1:size(block.events.feedbackTimes,2),ind2)) = block.events.interactiveOnTimes(setdiff(1:size(block.events.feedbackTimes,2),ind2));

     %   f_neg0 = moves;
        f_neg0 = block.events.feedbackTimes;
        eventTimes0 = f_neg0;
        eventTimes0 = applyCorrection(eventTimes0, bBlocktoTL);
        eventTimes0 = applyCorrection(eventTimes0, bTLtoMaster);
        
        
        
        
        
        depthBinSize = 80; % in units of the channel coordinates, in this case µm
        timeBinSize = 0.05; % seconds
        bslWin = [-5.7 -5.5]; % window in which to compute "baseline" rates for normalization
        psthType = 'norm'; % show the normalized version
        eventName = 'movement'; % for figure labeling
        window = [-1 1];
        
        
        
%         
%         myKsDirSorting = string(myKsDirSorting);
%         [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDirSorting);
%         depthBinSize = 80; % in units of the channel coordinates, in this case µm
%         timeBinSize = 0.01; % seconds
%         bslWin = [-0.7 -0.5]; % window in which to compute "baseline" rates for normalization
%         psthType = 'norm'; % show the normalized version
%         eventName = 'reward delivery'; % for figure labeling
%         window = [-1 1];
%         [timeBins, depthBins, allP, ~] = psthByDepth(spikeTimes, spikeDepths, ...
%             depthBinSize, timeBinSize, eventTimes0, window, bslWin);
%         
%         figure;
%         plotPSTHbyDepth(timeBins, depthBins, allP, eventName, psthType);
%         
%         
%         
        
        
        
        trialGroups = ones(size(block.events.feedbackValues));
        for i = 1:size(block.events.feedbackValues,2)
            if block.events.responseValues(i) ~=1
                trialGroups(i) = 2;
            end
        end
        trialGroups = trialGroups;
        
        
        %% by depth
        %trying to predict spiking of neurons at a certain depth using behavioural variables

        spCount = [];
        for index = startOfRange:endOfRange
            clear spikeReward
            index
            for i = 1:numel(sp(probe_number).cids)
                i
                clear spikeReward
            %    i
                neurons = sp(probe_number).cids(i);
               % sp(probe_number).templateYpos(neurons+1)
                st = sp(probe_number).st(ismember(sp(probe_number).clu,neurons));
%                eventTimes0
                
              
                for k = 1:size(eventTimes0,1)
                    
                    spikeReward(k, index) = numel(find(st>eventTimes0(k)+(index-50)*0.1 & st<eventTimes0(k)+(index-49)*0.1));
                    
                    
                end
% 
%                 q_left = [Q1session];
%                 q_right = [Q2session];
%                 corLeftI(i, index) = corr(spikeReward(:,index), q_left');
%                 corRightI(i, index) = corr(spikeReward(:,index), q_right');
                binnedSpikes = spikeReward;
                fb = [block.events.feedbackValues];
                fb = fb';
                rsp = [block.events.responseValues>0];
                rsp = rsp';
                pastfb = block.events.feedbackValues;
                pastfb = [0 pastfb(1:end-1)];
                pastfb = pastfb';
                
                action = [block.events.responseValues>0];
                pastAction = [0 action(1:end-1)]';
                WSLS = pastAction.*pastfb;
%                 q1old = [0.5 Q1session(1:end-1)];
%                 q2old = [0.5 Q2session(1:end-1)];
%                 Hold = [0.5 habits_session(1:end-1)];
                q1 = [Q1session]; %q left
                q2 = [Q2session]; %q right
                h = [habits_session];
                
                
                for l = 0:3
                    
                     shifted = circshift(binnedSpikes,l*10);
                    mdl1 = fitglm([fb'; action; q1; q2; h]', shifted(:, index));
                    dev1(l+1, index+101,i) =mdl1.Rsquared.Ordinary;
                    
                    mdl_withoutFeedback = fitglm([action; q1; q2; h]', shifted(:, index));
                %    dev_withoutFeedback(l+1, index+101,i) =mdl_withoutFeedback.Rsquared.Ordinary;
                    
                    mdl_withoutAction = fitglm([fb'; q1; q2; h]', shifted(:, index));
                %    dev_withoutAction(l+1, index+101,i) =mdl_withoutAction.Rsquared.Ordinary;
                    
                    mdl_withoutQ1 = fitglm([fb'; action; q2; h]', shifted(:, index));
                 %   dev_withoutQ1(l+1, index+101,i) =mdl_withoutQ1.Rsquared.Ordinary;
                    
                    mdl_withoutQ2 = fitglm([fb'; action; q1; h]', shifted(:, index));
                %    dev_withoutQ2(l+1, index+101,i) =mdl_withoutQ2.Rsquared.Ordinary;
                    
                    mdl_withoutH = fitglm([fb'; action; q1; q2]', shifted(:, index));
                 %   dev_withoutH(l+1, index+101,i) =mdl_withoutH.Rsquared.Ordinary;
                    
                    
                    CPD_feedback(l+1, index+101,i) = (mdl_withoutFeedback.SSE - mdl1.SSE)/mdl_withoutFeedback.SSE;
                    CPD_action(l+1, index+101,i) = (mdl_withoutAction.SSE - mdl1.SSE)/mdl_withoutAction.SSE;                    
                    CPD_Q1(l+1, index+101,i) = (mdl_withoutQ1.SSE - mdl1.SSE)/mdl_withoutQ1.SSE;
                    CPD_Q2(l+1, index+101,i) = (mdl_withoutQ2.SSE - mdl1.SSE)/mdl_withoutQ2.SSE;                    
                    CPD_H(l+1, index+101,i) = (mdl_withoutH.SSE - mdl1.SSE)/mdl_withoutH.SSE;
                    
                end
            end
        end
        
     xs = 1:31;
    fill_between_lines = @(X, Y1, Y2, C) fill([X, fliplr(X)], [Y1 fliplr(Y2)], C, 'facealpha', 0.2, 'edgecolor', 'none');
    h = figure; 
    
    color = 'blue';
    ys = mean(nanmean(CPD_feedback(2:4, 101+startOfRange:101+endOfRange,:),3));
    err = std(squeeze(nanmean(CPD_feedback(2:4, 101+startOfRange:101+endOfRange,:),3)));
    hold on; subplot(2,3,1); plot(nanmean(CPD_feedback(1,101+startOfRange:101+endOfRange,:),3)); %hold on; plot(nanmean(matChoice)); hold on; plot(nanmean(matQ1));  hold on; plot(nanmean(matQ2)); hold on; plot(nanmean(matH)); 
    hold on; fill_between_lines(xs, ys+err, ys-err, color)
    xticks([1,11,21, 31])
    xticklabels({'-2sec', '-1sec', 'response', '+1sec'})
    title('CPD reward')
    
    color = 'red';
   ys = mean(nanmean(CPD_action(2:4, 101+startOfRange:101+endOfRange,:),3));
    err = std(squeeze(nanmean(CPD_action(2:4, 101+startOfRange:101+endOfRange,:),3)));
    hold on; subplot(2,3,2); plot(nanmean(CPD_action(1,101+startOfRange:101+endOfRange,:),3)); %hold on; plot(nanmean(matChoice)); hold on; plot(nanmean(matQ1));  hold on; plot(nanmean(matQ2)); hold on; plot(nanmean(matH)); 
    hold on; fill_between_lines(xs, ys+err, ys-err, color)
    xticks([1,11,21, 31])
    xticklabels({'-2sec', '-1sec', 'response', '+1sec'})
    title('CPD choice')
    
    color = 'yellow';
    ys = mean(nanmean(CPD_Q1(2:4, 101+startOfRange:101+endOfRange,:),3));
    err = std(squeeze(nanmean(CPD_Q1(2:4, 101+startOfRange:101+endOfRange,:),3)));
    hold on; subplot(2,3,3); plot(nanmean(CPD_Q1(1,101+startOfRange:101+endOfRange,:),3)); %hold on; plot(nanmean(matChoice)); hold on; plot(nanmean(matQ1));  hold on; plot(nanmean(matQ2)); hold on; plot(nanmean(matH)); 
    hold on; fill_between_lines(xs, ys+err, ys-err, color)
    xticks([1,11,21, 31])
    xticklabels({'-2sec', '-1sec', 'response', '+1sec'})    
    title('CPD Q1')
    
    color = 'magenta';
    ys = mean(nanmean(CPD_Q2(2:4, 101+startOfRange:101+endOfRange,:),3));
    err = std(squeeze(nanmean(CPD_Q2(2:4, 101+startOfRange:101+endOfRange,:),3)));
    hold on; subplot(2,3,4); plot(nanmean(CPD_Q2(1,101+startOfRange:101+endOfRange,:),3)); %hold on; plot(nanmean(matChoice)); hold on; plot(nanmean(matQ1));  hold on; plot(nanmean(matQ2)); hold on; plot(nanmean(matH)); 
    hold on; fill_between_lines(xs, ys+err, ys-err, color)
    xticks([1,11,21, 31])
    xticklabels({'-2sec', '-1sec', 'response', '+1sec'})
    title('CPD Q2')
    
    color = 'green';
    ys = mean(nanmean(CPD_H(2:4, 101+startOfRange:101+endOfRange,:),3));
    err = std(squeeze(nanmean(CPD_H(2:4, 101+startOfRange:101+endOfRange,:),3)));
    hold on; subplot(2,3,5); plot(nanmean(CPD_H(1,101+startOfRange:101+endOfRange,:),3)); %hold on; plot(nanmean(matChoice)); hold on; plot(nanmean(matQ1));  hold on; plot(nanmean(matQ2)); hold on; plot(nanmean(matH)); 
    hold on; fill_between_lines(xs, ys+err, ys-err, color)
    xticks([1,11,21, 31])
    xticklabels({'-2sec', '-1sec', 'response', '+1sec'})
    title('CPD H')
    
%     legend('fb', 'choice', 'q left', 'q right', 'h')
%     xticks([1,11,21])
%     xticklabels({'-1sec', 'response', '+1sec'})
    [sortD, sortOr] = sort(sp(probe_number).templateYpos(sp(probe_number).cids+1), 'descend');
    savefig(h, fullfile('C:\Users\annaL\Documents\PhD\ephys_results\cpd_qvalBeep\', ...
        sprintf('%s_%s_%s_valueCor.fig', subject, date, sp(probe_number).name)), 'compact')
    
    
   

   



    h = figure; plot(nanmean(squeeze(dev1(:, 49+101, :)),2))
    title('cpd')
    xticks(1:4)
    xticklabels([0:10:30])
    title('R squared of a model predicting firing rates 100 ms:action based on Qleft only')
    xlabel('shift')


    savefig(h, fullfile('C:\Users\annaL\Documents\PhD\ephys_results\cpd_qval\', ...
        sprintf('%s_%s_%s_qleftModel.fig', subject, date, sp(probe_number).name)), 'compact')
    

     
    end
end

