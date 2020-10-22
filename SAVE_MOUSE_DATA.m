for sessionID =1
    
    clearvars -except sessionID;
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
%     db(8).subject = 'AL028';
%     db(8).date = '2019-10-29';%k2
    db(9).subject = 'AL026';
    db(9).date = '2019-11-01';%k1, k2
    db(10).subject = 'AL026';
    db(10).date = '2019-11-02';%k1, k2
    db(11).subject = 'AL026';
    db(11).date = '2019-11-03';%k1, k2
%     db(12).subject = 'AL028';
%     db(12).date = '2019-11-04';%k2
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
        mkdir(dataDir)
        bTLtoMaster = readNPY(fullfile(alignDir, ...
            sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{t})));
        
        %% trials
        trials.response_choice = block.events.responseValues;
        trials.feedbackType = block.events.feedbackValues;
        responseTimes = block.events.responseTimes; 
        responseTimes = applyCorrection(responseTimes, bBlocktoTL);
        responseTimes = applyCorrection(responseTimes, bTLtoMaster);
        trials.response_times = responseTimes;
        interactiveOnTimes =block.events.interactiveOnTimes(1:length(responseTimes)); 
        interactiveOnTimes = applyCorrection(interactiveOnTimes, bBlocktoTL);
        interactiveOnTimes = applyCorrection(interactiveOnTimes, bTLtoMaster);
        trials.goCue_times = interactiveOnTimes;
        feedbackTimes = block.events.feedbackTimes; 
        feedbackTimes = applyCorrection(feedbackTimes, bBlocktoTL);
        feedbackTimes = applyCorrection(feedbackTimes, bTLtoMaster);
        trials.feedback_times = feedbackTimes;
        trials.probLeft = 0.8*(block.events.correctResponseValues==-1) + 0.2*(block.events.correctResponseValues==1);
        trials.probRight = 0.8*(block.events.correctResponseValues==1) + 0.2*(block.events.correctResponseValues==-1);
        trials.probLeft = trials.probLeft(1:length(responseTimes));
        trials.probRight = trials.probRight(1:length(responseTimes));
        writeNPY(trials.response_choice, fullfile(dataDir, 'trials.response_choice.npy'))
        writeNPY(trials.feedbackType, fullfile(dataDir, 'trials.feedbackType.npy'))
        writeNPY(trials.response_times, fullfile(dataDir, 'trials.response_times.npy'))
        writeNPY(trials.goCue_times, fullfile(dataDir, 'trials.goCue_times.npy'))
        writeNPY(trials.feedback_times, fullfile(dataDir, 'trials.feedback_times.npy'))
        writeNPY(trials.probLeft, fullfile(dataDir, 'trials.probLeft.npy'))
        writeNPY(trials.probRight, fullfile(dataDir, 'trials.probRight.npy'))
        %% wheel
        rotEnc = double(Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder')));
        rotEnc(rotEnc>2^31) = rotEnc(rotEnc>2^31)-2^32;
        rotEncD = diff([0; rotEnc]);
        wheelVel = conv(rotEncD, myGaussWin(0.03, Fs), 'same');
        wheel.speed = wheelVel;
        wheel.position = rotEnc;
        tlTimes = Timeline.rawDAQTimestamps;
        tlTimes = applyCorrection(tlTimes, bTLtoMaster);
        et.wheel = tlTimes;
        alf.writeTimeseries(dataDir, 'wheel', et.wheel, [], []);
        wheel.timestamps = readNPY(fullfile(dataDir, 'wheel.timestamps.npy'));
        writeNPY(wheel.speed, fullfile(dataDir, 'wheel.speed.npy'))
        writeNPY(wheel.position, fullfile(dataDir, 'wheel.position.npy'))
        writeNPY(wheel.timestamps, fullfile(dataDir, 'wheel.timestamps.npy'))
        %% eye
        processedFile = fullfile(TLdir, 'eye_proc.mat');
        if ~exist (processedFile, 'file')
            fprintf(1, 'eye movie not processed\n');
        else
            load(processedFile);
            eye.area = proc.pupil.area;
            eye.xyPos = proc.pupil.com;
            %et.eye = tlTimes;
            %alf.writeTimeseries(dataDir, 'eye', et.eye, [], []);
            vidTs = load(fullfile(subjectsFolder,'1', 'eye_timestamps.mat'));
            eyeTs = vidTs.tVid;
            eyeTs = applyCorrection(eyeTs, bTLtoMaster);
            eye.timestamps(:,2) = eyeTs;
            eye.timestamps(:,1) = 1:length(eye.timestamps);
            writeNPY(eye.area, fullfile(dataDir, 'eye.area.npy'))
            writeNPY(eye.timestamps, fullfile(dataDir, 'eye.timestamps.npy'))
            writeNPY(eye.xyPos, fullfile(dataDir, 'eye.xyPos.npy'))
        end
        %% spikes
        sp = loadAllKsDir(subject, date);
        spikes.clusters = sp(t).spikeTemplates;
        spikes.times = sp(t).st;
       % s = spikes.times;
       % writeNPY(s, fullfile(dataDir, 'spikes.times.npy'))
        writeNPY(spikes.times, fullfile(dataDir, 'spikes.times.npy'))
        writeNPY(spikes.clusters, fullfile(dataDir, 'spikes.clusters.npy'))
        
        %% licks
        lick.deviation = double(Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name}, 'piezoLickDetector')));
        et.lick = tlTimes;
        alf.writeTimeseries(dataDir, 'lick', et.lick, [], []);
        lick.timestamps = readNPY(fullfile(dataDir, 'lick.timestamps.npy'));
        writeNPY(lick.deviation, fullfile(dataDir, 'lick.deviation.npy'))
        writeNPY(lick.timestamps, fullfile(dataDir, 'lick.timestamps.npy'))
        
        %% movements
        Fs = 1000;
        [moveOnsets, moveOffsets, moveAmps, vel,t, pos] = findWheelMoves3(block.inputs.wheelMMValues, block.inputs.wheelMMTimes, Fs, []);
        moveOnsets = applyCorrection(moveOnsets, bBlocktoTL);
        moveOnsets = applyCorrection(moveOnsets, bTLtoMaster);
        moveOffsets = applyCorrection(moveOffsets, bBlocktoTL);
        moveOffsets = applyCorrection(moveOffsets, bTLtoMaster);
        writeNPY(moveOnsets, fullfile(dataDir, 'moves.moveOnsets.npy'))
        writeNPY(moveOffsets, fullfile(dataDir, 'moves.moveOffsets.npy'))
        writeNPY(moveAmps, fullfile(dataDir, 'moves.moveAmps.npy'))
        
        %% Q, H
        %analysis_sides_multiple;
        %modeling_workflow_likelihoods;
        %cpd_qvalue; optional
        [Q1session, Q2session, habits_session] = SAVE_QVAL(sessionID);
        writeNPY(Q1session, fullfile(dataDir, 'Q1.npy'))
        writeNPY(Q2session, fullfile(dataDir, 'Q2.npy'))
        writeNPY(habits_session, fullfile(dataDir, 'h.npy'))
        CORR_VALUE;
    end
end
% 
% bins = zeros(length(trials.feedback_times), 250)
% for i = 1:length(trials.feedback_times)
%     for j = 1:250
%         bins(i,j) = trials.feedback_times(i) + (j-150)*0.01;
%     end
% end
% 
% eps = 0.01
% for i = 1:length(trials.feedback_times)
%     for j = 1:250
%         index = find(abs(eye.timestamps - trials.feedback_times(i, j))< eps,1);
%         eyeMat(i,j) = 
%     end
% end
