%for sessionID = [2:7,9:11,13:18]
for sessionID = 5

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
        
        fb = double([block.events.feedbackValues]);
            fb(find(fb==0)) = -1;
            fb = fb';
            rsp = double(block.events.responseValues);
            rsp = rsp';
            pastfb = double(block.events.feedbackValues);
            pastfb(find(pastfb==0)) = -1;
            pastfb = [1 pastfb(1:end-1)];
            pastfb = pastfb';
            
               %% spikes
        sp = loadAllKsDir(subject, date);
        spikes.clusters = sp(t).spikeTemplates;
        spikes.times = sp(t).st;
       % s = spikes.times;
       % writeNPY(s, fullfile(dataDir, 'spikes.times.npy'))
        st = readNPY(fullfile(dataDir, 'spikes.times.npy'));
        sc = readNPY(fullfile(dataDir, 'spikes.clusters.npy'));
        ft = readNPY(fullfile(dataDir, 'trials.feedback_times.npy'));
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
        for j = 1:length(allClu)
            currClu = allClu(j);
            currSt = st(find(sc==currClu));
            for i = 1:length(ft)
                st_bins(i, j) = numel(find(currSt<ft(i) & currSt>ft(i)-0.1)); % trial*cluster
            end
        end
        use_cluster = [];
        for j = 1:length(allClu)
            if sum(st_bins(:, j))>50
                use_cluster(end+1) = j;
            end
        end
        
         c_locationAll = readtable(fullfile(dataDir, 'clusters.brainLocation.csv'));
         c_location = c_locationAll(:,4);
         allLoc = c_location(use_cluster,1);
        l = load(fullfile('C:\Users\annaL\Documents\PhD\value_coding\regRsp', sprintf('%s_%s_%s_signifRegRsp.mat',subject, date, probeTag{1})), 'signifRegRsp')
        signifRegRsp = l.signifRegRsp;
        cl_num = use_cluster(find(signifRegRsp==1));
        crr = corr(st_bins(:, cl_num(:)), pastfb);
        
        onefb = randsample(find(pastfb==1),numel(find(pastfb~=1)));
        randfb1 = randsample(1:numel((pastfb~=1)),numel(find(pastfb~=1)));
        randfb2 = randsample(1:numel((pastfb~=1)),numel(find(pastfb~=1)));
        figure; histogram(10*st_bins(onefb, cl_num(5)), [0:10:160]);
        hold on; histogram(10*st_bins(find(pastfb~=1), cl_num(5)),  [0:10:160])
%          figure; histogram(st_bins(randfb1, cl_num(5)));
%         hold on; histogram(st_bins(randfb2, cl_num(5)))
        c_location(cl_num(5),1)
        xlabel('firing rate, spikes/s')
        ylabel('number of trials')
        legend('positive past feedback', 'negative past feedback')
        
        
%         sigRsp = c_location(find(signifRegRsp==1),1);
%         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\allLocsInter', sprintf('%s_%s_%s_allLoc.mat',subject, date, probeTag{1})), 'allLoc')
%         save(fullfile('C:\Users\annaL\Documents\PhD\value_coding\sigLocsInter', sprintf('%s_%s_%s_allLoc.mat',subject, date, probeTag{1})), 'sigRsp')
%         
    end
end
