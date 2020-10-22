%for sessionID = [1:7,9:11,13:18]
for sessionID = [1]
%for sessionID = [5]
    clearvars -except sessionID;
    % close all;
    
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
        pastrsp = double(block.events.responseValues);
        pastrsp = [1 pastrsp(1:end-1)];
        pastrsp = pastrsp';
        interaction = pastrsp.*pastfb;
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
        c_locationAll = readtable(fullfile(dataDir, 'clusters.brainLocation.csv'));
        c_location = c_locationAll(:,4);
        
        brain_groups = {
            ["RSP"]; % RSP
            ["FRP", "MOp", "MOs", "SSp", "SSs", "GU", "AUD", "VIS", "ACA", "PL", "ILA", "ORB", "AI", "PTLp", "TEa", "PERI", "ECT", "OLF", "MOB","AOB","AON", "DP", "PIR", "TT", "NLOT", "COA", "PAA", "TR"]; % non-visual cortex except RSP
            ["CL", "LD", "LGd", "LH", "LP", "MD", "MG", "PO", "PT", "RT", "SPF", "TH", "VAL", "VPL", "VPM", "Eth", "HY", "LGv", "MH", "PoT", "SubG"]; % thalamus, hypothalamus
            ["CA1", "CA2", "CA3", "DG", "SUB", "POST", "FC", "IG", "HPF"]; % hippocampal
            ["APN", "IC", "MB", "NB", "SAG", "PBG", "PAG", "NPC", "SNr", "VTA", "RR", "MRN", "PRT", "MRT", "NOT", "OP", "PPT", "CUN", "RPF"]; % midbrain except SC
            ["SC"]; % SC
            ["ACB", "CP", "FS", "GPe", "LSc", "LSr", "MS", "OT", "SNr", "SI", "STR", "BST", "NDB"];% basal ganglia
            ["LA", "BMA", "EP", "MEA", "CTXsp"]; % cortical subplate
            ["ZI", "FF"] %ZI
            ["IIn", "act", "alv", "chpl", "cing", "csc", "dhc", "fi", "fiber tracts", "och", "opt", "root"] %fiber tracts
            }
        
        time_period_beginning = [-0.7, -0.6, -0.5,-0.4,-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        time_period_end = [-0.6, -0.5, -0.4,-0.3,-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
        for region = [2] %brain group I'm interested in
            tic
            d = [];
            d2 = [];
            dev_nr_all = [];
            dev_r_all = [];
            dev_nr = 0;
            dev_r = 0;
            for tt = 1:length(time_period_beginning)
                neuron_in_BLA = [];
                for j = 1:length(allClu)
                    currClu = allClu(j);
                    currSt = st(find(sc==currClu));
                    for i = 1:length(ft)
                        st_bins(i, j) = numel(find(currSt<mt_before_ft(i)+time_period_end(tt) & currSt>mt_before_ft(i)+time_period_beginning(tt))); % trial*cluster
                    end
                end
                use_cluster = [];
                for j = 1:length(allClu)
                    if sum(st_bins(:, j))>50
                        use_cluster(end+1) = j;
                    end
                end
                allLoc = c_location(use_cluster,1);
                
                
                for i = 1:height(allLoc)
                    for j = 1:length(brain_groups{region})
                        if contains(allLoc{i,1},brain_groups{region}(j)) %& ~ismember(i, neuron_in_BLA)
                            if ismember(i, neuron_in_BLA)
                                i
                            end
                            neuron_in_BLA(end+1) = i;
                            break;
                        end
                    end
                end
                
                if numel(unique(neuron_in_BLA)) ~= numel((neuron_in_BLA))
                    "problem with regions"
                    break
                end
                if ~isempty(neuron_in_BLA)
                    %  random_as_neurons_in_BLA = randsample(height(allLoc),length(neuron_in_BLA));
                    spikes_rotated = [];
                    %spikesTrial = st_bins(:, use_cluster(neuron_in_BLA));%(:, intersect(neuron_in_BLA, use_cluster));
                    spikesTrial = st_bins(:,use_cluster(neuron_in_BLA));%(:, intersect(neuron_in_BLA, use_cluster));
                    for k = 1:length(ft)
                        spikes_rotated(:,:,k) = circshift(spikesTrial, k); % along second axis
                    end
                    
                    [mdl0_reverse, FitInfo] = lassoglm(spikesTrial,(pastfb+1)/2, 'binomial');
                    dev_nr = min(FitInfo.Deviance);
                    dev_r = [];
                    for k = 1:length(ft)/20-1
                        shifted = spikes_rotated(:,:,k*20);
                        [mdl0_reverse_rotated, FitInfo2] = lassoglm(shifted,(pastfb+1)/2, 'binomial');
                        dev_r(end+1) = min(FitInfo2.Deviance);
                    end
                    %             figure; histogram(dev_r, 15); hold on; plot([dev_nr,dev_nr], ylim)
                    %             title(['deviance of binomial model predicting past feedback from spikes in cortical subplate, session ', subject, date, probeTag{1}])
                    %             legend('circularly shifted', 'not shifted')
                    %
                else
                    break
                end
                dev_nr_all(tt) = dev_nr;
                dev_r_all(tt) = mean(dev_r);
                d(tt) = prctile(dev_r, 5);
                d2(tt) = prctile(dev_r, 95);
            end
            
            figure; plot(dev_nr_all, 'r', 'Linewidth',6)
            hold on; plot(dev_r_all, 'b', 'Linewidth',6)
            hold on; plot(d, 'b')
            hold on; plot(d2, 'b')
            xticks([1:length(d)])
            xticklabels(time_period_beginning+0.05)
            hold on; plot([7.5,7.5], ylim, 'k')
            xlabel('time from movement onset')
            ylabel('deviance')
            legend('real data', 'rotated data', '5th and 95th percentiles')
            title(['ability of model to predict previous trial feedback; data from ', num2str(region), ' ', subject, ' ', date, ' ', probeTag{1}])
            toc
             
            savefig(fullfile('C:\Users\annaL\Documents\PhD\value_coding\regions', sprintf('%s_%s_%s_%s_region.fig',subject, date, probeTag{1}, num2str(region))))
            
            
        end
       
        % figure; subplot(2,1,1); plot(dev_nr_all, 'r', 'Linewidth',6)
        
        
        
        %
        %         d = [];
        %         d2 = [];
        %         dev_nr_all = [];
        %         dev_r_all = [];
        %         time_period_beginning = [-0.7, -0.6, -0.5,-0.4,-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        %         time_period_end = [-0.6, -0.5, -0.4,-0.3,-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
        %         for tt = 1:length(time_period_beginning)
        %             for j = 1:length(allClu)
        %                 currClu = allClu(j);
        %                 currSt = st(find(sc==currClu));
        %                 for i = 1:length(ft)
        %                     st_bins(i, j) = numel(find(currSt<ft(i)+time_period_end(tt) & currSt>ft(i)+time_period_beginning(tt))); % trial*cluster
        %                 end
        %             end
        %             use_cluster = [];
        %             for j = 1:length(allClu)
        %                 if sum(st_bins(:, j))>50
        %                     use_cluster(end+1) = j;
        %                 end
        %             end
        %             allLoc = c_location(use_cluster,1);
        %             bg_counter2 = zeros(1, length(brain_groups));
        %             otherRegions2 = {};
        %             neuron_in_BLA = [];
        %             for i = 1:height(allLoc)
        %                 for k = 2 %brain group I'm interested in
        %                     for j = 1:length(brain_groups{k})
        %                         if contains(allLoc{i,1},brain_groups{k}(j))% & ~ismember(i, neuron_in_BLA)
        %                             neuron_in_BLA(end+1) = i;
        %                         end
        %                     end
        %                 end
        %             end
        %
        %             if ~isempty(neuron_in_BLA)
        %                 random_as_neurons_in_BLA = randsample(height(allLoc),length(neuron_in_BLA));
        %                 spikes_rotated = [];
        %                 %spikesTrial = st_bins(:, use_cluster(neuron_in_BLA));%(:, intersect(neuron_in_BLA, use_cluster));
        %                 spikesTrial = st_bins(:,use_cluster(neuron_in_BLA));%(:, intersect(neuron_in_BLA, use_cluster));
        %                 for k = 1:length(ft)
        %                     spikes_rotated(:,:,k) = circshift(spikesTrial, k); % along second axis
        %                 end
        %
        %                 mdl0_reverse = fitglm(spikesTrial,(pastfb+1)/2, 'linear', 'Distribution','binomial');
        %                 dev_nr = mdl0_reverse.Deviance;
        %                 dev_r = [];
        %                 for k = 1:length(ft)/10-1
        %                     shifted = spikes_rotated(:,:,k*10);
        %                     mdl0_reverse_rotated = fitglm(shifted,(pastfb+1)/2,'linear', 'Distribution','binomial');
        %                     dev_r(end+1) = mdl0_reverse_rotated.Deviance;
        %                 end
        %     %             figure; histogram(dev_r, 15); hold on; plot([dev_nr,dev_nr], ylim)
        %     %             title(['deviance of binomial model predicting past feedback from spikes in cortical subplate, session ', subject, date, probeTag{1}])
        %     %             legend('circularly shifted', 'not shifted')
        %
        %                 dev_nr_all(tt) = dev_nr;
        %                 dev_r_all(tt) = mean(dev_r);
        %                 d(tt) = prctile(dev_r, 5);
        %                 d2(tt) = prctile(dev_r, 95);
        %             end
        %         end
        %         subplot(2,1,2); plot(dev_nr_all, 'r', 'Linewidth',6)
        %         hold on; plot(dev_r_all, 'b', 'Linewidth',6)
        %         hold on; plot(d, 'b')
        %         hold on; plot(d2, 'b')
        %         xticks([1:length(d)])
        %         xticklabels(time_period_beginning+0.05)
        %         hold on; plot([7.5,7.5], ylim, 'k')
        %         xlabel('time from feedback')
        %         ylabel('deviance')
        %         legend('real data', 'rotated data', '5th and 95th percentiles')
        %         title(['ability of model to predict previous trial feedback; data from cortical subplate only', subject, date, probeTag{1}])
        %
        %
        %
        
        
        
        
        
    end
end
