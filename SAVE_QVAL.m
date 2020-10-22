function [Q1session, Q2session, habits_session] = SAVE_QVAL(sessionId)

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
    
    %     trials_before = cumsum(mousedata.trial_number_session(1:sess-1));
    %     starting_trial = trials_before(end)+1;
    %     final_trial = trials_before(end)+mousedata.trial_number_session(sess);
    for i = 1:length(mousedata.dates)
        DD(i) = strcmp([mousedata.dates{i}], date);
    end
    session = find(DD(:)==1);
    %% Define folders
    subjectsFolder = getRootDir(subject, date);
    alignDir = fullfile(subjectsFolder,  'alignments');
    codingPlotFolder = 'C:\Users\annaL\Documents\PhD\ephys_results\depth_graded';
    
    %% load the dataset
    Q1All = load(['q1_' subject '_all']);
    Q2All = load(['q2_' subject '_all']);
    habitsAll = load(['habits_' subject '_all']);
    %608 568
    Q1session = Q1All.q1(session)';
    Q2session = Q2All.q2(session)';
    habits_session = habitsAll.habits(session)';
    %writeNPY(Q1session, [subject '_' date '_Q1.npy'])
    %writeNPY(Q2session, [subject '_' date '_Q2.npy'])
    %writeNPY(habits_session, [subject '_' date '_h.npy'])
end
