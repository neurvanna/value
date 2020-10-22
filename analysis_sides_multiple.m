% first run analyse_sides_multiple
%then modeling_workflow_likelihoods
%then cpd_qvalue
close all
clear all
for num = 1
  %  close all
  %  clear all
    mouseName = 'AL040';

%    mouseName = mnAll((num-1)*5+1:(num-1)*5+5);
    %myfolder = '\\zserver.cortexlab.net\Data\Subjects'; %AL023 and further
   myfolder = '\\znas.cortexlab.net\Subjects'; %MW003 and before
    list = dir(fullfile(myfolder,mouseName));
    dirFlags = [list.isdir];
    % Extract only those that are directories.
    subFolders = list(dirFlags);
    subFolders = subFolders(3:end);
    expType = zeros(1, size(subFolders,1));
    previous_trials = 0;
    choices_all = [];
    feedback_all =[];
    date_all = [];
    correct_all = [];
    trials_before_switch_to_ephysRig = 0;
    trials_before_switch_to_sides3 = 0;
    trials_before_switch_to_sides2 = 0;
    rt_all = [];
    flag1 = 0;
    flag2 = 0;
    flag3 = 0;
    for j = size(subFolders,1)-size(list,1)+5:size(subFolders,1)-1
%  for j = 46:size(subFolders,1)
   %     for j =134:136
  %  for j = size(subFolders,1)-23:size(subFolders,1)-3
        date = subFolders(j).name;
        list2 = dir(fullfile(myfolder,mouseName, date));
        list2 = list2(3:end);
        ind = 0;
        for i = 1:size(list2)
            if ~isempty(str2num(list2(i).name))
                ind = ind+1;
            end
        end
        experNum = num2str(max(str2num([list2(ind).name]')));
  % experNum = num2str(2);
        filename = strcat(date, '_', experNum, '_', mouseName, '_', 'Block.mat');
        f = fullfile(myfolder,mouseName, date, experNum, filename)
        load(f);
        filename2 = strcat(date, '_', experNum, '_', mouseName, '_', 'parameters.mat');
        f2 = fullfile(myfolder,mouseName, date, experNum, filename2);
        load(f2);
        if (contains(block.expDef, 'sides'))
            if (contains(block.rigName, 'zgood') & flag1 ==0)
                trials_before_switch_to_ephysRig = previous_trials;
                flag1 = 1;
            end
            if (contains(block.expDef, 'sides2.m') & flag2 ==0)
                trials_before_switch_to_sides2 = previous_trials; 
                flag2 = 1;
            elseif (contains(block.expDef, 'sides3.m') & flag3 ==0)
                trials_before_switch_to_sides3 = previous_trials;
                flag3 = 1;
                firstSessionToConsider = j;
            end
            thisFolder = strcat(myfolder, '\', mouseName, '\', date, '\', experNum);
            choices_all = horzcat(choices_all,block.events.responseValues);
            feedback_all = horzcat(feedback_all,block.events.feedbackValues);
            correct_all = horzcat(correct_all, block.events.correctResponseValues(1:length(block.events.responseValues)));
            rt_all = horzcat(rt_all, block.events.responseTimes - block.events.interactiveOnTimes(1:length(block.events.responseTimes)));
            previous_trials = previous_trials+length(block.events.responseValues);
            dateAll = repelem({date}, length(block.events.responseValues));
            date_all = horzcat(date_all, dateAll);
            if flag3==1
%                 choices_allSession{j,:} = block.events.responseValues;
%                 feedback_allSession{j,:} = block.events.feedbackValues;
%                 correct_allSession{j, :} = block.events.correctResponseValues(1:length(block.events.responseValues));
                 dates_allSession{j - firstSessionToConsider +1} = date;
                trial_number_session(j - firstSessionToConsider + 1) = length(block.events.responseValues);
            end
    
        end
    end
    figure('Position', [300 300 1500 800]); yyaxis left;plot(movmean(correct_all==choices_all,500))
    ylim([0,1])
    yyaxis right; plot(movmean(rt_all,[500,0]));ylim([0,12]);
    hold on;plot([trials_before_switch_to_sides2 trials_before_switch_to_sides2], [ylim])
    hold on;plot([trials_before_switch_to_sides3 trials_before_switch_to_sides3], [ylim])
    hold on;plot([trials_before_switch_to_ephysRig trials_before_switch_to_ephysRig], [ylim])
    legend('accuracy', 'reaction time', 'switch of task', 'switch of task', 'switch to ephys rig')
    title(['accuracy and reaction time ', mouseName])
end



% save(['DATA_choices_allSession_' mouseName '.mat'], 'choices_allSession');
% save(['DATA_feedback_allSession_' mouseName '.mat'], 'feedback_allSession');
% save(['DATA_correct_allSession_' mouseName '.mat'], 'correct_allSession');
% save(['DATA_dates_allSession_' mouseName '.mat'], 'dates_allSession')








%after running analysis_sides_multiple
corresp = correct_all;
response = choices_all;
rewards = feedback_all;
probright = [];
probleft = [];
for i = 1:trials_before_switch_to_sides2
    if corresp(i)==1
        probleft(i)= 0;
        probright(i)= 1;
    else
        probleft(i)= 1;
        probright(i)= 0;
    end
end

for i = trials_before_switch_to_sides2+1:trials_before_switch_to_sides3
    if corresp(i)==1
        probleft(i)= 0.1;
        probright(i)= 0.9;
    else
        probleft(i)= 0.9;
        probright(i)= 0.1;
    end
end

for i = trials_before_switch_to_sides3+1:length(choices_all)
    if corresp(i)==1
        probleft(i)= 0.2;
        probright(i)= 0.8;
    else
        probleft(i)= 0.8;
        probright(i)= 0.2;
    end
end
        
tr_num = length(choices_all);

loaded = load('ratdata');
ratdata = loaded.ratdata;
probleft_all = [probleft]; 
probright_all = [probright];
response_all = [response]; 
rewards_all = [rewards];
probrigth_all = [probright]; 
trial_types=ratdata.trial_types  ;

if tr_num>length(trial_types)
    
    trial_types(length(trial_types)+1:tr_num) = 'f';
    trial_types(trial_types ~= 'f') ='f';
else
    trial_types = trial_types(1:tr_num);
    trial_types(trial_types ~= 'f') ='f';
end
ratdata.nTrials = tr_num;
ratdata.dates = date_all;
ratdata.sides = response_all;
ratdata.rewards = rewards_all;
ratdata.left_prob1 = probleft_all;
ratdata.right_prob1 = probrigth_all;
ratdata.right_prob1 = probright_all;
ratdata.ratname = mouseName;
ratdata.trial_types = trial_types;
ratdata2 = loaded.ratdata;

ratdata.sides = '';
for i = 1:length(response_all)
if response_all(i)==1
ratdata.sides(i) = 'r';
else
ratdata.sides(i) = 'l';
end
end
ratdata.sides = ratdata.sides';
ratdata.rewards = ratdata.rewards';
ratdata.left_prob1 = ratdata.left_prob1';
ratdata.right_prob1 = ratdata.right_prob1';
%bandits_glm(ratdata);

ratdata.dates = ratdata.dates(trials_before_switch_to_sides3+1:end);
ratdata.sides = ratdata.sides(trials_before_switch_to_sides3+1:end);
ratdata.rewards = ratdata.rewards(trials_before_switch_to_sides3+1:end);
ratdata.left_prob1 = ratdata.left_prob1(trials_before_switch_to_sides3+1:end);
ratdata.right_prob1 = ratdata.right_prob1(trials_before_switch_to_sides3+1:end);
ratdata.trial_types = ratdata.trial_types(trials_before_switch_to_sides3+1:end);
ratdata.nTrials = ratdata.nTrials - trials_before_switch_to_sides3;
bandits_glm(ratdata);

mousedata = ratdata;
mousedata.trial_number_session = trial_number_session;
mousedata.dates_allSession = dates_allSession;
save(['mousedata_' mouseName '.mat'], 'mousedata')

% fit_mouse = fit_model(mousedata, 'td_bias_habits')
% fit_tdbias = fit_model(mousedata, 'td_bias');
% 
% 
% nChoiceTrials = sum(mousedata.trial_types == 'f');
% normbic = @(ll, n_params) exp((ll - 0.5 * n_params * log(sum(nChoiceTrials))) / nChoiceTrials);
% 
% normbic_tdbias = normbic(fit_tdbias.log_likelihood, 4)
% normbic_tdbias_habits = normbic(fit_mouse.log_likelihood, 5)
% nombic_mouse = normbic(accur,npar)
% 
% 
% bandits_glm(mousedata, 5);
% modeldata = generate_simulated_data('td_bias_habits', fit_mouse.fit_params, mousedata)
% bandits_glm(modeldata, 5)
% 
% 
%  figure; scatter(ones([1,6]), normbic_tdbias_tdbias); hold on; scatter(2*ones([1,6]), normbic_tdbias_habits_habits);...
% hold on; scatter(3*ones([1,6]), normbic_td_ito);
% hold on; bar(1, mean(normbic_tdbias_tdbias), 'FaceColor', [0,0,1]); alpha(.1)
% hold on; bar(2, mean(normbic_tdbias_habits_habits), 'FaceColor', [1,0,0]); alpha(.1)
% hold on; bar(3, mean(normbic_td_ito), 'FaceColor', [ 0.9100 0.4100 0.1700]); alpha(.1)
% ylim([0.5, 1])
% 
% save_ratdata;