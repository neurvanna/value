%% Import some data
% first run analyse_sides_multiple
%then modeling_workflow_likelihoods
%then cpd_qvalue

clear all
%close all

 ratdata = load('mousedata_AL029.mat'); %I saved these data using function analyse_sides_multiple
 ratdata = ratdata.mousedata;


% UNCOMMENT THIS IS YOU WANT TO FIT GLM GENERATED DATA
% modeldata = generate_simulated_data('GLM', [0.05, 1, -0.8, 0.2, 0.8, 0.9, 2, 0.6, 0.4, 0.02], ratdata)
% nBack = 3;
% modeldata.rewards(modeldata.rewards==-1) = 0;
% regressors = construct_regressors(modeldata,nBack);
% choices = (modeldata.sides == 'r');
% [B,dev,stats] = glmfit(regressors,choices,'binomial');
% B
% 
% ratdata = modeldata;    

%% Split into train and test and fit habits agent
crossValLikelihoodAgent = [];
crossValLikelihoodAgentTrain = [];
crossValLikelihoodAgentHab = [];
crossValLikelihoodAgentTrainHab = [];

trnum = [];
trnumTrain = [];
nBack = 1;
likAgentTest = [];
likAgent = [];
for j = 1:10
    randSubset = [round(ratdata.nTrials*((j-1)/10))+1:round(ratdata.nTrials*((j)/10))];
    yTest = sort(randSubset);
    %        yTrain = [round(data.nTrials*0.2):round(data.nTrials*0.8)];
    yTrain = setdiff(1:ratdata.nTrials,yTest );
    
    mousedata_train.nTrials = length(yTrain);
    mousedata_train.sides = ratdata.sides(yTrain);
    mousedata_train.rewards = ratdata.rewards(yTrain);
    mousedata_train.left_prob1 = ratdata.left_prob1(yTrain);
    mousedata_train.right_prob1 = ratdata.right_prob1(yTrain);
    mousedata_train.trial_types = ratdata.trial_types(yTrain);
    
    fit_mouse_habits = fit_model(mousedata_train, 'td_bias_habits', 3)
%    fit_mouse_habitsF = fit_model(mousedata_train, 'td_bias_habitsF', 3)
%     fit_mouse_habits_fog = fit_model(mousedata_train, 'td_bias_habits_fog', 1)
%     %     fit_mouse = fit_model(mousedata_train, 'td_bias', 1)
%     %     fit_mouse_persev = fit_model(mousedata_train, 'td_bias_habits_persev', 1)
%     fit_mouse_ito = fit_model(mousedata_train, 'ito-doya', 1)
%     
    
    mousedata_test.nTrials = length(yTest);
    mousedata_test.sides = ratdata.sides(yTest);
    mousedata_test.rewards = ratdata.rewards(yTest);
    mousedata_test.left_prob1 = ratdata.left_prob1(yTest);
    mousedata_test.right_prob1 = ratdata.right_prob1(yTest);
    mousedata_test.trial_types = ratdata.trial_types(yTest);
    
    trnumTest(j) = mousedata_test.nTrials;
    trnumTrain(j) = mousedata_train.nTrials;
    [crossValLikelihoodAgentTrainHab(j), agentChoices, agentTrialConfidence, Q1all, Q2all, habitTraceAll] = banditLL_multiAgent(fit_mouse_habits.fit_params(1), fit_mouse_habits.fit_params(2),fit_mouse_habits.fit_params(3), fit_mouse_habits.fit_params(4), fit_mouse_habits.fit_params(5), mousedata_train);
    likAgent(j) = sum(agentChoices'==(mousedata_train.sides=='r'))/numel(agentChoices);
    
%     [crossValLikelihoodAgentTrainHabF(j), agentChoicesF, agentTrialConfidenceF] = banditLL_multiAgentF(fit_mouse_habitsF.fit_params(1), fit_mouse_habitsF.fit_params(2),fit_mouse_habitsF.fit_params(3), fit_mouse_habitsF.fit_params(4), fit_mouse_habitsF.fit_params(5), fit_mouse_habitsF.fit_params(6),mousedata_train);
%     likAgentF(j) = sum(agentChoicesF'==(mousedata_train.sides=='r'))/numel(agentChoicesF);
%     

    [crossValLikelihoodAgentHab(j),agentChoicesTest, ~, Q1allTest, Q2allTest, habitTraceAllTest] = banditLL_multiAgent(fit_mouse_habits.fit_params(1), fit_mouse_habits.fit_params(2),fit_mouse_habits.fit_params(3), fit_mouse_habits.fit_params(4), fit_mouse_habits.fit_params(5), mousedata_test);
    likAgentTest(j) = sum(agentChoicesTest'==(mousedata_test.sides=='r'))/numel(agentChoicesTest);
    


    
    for nBack = 0:20
        counterTrain = 0;
        counterTest = 0;
        GLMchoice = [];
        GLMchoiceTest = [];
        likTrain = [];
        likTest = [];
        resTrain = 0;
        resTest = 0;
        if nBack ~= 0
            regressors = construct_regressors(ratdata,nBack);
        else
            regressors = ones(size(ratdata.sides == 'r'));
        end
        choices = (ratdata.sides == 'r');
        [B,dev,stats] = glmfit(regressors(yTrain, :),choices(yTrain),'binomial');
        
        yfitTest = glmval(B,regressors(yTest, :),'logit');
        yfitTrain = glmval(B,regressors(yTrain, :),'logit');
        yRealTrain = choices(yTrain);
        yRealTest = choices(yTest);
        lll = [];
        
        for i = 1:length(yRealTrain)
            GLMchoice(end+1) = rand < yfitTrain(i);
            if yRealTrain(i)==1
                ll_trial = log(yfitTrain(i));
                lll(i) = ll_trial;
                likTrain(end+1) = yfitTrain(i);
            else
                ll_trial = log(1-yfitTrain(i));
                lll(i) = ll_trial;
                likTrain(end+1) = 1-yfitTrain(i);
            end
            resTrain = resTrain+ll_trial;
            counterTrain = counterTrain+1;
        end
        
        for i = 1:length(yRealTest)
            
            GLMchoiceTest(end+1) = rand < yfitTest(i);
            if yRealTest(i)==1
                ll_trial = log(yfitTest(i));
                likTest(end+1) = yfitTest(i);
            else
                ll_trial = log(1-yfitTest(i));
                likTest(end+1) = 1-yfitTest(i);
            end
            resTest = resTest+ll_trial;
            counterTest = counterTest+1;
            
        end
        
        

        
        
        
        
        likelihTest(nBack+1,j) = exp(resTest/counterTest);
        likelihTrain(nBack+1,j) = exp(resTrain/counterTrain);
        
        likGLM(j,nBack+1) = sum(GLMchoice'==(mousedata_train.sides=='r'))/ numel(GLMchoice);
        likGLMTest(j,nBack+1) = sum(GLMchoiceTest'==(mousedata_test.sides=='r'))/ numel(GLMchoiceTest);
        
      
    end
    rAgentrain(j)  =   exp(-crossValLikelihoodAgentTrainHab(j)./trnumTrain(j));
    rAgentrainTest(j)  =   exp(-crossValLikelihoodAgentHab(j)./trnumTest(j));
    Q1allMean(j,yTrain) = Q1all;
    Q2allMean(j,yTrain) = Q2all;   
    habitTraceAllMean(j,yTrain) = habitTraceAll;
    Q1allMean(j,yTest) = Q1allTest;
    Q2allMean(j,yTest) = Q2allTest;
    habitTraceAllMean(j,yTest) = habitTraceAllTest;
end


q1 = mean(Q1allMean);
q2 = mean(Q2allMean);
habits = mean(habitTraceAllMean);

figure; plot([0:20], mean(likelihTest,2))
hold on; plot([0:20], mean(likelihTrain,2))
hold on; plot(xlim, [mean(rAgentrain), mean(rAgentrain)])
hold on; plot(xlim, [mean(rAgentrainTest), mean(rAgentrainTest)])
title('fit to mouse')
legend('GLM test', 'GLM train', 'agent with no forgetting train', 'agent with no forgetting test', 'Location','southeast')
makepretty()


save(['q2_' ratdata.ratname '_all.mat'], 'q2')
save(['q1_' ratdata.ratname '_all.mat'], 'q1')
save(['habits_' ratdata.ratname '_all.mat'], 'habits')