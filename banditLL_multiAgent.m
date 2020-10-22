function [LL, agentchoice, ll_all, Q1all, Q2all, habitTraceAll] = banditLL_multiAgent(alphaRL,alphaHabits, betaRL, betaHabits,bias,ratdata)

% Returns the negative log-likelihood of a set of rat data from the
% context-free bandits task, given a set of model parameters

log_eps = log(1e-10); % Smallest probability that an action can have - this becomes important in the case of numerical problems with large betas

% Check that parameters are in bounds
if alphaRL < 0 || alphaRL > 1 ||alphaHabits < 0 || alphaHabits > 1
    LL = NaN;
    return
end
if ~exist('reset','var') % reset flag tells us whether to reset the values at the beginning of each session.  Default is "yes"
    reset = 1;
end

% We will consider forced-choice trials for learning, but not for choice
% Initialize the values of both sides at 0.5
Q = 0.5*ones(1,2);
Q1all = [];
Q2all = [];
habitTraceAll = [];
Qsim = 0.5;
habitTrace = 0;
LL = 0;
nTrials = ratdata.nTrials;
sideChosenbyAgent = [];
ll_all = [];
agentchoice = [];
for trial_i = 1:nTrials
    
    %% Update log-likelihood
    % UNCOMMENT THIS AND LINE 63-65 IF YOU WANT ASYMMETRIC AGENT
    Qeff = [Q(1)*betaRL + habitTrace*betaHabits + bias,...
        Q(2)*betaRL - habitTrace*betaHabits - bias];
    
%     Qeff = [Qsim*betaRL + habitTrace*betaHabits + bias,...
%         -Qsim*betaRL - habitTrace*betaHabits - bias];
    Q1all(end+1) = Q(1);
    Q2all(end+1) = Q(2);
    habitTraceAll(end+1) = habitTrace;
    logActionProbs = Qeff - logsumexp(Qeff,2);
    logActionProbs(logActionProbs==0) = log_eps;
    sideChosen = (ratdata.sides(trial_i) == 'r') + 1; % which side was chosen?
    ll_trial = logActionProbs(sideChosen); % What is the probability that the model would have chosen the same?
    
    ll_all(end+1) = ll_trial;
    agentchoice(end+1) = rand<exp(logActionProbs(2));
    
    if ~isnumeric(ll_trial) || isnan(ll_trial) || ll_trial >= 0 %|| ll_trial==log_eps
        error('Problem with log-likelihood');
    end
    
    LL = LL + ll_trial;
    
    
    %% Do the learning
    % Learning happens on all nonviolation trials
    if ratdata.sides(trial_i) == 'l' || ratdata.sides(trial_i) == 'r'
        % RL system
        sideChosen = double(ratdata.sides(trial_i) == 'r') + 1;
        if sideChosen==1
            sideNotChosen = 2;
        else
            sideNotChosen = 1;
        end
    %    UNCOMMENT THIS AND LINE 31-32 IF YOU WANT ASYMMETRIC AGENT
        PE = ratdata.rewards(trial_i) - Q(sideChosen);
        Q(sideChosen) = Q(sideChosen) + alphaRL*PE;
   %     Q(sideNotChosen) = Q(sideNotChosen) - alphaRL*PE;
%         PE = ratdata.rewards(trial_i) - Qsim;
%         Qsim = Qsim + alphaRL*PE;
   %     Q(-sideChosen+3) = Q(-sideChosen+3) - alphaRL*PE;
        % Habit System
        sideChosenForTrace = (ratdata.sides(trial_i) == 'l') - (ratdata.sides(trial_i) == 'r');  % Get a sidechosen metric that's +1 for left, -1 for right
        PEhabits = sideChosenForTrace - habitTrace;
        habitTrace = habitTrace + alphaHabits*PEhabits;
    end
end


LL = -LL;  % Later we'll use fmincon, so we need *negative* log-likelihood

end