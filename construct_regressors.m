function regressors = construct_regressors3(data,nBack)
% Constructs regressors are required by bandits_glm and bandits_glm_test
% This version generates regressors for choice, reward, and choice X reward
% interaction

if ~exist('nBack','var')
nBack = 5; % How many rewards into the past to look
end

regressor_ch = zeros(1,nBack);
regressor_cxr = zeros(1,nBack);
regressor_rew = zeros(1,nBack);


nChoices = length(data.sides);
regressors = zeros(nChoices,3*nBack);
for choice_i = 1:nChoices
   
   % Get context, side picked
   side = data.sides(choice_i);
   
   % Fill in the regressors
   regressors(choice_i,:) = [regressor_cxr,regressor_ch,regressor_rew]; 
   
   % update reward history vectors
   reward = data.rewards(choice_i);
   c = 0; % These letter variables will be one iff the choice_i'th letter is their letter. Set them to zero, so that if it's a violation trial and they don't get set, they'll be halfway between their set values of -1 and 1
   x = 0;
   r = 0;

   
   if side == 'l'  
       if reward % Left reward
           c = -1;
           x = -1;
           r = 1;

       else % Left omission
           c = -1;
           x = 1;
           r = -1;

       end
   elseif side == 'r'; 
       if reward % Right Reward
           c = 1;
           x = 1;
           r = 1;

       else % Right omission
           c = 1;
           x = -1;
           r = -1;

       end
   elseif side == 'v';
       % Do nothing - it's a violation trial
   else
       error('Problem with sides')
   end
   % Build the letter regressors
   regressor_ch = [c,regressor_ch(1:(end-1))];
   regressor_cxr = [x,regressor_cxr(1:(end-1))];
   regressor_rew = [r,regressor_rew(1:(end-1))];
   
   
end

end