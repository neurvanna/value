

mat = dir('C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ1full\*.mat'); 
for q = 1:length(mat) 
    signReg{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\sigLocsQ1full\' mat(q).name]);
end

mat = dir('C:\Users\annaL\Documents\PhD\value_coding\allLocsQ1full\*.mat'); 
for q =  1:length(mat) %1:length(mat) 
    allReg{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\allLocsQ1full\' mat(q).name]);
end

allSigRegions = [];
for j = 1:length(signReg)
    for k =1:height(signReg{1, j}.sigRsp)
        allSigRegions = [allSigRegions, signReg{1, j}.sigRsp{k,1}];
    end
end

allSnSRegions = [];
whichRec=[];
for j =  1:length(mat) %1:length(allReg)
    for k =1:height(allReg{1, j}.allLoc)
        allSnSRegions = [allSnSRegions, allReg{1, j}.allLoc{k,1}];
        whichRec(end+1) = j;
    end
end



% brain_groups = {
%     ["RSP"]; % RSP
%     ["FRP", "MOp", "MOs", "SSp", "SSs", "GU", "AUD", "VIS", "ACA", "PL", "ILA", "ORB", "AI", "PTLp", "TEa", "PERI", "ECT", "OLF", "MOB","AOB","AON", "DP", "PIR", "TT", "NLOT", "COA", "PAA", "TR"]; % non-visual cortex except RSP
%     ["CL", "LD", "LGd", "LH", "LP", "MD", "MG", "PO", "RT", "SPF", "TH", "VAL", "VPL", "VPM", "Eth", "HY", "LGv", "MH", "PoT", "SubG"]; % thalamus, hypothalamus
%     ["CA1", "CA2", "CA3", "DG", "SUB", "POST", "FC", "IG", "HPF"]; % hippocampal
%     ["APN", "IC", "MB", "NB", "SAG", "PBG", "PAG", "NPC", "SNr", "VTA", "RR", "MRN", "PRT", "MPT", "NOT", "OP", "PPT", "CUN", "RPF"]; % midbrain except SC
%     ["SC"]; % SC
%     ["ACB", "CP", "FS", "GPe", "LSc", "LSr", "MS", "OT", "SI", "STR", "BST", "NDB", "PAL"];% basal ganglia
%     ["LA", "BMA", "EP", "MEA", "CTXsp"]; % cortical subplate
%     ["ZI", "FF"] %ZI
%     ["IIn", "act", "alv", "chpl", "cing", "csc", "dhc", "fi", "fiber tracts", "och", "opt", "root"] %fiber tracts
%     }



        brain_groups = {
            ["FRP"];
            ["MOp"];
            ["MOs"];
            ["SSp", "SSs"]; 
            ["GU"];
            ["AUD"]; 
            ["VIS"]; 
            ["ACA"];
            ["PL"];
            ["ILA"];
            ["ORB"];
            ["AI"];
            ["PTLp"];
            ["TEa"];
            ["PERI"]; 
            ["ECT"];
            ["OLF", "MOB","AOB","AON", "DP", "PIR", "TT", "NLOT", "COA", "PAA", "TR"]; % non-visual cortex except RSP
            
            }


% 
% 
% 
% brain_groups = {
%     ["ACB"]; 
%     ["CP"]; 
%     ["FS"]; 
%     ["GPe"];
%     ["LSc"]; 
%     ["LSr"]; 
%     ["MS"]; 
%     ["OT"]; 
%     ["SI"]; 
%     ["STR"]; 
%     ["BST"]; 
%     ["NDB"];
%     ["PAL"];% basal ganglia
%     }
% 
% % 
% brain_groups = {
%     ["CA1"];
%     ["CA2"]; 
%     ["CA3"]; 
%     ["DG"];
%     ["SUB"]; 
%     ["POST"];
%     ["FC"];
%     ["IG"];
%     ["HPF"]; %hipp
% }

% brain_groups = {
%     ["BLA"];
%     ["BMA"];
%     ["EP"];
%     ["MEA"]; 
%     ["CTXsp"]; 
%     }

% brain_groups = {
%     ["APN"]; 
%     ["IC"];
%     ["MB"]; 
%     ["NB"]; 
%     ["SAG"]; 
%     ["PBG"]; 
%     ["PAG"]; 
%     ["NPC"]; 
%     ["SNr"]; 
%     ["VTA"]; 
%     ["RR"];
%     ["MRN"]; 
%     ["PRT"]; 
%     ["MPT"]; 
%     ["NOT"];
%     ["OP"];
%     ["PPT"]; 
%     ["CUN"]; 
%     ["RPF"];
% }
% 
% brain_groups = {
%     ["CL"];
%     ["LD"] ;
%     ["LGd"] ;
%     ["LH"];
%     ["LP"] ;
%     ["MD"] ;
%     ["MG"] ;
%     ["PO"] ;
%     ["RT"] ;
%     ["SPF"] ;
%     ["TH"];
%     ["VAL"] ;
%     ["VPL"] ;
%     ["VPM"] ;
%     ["Eth"] ;
%     ["HY"];
%     ["LGv"] ;
%     ["MH"];
%     ["PoT"] ;
%     ["SubG"];
% 
%     }




whichRecNeuron = [];
bg_counter2 = zeros(1, length(brain_groups));
otherRegions2 = {};
for i = 1:length(allSnSRegions)
    flag = 0;
    for k =8
        for j = 1:length(brain_groups{k})
            if contains(allSnSRegions{i},brain_groups{k}(j))
                bg_counter2(k) = bg_counter2(k)+1;
                flag = 1;
                whichRecNeuron(end+1) = whichRec(i);
            end
        end
    end
    if flag == 0
        otherRegions2{end+1} = allSnSRegions{i};
    end
end
mat(unique(whichRecNeuron)).name




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

whichRecNeuron = [];
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

% figure;
% X = categorical({'RSP', 'other ctx', 'thal+hyp', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate', 'ZI', 'fiber tracts', 'other'});
% X = reordercats(X,{'RSP', 'other ctx', 'thal+hyp', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate', 'ZI', 'fiber tracts', 'other'});
% 
X = categorical({'FRP','MOp','MOs','SSp', 'GU','AUD', 'VIS','ACA','PL','ILA','ORB','AI','PTLp','TEa','PERI', 'ECT','OLF', 'other'});
X = reordercats(X, {'FRP','MOp','MOs','SSp', 'GU','AUD', 'VIS','ACA','PL','ILA','ORB','AI','PTLp','TEa','PERI', 'ECT','OLF', 'other'});

% X = categorical({'CP','ACB','FS', 'OT','LS', 'SF','SH','AAA','BA','CEA','IA','MEA','GPe','GPi', 'SI','MA', 'MSC', 'TRS', 'BST', 'BAC', 'other'});
% X = reordercats(X, {'CP','ACB','FS', 'OT','LS', 'SF','SH','AAA','BA','CEA','IA','MEA','GPe','GPi', 'SI','MA', 'MSC', 'TRS', 'BST', 'BAC', 'other'});

% 
% X = categorical({'ACB', 'CP', 'FS', 'GPe', 'LSc', 'LSr', 'MS', 'OT', 'SI', 'STR', 'BST', 'NDB', 'PAL', 'other'});
% X = reordercats(X,{'ACB', 'CP', 'FS', 'GPe', 'LSc', 'LSr', 'MS', 'OT', 'SI', 'STR', 'BST', 'NDB', 'PAL', 'other'});
% 
% X = categorical({'CA1', 'CA2', 'CA3', 'DG', 'SUB', 'POST', 'FC', 'IG', 'HPF', 'other'}); % hippocampal})
% X = reordercats(X,{'CA1', 'CA2', 'CA3', 'DG', 'SUB', 'POST', 'FC', 'IG', 'HPF', 'other'});


% X = categorical({'SC'}); 
%  X = reordercats(X,{'SC'}); 


% X = categorical({'LA', 'BMA', 'EP', 'MEA', 'CTXsp', 'other'});
% X = reordercats(X,{'LA', 'BMA', 'EP', 'MEA', 'CTXsp', 'other'});

% X = categorical({'APN', 'IC', 'MB', 'NB', 'SAG', 'PBG', 'PAG', 'NPC', 'SNr', 'VTA', 'RR', 'MRN', 'PRT', 'MPT', 'NOT', 'OP', 'PPT', 'CUN', 'RPF', 'other'});
% X = reordercats(X,{'APN', 'IC', 'MB', 'NB', 'SAG', 'PBG', 'PAG', 'NPC', 'SNr', 'VTA', 'RR', 'MRN', 'PRT', 'MPT', 'NOT', 'OP', 'PPT', 'CUN', 'RPF', 'other'});
% 
% X = categorical({'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'RT', 'SPF', 'TH', 'VAL', 'VPL', 'VPM', 'Eth', 'HY', 'LGv', 'MH', 'PoT', 'SubG', 'other'});
% X = reordercats(X, {'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'RT', 'SPF', 'TH', 'VAL', 'VPL', 'VPM', 'Eth', 'HY', 'LGv', 'MH', 'PoT', 'SubG', 'other'});

Y = [bg_counter, numel(otherRegions)];
Y2 = [bg_counter2, numel(otherRegions2)];
bar(X,Y2)
hold on
bar(X,Y)
title('number of cells in each region')
legend('total', 'significantly coding Q1 at the time -100ms:0ms before movement')

% X = categorical({'RSP', 'other ctx', 'thal+hyp', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate', 'ZI', 'other'});
% X = reordercats(X,{'RSP', 'other ctx', 'thal+hyp', 'hipp', 'midbrain except SC', 'SC', 'basal ganglia', 'cortical subplate', 'ZI', 'other'});

X = categorical({'FRP','MOp','MOs','SSp', 'GU','AUD', 'VIS','ACA','PL','ILA','ORB','AI','PTLp','TEa','PERI', 'ECT','OLF'});
X = reordercats(X, {'FRP','MOp','MOs','SSp', 'GU','AUD', 'VIS','ACA','PL','ILA','ORB','AI','PTLp','TEa','PERI', 'ECT','OLF'});
% 
% % X = categorical({'CP','ACB','FS', 'OT','LS', 'SF','SH','AAA','BA','CEA','IA','MEA','GPe','GPi', 'SI','MA', 'MSC', 'TRS', 'BST', 'BAC'});
% X = reordercats(X, {'CP','ACB','FS', 'OT','LS', 'SF','SH','AAA','BA','CEA','IA','MEA','GPe','GPi', 'SI','MA', 'MSC', 'TRS', 'BST', 'BAC'});
% % 
% X = categorical({'ACB', 'CP', 'FS', 'GPe', 'LSc', 'LSr', 'MS', 'OT', 'SI', 'STR', 'BST', 'NDB', 'PAL'});
% X = reordercats(X,{'ACB', 'CP', 'FS', 'GPe', 'LSc', 'LSr', 'MS', 'OT', 'SI', 'STR', 'BST', 'NDB', 'PAL'});

% X = categorical({'CA1', 'CA2', 'CA3', 'DG', 'SUB', 'POST', 'FC', 'IG', 'HPF'}); % hippocampal})
% X = reordercats(X,{'CA1', 'CA2', 'CA3', 'DG', 'SUB', 'POST', 'FC', 'IG', 'HPF'});

% X = categorical({'SC'}); 
%  X = reordercats(X,{'SC'}); 

% 
% X = categorical({'LA', 'BMA', 'EP', 'MEA', 'CTXsp'});
% X = reordercats(X,{'LA', 'BMA', 'EP', 'MEA', 'CTXsp'});

% X = categorical({'APN', 'IC', 'MB', 'NB', 'SAG', 'PBG', 'PAG', 'NPC', 'SNr', 'VTA', 'RR', 'MRN', 'PRT', 'MPT', 'NOT', 'OP', 'PPT', 'CUN', 'RPF'});
% X = reordercats(X,{'APN', 'IC', 'MB', 'NB', 'SAG', 'PBG', 'PAG', 'NPC', 'SNr', 'VTA', 'RR', 'MRN', 'PRT', 'MPT', 'NOT', 'OP', 'PPT', 'CUN', 'RPF'});
% 
% 
% X = categorical({'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'RT', 'SPF', 'TH', 'VAL', 'VPL', 'VPM', 'Eth', 'HY', 'LGv', 'MH', 'PoT', 'SubG'});
% X = reordercats(X, {'CL', 'LD', 'LGd', 'LH', 'LP', 'MD', 'MG', 'PO', 'RT', 'SPF', 'TH', 'VAL', 'VPL', 'VPM', 'Eth', 'HY', 'LGv', 'MH', 'PoT', 'SubG'});


Y = bg_counter./bg_counter2;
figure; 
bar(X,Y)
hold on; plot(xlim, [0.05, 0.05], 'r')
title('percentages significant in each region')
