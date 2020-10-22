

clear signifPct
clear totalCellCount
sess = [1:34]
mat = dir('C:\Users\annaL\Documents\PhD\value_coding\signifPct\*.mat'); 
for q = sess%1:length(mat) 
    signifPct{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\signifPct\' mat(q).name]);
end
mat = dir('C:\Users\annaL\Documents\PhD\value_coding\totalCellCount\*.mat'); 
for q = sess%1:length(mat) 
    totalCellCount{q} = load(['C:\Users\annaL\Documents\PhD\value_coding\totalCellCount\' mat(q).name]);
end
sess=[1:34]
figure;hold on;
for i =sess%1:length(mat)
plot(signifPct{1, i}.signifPct  /totalCellCount{1, i}.totalCellCount  )
end
hold on; plot([0,100], [1,0], 'k')
xlabel('percentile of permuted data')
ylabel('percentage of cells better than Xth percentile of randomly permuted data')
title('coding of previous feedback')