function anovaRes = randANOVA(combined_dat,stmCnts)

combined_sh = combined_dat(randperm(length(combined_dat)),:); %
emoneu_sh = combined_sh(:,1);
RemFor_sh = combined_sh(:,2);
[anovaRes,~,~] = anovan(stmCnts,{emoneu_sh RemFor_sh},"Model","interaction","Varnames",["emo","mem"],"display","off");
% anovaRes = pa; % Save p-values in order: Emotion, Memory, Interaction