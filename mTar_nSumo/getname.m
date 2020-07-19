global numOFconstspecies allTargets numofSUMOontarget ;

dum1=num2cell(zeros(numOFconstspecies ,1));
for n = 1:numOFconstspecies
    dum1{n} = n; %value assignment to each cell.. dum{1}=1; dum{2}=2...
end

[npresumo,nsumo,nsenp,ne3] = deal(dum1{:});

%%
repeatingSPSnames=allTargets;

for ctrN = 1:max(size(repeatingSPSnames))
    ctr1=0;
    for ctr = n+1:n+numofSUMOontarget(ctrN)+1
        eval([repeatingSPSnames{ctrN} num2str(ctr1) '=ctr;']);
        ctr1= ctr1 + 1;
    end

    if  ctrN < 3
        ctr1=0;
        for ctr2 = ctr+1:ctr+numofSUMOontarget(ctrN)+1
            eval(['sumo' repeatingSPSnames{ctrN} num2str(ctr1) '=ctr2;']);
            ctr1= ctr1 + 1;
        end
        
         n=n+(2*(numofSUMOontarget(ctrN)+1));
    else
         n=n+numofSUMOontarget(ctrN)+1;
    end

   
end

clear ctr ctr1 ctr2 ctrN dum1 repeatingSPSnames n