function [acc,dprime] = IEMdepth_classifier(dat,reallabs,runlabs,classstr,voxelindsuse,subMean2)
% semi-generic function for classifying voxel activation patterns, use a
% leave-one-run-out cross validation scheme.

% dat is [trials x voxels]
% reallabs, runlabs are [trials x 1]
% runlabs stores the label of run for each trial (this is used to
% crossvalidate, leaving one run out and predicting based on the other
% runs)
% voxelindsuse is [nCrossVal x voxels], specify which voxels to use on each
% cross-validation fold.
% subMean2 is whether to subtract mean acros the voxel dimension
% classStr is which classifier to use

% output is the mean accuracy and d' over all crossvalidations

% MMH 9/6/17 
%%

% store a vector of the model predictions
predlabs=nan(length(runlabs),1);

unruns=unique(runlabs);
nruns=length(unruns);

if nargin<6 || isempty(voxelindsuse)
    voxelindsuse = boolean(ones(nruns,size(dat,2)));
end

for cv=1:nruns
    
    testinds=runlabs==unruns(cv);
    trninds=runlabs~=unruns(cv);
    
    %see if there's an unbalanced set
    un=unique(reallabs(trninds,:)); 
    numeach=zeros(length(un),1);
    for ii=1:length(un)
        numeach(ii)=sum(reallabs(trninds,:)==un(ii));
    end
    
    if any(numeach~=min(numeach))        
        error('training set not balanced')
    end

    %set is balanced.
    useinds=trninds;

    datuse = dat(:,voxelindsuse(cv,:));

    if subMean2
        datuse = datuse - repmat(mean(datuse,2),1,size(datuse,2));
    end

    %use it to predict on test set
    if strcmp(classstr,'svmtrain')   
        % note this should be calling a function in the libSVM package, not
        % a matlab built-in.
        obj=svmtrain(reallabs(useinds,:),datuse(useinds,:),'-t 0 -q');         
        predlabs(testinds,:)=svmpredict(reallabs(testinds),datuse(testinds,:),obj);
    elseif strcmp(classstr,'fitcdiscr')              
        obj=fitcdiscr(datuse(useinds,:),reallabs(useinds,:));
        predlabs(testinds,:)=predict(obj,datuse(testinds,:));
    end
end

acc=mean(predlabs==reallabs);
dprime = get_dprime(predlabs,reallabs,unique(reallabs));

end