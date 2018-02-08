function [acc,dprime] = IEMdepth_classifier(dat,reallabs,runlabs,classstr,nIter,voxelindsuse,subMean2)
% function for classifying voxel activation patterns 

% train and test are [trials x voxels]
% trainlabs, testlabs, runlabs are [trials x 1]
% runlabs stores the label of run for each trial (this is used to
% crossvalidate, leaving one run out and predicting based on the other
% runs)

% output is the mean accuracy over all crossvalidations

% MMH 9/6/17 updated
%%

% store a vector of the model predictions
predlabs=nan(length(runlabs),1);

unruns=unique(runlabs);
nruns=length(unruns);

if nargin<6 || isempty(voxelindsuse);
    voxelindsuse = boolean(ones(nruns,size(dat,2)));
end


% acceachrun=nan(nruns,1);
% ntesttrialseachrun=nan(nruns,1);

for cv=1:nruns

%     thesevoxelindsuse=voxelindsuse(cv,:);
    
    testinds=runlabs==unruns(cv);
    trninds=runlabs~=unruns(cv);
    
%     ntesttrialseachrun(cv)=sum(testinds);
    
    %see if there's an unbalanced set
    un=unique(reallabs(trninds,:)); 
    numeach=zeros(length(un),1);
    for ii=1:length(un)
        numeach(ii)=sum(reallabs(trninds,:)==un(ii));
    end
    
    if any(numeach~=min(numeach))        
        error('training set not balanced')
    end
        
%         
%         fprintf('\nbalancing training set...\n');
% 
%         theseaccs=nan(nIter,1);
%         parfor xx=1:nIter
%                     
%             thesepredlabs=[];
%             %get one possible balanced training set
%             useinds=[];
%             for ii=1:length(un)
%                 theseinds=find(trninds & reallabs==un(ii));
%                 if numeach(ii)>min(numeach)
%                     %this is a larger set, take a random set of these trials
%                     randinds=randperm(length(theseinds));
%                     useinds=[useinds;theseinds(randinds(1:min(numeach)))];
%                 else
%                     %this is the smallest set
%                     useinds=[useinds;theseinds];
%                 end
%             end
%             
%             %use it to predict on test set (whole test set used every time)
%             if strcmp(classstr,'svmtrain');    
%                 obj=svmtrain(reallabs(useinds,:),dat(useinds,thesevoxelindsuse),'-t 0 -q');         
%                 thesepredlabs=svmpredict(reallabs(testinds),dat(testinds,thesevoxelindsuse),obj);
%             elseif strcmp(classstr,'fitcdiscr');                 
%                 obj=fitcdiscr(dat(useinds,thesevoxelindsuse),reallabs(useinds,:));
%                 thesepredlabs=predict(obj,dat(testinds,thesevoxelindsuse));
%             end
%             
%             theseaccs(xx) = mean(thesepredlabs==reallabs(testinds,:));
%             
%         end
% 
%         acceachrun(cv)=mean(theseaccs);
%         
%         
%     else
        %set already balanced!
        useinds=trninds;
        
        datuse = dat(:,voxelindsuse(cv,:));
        
        if subMean2
            datuse = datuse - repmat(mean(datuse,2),1,size(datuse,2));
        end
        
        %use it to predict on test set
        if strcmp(classstr,'svmtrain');    
            obj=svmtrain(reallabs(useinds,:),datuse(useinds,:),'-t 0 -q');         
            predlabs(testinds,:)=svmpredict(reallabs(testinds),datuse(testinds,:),obj);
        elseif strcmp(classstr,'fitcdiscr');                 
            obj=fitcdiscr(datuse(useinds,:),reallabs(useinds,:));
            predlabs(testinds,:)=predict(obj,datuse(testinds,:));
%         elseif strcmp(classstr, 'svmtrain_gpu')
%             predlabs(testinds,:) = run_SVM_gpu(dat(useinds,voxelindsuse(cv,:)),dat(testinds,voxelindsuse(cv,:)), reallabs(useinds,:), reallabs(testinds,:),'-t 0 q');
        end

%         acceachrun(cv) = mean(predlabs(testinds,:)==reallabs(testinds,:));

    
end

%get the mean accuracy over crossvals (weight according to ntrials in test
%set, which is usually the same for all crossvals)
% acc=sum(acceachrun.*ntesttrialseachrun)/sum(ntesttrialseachrun);

acc=mean(predlabs==reallabs);
dprime = get_dprime(predlabs,reallabs,unique(reallabs));

end