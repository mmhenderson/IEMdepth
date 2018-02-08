%% run a bunch of classifier scripts in sequence, goes overnight

subj = {'AI','AP','BB','BC','BD','BJ','BM','BN','BO'};

% subj = {'BB','BC','BD','BJ','BM','BN','BO'};
VOIs={'V1','V2','V3','V4','V3A','V3B','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

nIter=1000;
subMean2=0;
parpool(8);

nVox2Use = [];

IEMdepth_classify_XZ2(subj,VOIs,nVox2Use,nIter,subMean2);

IEMdepth_classify_Z6(subj,VOIs,nVox2Use,nIter,subMean2);

nVox2Use = 50;

IEMdepth_classify_XZ2(subj,VOIs,nVox2Use,nIter,subMean2);

IEMdepth_classify_Z6(subj,VOIs,nVox2Use,nIter,subMean2);
