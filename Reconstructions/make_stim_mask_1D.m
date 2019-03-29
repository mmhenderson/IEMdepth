% make_stim_mask.m
%
% Creates a mask of a round stimulus for each stimulus (trial), which is an
% element of stimX, stimY, each is stim_size (either a single value or row
% the same length as stimX and stimY), evaluated at [xx yy], which is
% reshaped into columns, so return is n_stimpts x length(stimX)




function stim_mask = make_stim_mask_1D(stimZ,stim_size,zz)


% TODO: check inputs

if length(stim_size)==1
    stim_size = stim_size * ones(length(stimZ),1);
end

stim_mask = zeros(length(zz),length(stimZ));

for ii = 1:length(stimZ)
    
    rr = abs(zz-stimZ(ii));
    stim_mask(rr <= stim_size(ii),ii) = 1;
    
end



return