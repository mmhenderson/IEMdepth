function basis = cos_basis_1D(centers,eval_at,size)
% pp(1) is mu, pp(2) is size, pp(3) is size, pow is 7
% eval_at is 1D vectors
% ff is 

basis=nan(length(eval_at),length(centers));

for cc=1:length(centers)
    
    myr = abs(eval_at-centers(cc));
    ff = ((0.5*(1 + cos(myr*pi/size) )).^7) .* (myr<=size)  ;
    
    basis(:,cc)=ff;
    
end

return