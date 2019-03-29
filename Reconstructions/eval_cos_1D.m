function ff = eval_cos_1D(eval_at,pp,varargin)
% pp(1) is mu, pp(2) is size, pow is 7

if nargin == 2
    pow = 7;
else
    pow = varargin{1};
end
myr = abs(eval_at(:)-pp(1));

ff = ((0.5*(1 + cos(myr*pi/pp(2)) )).^pow) .* (myr<=pp(2))  ;
% ff = ((0.5*(1 + cos(myr*pi/pp(2)) )).^pow);

return