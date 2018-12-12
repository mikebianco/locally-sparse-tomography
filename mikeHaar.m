function [overcomplete_haar ] = mikeHaar(sigDim,maxOrder)
% [overcomplete_haar ] = mikeHaar(sigDim,maxOrder)
% overcomplete haar dictionary
% 8/16/17
% INPUTS:
% sigDim = dimension of signal (must be power of 2)
% maxOrder = max level of haar wavelet
% OUTPUT:
% overcomplete_haar

% checking sigDim
if log2(sigDim) ~= floor(log2(sigDim))
    disp('signal dimension not power of 2!')
    return;
end
x = linspace(0,1,sigDim);

dict = [];
dict(:,1)=ones(length(x),1);

count = 1;
for j = 0:maxOrder
    arg=(2^j*x);
    if j == 0
        lim = length(x)/2;
    else
        lim = length(x);
    end
    for l = 0:lim-1
        shift= l;
        arg_sh = [arg(end-shift+1:end),arg(1:end-shift)];
        psi = double(and(arg_sh>=0,arg_sh<0.5)-and(arg_sh>0.5,arg_sh<=1));
        count = count+1;
        dict(:,count)=psi';
    end
end

% normalizing columns
dict = dict*sqrt(diag(1./sum(abs(dict))));

overcomplete_haar = kron(dict,dict);


end

