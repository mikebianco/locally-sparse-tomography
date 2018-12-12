function [X,error] =OMP_N(A,Y,K)
% [x,error] = OMP_N(A,y,K)
% Finds columns in A that best describe data using OMP
% ****ASSUMES COLUMN NORMALIZED A******
% ---- Inputs
% A = dictionary
% y = observations
% K = sparsity
% ---- Outputs
% x = coefficient vector
% error = l2-squared error for reconstruction
%
% Mike Bianco 3/6/16


% initializing variables
[~,nA] = size(A); % number of columns in A
[~,nY] = size(Y); % number of data vectors

X = zeros(nA,nY);

for n = 1:nY
    y = Y(:,n);
    x = zeros(nA,1);
    ai = [];
    for m = 1:K
        r = y-A*x; % calculating residual vector
        aProj = abs(A'*r); % finding projections
        aProj(ai) = -1; %making sure no index is used twice (impose 'orthogonality')
        [~,I] = max(aProj); % max projection corresponds to optimal column
        ai = [ai,I]; % accumulating indices
        xp = pinv(A(:,ai))*y; % updating coefficients
        x(ai)=xp;
    end
    X(:,n)=x;
end
R = Y-A*X; % calculating residual matrix

error = mean(mean(abs(R)));
