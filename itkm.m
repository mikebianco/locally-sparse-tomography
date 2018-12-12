function dico=itkm(data,K,S,maxit,dinit)

% syntax: dico=itkm(data,K,S,maxit,dinit)
%
% Iterative Thresholding& K signed Means 
% dictionary learning algorithm as described in
% 'Local Identification of Overcomplete Dictionaries'
% arXiv:1401.6354
%
% input:
% data... d x N matrix containing the training signals as its columns
% K... number of dictionary atoms/dictionary size - default d
% S... desired/estimated sparsity level of the signals - default 1
% maxit... number of iterations - default 1000
% dinit... initialisation, d x K unit norm column matrix - default random 
%
% output:
% dico... d x K dictionary 
%
% Karin Schnass 24.01.14


%%%% preparations
if(nargin < 1)
    disp('syntax: dico=itkm(data,K,S,maxit,dinit)');
    dico=[];
    return;
end

[d,N]=size(data);

if(nargin < 2)
    K=d;
end

if (N < K+1)
    disp('less training signals than atoms => trivial solution');
    dico=data;
    return;
end

if(nargin < 5) 
    dinit = randn(d,K); 
    scale = sum(dinit.*dinit);
    dinit=dinit*diag(1./sqrt(scale));  	
end

if(nargin < 4)
    maxit = 1000;
end

if(nargin < 3)
    S=1;
end

if size(dinit)~=[d,K]
    disp('initialisation does not match dictionary size - random initialisation used');	
    dinit = randn(d,K); 
    scale = sum(dinit.*dinit);
    dinit=dinit*diag(1./sqrt(scale));
end
	

%%%% algorithm

dold=dinit;

for it=1: maxit
    
    ip=dold'*data;
    absip=abs(ip);
    signip=sign(ip);
    [sortip,I] = sort(absip,1,'descend');
    dnew=zeros(d,K);
    for n=1:N
        dnew(:,I(1:S,n))=dnew(:,I(1:S,n))+ data(:,n)*signip(I(1:S,n),n)';
    end
    scale=sum(dnew.*dnew);
    nonzero=find(scale > 0.001);
    dnew(:,nonzero)=dnew(:,nonzero)*diag(1./sqrt(scale(nonzero)));
    dold(:,nonzero)=dnew(:,nonzero);
end

dico=dold;
