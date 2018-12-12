function [ D_im] = plotDict2D_2(D,gap,backShade,dims)
% [D_im] = plotDict2D_2(D,gap,backShade,dims)
% function to make 2D dictionary image
% D = dictionary
% nib = number of pixes per block side (square)
% D_im = 2D dictionary image (sorted by variance, contrast normalized)
%
% rev 2, 8/13/17: added bigger gaps bw dictionary entries and adjusted
% background color, also sorting by variance after stretching the grayscale
% values of the individual entries
%
% % D = Q;
% % nargin = 2;
% % gap = 3;

if nargin < 3
    backShade = 0.5;
elseif nargin < 4
    dims = ceil(sqrt(size(D,2)));
end

nib = sqrt(size(D,1));


D_im = backShade*ones(dims*(nib)+(dims+1)*gap); % dims+1 is border

Dcn = [];
% contrast normalization
for mm = 1:size(D,2)
    Dim0 = D(:,mm);
    
    % contrast normalization
    Dim0 = Dim0-min(Dim0);
    Dim0 = Dim0./max(Dim0);
    
    % accounting for a constant block (if there is one)
    if isnan(Dim0(1))
        Dim0 = ones(nib^2,1);
    end
    
    Dcn(:,mm)=Dim0;
end


count00 = 0;
[~,ii]=sort(std(Dcn),'descend');
ii=1:length(Dcn);
Dprint = Dcn(:,ii);
for m = 1:dims
    for n = 1:dims
        if count00 < size(D,2)
            count00 = count00+1;
            Dim0 = Dprint(:,count00);   
            Dim0 = reshape(Dim0,[nib,nib]);
            
            D_im(m*(nib+gap)-nib+1:m*(nib+gap),...
                n*(nib+gap)-nib+1:n*(nib+gap))= Dim0;
        else
            continue
        end
    end
    if count00 >= size(D,2)
        continue
    end
end

% % figure(999)
% % clf;
% % imagesc(D_im)
% % colormap gray
end

