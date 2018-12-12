function [ percZero ] = patchSamp(A,blocks)
% [ percZero ] = patchSamp(A,blocks)
% code to quantify sampling of pixels in each block (patch) by rays
% made from defineBlocks.m code
% 9/24/17 Mike Bianco
% ----------------
% INPUTS:
% A = tomography matrix
% blocks = indices of blocks (patches)
% OUTPUT:
% percZero = percentage of pixels in blocks (patches) which aren't sampled
% by rays
%

[nPatch,nblocks] = size(blocks);
Bpass = []; % rays passing thru the pixels in each block (looking for pixels with no sampling)
A = full(A);
for k = 1:nblocks
    Bpass(:,k) = sum(A(:,blocks(:,k)))';
end

Bpass_zero = Bpass==0;
zeroCheck = sum(Bpass_zero);
percZero = zeroCheck/nPatch;

% % % looking at distribution of percentage zero crossings
% % figure;
% % plot(percZero)

end

