function [ patches ] = getPatches(W1,W2,nib)
% [ patches ] = getPatches(W1,W2,nib)
% code to get all patches from image (with wrap-around)
% made from getBlocks.m code (made 9/24/17)
% 4/2/18 Mike Bianco
% ----------------
% INPUTS:
% W1 = slowness map size (vertical pixels)
% W2 = slowness map size (horizontal pixels)
% nib = number of pixels on block (patch) side
% stride = number of pixels from one block to the next (stride=1 is all possible blocks)
% OUTPUT:
% blocks = indices of patches
%
% W1=20; W2=20;
% nib=10;

inds=1:W1*W2;
im=reshape(inds,W1,W2);

% breaking im into sub-ims
imTop=im(1:nib-1,:);
imLeft=im(:,1:nib-1);
imUL=imTop(1:nib-1,1:nib-1);

% combining into wrap-around indices
imInds=[im,imLeft;imTop,imUL];

patches=im2col(imInds,[nib nib]);

% pixel_sum = zeros(W1,W2);
% for m = 1:length(patches)
%     pixel_sum(patches(:,m)) = pixel_sum(patches(:,m))+ones(nib^2,1);
% end


end

