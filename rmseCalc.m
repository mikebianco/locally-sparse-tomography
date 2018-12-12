function [ rmse ] = rmseCalc(sEst,sTrue,validBound)
% function to calculate rmse for LST paper
% Bianco 4/10/18
%  

vb=validBound(:);
nValid=sum(vb);

se=(sEst(:)-sTrue(:)).^2.*vb;
rmse=sqrt(sum(se)/nValid);

end

