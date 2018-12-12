function [ds] = conventional_tomo(inputs)
%============================================================
% [ds] = conventional_tomo(inputs)
% 
% Implementation of conventional tomography based on Bayesian MAP
% (least-squares) estimator from Rodgers 2000 (please see below for
% citation).
%
% Inputs:
% inputs.eta:             \eta regulariation parameter
% inputs.L:               smoothness length scale L
% inputs.tomoMatrix:      tomography matrix
% inputs.refSlowness:     reference slowness
% inputs.travelTime:      travel times
% inputs.validBounds:     valid boundaries for TV-regularized inversion
% inputs.normNoise:       Euclidian norm of noise vector
% inputs.sTrue:           true slowness
% inputs.noiseRealiz:     noise realization number (for display purposes
% only)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The method implemented here is the same as described in:
% 1. M.J. Bianco and P. Gerstoft, "Travel time tomography with adaptive 
% dictionaries," IEEE Trans. on Computational Imaging, Vol. 4, No. 4, 2018.
% 
% Implemented here is a Bayesian MAP (least-squares) estimator from Rodgers
% 2000:
% 2. C.D. Rodgers, Inverse methods for atmospheric sounding: theory and
% practice, World Scientific 2000.
% 
% If you find this code useful in your research, please cite the above
% references (1-2).
%
% LST version 1.0
% Michael J. Bianco, December 11th, 2018.
% email: mbianco@ucsd.edu
%============================================================

% code to simulate conventional tomography
% Mike Bianco 
% 11/5/17
%
% experimenting with conv. tomo. inversions with noise
%

eta=inputs.eta;
L=inputs.L;
gg=inputs.noiseRealiz;
A =        inputs.tomoMatrix;
sRef=      inputs.refSlowness;
Tarr =     inputs.travelTime;
sTrue=inputs.sTrue;
[W1,W2]=size(sTrue);

[xxc,yyc]=meshgrid(1:W1,1:W2);
npix = W1*W2;


% precalculating slowness priors (inverting pixel covariance)
disp(['Conventional: Realization #',num2str(gg),', Inverting covariance matrix'])
Sig_L=zeros(npix);
for ii=1:npix
    %%Determine the coordinates of the neighboor
    distc= sqrt((xxc-xxc(ii)).^2+(yyc-yyc(ii)).^2);
    distc = distc(:)';
    Sig_L(ii,:)=exp(-distc/L);
end
invSig_L = Sig_L\eye(npix);

%% inverting for slowness
Tref = A*(sRef*ones(npix,1)); % reference travel time

% travel time perturbation
dT=Tarr-Tref;

% inverting for slowness
disp(['Conventional: Realization #',num2str(gg),', Inverting for slowness'])
G = A'*A+eta*invSig_L;
ds = G\A'*dT;

if inputs.plots==true
    figure;
    imagesc(reshape(ds,size(sTrue)).*inputs.validBounds+sRef,inputs.lims)
    xlabel('Range (km)')
    ylabel('Range (km)')
    h=colorbar;
    xlabel(h,'Slowness (s/km)')
    title('Conventional inversion example (Bayesian MAP)','fontsize',16,'interpreter','latex')
end

end
