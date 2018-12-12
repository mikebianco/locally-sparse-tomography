%============================================================
% Locally sparse travel time tomography (LST): algorithm which uses sparse
% modeling and dictionary learning to estimate a 2D wavefield slowness map 
% based on wave travel times from an array of sensors.
%
% This is the 'main' (running) file implementing a synthetic test of the 
% LST travel time tomography algorithm, as well as conventional and total 
% variation (TV)-regularized tomography. The methods implemented here are 
% the same as described in:
% 1. M.J. Bianco and P. Gerstoft, "Travel time tomography with adaptive 
% dictionaries," IEEE Trans. on Computational Imaging, Vol. 4, No. 4, 2018.
% 
% Implemented here is also the iterative thresholding and K-means (ITKM)
% algorithm which is developed in:
% 2. K. Schnass, "Local identification of overcomplete dictionaries," J.
% Machine Learning Research, Vol. 16, 2015.
%
% Further implemented is the LSQR sparse least squares algorithm and the 
% orthogonal matching pursuit (OMP) algorithm which are developed,
% respectively, in the following references:
% 3. C.C. Paige and M.A. Saunders, "LSQR: And algorithm for sparse linear
% equations and sparse least squares," ACM Trans. on Mathematical Software,
% vol., no. 1, 1982.
% 4. Y.C. Pati, R. Rezaiifar, and P.S. Krishnaprasad, "Orthogonal matching
% pursuit: Recursive function approximation with applications to wavelet
% decomposition," in Proc. IEEE 27th Annual Asilomar Conf. on Signals,
% Systems and Computers, 1993.
% 
% If you find this code useful in your research, please cite the above
% references (1-4).
%
% LST version 1.0
% Michael J. Bianco, December 11th, 2018.
% email: mbianco@ucsd.edu
%============================================================


clc
clear all
% close all
set(0,'DefaultAxesFontSize',14)
set(0,'defaultfigurecolor',[1 1 1])

map='ch';
sTrue=slownessMap3(map); % choosing slowness map, 'ch'=checkerboard, 'sd'=smooth-discontinuous
[W1,W2] = size(sTrue);

stations_rays % code to configure sensor array, and to calculate travel times assuming straight-rays

% change noise fraction from zero to simulate noise realizations
noiseFrac = 0; %0.02; % noise STD as fraction of mean value of travel time (0=noise free case)
if noiseFrac == 0
    nRuns=1;
else
    nRuns=10; % number of noise realizations (calls all codes nRuns times)
end

outputs={}; % cell array for results

for nn=1:nRuns
rngSeed=nn;
% adding noise to travel time observations
stdNoise = mean(Tarr)*noiseFrac;
rng(rngSeed) % fixing random seed
noise = stdNoise*randn(length(Tarr),1);
Tarr_n  = Tarr+noise;

% estimating referense slowness from travel time observations
Asum = sum(A,2);
invAsum = pinv(Asum);
sRef = invAsum*Tarr_n;

%% Running LST
% INPUTS
in_lst=[];
in_lst.lam2 = 0;             % reg. param 2

if noiseFrac==0
    in_lst.lam1 = 0;           % reg param 1
    if strcmp(map,'ch')
        Tcoeff_dlearn=1;       % number of non-zero (sparse) coefficients
        Tcoeff_nolearn=5;
    elseif strcmp(map,'sd')
        Tcoeff_dlearn=1;
        Tcoeff_nolearn=2;
    end
    
    in_lst.plots = true;       % option to plot results
else
    
    if strcmp(map,'ch')
        in_lst.lam1 = 2;          
        Tcoeff_dlearn=2;
        Tcoeff_nolearn=5;
    elseif strcmp(map,'sd')
        in_lst.lam1 = 10;           
        Tcoeff_dlearn=2;
        Tcoeff_nolearn=2;
    end
    
    in_lst.plots = false;
end

in_lst.solIter = 2;        % number of iterations of LST algorithm (set to 100 in Bianco and Gerstoft 2018 LST)
in_lst.itkmIter = 50;        % # itkm iterations
in_lst.percZeroThresh = 0.1; % threshold on the allowable percentage of unsampled pixels in patch
in_lst.noiseRealiz = 1;

in_lst.rngSeed=rngSeed; % random seed for dictionary initialization
in_lst.tomoMatrix = A;
in_lst.refSlowness = sRef;
in_lst.travelTime = Tarr_n;
in_lst.validBounds=vb2;
in_lst.normNoise = norm(noise);
in_lst.sTrue = sTrue;
in_lst.lims=[0.3 0.5];

% LST with dictionary learning
in_lst.dictType = 'Learned';     % choices are 'haar','dct','learned'
in_lst.nD = 150;             % number of dictionary atoms
in_lst.Tcoeff = Tcoeff_dlearn; % number of non-zero coefficients
in_lst.nib = 10;              % sqrt(n), number of pixels on side of block (patch)
in_lst.figNo=2; % figure number
[ss_dlearn,~,~,D] = LST(in_lst);
outputs{1,1}='LST Dict. Learn.';
ss_dlearn_save(:,nn)=ss_dlearn;
outputs{1,2}=ss_dlearn_save;
outputs{1,3}=vb2(:);


% LST with Haar wavelets
in_lst.dictType = 'Haar';     % choices are 'haar','dct','learned'
in_lst.nD = 169;             % number of dictionary atoms
in_lst.nib = 8;              % sqrt(n), number of pixels on side of block (patch)
in_lst.Tcoeff = Tcoeff_nolearn; % number of non-zero coefficients
in_lst.figNo=3; % figure number
[ss_haar,~,~,D] = LST(in_lst);
outputs{2,1}='LST Haar';
ss_haar_save(:,nn)=ss_haar;
outputs{2,2}=ss_haar_save;
outputs{2,3}=vb2(:);


% LST with DCT
in_lst.dictType = 'DCT';     % choices are 'haar','dct','learned'
in_lst.nD = 169;             % number of dictionary atoms
in_lst.nib = 8;              % sqrt(n), number of pixels on side of block (patch)
in_lst.Tcoeff = Tcoeff_nolearn; % number of non-zero coefficients
in_lst.figNo=4; % figure number
[ss_dct,~,~,D] = LST(in_lst);
outputs{3,1}='LST DCT';
ss_dct_save(:,nn)=ss_dct;
outputs{3,2}=ss_dct_save;
outputs{3,3}=vb2(:);


%% Running TV-regularized tomography
% INPUTS
in_tv=[];
if noiseFrac==0
    in_tv.lam1=0;
    in_tv.lamTV=0.01;
    in_tv.plots=true;
else
    in_tv.lam1=5;
    in_tv.lamTV=0.02;
    in_tv.plots=false;
end

in_tv.solIter = 2;        % number of iterations of TV-tomo algorithm
in_tv.tomoMatrix = in_lst.tomoMatrix;
in_tv.refSlowness = in_lst.refSlowness;
in_tv.travelTime = in_lst.travelTime;
in_tv.validBounds=vb;
in_tv.normNoise = in_lst.normNoise;
in_tv.sTrue = in_lst.sTrue;
in_tv.lims=in_lst.lims;
in_tv.noiseRealiz=in_lst.noiseRealiz;

% Calling TV code
in_tv.figNo=5; % figure number
s_TV=TV_tomo(in_tv);
s_TV_save(:,nn)=s_TV;
outputs{4,1}='Total variation (TV)';
outputs{4,2}=s_TV_save;
outputs{4,3}=vb(:);




%% Running conventional tomography code
%  (damped least squares with non-diagonal pixel covariance)

% INPUTS
in_conv=[];

if noiseFrac==0
    in_conv.eta = 0.1;       % conventional \eta regularization parameter
    in_conv.L=10;            % smoothness length scale
    in_conv.plots = true;
else
    in_conv.eta = 10;       
    in_conv.L=6;
    in_conv.plots = false;
end

in_conv.tomoMatrix = in_lst.tomoMatrix;
in_conv.refSlowness = in_lst.refSlowness;
in_conv.travelTime = in_lst.travelTime;
in_conv.sTrue = in_lst.sTrue;
in_conv.lims=in_lst.lims;
in_conv.validBounds=in_lst.validBounds;
in_conv.noiseRealiz=in_lst.noiseRealiz;


% Calling conventional tomography code
s_conv=conventional_tomo(in_conv);
s_conv_save(:,nn)=s_conv;
outputs{5,1}='Conventional';
outputs{5,2}=s_conv_save;
outputs{5,3}=vb(:);


end

%% plotting results for comparison
nCases=size(outputs,1);

figure;
                        %  edge gap   vertical horizontal
ha=tight_subplot(nCases,3,[0.04 0.09],[.1 .1],[.12 .05]);


sliceLoc=48; % horizontal slice location
sliceTrue=sTrue(sliceLoc,:);
x=1:length(sliceTrue);

for nn=1:nCases
    
    ds=outputs{nn,2};
    valb=outputs{nn,3};
    nds=size(ds,2);
    sEst=ds.*repmat(valb,[1,nds])+sRef;

    sEst=reshape(sEst,[size(sTrue),nds]);
    sMean=mean(sEst,3);
    sStd=std(sEst,0,3);
    
    sliceEst=sEst(sliceLoc,:,1); % arbitrarily choosing the first noise realization
    sliceMean=sMean(sliceLoc,:);
    sliceStd=sStd(sliceLoc,:);
    rmse=rmseCalc(sMean,sTrue,vb)*1000;    
    
    axes(ha(nn*3-2))
    a = ha(1);
    pos = get(a,'Position');
    imagesc(sMean,[0.3 0.5])
    if nn==1
        h=colorbar('northoutside');
        set(a,'Position',pos) % fixing image size from colorbar addition
        xlabel(h,'Slowness (s/km)')
    end
    xlabel('Range (km)')
    ylabel('Range (km)')
    set(gca,'Xtick',[1,20:20:100])
    set(gca,'Ytick',[20:20:100])
    hold on
    plot(x,sliceLoc*ones(size(x)),'k')
    hold off
    text(2,7,[outputs{nn,1},', 2D'],'fontsize',14)
    
    % slice of slowness
    axes(ha(nn*3-1))
    plot(x,sliceTrue,'k',x,sliceMean,'r')
    hold on
    y1 = sliceStd+sliceMean;
    y2 = sliceMean-sliceStd;
    X=[x,fliplr(x)];               
    Y=[y1,fliplr(y2)];
    h=fill(X,Y,'b');
    set(h,'EdgeColor','None')
    alpha(0.3)
    hold off
    ylim([0.25 0.55])
    ylabel('Slowness (s/km)')
    xlabel('Range (km)')
    text(2,8,outputs{nn,1},'fontsize',14)
    set(gca,'Xtick',[1,20:20:100])
    set(gca,'Ytick',[.3:.1:.5])
    text(2,.53,[outputs{nn,1},', 1D'],'fontsize',14)
    text(5,0.27,['RMSE=',num2str(rmse,'%.2f'),' ms/km'],'fontsize',14)
    if nn==1
        legend({'True','Est. mean','Est. 1STD'},'Position',[0.53 0.905 0.1 0.05])
    end
    
    % slice of slowness
    axes(ha(nn*3))
    plot(x,sliceTrue,'k',x,sliceEst,'r')
    ylim([0.25 0.55])
    ylabel('Slowness (s/km)')
    xlabel('Range (km)')
    text(2,8,outputs{nn,1},'fontsize',14)
    set(gca,'Xtick',[1,20:20:100])
    set(gca,'Ytick',[.3:.1:.5])
    text(2,.53,[outputs{nn,1},', 1D'],'fontsize',14)
    text(5,0.27,['RMSE=',num2str(rmse,'%.2f'),' ms/km'],'fontsize',14)
    if nn==1
        legend({'True','Est. one realization'},'Position',[0.805 0.905 0.1 0.05])
    end
    

    
end

% adding overall title
bigTitle=['Comparison of LST, TV, and conventional tomography, ',...
    '$\sigma_\epsilon=$ ',num2str(noiseFrac),'$\bar{t}$'];
axes(ha(1))
text(-20,-60,bigTitle,'fontsize',16,'interpreter','latex')

