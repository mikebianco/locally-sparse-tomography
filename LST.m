function [ss,sg, normError,D] = LST(inputs)
%============================================================
% [ss,sg, normError,D] = LST(inputs)
% 
% Locally sparse travel time tomography (LST): algorithm which uses sparse
% modeling and dictionary learning to estimate 2D wave slowness map 
% based on wave travel times from an array of sensors.
%
% Inputs:
% inputs.lam1:            regularization param. 1
% inputs.lam2:            regularization param. 2
% in_lst.Tcoeff:          number of non-zero (sparse) coefficients
% inputs.solIter:         number of iterations of lst algorithm
% in_lst.itkmIter:        number of itkm iterations
% in_lst.percZeroThresh:  threshold on the allowable fraction of unsampled
% pixels in patches
% in_lst.rngSeed:         random seed for dictionary initialization
% in_lst.tomoMatrix:      tomography matrix
% in_lst.refSlowness:     reference slowness
% in_lst.travelTime:      travel times
% in_lst.validBounds:     valid boundary for LST inversion
% in_lst.normNoise:       Euclidian norm of noise vector
% in_lst.sTrue:           true slowness
% in_lst.lims:            slowness map image range
% in_lst.dictType:        generic dictionary of dictionary learning,
% choices are 'Haar','DCT', or 'Learned'
% in_lst.nD:              number of dictionary atoms (for learned
% dictionary)
% in_lst.nib:             number of pixels on side of patch
% in_lst.figNo:           figure number
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The method implemented here are the same as described in:
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



A =        inputs.tomoMatrix;
sRef=      inputs.refSlowness;
Tarr =     inputs.travelTime;
percZeroThresh = inputs.percZeroThresh;
vb2 =      inputs.validBounds;
lam1=      inputs.lam1;
lam2=      inputs.lam2;
T =        inputs.Tcoeff;
normNoise = inputs.normNoise;
sTrue = inputs.sTrue;
nib=inputs.nib;
dictType=inputs.dictType;
nD=inputs.nD;

[W1,W2]=size(sTrue);
patches=getPatches(W1,W2,nib); % calculating image patch indices
percZero = patchSamp(A,patches); % percentage ray coverage in patches

% defining (or initializing) dictionary
if strcmp(dictType,'Haar')
    %calculating haar wavelets
    D0=mikeHaar(nib,1);
    dictLearning=false;
elseif strcmp(dictType,'DCT')
    % calculating dct
    Pn=ceil(sqrt(nD));
    DCT=zeros(nib,Pn);
    for k=0:1:Pn-1
        V=cos([0:1:nib-1]'*k*pi/Pn);
        if k>0, V=V-mean(V); end
        DCT(:,k+1)=V/norm(V);
    end
    DCT=kron(DCT,DCT);
    D0 = DCT;
    dictLearning=false;
elseif strcmp(dictType,'Learned')
    rng(inputs.rngSeed)
    Drand = randn(nib^2,nD);
    Drand = Drand*diag(sqrt(1./diag(Drand'*Drand)));
    D0=Drand;
    dictLearning=true;
else
    disp('Invalid dictionary choice!')
    %     return;
end



ss=zeros(length(vb2(:)),1); % initializating sparse slowness to be 0

[npp,npatches] = size(patches); % npp=number of patches per pixel (which is constant with wrap-around)

normError = norm(Tarr-A*(ss+sRef)); % initial error (of reference solution)

[nrays,npix] = size(A);
for nn = 1:inputs.solIter
    % global slowness
    
    dt = Tarr-A*(ss+sRef);
    A = sparse(A);
    ds = lsqrSOL(nrays,npix,A,dt,lam1,1.0e-6,1.0e-6,1.0e+5,1e3,0);
    sg = ds+ss;
    
    % patch slowness
    Y = sg(patches);
    meanY = mean(Y);
    meanYmtx = repmat(meanY,[size(Y,1),1]);
    Yc = Y-meanYmtx; % centering patch slownesses
    
    % dictionary learning or solving with predefined dictionary
    if dictLearning == true
        Yl = Yc(:,percZero<=percZeroThresh); % Ylearning, exluding patches which are undersampled
        D = itkm(Yl,size(D0,2),T,inputs.itkmIter,D0);
        D0 = D;
    else
        D = D0;
    end
    
    disp(['LST: Realization #',num2str(inputs.rngSeed),', Iteration #',num2str(nn),', ',inputs.dictType,' dictionary'])
    
    [X,~]=OMP_N(D,Yc,T);
    ss_b=D*X+meanYmtx; % adding-back mean values of patches
    
    %     ss_b=ss_b.*valBlocks; % not including pixels outside of valid region
    
    % updating reference slowness (averaging patch solutions)
    ss_p_sum = zeros(size(ss));
    for m = 1:npatches
        ss_p_sum(patches(:,m))=ss_p_sum(patches(:,m))+ss_b(:,m);
    end
    
    ss_f = (lam2*sg+ss_p_sum)/(lam2+npp);
    
    valb = double(vb2(:));
    ss=ss_f.*valb;
    
    Tref = A*(ss+sRef); % new reference travel time for global estimate
    
    % calculating errors
    normError_new = norm(Tref-Tarr); % error in travel time
    normError = [normError,normError_new];
    
    if inputs.plots ==true
        [nrow,ncol]=size(sTrue);
        
        sg_plot = reshape(sg+sRef,nrow,ncol);
        ss_plot = reshape(ss+sRef,nrow,ncol);
        
        figure(inputs.figNo)
        clf;
        subplot(3,2,1)
        imagesc(sg_plot,inputs.lims)
        colormap(gca,'default')
        title('$\widehat{\mathbf{s}}_\mathrm{g}$','interpreter','latex','fontsize',16)
        
        
        bigTitle=['LST inversion example (Bianco and Gerstoft 2018, IEEE TCI), '...
            dictType,' dictionary'];
        text(0,-20,bigTitle,'fontsize',15,'interpreter','latex')
        
        h= colorbar;
        ylabel(h,'Slowness (s/km)')
        xlabel('Range (km)')
        ylabel('Range (km)')
        
        subplot(3,2,2)
        colormap(gca,'default')
        imagesc(ss_plot,inputs.lims)
        h= colorbar;
        ylabel(h,'Slowness (s/km)')
        xlabel('Range (km)')
        ylabel('Range (km)')
        title('$\widehat{\mathbf{s}}_\mathrm{s}$','interpreter','latex','fontsize',16)
        
        
        % plotting slices
        % horizontal slice
        sliceLoc = 48;
        
        uSlice = ss_plot(sliceLoc,:);
        sTrue_slice = sTrue(sliceLoc,:);
        
        subplot(3,2,3)
        plot(1:100,sTrue_slice,'k',1:100,uSlice,'r')
        legend('True','Estimated')
        ylabel('Slowness (s/km)')
        xlabel('Range (km)')
        title('$\widehat{\mathbf{s}}_\mathrm{s}$, horizontal slice',...
            'interpreter','latex','fontsize',16)
        
        % vertical slice
        sliceLoc = 30;
        uSlice = ss_plot(:,sliceLoc);
        sTrue_slice = sTrue(:,sliceLoc);
        
        subplot(3,2,4)
        plot(1:100,sTrue_slice,'k',1:100,uSlice,'r')
        legend('True','Estimated')
        ylabel('Slowness (s/km)')
        xlabel('Range (km)')
        view(-90,90)
        title('$\widehat{\mathbf{s}}_\mathrm{s}$, vertical slice',...
            'interpreter','latex','fontsize',16)
        
        
        subplot(3,2,5)
        plot(0:nn,normError,'b-o',0:nn,ones(1,nn+1)*normNoise,'r')
        ylabel('Travel time error norm (s)')
        xlabel('iteration #')
        legend('error','noise')
        %         ylim([0 0.1])
        Dim = plotDict2D_2(D,3,0.45);
        title('Travel time error vs. iteration',...
            'interpreter','latex','fontsize',16)
        
        subplot(3,2,6)
        imagesc(Dim);
        colormap(gca,'gray')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        title([dictType,' dictionary'],...
            'interpreter','latex','fontsize',16)
        drawnow
        
    end
    
end

