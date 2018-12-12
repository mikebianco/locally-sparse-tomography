function [sTV] = TV_tomo(inputs)
%============================================================
% [ss,sg, normError,D] = TV_tomo(inputs)
% 
% Implementation of total variation (TV)-regularized tomography, based on
% MTV regularization from Lin et al. 2015. and Chambolle 2004 (please see
% below for citations).
%
% Inputs:
% inputs.lam1:            regularization param. 1
% inputs.lamTV:           TV regularization parameter
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
% The methods implemented here are the same as described in:
% 1. M.J. Bianco and P. Gerstoft, "Travel time tomography with adaptive 
% dictionaries," IEEE Trans. on Computational Imaging, Vol. 4, No. 4, 2018.
% 
% Implemented here is a TV-regulared tomography method based on MTV
% regularization, developed in:
% 2. Y. Lin and L. Huang, "Quantifying subsurface geophysical properties
% changes using double-difference seismic-waveform inversion with a
% modificed total-variation regularization scheme," Geophysical Journal
% Int., vol. 203, no. 2, 2015.
%
% This implementation of MTV is uses TV image denoising code from
% Chambolle 2004:
% 3. A. Chambolle, "An algorithm for total variation minimization and
% applications," J. Mathematical Imaging and Vision, vol. 20, no. 1--2,
% 2004.
% 
% Further implemented is the LSQR sparse least squares algorithm, developed
% in the reference:
% 4. C.C. Paige and M.A. Saunders, "LSQR: And algorithm for sparse linear
% equations and sparse least squares," ACM Trans. on Mathematical Software,
% vol., no. 1, 1982.
% 
% If you find this code useful in your research, please cite the above
% references (1-4).
%
% LST version 1.0
% Michael J. Bianco, December 11th, 2018.
% email: mbianco@ucsd.edu
%============================================================

lam1=inputs.lam1;
lamTV=inputs.lamTV;
A=inputs.tomoMatrix;
s0=inputs.refSlowness;
Tarr=inputs.travelTime;
vb=inputs.validBounds;
normNoise=inputs.normNoise;
sTrue=inputs.sTrue;
gg=inputs.noiseRealiz;


sTV=zeros(size(sTrue(:))); % initializating reference to be 0


normError = norm(Tarr-A*s0); % initial error (of reference solution)
A = sparse(A);

[nrays,npix] = size(A);
for nn = 1:inputs.solIter
    % global slowness
    
    disp(['TV: Realization #',num2str(gg),', Iteration #',num2str(nn)])
    
    dt = Tarr-A*(sTV+s0);
    
    ds = lsqrSOL(nrays,npix,A,dt,lam1,1.0e-6,1.0e-6,1.0e+5,1e3,0);
    sg = ds+sTV;
    
    valb = double(vb(:));
    
    n=length(sTrue);
    f=sg;
    alpha = 0.25; % optimal per chambolle 2004
    NIT=5000;
    GapTol=1e-2;
    verbose=0;%true;
    lamCham=2/lamTV; % per Zhu 2010
    [u, ~,~, ~,~, ~] ...
        = TV_Chambolle(zeros(n),zeros(n),reshape(f,size(sTrue)),lamCham,alpha,NIT,GapTol,verbose);
    
    %   ss=u(:);
    
    sTV=u(:).*valb;
    
    s=sTV+s0;
    
    Tref = A*(s); % new reference travel time for global estimate
    
    
    % calculating errors
    normError_new = norm(Tref-Tarr); % error in travel time
    normError = [normError,normError_new];
    
    if inputs.plots ==true
        [nrow,ncol]=size(sTrue);
        
        Splot = reshape(sg+s0,nrow,ncol);
        uplot = reshape(s,nrow,ncol);
        
        
        figure(inputs.figNo)
        clf;
        subplot(3,2,1)
        imagesc(Splot,inputs.lims)
        colormap(gca,'default')
        title('$\widehat{\mathbf{s}}_\mathrm{g}$','interpreter','latex','fontsize',16)
        
        
        bigTitle='Total variation (TV) inversion example';
        text(0,-20,bigTitle,'fontsize',15,'interpreter','latex')
        
        h= colorbar;
        ylabel(h,'Slowness (s/km)')
        xlabel('Range (km)')
        ylabel('Range (km)')
        
        subplot(3,2,2)
        colormap(gca,'default')
        imagesc(uplot,inputs.lims)
        h= colorbar;
        ylabel(h,'Slowness (s/km)')
        xlabel('Range (km)')
        ylabel('Range (km)')
        title('$\widehat{\mathbf{s}}_\mathrm{TV}$','interpreter','latex','fontsize',16)
        
        % plotting slices
        % horizontal slice
        sliceLoc = 48;
        
        uSlice = uplot(sliceLoc,:);
        sTrue_slice = sTrue(sliceLoc,:);
        
        subplot(3,2,3)
        plot(1:100,sTrue_slice,'k',1:100,uSlice,'r')
        legend('True','Estimated')
        ylabel('Slowness (s/km)')
        xlabel('Range (km)')
        title('$\widehat{\mathbf{s}}_\mathrm{TV}$, horizontal slice',...
            'interpreter','latex','fontsize',16)
        
        % vertical slice
        sliceLoc = 30;
        uSlice = uplot(:,sliceLoc);
        sTrue_slice = sTrue(:,sliceLoc);
        
        subplot(3,2,4)
        plot(1:100,sTrue_slice,'k',1:100,uSlice,'r')
        legend('True','Estimated')
        ylabel('Slowness (s/km)')
        xlabel('Range (km)')
        title('$\widehat{\mathbf{s}}_\mathrm{TV}$, vertical slice',...
            'interpreter','latex','fontsize',16)
        
        
        subplot(3,2,5)
        plot(0:nn,normError,'b-o',0:nn,ones(1,nn+1)*normNoise,'r')
        ylabel('Travel time error norm (s)')
        xlabel('iteration #')
        legend('error','noise')
        title('Travel time error vs. iteration',...
            'interpreter','latex','fontsize',16)
        drawnow
        
    end
    
end
end
