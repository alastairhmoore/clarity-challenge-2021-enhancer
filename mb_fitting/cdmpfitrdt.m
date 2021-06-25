function [rmft,sdft,doft,erft]=cdmpfitrdt(ry,c,vq,dift,cft)
% Fit a multicomponent model to a microphone-array covariance matrix with a priori direction
%  Inputs: ry(nm^2,nt,nf)   Smoothed real-vectorized observed covariance matrices
%          c(nm^2,nc,nf)    Real-vectorized covariance matrices for constant fixed components
%          vq(nm^2,nq,nf)   Real-vectorized covariance matrices for directional component [use [] or omit if unwanted]
%          dift(nt,nf)      Specify the direction index a priori [use [] or omit if unwanted]
%          cft(nm^2,nd,nt,nf)  Real-vectorized covariance matrices for time-varying fixed component [use [] or omit if unwanted]
%
% Outputs: rmft(nm^2,nt,nf) Model fit to ry
%          sdft(nb,nt,nf)   Component weights: directional component, time-varying and then constant fixed components
%          doft(nt,nf)      Index of directional component in vq array (will equal dift)
%          erft(nt,nf)      Squared Frobenius norm of fitting error
%
%    where nm = number of microphone channels
%          nt = number of time frames
%          nf = number of frequency bins
%          nc = number of constant fixed components
%          nd = number of time-varying fixed components
%          nv = number of directional components (0 or 1)
%          nb = nv+nd+nc number of components in the model
%          nq = number of search directions
%
% Versions:
%
% 2020_0619 Initial version.
% 2020_0718 Added cft input
%
%
[nm2,nt,nf]=size(ry);       % extract input dimensions
if nargin<5
    cft=[];
    if nargin<4
        dift=[];
        if nargin<3
            vq=zeros(nm2,0,nf);
        end
    end
end
nc=size(c,2);           % # of constant fixed components
nd=size(cft,2);         % # of time-varying fixed components
nq=size(vq,2);      	% # of directions
nb=(nq>0)+nd+nc;         	% # of total basis-components
if nd>0 || ~isempty(dift)
    % calculate index offset corresponding to the third subscript of vq
    % since dift and/or cft is different in each time frame, we call cdmpfitr with only one time frame but nt*nf frequency bins
    if isempty(dift)
        [rmft,sdft,doft,erft]=cdmpfitras(reshape(ry,nm2,1,nt*nf),cat(2,reshape(cft,[nm2 nd nt*nf]),reshape(repmat(c,[1 nt 1]),[nm2 nc nt*nf])),reshape(repmat(vq,[1 nt 1]),[nm2 nq nt*nf]));
        doft=reshape(doft,nt,nf);
    else
        [rmft,sdft,x,erft]=cdmpfitras(reshape(ry,nm2,1,nt*nf),cat(2,reshape(cft,[nm2 nd nt*nf]),reshape(repmat(c,[1 nt 1]),[nm2 nc nt*nf])),reshape(vq(:,dift+repmat(nq*(0:nf-1),nt,1)),nm2,1,nt*nf));
        doft=dift;                  % force the selected direction to the a priori value
    end
    rmft=reshape(rmft,nm2,nt,nf);
    sdft=reshape(sdft,nb,nt,nf);
    erft=reshape(erft,nt,nf);
else
    [rmft,sdft,doft,erft]=cdmpfitras(ry,c,vq);
end

