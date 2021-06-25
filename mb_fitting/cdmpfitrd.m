function [rmft,sdft,doft,erft]=cdmpfitrd(ry,c,vq,dift)
% Fit a multicomponent model to a microphone-array covariance matrix with a priori direction
%  Inputs: ry(nm^2,nt,nf)   Smoothed real-vectorized observed covariance matrices
%          c(nm^2,nb-1,nf)  Real-vectorized covariance matrices for fixed components
%          vq(nm^2,nq,nf)   Real-vectorized covariance matrices for directional component
%          dift(nt,nf)      Specify the direction index a priori
%
% Outputs: rmft(nm^2,nt,nf) Model fit to ry
%          sdft(nb,nt,nf)   Component weights: directional component first followed by fixed components
%          doft(nt,nf)      Index of directional component in vq array (will equal dift if specified)
%          erft(nt,nf)      Squared Frobenius norm of fitting error
%
%    where nm = number of microphone channels
%          nt = number of time frames
%          nf = number of frequency bins
%          nb = number of components in the model
%          nq = number of search directions
%
% This implements the fitting stage of the covariance matrix model described in [1].
% The input arguments are all in the form of real-vectorized covariance matrices; given a
% complex nm x nm Hermitian matrix, X, the real-vectorized form is Y(:) where
%             { sqrt(2)*real(X(i,j))   for i<j
%    Y(i,j) = { X(i,j)                 for i=j
%             { sqrt(2)*imag(X(i,j))   for i>j
% The functions cdmpcovr() and cdmpcovs() create real-vectorized vectors
% from covariance matrices and STFT-domain signals respectively.
%
% Versions:
%
% 2020_0619 Initial version.
%
%
[nm2,nt,nf]=size(ry);       % extract input dimensions
nb=size(c,2)+1;             % # total basis-components
nq=size(vq,2);              % number of directions
% create space for output arrays
rmft=zeros(nm2,nt,nf);      % Model fit to ry
sdft=zeros(nb,nt,nf);       %   Component weights: directional component first followed by fixed components
erft=zeros(nt,nf);          %    Squared Frobenius norm of fitting error
% calculate index offset corresponding to the third subscript of vq
% since dift is different in each time frame, we call cdmpfitr with only one time frame but nt*nf frequency bins
% [rmft,sdft,x,erft]=cdmpfitr(reshape(ry,nm2,1,nt*nf),reshape(repmat(c,[1 nt 1]),[nm2 nb-1 1 nt*nf]),reshape(vq(:,dift+repmat(nq*(0:nf-1),nt,1)),nm2,1,nt*nf));
[rmft,sdft,x,erft]=cdmpfitr(reshape(ry,nm2,1,nt*nf),reshape(repmat(c,[1 nt 1]),[nm2 nb-1 nt*nf]),reshape(vq(:,dift+repmat(nq*(0:nf-1),nt,1)),nm2,1,nt*nf));

rmft=reshape(rmft,nm2,nt,nf);
sdft=reshape(sdft,nb,nt,nf);
doft=dift;                  % force the selected direction to the a priori value
erft=reshape(erft,nt,nf);

