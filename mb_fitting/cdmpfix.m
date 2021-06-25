function [msku,dch,qph]=cdmpfix(ci,vqi)
% precalculate different combinations of free variables
%  Inputs: ci(nm^2,n)  Real-vectorized covariance matrices for fixed components
%          vqi(nm^2,nq) Directional basis component for each of nq directions
%
% Outputs: msku(2^n,n) Mask of unconstrained variables
%          dch{2^n,1}    cell array of (ci'*ci)\ci' for unconstrained columns of ci
%                        dch{k} is of dimension (u,nm^2) where u(k) is number of freee variables in state k
n=size(ci,2); % number of fixed variables
nq=size(vqi,2); % number of search directions
k=pow2(n); % number of variable combinations
msku=false(k,n);
dch=cell(k,1);
qph=cell(k,nq);
bmweights=pow2(1-n:0); % weights to convert bitmap to mask
warning('off'); % at DC especially, (vcui'*vcui) will be singuar
for i=1:k
    bmi=i-1; % bitmap of uncnnstrained variables
    mskui=rem(floor(bmi*bmweights),2)>0; % convert bitmap to mask
    msku(i,:)=mskui;
    cui=ci(:,mskui);
    dch{i}=(cui'*cui)\cui';
    vcui=[vqi(:,1) cui];
    for j=1:nq
        vcui(:,1)=vqi(:,j);
        qph{i,j}=(vcui'*vcui)\vcui';
    end
end
warning('on');