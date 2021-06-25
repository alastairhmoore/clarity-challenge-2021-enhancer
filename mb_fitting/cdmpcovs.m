function y=cdmpcovs(x,p)
% Calculated real-valued vectorized covariance matrices from a multichannel STFT-domain signal
% Inputs:   x(nt,nf,nm)   multichannel STFT-domain signal with nt frames, nf frequency bins and nm channels
%           p             smoothing factor or 0 for no smoothing [default: 0]
%                         p=exp(-T/tau) where T=nhop/fs is frame hop and tau is smoothing time constant
% Outputs:  y(nm^2,nt,nf)   Real-valued vectorized matrices (one per time-frequency cell)
%
% For a Hermition input X, the output vector is y=Y(:) where
%             { sqrt(2)*real(X(i,j))   for i<j
%    Y(i,j) = { X(i,j)                 for i=j
%             { sqrt(2)*imag(X(i,j))   for i>j
persistent m k a b c d e f g h q
[nt,nf,nm]=size(x);
if isempty(m) || nm~=m  % calculate subscript indices if necessary
    m=nm;
    k=m^2;
    cvrow=1+mod(0:k-1,m)'; % covariance matrix element row indices
    cvcol=1+floor((0:k-1)/m)'; % covariance matrix element column indices
    a=find(cvrow==cvcol); % diagonal elements of vectorized matrix
    b=find(cvrow<cvcol); % strictly upper triangular elements of vectorized matrix
    c=find(cvrow>cvcol);  % strictly lower triangular elements of vectorized matrix
    lti=cvrow>=cvcol; % lower triangle mask
    g=cvrow(lti); % row indices of lower triangle
    h=cvcol(lti);  % column indices of lower triangle
    d=find(g==h);  % diagonal elements of vectorized lower triangle
    f=find(g>h);   % strictly lower triangle elements of vectorized lower triangle
    [ur,ix]=sort(g(f));
    e=f(ix);
    q=sqrt(2);
end
v=x(:,:,g).*conj(x(:,:,h)); % form products for lower triangular part of covariance matrix
w=zeros(nt,nf,k); % space for unfiltered version
w(:,:,a)=real(v(:,:,d));
w(:,:,b)=q*real(v(:,:,e));
w(:,:,c)=q*imag(v(:,:,f));
if nargin>1
    y=permute(reshape(filter(1,[1 -p],reshape(w,nt,nf*k)),[nt,nf,k]),[3 1 2]); % filter the observed covariance matrix
else
    y=permute(w,[3,1,2]);
end
