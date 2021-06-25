function y=cdmpcovc(x)
% Convert the real-valued vectorized form of a covariance matrix to complex hermitian form
% Inputs:   x(n^2,...)   Real-valued vectorized matrix array
% Outputs:  y(n,n,...)   Hermitian matrix array
%
% For a Hermition input X, the output vector is y=Y(:) where
%             { sqrt(2)*real(X(i,j))   for i<j
%    Y(i,j) = { X(i,j)                 for i=j
%             { sqrt(2)*imag(X(i,j))   for i>j
%
% Ix x is a column vector of length n^2, then y is a square nxn matrix.
% Otherwise y has one more dimension than x.
%
% The inverse of this routine is cdmpcovr
%
persistent m k a b c p q
s=size(x);
n=s(1); % number of channels squared
if isempty(k) || n~=k % calculate subscript indices if necessary
    k=n;
    m=sqrt(k);
    if k~=m^2
        error('First dimension must be of size n^2')
    end
    q=sqrt(0.5);
    p=q*1i;
    cvrow=1+mod(0:k-1,m)'; % covariance matrix element row indices
    cvcol=1+floor((0:k-1)/m)'; % covariance matrix element column indices
    a=find(cvrow==cvcol);  % list of diagonal elements
    c=find(cvrow>cvcol); % list of lower triangle elements
    b=cvcol(c)+m*(cvrow(c)-1); % indices of transposed lower triangle
end
v=reshape(x,k,[]); %
y=zeros(size(v)); % create space for the output
y(a,:)=v(a,:);
y(c,:)=q*v(b,:)+p*v(c,:);
y(b,:)=conj(y(c,:));
if numel(y)>k
    y=reshape(y,[m m s(2:end)]);
else
    y=reshape(y,[m m]);
end