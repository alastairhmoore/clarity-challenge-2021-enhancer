function y=cdmpcovr(x)
% Convert a complex hermitian covariance matrix to its real-valued vectorized form
% Inputs:   x(n,n,...)   Hermitian matrix array
% Outputs:  y(n^2,...)   Real-valued vectorized matrices
%
% For a Hermition input X, the output vector is y=Y(:) where
%             { sqrt(2)*real(X(i,j))   for i<j
%    Y(i,j) = { X(i,j)                 for i=j
%             { sqrt(2)*imag(X(i,j))   for i>j
persistent m k a b c q
s=size(x);
n=s(1); % number of channels (must have s(2)=n also)
if isempty(m) || n~=m % calculate subscript indices if necessary
    q=sqrt(2);
    m=n;
    k=m^2;
    cvrow=1+mod(0:k-1,m)'; % covariance matrix element row indices
    cvcol=1+floor((0:k-1)/m)'; % covariance matrix element column indices
    a=cvrow==cvcol;
    b=cvrow<cvcol;
    c=cvrow>cvcol;
end
v=reshape(x,k,[]); %
y=zeros(size(v)); % create space for the output
y(a,:)=real(v(a,:));
y(b,:)=q*real(v(b,:));
y(c,:)=q*imag(v(c,:));
if length(s)>3
    y=reshape(y,[k s(3:end)]);
end