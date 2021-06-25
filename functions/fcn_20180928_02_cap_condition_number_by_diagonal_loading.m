function[R] = fcn_20180928_02_cap_condition_number_by_diagonal_loading(R,threshold)

n = size(R,1);
cx0=cond(R);    % save original condition number
%
% now correct the condition number
%
if cx0>threshold      % correct if too high
    e=eig(R);   % find eigenvalues
    R=R+(max(e)-threshold*min(e))/(threshold-1)*eye(n); % add a multiple of the identity
end
cx=cond(R); % this will always be <= thr