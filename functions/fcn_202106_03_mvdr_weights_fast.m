function[w,issingular] = fcn_202106_03_mvdr_weights_fast(R,d,debug_level)
%R is spatial covariance of stft of received signal [nChans,nChans,nFreq,nFrames]
%d is steering vector for each TF bin [nChans,nFreq,nFrames]
%debug_level controls how to handle singular or ill-conditioned matrices
% 0: ignore
% 1: print warning to console (as standard)
% 2: catch and allow debugging
% 3: throw error

% define the msg id for singular matrix inversion
MSGID_NEAR_SINGULAR = 'MATLAB:nearlySingularMatrix';
MSGID_SINGULAR = 'MATLAB:singularMatrix';
if 0
  % handy to have these lines for debugging
  warning('on', MSGID_NEAR_SINGULAR);
  warning('on', MSGID_SINGULAR);
end

if nargout > 1 || debug_level~=0
    error('No chance for niceness with this fast code. Try using fcn_20170222_01_mvdr_weights')
end

% don't care about singular values - suppress warnings
prev_warn_state(1) = warning('off', MSGID_NEAR_SINGULAR);
prev_warn_state(2) = warning('off', MSGID_SINGULAR);


%preallocate
[nChans,nChans2,nFreq,nFrames] = size(R);
w = zeros(nChans,nFreq,nFrames);

%for clarity and simplicity, loop over time and frequency indices
for ifreq = 1:nFreq
    for iframe = 1:nFrames
        % ---------
        invRd = R(:,:,ifreq,iframe) \ d(:,ifreq,iframe);
        w(:,ifreq,iframe) = invRd / (d(:,ifreq,iframe)' * invRd);
        % ---------
    end
end

warning(prev_warn_state)
