function[out] = convert_in_params_to_stft_params(in_params)
% wrap up everything to do with frame length conversion
% return struct containing the required inputs to stft function
% namely:
%
%        win_anal
%   sig_frame_inc
%            nfft
%              fs

required_fields = {'fs','frame_duration','frame_overlap_frac'};
if any(~ismember(required_fields,fieldnames(in_params)))
    error('in_params must contain all required fields:\n%s',...
        sprintf('\t%s\n',string(required_fields.')))
end

% setup
fs = in_params.fs;

% frame lengths come directly from params
sig_frame_len = round(in_params.frame_duration * fs);
sig_frame_inc = round((1-in_params.frame_overlap_frac) * in_params.frame_duration * fs);

% using naive frequency domain processing - fft is same size as frame
nfft = sig_frame_len;
w = sqrt(hamming(sig_frame_len,'periodic'));
% w=w/sqrt(sum(w(1:nh:nw).^2)*nh*nw);
win_anal{1} = w ./ sqrt(sum(w(1:sig_frame_inc:sig_frame_len).^2 * sig_frame_len * sig_frame_inc));
% w=w*sqrt(nh*nw/sum(w(1:nh:nw).^2)); % scale to give perfect reconstruction
win_anal{2} = w .* sqrt(sig_frame_len * sig_frame_inc / sum(w(1:sig_frame_inc:sig_frame_len).^2)); %

out.fs = fs;
out.nfft = nfft;
out.sig_frame_inc = sig_frame_inc;
out.win_anal = win_anal;