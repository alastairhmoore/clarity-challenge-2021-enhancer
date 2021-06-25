function[frame_end_times] = get_frame_end_times(pm)
% stft_v2 gives us the center of each frame using
%  pm.t = (pm.fr_st - 1 - pm.pre_pad_len + (pm.n_win_fwd+1)/2) ./ pm.fs;
% where
% pm.fr_st is 1-based index of first sample in each frame and the samples have
% already been pre-padded

frame_end_times = (pm.fr_st - 1 - pm.pre_pad_len + (pm.n_win_fwd+1)) ./ pm.fs;
