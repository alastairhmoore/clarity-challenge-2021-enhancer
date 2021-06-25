function[] = mvdr_clarity_challenge_cec1(...
    in_wav_file_stub ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % defines the HL compensation
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir ...            % path to folder where temporary files should be written - normally empty but may be useful for debugging
    )
% mvdr_clarity_challenge_cec1 implements an mvdr beamformer where the 
% PSD/covariance matrix is obtained from a noise-only reference signal.
% Given an HRIR set, the direction of arrival and particular head within
% the HRIR is estimated/selected.
%
% Data is expected to saved in the format used by the clarity challenge
% where microphone signals are stored in stereo pairs. File paths are
% therefore actually stubs which require the CH channel number suffix to be
% appended.
%
% Usage:
%  mvdr_clarity_challenge_cec1(...
%      in_wav_file_path, out_wav_file_path, ht_file_path, listener_characteristics, ...
%      in_params, oracle_data, saved_data_dir, temp_data_dir)
%
% Outputs:
%  None
%
% Inputs:
%  in_wav_file_stub:
%      received signal
%  out_wav_file_path:
%      path where enhanced signal(s) should be written
%  ht_file_path:
%      head tracker signal
%  listener_characteristics:
%      struct containing field 'openmha_gains' as returned by 
%      fcn_20210615_04_get_openmha_correction.m
%  in_params:
%      struct
%  oracle_data:
%      sruct containing information which would not be available in a
%      real-world application but may be useful in development
%  saved_data_dir:
%      folder for storing intermediate data for exchange between modules or
%      supplementary results which should not be included in metrics.
%      May be empty if test bench dictates that this data should not be retained
%      after the experiment, in which case intermediate data should be
%      saved to temporary storage (see tempname).
%  temp_data_dir:
%      folder for storing data which is not required once processing is
%      complete. Normally empty, in which case use standard temporary directory
%      (see tempname) but useful for debugging.
%
% Alastair Moore, June 2021


%% validate the input

% final_data_dir
final_data_dir = fileparts(out_wav_file_path);
check_output_dir_exists(final_data_dir)

% ht_file_path
% if ~exist(ht_file_path,'file'), error('Couldn''t find %s',ht_file_path), end

% oracle_data (struct)
% check the required field(s) are present
required_oracle_data = {...
%     'target_azimuth_degrees'... % direction of target (used for metrics)
%     'hrir_id'...                % file containing the ground truth hrtfs
    };
for ireq = 1:length(required_oracle_data)
    if ~isfield(oracle_data,required_oracle_data{ireq})
        error('oracle_data is missing %s field',required_oracle_data{ireq})
    end
end

% saved_data_dir (string)
if nargin<7 || isempty(saved_data_dir)
    error('saved_data_dir should not be empty - want to retain the data')
end
check_output_dir_exists(saved_data_dir);

% temp_data_dir (string)
% if nargin<8 || isempty(temp_data_dir)
%     temp_data_dir = tempname;
% end
% check_output_dir_exists(temp_data_dir);

%% defualt parameters
params.fs = [];
params.c = soundspeed();
params.stft_params = [];
params.max_condition_number = 1000;
params.mvdr_fcn = @fcn_202106_03_mvdr_weights_fast;
params.mvdr_debug_level = 0;%


params.binaural_mode = 'binaural';
allowed_values.binaural_mode = {'binaural','bilateral'};
params.clarity_array_type = 'bte3';
allowed_values.clarity_array_type = {'bte1','bte2','bte3'};
params.rtf_ref = 'array'; % or 'in-ear'
params.hrir_database_root = [];
params.noise_only_duration = 2;
params.noise_cov_update_times = [0:0.1:2]';
params.max_lookahead_secs = 0.005;
params.linear_scale_factor = 1;
params.target_platform = 'mbstoi';

params.steering_vectors = [];
params.est_doa_hrir_path = [];

%% update with input parameters
if ~isempty(in_params)
    params = override_valid_fields(params,in_params,allowed_values);
end

if strcmp(params.binaural_mode, 'bilateral') && ...
        strcmp(params.clarity_array_type, 'bte1')
    error('Cannot do beamforming - only 1 channel per filter!')
end


%% setup our microphone array
% shoehorn the clarity data into a format compatibile with our normal 
% ElobesMicArray class
ema = fake_ema_for_clarity(params.clarity_array_type);
switch params.rtf_ref
    case 'array'
        % do nothing - use ema.refChanLeft and ema.refChanRight as is
    case 'in-ear'
        % slight hack - we will append the in-ear pair to the end of the
        % array
        ema.refChanLeft = ema.channelsLeft(end) + 2;
        ema.refChanRight = ema.channelsRight(end) + 2;
    otherwise
        error('Unknown option for params.rtf_ref')
end

%% read in data and do pre-processing
%fprintf('Loading wav file...\n')
params.fs = params.stft_params.fs;
[y, fs_in] = load_clarity_multichannel_wav(in_wav_file_stub,params.clarity_array_type); %[nSamples,nChans]
y = resample(y,params.fs,fs_in);

% need this for fft sizes
stft_params = params.stft_params;

%% steering vectors
if isempty(params.steering_vectors)
    warning('Calculating steering vectors - passing these in will speed up batch processing')
    sv_params.hrir_database_root = params.hrir_database_root;
    sv_params.stft_params = stft_params;
    sv_params.fs = [];
    sv_params.clarity_array_type = params.clarity_array_type;
    sv_params.binaural_mode = params.binaural_mode;
    sv_params.look_directions_deg = [30:-7.5:-30];
    sv_params.rtf_ref = params.rtf_ref;
    sv = fcn_20210615_02_precompute_steering_vectors(sv_params);
else
    sv = params.steering_vectors;
end
hrir_names = fieldnames(sv);


% convert input signal
[Y,~,pm] = stft_v2('fwd', y, ...
    stft_params.win_anal, ...
    stft_params.sig_frame_inc,...
    stft_params.nfft,...
    stft_params.fs);
[nFreq,nChan,nFrames] = size(Y);

% instantaneous outer product - used for covariance matrix calculation
Ry = bsxfun(@times,permute(Y,[2 4 1 3]),conj(permute(Y,[4 2 1 3]))); %[nChan nChan, nFreq,nFrames]

% get HL compensation gains
fft_gains = fcn_20210615_05_HL_linear_compensation(pm.fs,pm.f,...
    listener_characteristics.openmha_gains); % struct with elements: [nFreq,1]


%% load in the precomputed estimates of doa and hrir_id
%  update_time: time when estimate was made [seconds]
%      hrir_id: estimated name of hrir (str);
%   est_az_deg: estimated target azimuth [degrees]
if ~isempty(params.est_doa_hrir_path)
    in = load(params.est_doa_hrir_path,'doa_hrir_estimates');
else
    % first make the prerequisite datastructure
    dp1.fs = params.fs;
    dp1.clarity_array_type = params.clarity_array_type;
    dp1.doa_estimation_frame_duration = 0.04998; %2204 samples at 44100
    dp1.hrir_database_root = params.hrir_database_root;
    doa_pre_dat = fcn_20210614_01_precompute_doa_estimation_steering_vectors(...
        dp1);

    % do the estimation
    dp2.fs = params.fs;
    dp2.doa_estimation_f_max = 16000;
    in.doa_hrir_estimates = fcn_20210614_03_estimate_doa_and_hrir(y,...
            dp2, doa_pre_dat);
        
end        
% get the steering vector required for each estimate

% metadata
% noise_updates:            one entry per update - only during noise-only period
% doa_updates:              one entry per update - only during target activity
% bf_steering_vectors:      unique set of vectors
% bf_updates:               defines the switching between beamformers

% define default first
bf_steering_vectors.hrir_id{1} = hrir_names{1}; % just pick the first one
bf_steering_vectors.look_az_deg(1) = 0;         % assume straight ahead
bf_steering_vectors.d_l_t{1} = local_extract_sv(sv, hrir_names{1}, 0,'left');
bf_steering_vectors.d_r_t{1} = local_extract_sv(sv, hrir_names{1}, 0,'right');

nDOA_hrir_estimates = length(in.doa_hrir_estimates);
doa_updates.update_time = zeros(nDOA_hrir_estimates,1);
doa_updates.sv_id = zeros(nDOA_hrir_estimates,1);
for iupdate = 1:nDOA_hrir_estimates
    doa_updates.update_time(iupdate) = in.doa_hrir_estimates(iupdate).update_time;
    
    % get the id
    hrir_id = in.doa_hrir_estimates(iupdate).hrir_id;
    look_az_deg = in.doa_hrir_estimates(iupdate).est_az_deg;
    % check to see if it already exists
    i_sv = find(strcmp(hrir_id, bf_steering_vectors.hrir_id) ...
                & look_az_deg == bf_steering_vectors.look_az_deg);
    if isempty(i_sv)
        % add it
        i_sv = length(bf_steering_vectors.hrir_id) + 1;
        bf_steering_vectors.hrir_id{i_sv} = hrir_id;
        bf_steering_vectors.look_az_deg(i_sv) = look_az_deg;
        bf_steering_vectors.d_l_t{i_sv} = local_extract_sv(sv, hrir_id, look_az_deg,'left');
        bf_steering_vectors.d_r_t{i_sv} = local_extract_sv(sv, hrir_id, look_az_deg,'right');
        doa_updates.sv_id(iupdate) = i_sv;
    elseif numel(i_sv)==1
        doa_updates.sv_id(iupdate) = i_sv;
    else
        error('Expected 0 or 1 matching entries');
    end
end


%% sequential estimates of NCM using same stft_params as the beamforming
if any(params.noise_cov_update_times > params.noise_only_duration)
    warning('Noise cov is being updated after noise only duration')
end


frame_end_times = get_frame_end_times(pm);

for iupdate = 1:length(params.noise_cov_update_times)
    update_time = params.noise_cov_update_times(iupdate);
    num_used_frames = sum(frame_end_times + params.max_lookahead_secs < update_time);
    if num_used_frames < nChan
        % matrix will be rank defficient, so default to isotropic
        % assumption
        % for computation speed - use spatially white assumption instead
        noise_updates(iupdate).Rv_est = repmat(eye(nChan),[1 1 nFreq nFrames]); %model.Riso;
    else
        noise_updates(iupdate).Rv_est = mean(...
            Ry(:,:,:,frame_end_times + params.max_lookahead_secs < update_time),...
            4);
    end
    noise_updates(iupdate).update_time = update_time;
end

%% create beamformer for each doa and noise estimate
% channel selection depends on binaural mode
switch params.binaural_mode
    case 'binaural' % uses all nChan channels from bte array    
        idc_binaural = 1:nChan;
        idc_left = idc_binaural;
        idc_right = idc_binaural;
    case 'bilateral' % uses only the channels from corresponding side
        idc_left = ema.channelsLeft;
        idc_right = ema.channelsRight;
end


for i_sv = 1:length(bf_steering_vectors.d_l_t)
    d_l_t = bf_steering_vectors.d_l_t{i_sv}; % [nChan,nFreq]
    d_r_t = bf_steering_vectors.d_r_t{i_sv}; % [nChan,nFreq]

    for iupdate = 1:length(noise_updates)
        % get the NCM to use
        Rv_est = noise_updates(iupdate).Rv_est;
               
        %% loop over freq to calculate weights
        % just use one switch statement for ease of reading
        w_l = zeros(nFreq,nChan);
        w_r = zeros(nFreq,nChan);

        switch params.binaural_mode
            case 'binaural'
                % uses all nChan channels from bte array
                for ifreq = 1:nFreq
                    % impose diagonal loading if necessary
                    R = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
                        Rv_est(idc_binaural,idc_binaural,ifreq), params.max_condition_number);

                    % compute beamformer weights
                    w_l(ifreq,idc_binaural) = conj(params.mvdr_fcn(R, d_l_t(:,ifreq), params.mvdr_debug_level).'); %
                    w_r(ifreq,idc_binaural) = conj(params.mvdr_fcn(R, d_r_t(:,ifreq), params.mvdr_debug_level).'); %
                end

            case 'bilateral' % uses only the channels from corresponding side
                for ifreq = 1:nFreq
                    % impose diagonal loading if necessary
                    R_l = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
                        Rv_est(idc_left,idc_left,ifreq), params.max_condition_number);
                    R_r = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
                        Rv_est(idc_right,idc_right,ifreq), params.max_condition_number);

                    % compute beamformer weights
                    w_l(ifreq,idc_left) = conj(params.mvdr_fcn(R_l, d_l_t(:,ifreq), params.mvdr_debug_level).'); %  
                    w_r(ifreq,idc_right) = conj(params.mvdr_fcn(R_r, d_r_t(:,ifreq), params.mvdr_debug_level).'); %
                end
        end

        % apply directly in frequency domain then obtain time domain signals
        Z(:,1,:) = sum(bsxfun(@times,w_l,Y),2); % broadcast over time
        Z(:,2,:) = sum(bsxfun(@times,w_r,Y),2); % broadcast over time
        
        
        % HL compensation
        Z(:,1,:) = bsxfun(@times,Z(:,1,:),fft_gains.l(:));
        Z(:,2,:) = bsxfun(@times,Z(:,2,:),fft_gains.r(:));

        % collect signals together
        out_cell{iupdate,i_sv} = stft_v2('inv',Z,pm);
        
    end
end

%% all the beamforming has now been done. Now consolidate beamformer selection
% into bf_updates structure

% during noise-only we cannot know where the source is so use the defualt
i_sv = 1;
for iupdate = 1:length(noise_updates)
    bf_updates(iupdate).update_time = noise_updates(iupdate).update_time;
    bf_updates(iupdate).noise_select = iupdate;
    bf_updates(iupdate).sv_select = i_sv;  
end

% once target is active, freeze the noise estimate, update the noisy
% signal estimate
for idoa_update = 1:length(doa_updates.update_time)
    iupdate = iupdate+1;
    update_time = doa_updates.update_time(idoa_update);
    bf_updates(iupdate).update_time = update_time;
    bf_updates(iupdate).noise_select = length(noise_updates); % keep using last one
    bf_updates(iupdate).sv_select = doa_updates.sv_id(idoa_update);  
end


%%
[nNoiseUpdates,nSteeringVectors] = size(out_cell);
nUpdates = length(bf_updates);
nSamples = size(y,1);
bf_xfade_gain = zeros(nSamples,nNoiseUpdates,nSteeringVectors);

% initialise values to those of first update
sv_select = bf_updates(1).sv_select;
noise_select = bf_updates(1).noise_select;

for iupdate=1:nUpdates
    istart = max(1,floor(params.fs * bf_updates(iupdate).update_time));
    
    if iupdate<nUpdates
        iend = floor(params.fs * bf_updates(iupdate+1).update_time)-1;
    else
        %use old value of seg_len
        iend = min(istart+seg_len-1,nSamples);
    end
    seg_len = length(istart:iend);

    % do fade out old beam
    bf_xfade_gain(istart:iend, noise_select, sv_select) = ...
        bf_xfade_gain(istart:iend, noise_select, sv_select) + linspace(1,0,seg_len).';
    
    % update state
    sv_select = bf_updates(iupdate).sv_select;
    noise_select = bf_updates(iupdate).noise_select;
    

    % fade in new
    bf_xfade_gain(istart:iend, noise_select, sv_select) = ...
        bf_xfade_gain(istart:iend, noise_select, sv_select) + linspace(0,1,seg_len).';
end

% set constant after fade in using existing beam
bf_xfade_gain((iend+1):end, noise_select, sv_select) = 1;

%% final output is weighted sum of beams
total_gain = zeros(nSamples,1);
out = zeros(nSamples,2);
for iel = 1:numel(out_cell)
    total_gain = total_gain + bf_xfade_gain(:,iel);
    out = out + bsxfun(@times,bf_xfade_gain(:,iel), out_cell{iel});
end
if any(total_gain-1 > eps)
    warning('total_gain must be 1')
   keyboard
end

%% compute agc
% - use the input signal to determine the gain per frame
agc_frame_gains = fcn_20210615_01_agc(Y,pm);
% - convert to gain per sample
agc_gain = fcn_20210615_03_frame_gains_to_sample_gains(size(y,1),...
                                                       agc_frame_gains,...
                                                       pm);
% - apply gain to output
out = bsxfun(@times,agc_gain,out);                                                   
                                                   

%% write out
clarity_writewav_wrapper( params.linear_scale_factor * out, ...
                          params.fs, ...
                          out_wav_file_path, ...
                          params.target_platform);
 

                      
function[d_t] = local_extract_sv(sv, hrir_id, look_az_deg, ear)

this_hrir = sv.(hrir_id);
[~, i_doa] = find(look_az_deg == this_hrir.az_deg);
if numel(i_doa)~=1, error('Expect a unique doa'), end
switch ear
    case 'left'
        d_t = this_hrir.d_l_t(:,:,i_doa);
    case 'right'
        d_t = this_hrir.d_r_t(:,:,i_doa);
end
