function[sv] = fcn_20210615_02_precompute_steering_vectors(in_params)

params.hrir_database_root = [];
params.stft_params = [];
params.clarity_array_type = 'bte3';
params.binaural_mode = 'binaural';
params.rtf_ref = 'array';
params.look_directions_deg = [30:-7.5:-30];

params = override_valid_fields(params,in_params);

%% get a list of hrir names
files=dir(fullfile(params.hrir_database_root,'HRIRs_MAT','*.mat'));
hrir_names=cell(0,1);
for i=1:length(files)
    id=find(files(i).name=='-',1,'first');
    hrir_names{end+1}=files(i).name(1:id-1);
end
hrir_names=unique(hrir_names)';


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

% channel selection depends on binaural mode
switch params.binaural_mode
    case 'binaural' % uses all nChan channels from bte array    
        idc_binaural = sort([ema.channelsLeft(:);ema.channelsRight(:)]);
        idc_left = idc_binaural;
        idc_right = idc_binaural;
    case 'bilateral' % uses only the channels from corresponding side
        idc_left = ema.channelsLeft;
        idc_right = ema.channelsRight;
end

% need stft params so we can do the right size fft
stft_params = params.stft_params;


%% compute the steering vectors
nHRIR = length(hrir_names);
nLook = length(params.look_directions_deg);
sv = struct();
for i_hrir = 1:nHRIR
    hrir_id = hrir_names{i_hrir};
    clear h
    for idoa = 1:nLook
        [htmp, fs_in] = load_clarity_hrir(hrir_id, ...
                                       params.clarity_array_type,...
                                       params.look_directions_deg(idoa)); %[nSamples,nChans]
        h(:,:,idoa) = resample(htmp, stft_params.fs, fs_in);
    end
    H = v_rfft(h,stft_params.nfft,1);    

    %% rtf for target
    sv.(hrir_id).d_l_t = permute(bsxfun(@rdivide,H(:,idc_left,:),H(:,ema.refChanLeft,:)),[2 1 3]); % [nChan,nFreq,nLook]
    sv.(hrir_id).d_r_t = permute(bsxfun(@rdivide,H(:,idc_right,:),H(:,ema.refChanRight,:)),[2 1 3]); % [nChan,nFreq,nLook]
    sv.(hrir_id).az_deg = params.look_directions_deg;
end
