function[out] = fcn_20210329_ahm_01_get_compact_model_bases_for_clarity(in_params)

required_fields = {'fs','stft_params','hrir_id','clarity_array_type',...
    'doa_circ',};
if any(~ismember(required_fields,fieldnames(in_params)))
    error('in_params must contain all required fields:\n%s',...
        sprintf('\t%s\n',string(required_fields.')))
end


for idoa=length(in_params.doa_circ.az_deg):-1:1
    [h_in, fs_in] = load_clarity_hrir(in_params.hrir_id,...
        in_params.clarity_array_type,...
        in_params.doa_circ.az_deg(idoa));
    h_a(:,:,idoa) = resample(h_in,in_params.fs,fs_in);
end

H_a = rfft(h_a,in_params.stft_params.nfft,1);
[nFreq,nChan,nDOA] = size(H_a);


% covariance for each direction
out.Rcirc = bsxfun(@times,permute(H_a,[2 4 1 3]),conj(permute(H_a,[4 2 1 3]))); %[nChan nChan nFreq nDOA]

% covariance assuming same power from each direction
out.Riso = mean(out.Rcirc,4); %[nChan nChan nFreq]

% covariance assuming same power in each mic
out.Rwhite = repmat(eye(nChan),1,1,nFreq); %[nChan nChan nFreq]


