function[out] = fcn_20210614_03_estimate_doa_and_hrir(y,in_params,doa_struct)

params.fs = [];
params.noise_only_duration = 2;
params.doa_estimation_update_times = [2.1:0.1:2.5];
params.doa_estimation_f_min = 0;
params.doa_estimation_f_max = 5000;
params.max_lookahead_secs = 0.005;

params = override_valid_fields(params,in_params);


% iherit settings from the precomputed doa_struct
stft_params = doa_struct.stft_params;

if params.fs~=stft_params.fs
    error('should run everything at the same sample rate')
end
fs = stft_params.fs;
az_deg_repeated = doa_struct.az_deg_repeated;


%% covert signal to stft domain
[Y,x_tail_anal,pm] = stft_v2('fwd', y, ...
    stft_params.win_anal, ...
    stft_params.sig_frame_inc,...
    stft_params.nfft,...
    stft_params.fs);
[nFreq,nChan,nFrames] = size(Y);

% important for ensuring causality
frame_end_times = get_frame_end_times(pm);

% instantaneous outer product of received signal
Ry = bsxfun(@times,permute(Y,[2 4 1 3]),conj(permute(Y,[4 2 1 3]))); %[nChan nChan, nFreq,nFrames]


% noise only covariance
Rv_est = mean(...
            Ry(:,:,:,frame_end_times + params.max_lookahead_secs < params.noise_only_duration),...
            4);
r_v = permute(cdmpcovr(Rv_est),[1 3 2]);

% loop over updated times
for iupdate = 1:length(params.doa_estimation_update_times)
    update_time = params.doa_estimation_update_times(iupdate);
    

    Ry_est = mean(Ry(:,:,:,frame_end_times + params.max_lookahead_secs < update_time),...
                  4);
    r_y = permute(cdmpcovr(Ry_est),[1 3 2]);

    % estimate the model parameters and the reconstructed NCM
    [rmodel,pow_component,imin,err]=cdmpfitras(r_y, r_v, doa_struct.r_circ);
    
    %%
    idc_valid = pm.f > 0 & pm.f < min(fs/2, 20000) & ...
                squeeze(pow_component(1,1,:))>0;
    idc_selected = pm.f > params.doa_estimation_f_min & ...
                   pm.f < params.doa_estimation_f_max & ...
                   az_deg_repeated(imin) >= -30 & ...
                   az_deg_repeated(imin) <=  30 & ...
                   idc_valid;

    % how many frequencies are assigned to each DOA
    [counts, doa_categories] = histcounts(categorical(az_deg_repeated(imin(idc_selected))));
    [~, imax] = max(counts);
    est_az_deg = str2double(doa_categories{imax});

    [doa, idoa] = find_nearest(est_az_deg, doa_struct.az_deg);
    d_az=abs(doa_struct.az_deg(2)-doa_struct.az_deg(1));
    %[rows cols]=ind2sub([nAz nRIR],imin(find(az_deg_repeated(imin)==doa))); % row: doa - col: rir
    [rows, cols]=ind2sub([doa_struct.nAz doa_struct.nHRIR],imin(find(abs(az_deg_repeated(imin)-doa)<=d_az))); % row: doa - col: rir
    
    
    [rir_count, rir_edges] = histcounts(cols,[0:doa_struct.nHRIR]+.5);
    
    %[~,est_sv_id]=max(rir_count);
    est_hrir_id=find(rir_count==max(rir_count));
    est_doa_id=idoa;
    
    out(iupdate).update_time = update_time;
    out(iupdate).hrir_id = doa_struct.hrir_names{est_hrir_id(1)};
    out(iupdate).est_az_deg = est_az_deg;
end
