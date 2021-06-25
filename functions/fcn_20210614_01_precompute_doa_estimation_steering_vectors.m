function[out] = fcn_20210614_01_precompute_doa_estimation_steering_vectors(...
    in_params)

params.fs = [];
params.az_deg = (-180:7.5:(180-7.5)).';
params.hrir_database_root = '';
params.clarity_array_type = 'bte3';
allowed_values.clarity_array_type = {'bte1','bte2','bte3'};
params.doa_estimation_frame_duration = 0.032;
params.doa_estimation_frame_overlap_frac = 0.5;


params = override_valid_fields(params,in_params,allowed_values);



%% get a list of hrir names
files=dir(fullfile(params.hrir_database_root,'HRIRs_MAT','*.mat'));
hrir_names=cell(0,1);
for i=1:length(files)
    id=find(files(i).name=='-',1,'first');
    hrir_names{end+1}=files(i).name(1:id-1);
end
hrir_names=unique(hrir_names)';


% need stft params so we can do the right size fft
params.frame_duration = params.doa_estimation_frame_duration;
params.frame_overlap_frac = params.doa_estimation_frame_overlap_frac;
stft_params = convert_in_params_to_stft_params(params);


model_params.fs = params.fs;
model_params.stft_params = stft_params;
model_params.clarity_array_type = params.clarity_array_type;
model_params.doa_circ.az_deg = params.az_deg(:); 
model_params.hrir_path = params.hrir_database_root;


% hrirs setup
nAz = length(params.az_deg);
nHRIR= length(hrir_names);
r_circ=[];%zeros(nChan^2,0,nFreq);

for i_hrir=1:nHRIR   
    % update with new hrir
    model_params.hrir_id = hrir_names{i_hrir};
    model = fcn_20210329_ahm_01_get_compact_model_bases_for_clarity(model_params);

    % convert bases to REAL vectorised versions with right dimensions 
    this_r = permute(cdmpcovr(model.Rcirc),[1 3 2]);     %[nChan^2, nDOA, nFreq]    
    r_circ = cat(2,r_circ,this_r);  
end

out.r_circ = r_circ;
out.az_deg_repeated = reshape(repmat(model_params.doa_circ.az_deg(:),...
                                     1,nHRIR),...
                              [],1); % nDOA*nRIR x 1
out.nAz = nAz;
out.nHRIR = nHRIR;
out.stft_params = stft_params;
out.az_deg = params.az_deg;
out.hrir_names = hrir_names;
