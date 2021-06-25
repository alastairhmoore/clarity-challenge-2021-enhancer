% example.m
%
% Demonstrates the use of the ELO-SPHERES consortium enhancer
%
% To make this repository self-sufficient a few intermediate files which
% were generated using tools provided by the cec1 software are included
% by way of example.

if ~exist('hrir_dir','var')
    fprintf('hrir_dir should specify the directory of the HRIRs\n')
    input('press return to find it now')
    hrir_dir = uigetdir('Path to HRIRs');
end


%% data configuration
clarity_data_dir = './clarity_data'; % point this to your 
dataset = 'dummy';
scene_id = 'S08501';
listener_id = 'L0000';


%% experiment parameters
% give the experiment a name
EXP_ID = 'example';

% path of the openmha cfg file directory, relative to dataset root
OPENMHA_REL_DIR = 'saved_cfgs'; % generated separately

% assumed level of input signal - used to read the HL correction gain table
NOMINAL_SIGNAL_LEVEL_DB_SPL = 65;

% specify target - controls wav file format for saving
TARGET_PLATFORM = 'listen@home';

% account for the differce in calibration level between metric and listener
% evaluation
OFFSET_GAIN_DB = 20;



%% Setup the inputs/parameter structures required by beamformer function
% frame/stft settings
p.fs = 44100;
p.frame_duration = 0.004989; %200 samples at 44100
p.frame_overlap_frac = 0.5;
in_params.stft_params = convert_in_params_to_stft_params(p);

% required for loading steering vectors
in_params.hrir_database_root = hrir_dir;

% set the ouput format/scaling
in_params.target_platform = TARGET_PLATFORM;
in_params.linear_scale_factor = 10^(OFFSET_GAIN_DB/20);



% input audio
scene_dir = fullfile(clarity_data_dir, dataset, 'scenes');
% output audio
enhanced_dir = fullfile(clarity_data_dir,dataset,sprintf('enhanced_%s',EXP_ID));


% path of input/output signals
% - clarity uses stereo pairs so just pass the stub and rely on mvdr
%   function to call the custom loading function
in_wav_file_stub = fullfile(scene_dir,sprintf('%s_mixed',scene_id));
out_wav_file_path = fullfile(enhanced_dir,...
    sprintf('%s_%s_HA-output.wav',scene_id,listener_id));

% path of headtracker file - not used in clarity    
ht_file_path = '';

% listener-specific configuration
% - read from openmha cfg file
listener_characteristics.openmha_gains = ...
    fcn_20210615_04_get_openmha_correction(...
        fullfile(clarity_data_dir,dataset,OPENMHA_REL_DIR,...
        sprintf('%s_openmha.cfg',listener_id)),...
        NOMINAL_SIGNAL_LEVEL_DB_SPL); % gains at specified frequencies

% pass any oracle data
% - not used for final evaluation
oracle_data = struct();     % empty for eval

% additional data can be saved
% - here we just save it in the same place as the output wavs
saved_data_dir = enhanced_dir;
temp_data_dir = enhanced_dir;



%% finally call the mvdr function which does the processing
mvdr_clarity_challenge_cec1(...
    in_wav_file_stub ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % corrections to apply
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir ...            % path to folder where temporary files should be written - normally empty but may be useful for debugging
    )
        
        
