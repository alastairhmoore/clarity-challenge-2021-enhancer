function[out_path] = clarity_writewav_wrapper(z,fs,out_path,platform)

% drop-in replacement for elobes_writewav_wrapper
% - hardcoded switch for output target
% - assumes input is using 0 dB FS = 120 dB

if nargin < 4 || isempty(platform)
    platform = 'mbstoi';
end

allowed_values = {'mbstoi','listen@home'};
if ~ismember(platform, allowed_values)
    error('platform was not one of allowed values')
end


DO_CLIP = 1;
    

switch platform
    case 'mbstoi'
        % forces writewav to scale the output and save as 32 bit floating point 
        % values. Use readwave with 'g' option to recover the absolute values
        mode_str = 'gesvL'; 

        % scale the values so that this value maps to +/-1. Since floating point is 
        % used there will be no clipping if this value is exceeded.
        ref_amplitude = 1;
        clip_vals = [-1 1];
        
        if fs~=44100
            error('Expected input to be at 44.1kHz sample rate')
        end
    case 'listen@home'
        % 16bit at 32kHz
        fs_out = 32000;
        mode_str = 'gesc16';
        ref_amplitude = 1;
        clip_vals = [-2^15, (2^15)-1]./ (2^15);
        z = resample(z,fs_out,fs);
        fs = fs_out;
end
        

% ensure the folder structure is in place
[outdir,~] = fileparts(out_path);
check_output_dir_exists(outdir)

try
    if DO_CLIP
        if 1
            nClipped = sum(z(:)<clip_vals(1)) + sum(z(:)>clip_vals(2));
            nSamples = numel(z);
            if nClipped~=0
                [~, fn] = fileparts(out_path);
                fprintf('%s: %d of %d (%2.3f%%) samples clipped\n',...
                    fn, nClipped, nSamples, 100*nClipped/nSamples)
            end
        end
        z(z>clip_vals(2)) = clip_vals(2);
        z(z<clip_vals(1)) = clip_vals(1);
    end
    v_writewav(z,fs,out_path,mode_str,[],[],ref_amplitude);
    %fprintf('\n\nProcessed audio saved to %s\n',out_path);
catch
    error('Error saving file to %s',out_path);
end


