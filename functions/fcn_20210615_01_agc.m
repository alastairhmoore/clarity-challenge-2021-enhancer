function[per_frame_gain,Y] = fcn_20210615_01_agc(Y,pm,in_params)

params.target_level_dbfs = -55; %(0 dBFS = 120 dB SPL -> -55 dBFS = 120-55 = 65 dB SPL)
params.headroom_db = 6;         % after initial adaptation, threshold is this much higher
params.level_estimation_time_constant = 0.200;
params.init_time = 0.1;         % average power over this duration at start 
params.clarity_ramp_duration = 0.5; % time over which signals are faded in 

params.do_plots = 0;

if nargin > 2 && ~isempty(in_params)
    params = override_valid_fields(params,in_params);
end

[nFreq,nChan,nFrames] = size(Y);


fade_in_t = pm.t;
fade_in_t(fade_in_t>params.clarity_ramp_duration) = [];
fade_in_gain = cos(2*pi*(fade_in_t + params.clarity_ramp_duration) ...
                   ./(2*params.clarity_ramp_duration));
fade_in_gain = (fade_in_gain + 1)./2;
clarity_inv_ramp = ones(size(pm.t));
clarity_inv_ramp(1:length(fade_in_gain)) = 1./fade_in_gain;
clarity_inv_ramp(fade_in_t<params.init_time ) = 1;



%%


nInitFrames = sum(pm.t<params.init_time);
smoothed_power = mean(frame_power(Y(:,1:2,1:nInitFrames),pm),1); % [1 2]
per_frame_gain = ones(nFrames,1); % first nInitFrames are untouched, others are overriden

if params.do_plots
    power_db_history = zeros(nFrames,nChan);
    smoothed_power_db_history = zeros(nFrames,nChan);
    power_db_history(1:nInitFrames,:) = 10*log10(frame_power(Y(:,1:2,1:nInitFrames),pm));
    smoothed_power_db_history(1:nInitFrames,:) = repmat(10*log10(smoothed_power),nInitFrames,1);
end

tinc = diff(pm.t(1:2));
ax = exp(-tinc./params.level_estimation_time_constant);
agc_gain_lin = 1;
for iframe = (nInitFrames+1):nFrames
    this_frame_power = frame_power(agc_gain_lin * clarity_inv_ramp(iframe) * Y(:,1:2,iframe),pm);
    smoothed_power = ax * smoothed_power + (1-ax) * this_frame_power;
    max_smoothed_power_db = 10*log10(max(smoothed_power));
    if (pm.t(iframe) < 1) && (max_smoothed_power_db > params.target_level_dbfs)
        gain_step = 10^((params.target_level_dbfs-max_smoothed_power_db)/20);
        agc_gain_lin = gain_step * agc_gain_lin;
        smoothed_power = gain_step.^2 * smoothed_power;
    elseif (max_smoothed_power_db > (params.target_level_dbfs + params.headroom_db))
        gain_step = 10^((params.target_level_dbfs-max_smoothed_power_db)/20);
        agc_gain_lin = gain_step * agc_gain_lin;
        smoothed_power = gain_step.^2 * smoothed_power;
    end

    % accumulation
    per_frame_gain(iframe) = agc_gain_lin;
    if params.do_plots
        power_db_history(iframe,:) = 10*log10(this_frame_power);
        smoothed_power_db_history(iframe,:) = 10*log10(smoothed_power);
    end
end

if nargout>1 || params.do_plots
    Y = bsxfun(@times,Y,permute(per_frame_gain,[3 2 1]));
end

if params.do_plots
    out = stft_v2('inv',Y,pm);   
    
    my_ramp = ones(size(y,1),1);
    idc_ramp = (round(params.init_time*fs)+1):round(params.clarity_ramp_duration*fs);
    my_ramp(1:idc_ramp(1)-1)=0;
    my_ramp(idc_ramp) = linspace(0,1,length(idc_ramp));
 
    out = bsxfun(@times, out, my_ramp);
    rms_out = 10*log10(mean(out.^2,1));

    figure;
    % % idc_frames = (nInitFrames+1):nFrames;
    % plot(pm.t(idc_frames),smoothed_power_db_history(idc_frames,:))
    subplot(2,1,1)
    plot(pm.t,10*log10(frame_power(Y(:,1:2,:),pm)))
    hold all;
    plot(pm.t,power_db_history)
    plot(pm.t,smoothed_power_db_history)
    plot(pm.t([1 end]),params.target_level_dbfs([1 1]),'r--');
    plot(pm.t([1 end]),rms_out([1 1]));
    ylabel('Signal power [dB]')
    xlabel('Time [s]')
    subplot(2,1,2)
    plot(pm.t,20*log10(per_frame_gain))
    ylabel('AGC Gain [dB]')
    xlabel('Time [s]')
    linkaxes(get(gcf,'children'),'x')
end