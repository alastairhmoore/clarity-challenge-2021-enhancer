function[g] = fcn_20210615_03_frame_gains_to_sample_gains(nSamples,frame_gains,pm)

params.clarity_ramp_duration = 0.5;
params.mute_duration = 0.2;
fs = pm.fs;
time_constant = 0.01;

% base all updates on the end of the frame to absolutely guarantee
% causality
s_end = round(fs * get_frame_end_times(pm));


% loop over frames
g = ones(nSamples,1);
s_inc = (1:pm.inc).';
tinc = 1/fs;
ax = exp(-tinc./time_constant);

for iframe = 1:length(frame_gains)
    idc = s_end(iframe) + s_inc;        % index into the following frame
    idc(idc > nSamples) = []; % lop off end of last frame
    for ii = 1:length(idc)
        g(idc(ii)) = ax * g(idc(ii)-1) + (1-ax) * frame_gains(iframe);
    end
end
    

% finally apply a fade in to avoid intial glitching
my_ramp = ones(nSamples,1);
idc_ramp = (round(params.mute_duration*fs)+1):round(params.clarity_ramp_duration*fs);
my_ramp(1:idc_ramp(1)-1)=0;
my_ramp(idc_ramp) = linspace(0,1,length(idc_ramp));
g = g .* my_ramp;

if 0
	figure;plot((1:nSamples)./fs,g);hold all;plot(pm.t, frame_gains);
end