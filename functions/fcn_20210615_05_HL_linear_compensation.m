function[fft_gains] = fcn_20210615_05_HL_linear_compensation(fs,fft_freqs,openmha_correction)


fft_gains.l = 10.^(freq_interp_sh(openmha_correction.f,...
    openmha_correction.left_gain_db,...
    fft_freqs) ./20);

fft_gains.r = 10.^(freq_interp_sh(openmha_correction.f,...
    openmha_correction.right_gain_db,...
    fft_freqs) ./20);

% figure;
% plot(openmha_correction.f,openmha_correction.left_gain_db);
% hold all
% plot(openmha_correction.f,openmha_correction.right_gain_db);
% plot(fft_freqs,20*log10(fft_gains.l));
% plot(fft_freqs,20*log10(fft_gains.r));

function y = freq_interp_sh( f_in, y_in, f )
% FREQ_INTERP - linear interpolation on logarithmic frequency scale
% with sample and hold on edges
%
% Usage:
% y = freq_interp_sh( f_in, y_in, f );

% This file is part of the HörTech Open Master Hearing Aid (openMHA)
% Copyright © 2011 2013 2017 HörTech gGmbH
%
% openMHA is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% openMHA is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License, version 3 for more details.
%
% You should have received a copy of the GNU Affero General Public License, 
% version 3 along with openMHA.  If not, see <http://www.gnu.org/licenses/>.


  ;
  f_in = max(eps,f_in);
  y = interp1(log([0.5*f_in(1);f_in(:);2*f_in(end)]),...
              [y_in(1);y_in(:);y_in(end)],...
              log(f),'linear','extrap');

