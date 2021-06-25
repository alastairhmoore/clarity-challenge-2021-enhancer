% makes sure we are only using built-in functions and those exposed at the
% following paths
fprintf('---- WARNING ----\n')
fprintf('Running this script will clear your workspace and reset your path\n')
input('Press return to continue or ctrl+c to cancel\n')

restoredefaultpath

% https://github.com/ImperialCollegeLondon/sap-voicebox
addpath('submodules/sap-voicebox/voicebox')

% https://github.com/ImperialCollegeLondon/sap-elobes-utilities
addpath('submodules/sap-elobes-utilities/')


addpath('./functions')
addpath('./mb_fitting')