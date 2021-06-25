function[gt] = fcn_20210615_04_get_openmha_correction(filepath, input_level_db_spl)

identifier='mha.mhachain.overlapadd.smoothgains_bridge.mhachain.fftfilterbank.f';
gt.f = ExtractMatrixFromText(filepath,identifier);

identifier='mha.mhachain.overlapadd.smoothgains_bridge.mhachain.dc.gtdata';
data = ExtractMatrixFromText(filepath,identifier);

identifier='mha.mhachain.overlapadd.smoothgains_bridge.mhachain.dc.gtmin';
min = ExtractMatrixFromText(filepath,identifier);

identifier='mha.mhachain.overlapadd.smoothgains_bridge.mhachain.dc.gtstep';
step_db = ExtractMatrixFromText(filepath,identifier);

nBands = length(gt.f);
[nRows,nLevels] = size(data);
if nRows ~= 2*nBands
    error('Shape is wrong')
end
table_levels = min + step_db * (0:nLevels-1);
[~,ilevel] = find_nearest(input_level_db_spl,table_levels);

idc_left = 1:nBands;
idc_right = nBands + (1:nBands);
gt.left_gain_db = data(idc_left,ilevel);
gt.right_gain_db = data(idc_right,ilevel);



%% Sina made this
function out = ExtractMatrixFromText(textpath,identifier)
%UNTITLED4 Summary of this function goes here
%   INPUT:
%   textpath            (string) full path to the text file
%   identifier          (string) name of the matrix in the text file

cfg  = fopen(textpath);
while( ~feof(cfg) )
    line = fgetl(cfg);
    if contains(line,[identifier ' ='])
        eval([line ';']);
    end
end
out=eval(identifier);
fclose(cfg);

