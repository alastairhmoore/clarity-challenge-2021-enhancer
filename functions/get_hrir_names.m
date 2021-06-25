%hrir_path='/Users/sinahafezi/Imperial College London/ELO-SPHERES - Clarity challenge/cec1/clarity_data/hrir/';
files=dir(fullfile(hrir_path,'HRIRs_MAT','*.mat'));
hrir_names=cell(0,1);
for i=1:length(files)
id=find(files(i).name=='-',1,'first');
hrir_names{end+1}=files(i).name(1:id-1);
end
hrir_names=unique(hrir_names)'