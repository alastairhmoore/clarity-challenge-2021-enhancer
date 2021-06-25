function[ir,fs] = load_clarity_hrir(id,type,azimuth_deg,database_root)

AZ_DEG_LIST = -180:7.5:180;

azimuth_deg = -azimuth_deg; % database labels negative azimuth to listeners left

if numel(azimuth_deg)~=1
    error('Only one azimuth')
end

if nargin<4 || isempty(database_root)
    database_root = '/Volumes/T7/clarity_CEC1_data/clarity_data/hrir';
end
if ~exist(database_root,'dir')
    error('Couln''t find directory at %s',database_root)
end

switch type
    case 'in-ear'
        filename = fullfile(database_root,'HRIRs_MAT',sprintf('%s-ED.mat',id));
        chan_select = [1 2];        
    case 'bte1'   
        filename = fullfile(database_root,'HRIRs_MAT',sprintf('%s-BTE_MultiCh.mat',id));
        chan_select = [1 2];
    case 'bte2'   
        filename = fullfile(database_root,'HRIRs_MAT',sprintf('%s-BTE_MultiCh.mat',id));
        chan_select = [1 2 5 6];
    case 'bte3'   
        filename = fullfile(database_root,'HRIRs_MAT',sprintf('%s-BTE_MultiCh.mat',id));
        chan_select = [1:6];
    otherwise
        error('unknown type')
end

Dat = load(filename);
nearest_az = find_nearest(azimuth_deg,AZ_DEG_LIST);
assert(abs(nearest_az-azimuth_deg)<0.5);
% account for periodicity
if nearest_az==180
    nearest_az=-180;
end
doa_select = Dat.M_directions(1,:)==nearest_az & Dat.M_directions(2,:)==0;
ir = permute(Dat.M_data(:,doa_select,chan_select),[1,3,2]);
fs = Dat.srate;
