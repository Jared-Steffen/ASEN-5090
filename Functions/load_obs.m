datadir = '../HW3_DATA/';
f = 'NIST00USA_R_20262310000_01D_30S_MO'
eval(['!' datadir 'crx2rnx '  datadir f '.crx'])
data_o =rinexread([datadir f '.rnx']);
save([datadir 'data_o.mat'],'data_o')
load([datadir 'data_o.mat'])