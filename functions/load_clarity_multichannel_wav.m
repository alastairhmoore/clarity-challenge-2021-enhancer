function[x,fs] = load_clarity_multichannel_wav(stub,type)

switch type
    case 'in-ear'
        suffix_list = {'CH0'};
    case 'bte1'
        suffix_list = {'CH1'};
    case 'bte2'
        suffix_list = {'CH1','CH3'};
    case 'bte3'
        suffix_list = {'CH1','CH2','CH3'};
    otherwise
        error('Unknown type')
end

x = [];
for ipair = 1:length(suffix_list)
    filepath = [stub '_'  suffix_list{ipair} '.wav'];
    [x_in,fs] = v_readwav(filepath);
    x = [x, x_in];
end

