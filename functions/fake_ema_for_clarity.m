function[ema] = fake_ema_for_clarity(type)

switch type
    case 'in-ear'
        error('in-ear isn''t a valid array type')
    case 'bte1'       
        ema.channelsLeft = [1];
        ema.channelsRight = [2];
    case 'bte2'
        ema.channelsLeft = [1 3];
        ema.channelsRight = [2 4];
    case 'bte3'
        ema.channelsLeft = [1 3 5];
        ema.channelsRight = [2 4 6];
    otherwise
        error('Unknown type')
end

ema.refChan = 1;
ema.isBinaural = 1;
ema.refChanLeft = 1;
ema.refChanRight = 2;
