function ret=greedy_red_syn_th(Am,Su,ii,q,eps,measure)

switch measure
    case 'iir'
        ret=iir_syn_red_th(Am,Su,ii,q,eps);
    case 'rsi'
        ret=rsi_syn_red_th(Am,Su,ii,q,eps);
    case 'oir'
        ret=oir_syn_red_th(Am,Su,ii,q,eps);
    otherwise
        warning('Unrecognized measures. Please choose between iir, rsi or oir.');
end

end