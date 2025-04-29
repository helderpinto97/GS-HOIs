function ret=greedy_red_syn_est(Y,p,ii,q,numsurr,alpha,measure)

switch measure
    case 'iir'
        ret=iir_syn_red_est(Y,p,ii,q,numsurr,alpha);
    case 'rsi'
        ret=rsi_syn_red_est(Y,p,ii,q,numsurr,alpha);
    case 'oir'
        ret=oir_syn_red_est(Y,p,ii,q,numsurr,alpha);
    otherwise
        warning('Unrecognized measures. Please choose between iir, rsi or oir.');
end

end