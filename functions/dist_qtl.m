function dist=dist_qtl(y)
    qtl=quantile(y,[0.95 0.05]);
    dist=[qtl(1,:)-median(y); median(y)-qtl(2,:)];
end