%% generate surrogate time series with sample shuffling (destroys any temporal correlation, may maintain static cross-correlations) 

% x: multivariate time series, N datapoints x M time series
% typeshuf: for multivariate data, apply the same ('same') or independent ('ind') random permutation to the columns of x

function [xs,p]=surr_shuf(x,typeshuf)

if nargin<2, typeshuf='same'; end

[N,M]=size(x);
if M==1, typeshuf='same'; end

xs=zeros(N,M); p=xs;
switch typeshuf
    case 'same'
        p=randperm(N)';
        xs=x(p,:);
    case 'ind'
        for m=1:M
            p(:,m)=randperm(N)';
            xs(:,m)=x(p(:,m),m);
        end
end

end

% xs2=zeros(sx(1),1);
% for k = 1:sx(1)
% 	xs2(k)=x(p(k));
% end
% plot(xs,'k'); hold on; plot(xs2,'r--');