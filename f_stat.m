function [F] = f_stat(data,window,overlap)

% function F_STAT(DATA,WINDOW)
%
% INPUTS:
% data - pre-time-lagged data matrix, format = data(samples,channels)
% window size [samples] (optional)
%
% OUTPUTS:
% F - Fisher Statistic time series
%
% Calculates Fisher Statisric time series from all input channels over the
% specified time window

% initialize sizes
L = max(size(data));
J = min(size(data));

% condition data
data = data - mean(data);

if nargin < 2
    window = L;
    overlap = 0;
end

ii = 1;
for n = 1:window-1:L-window+1
F(n:n+window-1) = (J-1)/J .* sum(sum(data(n:n+window-1,:).^2)) ./ sum(sum((data(n:n+window-1,:)-repmat((1./J.*sum(data(n:n+window-1,:),2)),1,J)).^2));
ii = 1+ii;
end
% for plotting:
% plot([window/2:window:L-window/2]/200/60,F,'.k')


% figure(1222);
% subplot(221)
% plot([0:L-1]/200/60/60,data)
% xlabel('Hours')
% subplot(223)
% plot([window/2:window:L-window/2]/200/60,F,'.k')


% good = find(F > 1.64);
% bad = find(F <= 1.64);
% plot([good+window/2]/200/60/60,F(good),'.g')
% hold on,
% plot([bad+window/2]/200/60/60,F(bad),'.r')
% set(gca,'yscale','log')


