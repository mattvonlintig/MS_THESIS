function [cf, e_idx, c_idx, ic] = MF_PickerV2(F,CF)
%  function MF_PICKER automates event picking using xcorr scores and sta/lta
%  INPUTS:
%  
%  stacka - single channel or stack
%  xcn - normailized xc-scores
%  F - F-stat vector
%  CF - F-stat cutoff [low, high]
%  
%  
%  OUTPUTS:
%  
%  idx - indices of detections
%  val - associated peak value
tic

e_idx = cell(numel(CF),1);
c_idx = cell(numel(CF),1);
% for ll = 1:numel(CF)
ll = 1;

[idx1] = find(log(F) > CF(ll)-1);% 1.5 
[idx2] = find(log(F) < CF(ll)+1); % 0.85
ax = find(diff(fliplr(idx1)) ~= -1);
ay = find(diff(idx1) ~= 1);
ay = sort([ay ax]);
axx = find(diff(fliplr(idx2)) ~= -1);
ayy = find(diff(idx2) ~= 1);
ayy = sort([ayy axx]);
idx2 = sort([idx2(ayy(ayy>1)) idx1(ay(ay>1))]);
idxb = idx2; % closures
icc = [];


% use the percentile here as a tolerance for closedness:
% if the X-th %ile is too high, the sample period is removed
ii = numel(idxb)-1;
while ii >= 2
    ii = ii-1;
    if (prctile(log(F(idxb(ii):idx2(ii+1))),90) > CF(ll)+1 ) || idx2(ii+1)-idxb(ii) < 60*200 %|| rms(stack(idxb(ii):idx2(ii+1))) > 0.008
        idxb(ii) = [];
        idx2(ii+1) = [];
    end
end

% use the percentile here as a tolerance for openness:
% if the X-th %ile is too low, the sample period is removed
ii = 0;
while ii < numel(idx2)-1
    ii = ii+1;
    if (prctile(log(F(idx2(ii):idxb(ii))),80) < CF(ll)-1 ) || idxb(ii)-idx2(ii) < 20*200 %|| rms(stack(idxb(ii):idx2(ii+1))) > 0.008
        idxb(ii) = [];
        idx2(ii) = [];
    end
end

% use the percentile here as a tolerance for closedness:
% if the X-th %ile is too high the sample period is removed
ii = 0;
while ii < numel(idx2)-1
    ii = ii+1;
    if (prctile(log(F(idxb(ii):idx2(ii+1))),70) > CF(ll)+1 ) 
        idxb(ii) = [];
        idx2(ii+1) = [];
    end
end

% use the percentile here as a tolerance for openness:
% if the X-th %ile is too low, the sample period is removed
ii = 0;
while ii < numel(idx2)-1
    ii = ii+1;
    if (prctile(log(F(idx2(ii):idxb(ii))),99) < CF(ll)-1 ) || idxb(ii)-idx2(ii) < 20*200 %|| rms(stack(idxb(ii):idx2(ii+1))) > 0.008
        idxb(ii) = [];
        idx2(ii) = [];
    end
end

% use the percentile here as a tolerance for closedness:
% if the X-th %ile is too high, the sample period is removed
ii = numel(idxb)-1;
while ii >= 2
    ii = ii-1;
    if (prctile(log(F(idxb(ii):idx2(ii+1))),90) > CF(ll)+1 ) || idx2(ii+1)-idxb(ii) < 60*200 %|| rms(stack(idxb(ii):idx2(ii+1))) > 0.008
        idxb(ii) = [];
        idx2(ii+1) = [];
    end
end

ii = 0;
while ii < numel(idx2)-1
    ii = ii+1;
    if (prctile(log(F(idx2(ii):idxb(ii))),99) < CF(ll)-1 ) || idxb(ii)-idx2(ii) < 20*200 %|| rms(stack(idxb(ii):idx2(ii+1))) > 0.008
        idxb(ii) = [];
        idx2(ii) = [];
    end
end

ax = find(idxb - idx2 == 0);
idxb(ax) = [];
idx2(ax) = [];

% index correlated noise intervals

% need to change this so that explosive intervals are not counted in the
% i^th count
ii = numel(idxb);
while ii >= 2
    ii = ii-1;
    if prctile(log(F(idxb(ii):idx2(ii+1))),25) > 0.75 && mean(log(F(idxb(ii):idx2(ii+1)))) > 1% prctile(log(F(idxb(ii):idx2(ii+1))),25) > 0.5
        icc(ii) = 1;
    end
    if (prctile(log(F(idx2(ii):idxb(ii))),95) < 3.5 && mean(log(F(idx2(ii):idxb(ii)))) < 3) %|| max(stack(idx2(ii):idxb(ii))) < 0.99
        iee(ii) = 1;
    end
end
ic{1} = icc(2:end);
ic{2} = iee(2:end);




e_idx{ll} = idx2(2:end)/200/60;
c_idx{ll} = idxb(2:end)/200/60;


cf(ll) = CF(ll);
% end


toc
end




