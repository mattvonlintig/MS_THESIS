function [ P, Amp , Frame3 ] = f_adaptV2( F , c , stack )
% F_ADAPT - adaptive window F detector
%   Analysis of F-statstic distribution using a forced background pdf
%   method with applications towards real-time volcano monitoring using
%   infrasound microphone array
%
%   Inputs:
% 
%       F : F-statistic vector
%        c: Confidence bound cutoff (recommended 
%    stack: (Infrasound) Data time series.. can be single channel
%
%  Outputs:
%
%        P: Probability values for each small window interval, based on
%           updated background pdf
% 
% AUTHOR: 
% MATTHEW R VON LINTIG
% Boise State University
% 2/08/2018
%

tic
% c = .9995;% confidence threshold
sps = 200;% sample rate        %%%%%%%%%%%%%%%%%%
bw = sps*60*60*3; % big window %% FOR PLOTTING %%
sw = sps*60*0.5;% small window   %%%%%%%%%%%%%%%%%%
npts = floor(numel(F)/bw);
mpts = bw/sw;
kk = 0;
event_num = 0; % for counting event windows
p_obs = 0;
idx_expl = [];
M = 4;
f = reshape(F(1:npts*bw),[],npts);
s = reshape(stack(1:npts*bw),[],npts);
x = [0.1:0.01:600];% x-vector for cdf calculation
f_cdf = fcdf(x,2*2*2.5*10,2*2.5);% compute theoretical background cdf
% [NF,EP] = histcounts(F,x,'normalization','probability');
% CDFo = cumsum(NF)/sum(NF);
% position vectors for subfigs: [left,bottom,width,height] between 0-1
p1 = [0.1 0.5 0.8 0.45];
p2 = [0.1 0.1 0.35 0.35];
p3 = [0.55 0.1 0.35 0.35];
figure(97);cla;
% for ii = [1,3:24,27:npts]
%         subplot('position',p1)
%         plot([0:bw-1]/200/60,s(:,ii)-5*ii,'-k')
%         hold on;
%         ylim([-275 15])
%         title('7 Days of Vent-sourced Infrasound Activity at Sakurajima')
%         xlabel('Minutes')
% end
    [~,drift_idx1] = min(abs(f_cdf - c));% index confidence threshold
    idx = zeros(npts,mpts);
    bf = fpdf(x,2*2*2.5*10,2*2.5);% build apriori background pdf
    gf = [];

for ii = [1:15]%[3:24,27:npts] % data gaps at 2,25&26 for f2
        ff = reshape(f(:,ii),[],mpts);
    for jj = 1:mpts
%%%%%%%%% Updating Background pdf %%%%%%%%%%%%
% [NF,~] = histcounts(F(1:ii*bw),x,'normalization','probability');
% CDFo = cumsum(NF)/sum(NF);
[NF,~] = histcounts(bf,x,'normalization','probability');
CDFb = cumsum(NF)/sum(NF);

%%%%%%%%% Calculate Drift Correlction %%%%%%%%%%
[~,drift_idx2] = min(abs(CDFb - c));
drift_corr = x(drift_idx1)/x(drift_idx2);% scales observed background pdf into theoretical pdf

%%%%%%%% PLOTTING %%%%%%%%
%         clf;
% if ii == npts
        subplot('position',p1)
        cla;
        plot([0:bw-1]/200/60,s(:,ii)-5*ii,'-k')
        hold on;
        plot([jj/2 jj/2], [-5*ii-10 -5*ii+10],'-y','linewidth',3)
        ylim([-ii*5-10 -ii*5+10])
        title('Detecting Vent-sourced Infrasound Activity at Sakurajima')
%         xlabel('Minutes')
text(82,-5*ii-9,'\fontsize{18}{Minutes}')
if event_num > 0
text(82,-5*ii+8,['\fontsize{18}{RSMR = }', num2str(sum(p_obs(end-kk:end)))])
end
        xlim([0 180]);
         set(gca,'ytick',[-5*ii-10 -5*ii+10]);
         set(gca,'yticklabel',[-10 10]);
        subplot('position',p2)
        cla;
        histogram(log(bf),[-1:0.25:6],'normalization','probability','facealpha',0.9,'facecolor','r')
        hold on;
        histogram(log(gf),[-1:0.25:6],'normalization','probability','facealpha',0.9,'facecolor','g')
%         title('Categorized F Stat Probability Distribution')
        legend('Observed Background pdf','Vent-signal pdf')
        xlabel('log(F)')
        subplot('position',p3)
        cla;
        plot(x,f_cdf,'-r','linewidth',2);
        hold on;
    plot(x(1:end-1),CDFb,'linewidth',2);
%     title('Theoretical and Empirical Background F CDFs')
    xlabel('F')
    legend('Theoretical cdf','Observed cdf','location','southeast');
%     legend([h1, h2, hg, hb],{'Theoretical Background pdf','Empirical Background pdf',...
%         '> 99.95% Confidence','< 99.95% Confidence'},'Location','southeast')
    set(gca,'xscale','log')
        hold on;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ff = reshape(f(:,ii),[],mpts);
%     for jj = 1:mpts
%%%%%%%%% Index max F value w/ Drift Correction applied %%%%%%%%%%
[~,fii] = min(abs(max(ff(:,jj))/drift_corr - x(1:end-1)));
fyy = CDFb(fii);% compute confidence from updated background pdf


%         if fyy >= c % it's good!
%             gf = [gf ff(:,jj)'];
% 
%         else
%         bf = [bf ff(:,jj)'];% background F
% 
%         end

if fyy >= c % it's good!
%             plot([(jj-1)*sw+1:(jj*sw)]/200/60,s((jj-1)*sw+1:(jj*sw),ii)-5*ii,'-','color',[0.2 0.58 0.2])
            gf = [gf ff(:,jj)'];
            if kk ~= 0 % and it's still good!
            event_num = event_num+1;% only indexing new events
            L(event_num) = kk;% length of quiescence prior to event
            mAmp(event_num) = max(abs(s((jj-1)*sw+1:(jj*sw),ii)));% max amplitude of event
            if mAmp(event_num) > M 
                idx_expl(event_num) = 1;
            else
                idx_expl(event_num) = 0;
            end
            %elseif event_num > 1 && L(event_num-1) == 1
%                 L(event_num-1) = L(event_num-1)+L(event_num);
%                 L(event_num) = [];
%                 idx_expl(event_num) = [];
%                 event_num = event_num-1;
            end
            kk = 0;
            % index max amplitude over entire event lenght
            mAmp(event_num) = max(mAmp(event_num),max(abs(s((jj-1)*sw+1:(jj*sw),ii))));
            if mAmp(event_num) > M 
                idx_expl(event_num) = 1;
            else
              idx_expl(event_num) = 0;
            end
        
        else
            bf = [bf ff(:,jj)'];% background F
%             plot([(jj-1)*sw+1:(jj*sw)]/200/60,s((jj-1)*sw+1:(jj*sw),ii)-5*ii,'-','color',[0.88 0.2 0.2])

            kk = kk+1;% counting consecutive 15s closure intervals
        end

idx(ii,jj) = (fyy);% save all interval probability values
Amp(ii,jj) = max(abs(s((jj-1)*sw+1:(jj*sw),ii)));

% keep running tab of closure lengths and resulting activity

if numel(idx_expl) > 0 
xL = L(find(idx_expl == 1));% pre-explosion repose interval lengths
nL = L(find(idx_expl == 0));% non-expl repose lengths
% px = prctile(xL,[10, 25, 50, 75, 95]);
% pn = prctile(nL,[10, 25, 50, 75, 95]);
icdf_x = [0.0 diff(sort(xL)/max(xL))];
% icdf_n = 1 - sort(nL)/max(nL);
icdf_n = [0.0 diff(sort(nL)/max(nL))];
[~,Ix] = min(abs(kk - sort(xL))); 
[~,In] = min(abs(kk - sort(nL)));
% [~,Ix50] = min(abs(kk - sort(xL))); 
% [~,In50] = min(abs(kk - sort(nL)));
Ix = max(Ix);
In = max(In);
if numel(Ix) == 0
    Ix = 1;
    icdf_x = 0;
end
if numel(In) == 0
    In = 1;
    icdf_n = 0;
end
if isempty(nL)
    nL = 0;
end
% keyboard;
if isempty(xL)
    xL = 0;
    if exist('idx_expl') == 0
    idx_expl = 0;
    end
end

err1 = (abs(median(nL)-L(event_num)))/(numel(nL));
err2 = (abs(median(xL)-L(event_num)))/(numel(xL));
% p_obs(ii,jj) = (icdf_x(Ix)/icdf_n(In))^2;
p_obs(ii,jj) = (err1/err2)^2;% relative probability of explosion
mp(event_num) = sum(p_obs(end-L(event_num):end));

if mp(event_num) == Inf
    mp(event_num) = 1000;
elseif mp(event_num) == 0
    mp(event_num) = 0.01;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% idx(ii,jj) = (fyy);% save all interval probability values
% Amp(ii,jj) = max(abs(s((jj-1)*sw+1:(jj*sw),ii)));
% [~,~,butt(ii,jj)] = ginput(1);
    drawnow;
%     Frame3(ii,jj) = getframe(gcf);
    end

%    pause(0.5);% for cummulative plotting
end
    P = idx;% return interval probabilities
toc
end