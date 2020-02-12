function [ L , az  ] = f_detect( c , w , data )
% F_ADAPT - adaptive window F detector
%   (Infrasound) Signal detection using adaptive F statistic thresholding 
%
%   Inputs:
% 
%       w : window size for F statistic calculation..........
%       c : Confidence bound cutoff (0.95, 0.9, 0.75... etc.) 
%    data : Multichannel (Infrasound) time series data.......
%
%  Outputs:
%
%       L : Length of each detected event
%       R : Repose length between detected events
%
% AUTHOR: 
% MATTHEW R VON LINTIG
% Boise State University
% 6/21/2018
%

tic
% c = .9995;% confidence threshold
F = f_stat(data,w);
sps = 200;% sample rate        %%%%%%%%%%%%%%%%%%
bw = sps*60*60*3; % big window %% FOR PLOTTING %%
sw = sps*60*1;% small window   %%%%%%%%%%%%%%%%%%
npts = floor(numel(F)/bw);
mpts = bw/sw;
f = reshape(F(1:npts*bw),[],npts);
s = reshape(data(1:npts*bw),[],npts);
x = [0.1:0.01:600];% x-vector for cdf calculation
f_cdf = fcdf(x,2*2*2.5*10,2*2.5);% compute theoretical background cdf
% [NF,EP] = histcounts(F,x,'normalization','probability');
% CDFo = cumsum(NF)/sum(NF);
% position vectors for subfigs: [left,bottom,width,height] between 0-1
% p1 = [0.1 0.5 0.8 0.45];
% p2 = [0.1 0.1 0.35 0.35];
% p3 = [0.55 0.1 0.35 0.35];
% figure(97);cla;
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
    
for ii = [3:24,27:npts-2] % data gaps at 2,25&26 for f2
%%%%%%%%% Updating Background pdf %%%%%%%%%%%%
% [NF,~] = histcounts(F(1:ii*bw),x,'normalization','probability');
% CDFo = cumsum(NF)/sum(NF);
[NF,~] = histcounts(bf,x,'normalization','probability');
CDFb = cumsum(NF)/sum(NF);

%%%%%%%%% Calculate Drift Correction %%%%%%%%%%
[~,drift_idx2] = min(abs(CDFb - c));
drift_corr = x(drift_idx1)/x(drift_idx2);% scales observed background pdf into theoretical pdf

%%%%%%%% PLOTTING %%%%%%%%
%         clf;
% if ii == npts
%         subplot('position',p1)
% %         plot([0:bw-1]/200/60,s(:,ii)-5*ii,'-k')
%         ylim([-ii*5-5 -5])
%         title('7 Days of Vent-sourced Infrasound Activity at Sakurajima')
%         xlabel('Minutes')
%         xlim([0 180]);
%          set(gca,'ytick',[]);
%         subplot('position',p2)
%         cla;
%         histogram(log(bf),[-1:0.25:6],'normalization','probability','facealpha',0.9,'facecolor','r')
%         hold on;
%         histogram(log(gf),[-1:0.25:6],'normalization','probability','facealpha',0.9,'facecolor','g')
% %         title('Categorized F Stat Probability Distribution')
%         legend('Observed Background pdf','Vent-signal pdf')
%         xlabel('log(F)')
%         subplot('position',p3)
%         cla;
%         plot(x,f_cdf,'-r','linewidth',2);
%         hold on;
%     plot(x(1:end-1),CDFb,'linewidth',2);
% %     title('Theoretical and Empirical Background F CDFs')
%     xlabel('F')
%     legend('Theoretical Background pdf','Empirical Background pdf','location','southeast');
% %     legend([h1, h2, hg, hb],{'Theoretical Background pdf','Empirical Background pdf',...
% %         '> 99.95% Confidence','< 99.95% Confidence'},'Location','southeast')
%     set(gca,'xscale','log')
%         hold on;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%
    ff = reshape(f(:,ii),[],mpts);
    for jj = 1:mpts
%%%%%%%%% Index max F value w/ Drift Correction applied %%%%%%%%%%
[~,fii] = min(abs(max(ff(:,jj))/drift_corr - x(1:end)));
fyy = CDFb(fii);% compute confidence from updated background pdf

%%%%%%%% PLOTTING %%%%%%%%%%
%         subplot('position',p1)
%         hold on;
% subplot('position',p3)
%         cla;
%         plot(x,f_cdf,'-r','linewidth',2);
%         hold on;
%     plot(x(1:end-1),CDFb,'linewidth',2);
%     title('Theoretical and Empirical Background F CDFs')
%     xlabel('F Score')
%     set(gca,'xscale','log')
%         hold on;
        if fyy >= c % it's good!
%             subplot('position',p1),hold on;
%             plot([(jj-1)*sw+1:(jj*sw)]/200/60,s((jj-1)*sw+1:(jj*sw),ii)-5*ii,'-','color',[0.1 0.77 0.2])
            gf = [gf ff(:,jj)'];
%             xlim([0 180]);
%             set(gca,'ytick',[]);
% %             if ii == npts-1
%         subplot('position',p3)
% %         subplot('position',p3)
%         cla;
%         plot(x,f_cdf,'-r','linewidth',2);
%         hold on;
%     plot(x(1:end-1),CDFb,'linewidth',2);
% %     title('Theoretical and Empirical Background F CDFs')
%     xlabel('F')
%     legend([h1, h2, hg, hb],{'Theoretical Background pdf','Empirical Background pdf',...
%         '> 99.95% Confidence','< 99.95% Confidence'},'Location','southeast')
%     set(gca,'xscale','log')
% %         hold on;
%         hold on;
%         plot(max(ff(:,jj))/drift_corr,fyy,'pk','markerfacecolor','g','markersize',8);
%         legend('Theoretical Background pdf','Empirical Background pdf','location','southeast');
%             end
        %             pt(k) = ii*jj;
        else
        bf = [bf ff(:,jj)'];% background F
%             subplot('position',p1),hold on;
%             plot([(jj-1)*sw+1:(jj*sw)]/200/60,s((jj-1)*sw+1:(jj*sw),ii)-5*ii,'-','color',[1 0 0.2])
%              set(gca,'ytick',[]);
% %         if ii == npts-2
%             subplot('position',p3)
%             subplot('position',p3)
%         cla;
%         plot(x,f_cdf,'-r','linewidth',2);
%         hold on;
%     plot(x(1:end-1),CDFb,'linewidth',2);
% %     title('Theoretical and Empirical Background F CDFs')
%     xlabel('F')
% %     legend([h1, h2, hg, hb],{'Theoretical Background pdf','Empirical Background pdf',...
% %         '> 99.95% Confidence','< 99.95% Confidence'},'Location','southeast')
%     set(gca,'xscale','log')
%         hold on;
% %         hold on;
%         plot(max(ff(:,jj))/drift_corr,fyy,'sk','markerfacecolor','r','markersize',8);
%         legend('Theoretical Background pdf','Empirical Background pdf','location','southeast');
%         end
        end
%         subplot('position',p2)
%         cla;
%         histogram(log(bf),[-1:0.25:6],'normalization','probability','facealpha',0.9,'facecolor','r')
%         hold on;
%         histogram(log(gf),[-1:0.25:6],'normalization','probability','facealpha',0.9,'facecolor','g')
% %         title('Categorized F Stat Probability Distribution')
%         legend('Observed Background pdf','Vent-signal pdf')
        
%                 pause(0.05);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx(ii,jj) = (fyy);% save all interval probability values
% Amp(ii,jj) = max(abs(s((jj-1)*sw+1:(jj*sw),ii)));
% [~,~,butt(ii,jj)] = ginput(1);
%     drawnow;
%     Frame3(jj) = getframe(gcf);
    end

%    pause(0.5);% for cummulative plotting
end
    P = idx;% return interval probabilities
toc
end