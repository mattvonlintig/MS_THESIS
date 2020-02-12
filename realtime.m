%% P detections for real time Application
f_comp = f(2*200*60*60:end);
% f_comp(26*200*60*60+11*60*200:26*200*60*60+14*60*200) = [];
stack_comp = stack(2*200*60*60:end);
% stack_comp(26*200*60*60+11*60*200:26*200*60*60+14*60*200) = [];
%% video of detections
[ ~, ~ , Frame3 ] = f_adaptV2( f_comp(200*60*60*57.5:end) , 0.995 , stack_comp(200*60*60*57.5:end) );
isWriteVideo = 1;
        if isWriteVideo
            v = VideoWriter('F_detect.avi');
            v.FrameRate = 64;
            open(v);
            for ii = 1:15
                for jj = 1:4:360
                        writeVideo(v,Frame3(ii,jj));
                end
            end
            close(v);
        end

% implay('F_detect.avi');

%%
figure(29),clf;
for ii = 0:15
[ L, mAmp, p_obs, Amp, P, mp, idx_expl , drift_corr ] = fp_test( [f_comp] , 0.99999-0.00001*ii^3 , [stack_comp] , ii);
figure(33);
% subplot(212)
% cla;
% imagesc(P)
% subplot(211)
plot(drift_corr(:))
set(gca,'xscale','log')
% ylim([0 40])
hold on;
drawnow;
end

[ L, mAmp, p_obs, Amp, P, mp, idx_expl ] = fp_L( f_comp(end/4:end/2) , 0.9995 , stack_comp(end/4:end/2) );
[ L1, mAmp1, p_obs1, Amp1, P1, mp1, idx_ex1 ] = fp_adapt( f_comp , 0.9995 , stack_comp );

figure(24),clf;
scatter(L/12,mp,20,[0 1 1].*[1:numel(mp)]'./numel(mp),'+')
hold on;
scatter(L(find(idx_expl==1))/12,mp(find(idx_expl==1)),15+mAmp(find(idx_expl==1))/2,find(idx_expl==1)'.*[1 0 0]./numel(mp),'filled')
set(gca,'yscale','log')
set(gca,'xscale','log')

xdot = logspace(log(min(mp(mp~=0))),log(max(mp)),(log(max(mp))+1)*100);
    
ff = [];
TP = [];
for ii = xdot%0.3374
    ff = [ff,numel(find(idx_ex==0 & mp>ii))/(numel(find(mp>ii)))];
    TP = [TP, numel(find(idx_ex==1 & mp>ii))];
end

gg = [];
TN = [];
for ii = xdot%0.3374
    gg = [gg,numel(find(idx_ex==1 & mp>ii))/numel(find(idx_ex==1))];
    TN = [TN, numel(find(idx_ex==0 & mp<ii))];
end
% figure(2);plot(0:max(mp)/numel(mp)/50:max(mp),ff,0:max(mp)/numel(mp)/50:max(mp),gg)
% legend('N_x/N_D','N_x/N_{Tx}')

figure(4);clf;
hdot1 = plot([0 1],[0 1],':b','linewidth',2);
hold on;
for ii = 1:numel(gg)
plot([ff(ii), (ff(ii)+gg(ii))/2],[gg(ii), (ff(ii)+gg(ii))/2],'-k')
% L1(ii) = sqrt((gg(ii)-(1-ff(ii)+gg(ii)/2))^2 + (1-ff(ii) - (1-ff(ii)+gg(ii))/2)^2);
end
plot(ff,gg,'--k','linewidth',1.5)
hdot = plot(ff,gg,'ok','markerfacecolor','r','markersize',9);
hold on;
xlabel('P_{FP}')
ylabel('P_{TP}')
title('Real-time Explosion Forecasting ROC Curve from Sakurajima Data')
legend([hdot1 hdot],'1:1 Reference Line','Detection Statistics','Location','Northwest')
% title('Real-time Explosion Forecasting ROC Curve')
% maxp = max(

ROC = trapz(ff,gg);
hdot2 = plot(TN,TP,'ok','markerfacecolor','g','markersize',9);
find(gg./ff==max(gg./ff))

%% Poisson time interval statistics
jj = 0;
kk = 0;
for ii = logspace(0,3,9)
    jj = jj+1;
    for kk = 1:numel(L)-10
        kk = kk+1;
    npstat = cumsum(L(kk:end)/2) - ii;
    npstat(npstat > 0) = [];
    npval(jj,kk) = numel(find(npstat > -ii & npstat ~= 0));
    end
end
    figure(),plot(npval(9,:),'-')
    
    for ii = 1:9
        [pd1,pci1] = poissfit(npval(ii,:));
        pcdf = poisscdf(0:1:max(npval(ii,:)),pd1);
figure(3); histogram(npval(ii,:),35,'normalization','cdf')
hold on;
plot(0:1:max(npval(ii,:)),pcdf,'-')
    end
    set(gca,'xscale','log')
    legend('')
%%
% hmmm, fix the linear calculator
for ii = 1:numel(gg)
% plot([1-ff(ii), (1-ff(ii)+gg(ii))/2],[gg(ii), (1-ff(ii)+gg(ii))/2],'-r')
L1(ii) = sqrt((gg(ii)-(ff(ii)+gg(ii))/2)^2 + (ff(ii) - (ff(ii)+gg(ii))/2)^2);
end
L1 = L1';
pc = xdot(find(L1==max(L1)));


%% movie of detections

figure(5),clf;
for ii = 1:numel(mp)
scatter(L(ii)/2,mp(ii),20,[0 1 1].*[1:numel(mp(ii))]'./numel(mp),'+')
hold on;
if idx_ex(ii) == 1
scatter(L(ii)/2,mp(ii),15+mAmp(ii)/4,ii.*[1 0 0]./numel(mp),'filled')
end
set(gca,'yscale','log')
set(gca,'xscale','log')
pause(0.01);
end



ii = 100;
while ii <= numel(mp) 
figure(5),cla;
scatter(L(99+find(idx_ex(100:ii) == 0)),mp(99+find(idx_ex(100:ii)==0)),30,[0 1 1].*[1:numel(find(idx_ex(100:ii)==0))]'./(numel(mp)-100),'+')
hold on;
if ii >= min(find(idx_ex == 1))
scatter(L(99+find(idx_ex(100:ii)==1)),mp(99+find(idx_ex(100:ii)==1)),25+mAmp(99+find(idx_ex(100:ii)==1))/4,find(idx_ex(100:ii)==1)'.*[1 0 0]./(numel(mp)-100),'filled')
% colorbar(find(idx_ex(1:ii)==1)'.*[1 0 0]./numel(mp))
end
title({'Real-time Explosion Forecasting','using Relative Squared Median Residual'})
xlabel('Current Closure Length [min]')
ylabel('RSMR')
set(gca,'yscale','log')
set(gca,'xscale','log')
legend('Non-explosion Events','Explosion Events','Location','northwest')
grid on;
Frame2(ii-99) = getframe(gcf);
ii = ii+1;
end
     

     isWriteVideo = 1;
        if isWriteVideo
            v = VideoWriter('SG_movie.avi');
            v.FrameRate = 16;
            open(v);
            for ii = 1:357
%                 for jj = -10:ii+12
%                     if jj>ii
                        writeVideo(v,Frame(ii));
%                     elseif jj<2
%                         writeVideo(v,Frame(ii));
%                     else
%                     writeVideo(v,Frame(jj));
%                     end
%                 end
            end
            close(v);
        end

figure();
implay('median_detection.avi',100);