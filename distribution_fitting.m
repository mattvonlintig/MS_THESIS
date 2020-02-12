

% contact


%% break the data apart into subsets
% test each set to assess confidence
tic
clear A B MU CI1 pci2;
jmax = 10;
figure(501),clf;
% for jj = 1
% jj = 1;
% a = contact(2:2:end,4)/200/60;% for repose intervals
% a = cll;
a = L/12;
% a = contact(2:2:end,3);% for event magnitudes
b = a(a>0);

% chi-square test
clear r;
[N1,E1] = histcounts((b),logspace(-2,6,250),'normalization','probability');% binned data
parfor ii = 1:10000
    r= randsample((b),250,'true');
    [pd,pci] = expfit(r);
    MU(ii) = pd;
    CI3(:,ii) = pci;
    [pd2,pci2] = gamfit(r,0.1);
    A(ii) = pd2(1);
    B(ii) = pd2(2);
%     R(:,ii) = r;
    CI2(:,ii) = pci2(:);
%     [pd1,pci1] = poissfit(r);
%     LAMBDA(ii) = pd1;
%     CI1(:,ii) = pci1;
    [pd1,pci1] = wblfit(r,0.001);
    AA(ii) = pd1(1);
    BB(ii) = pd1(2);
    CI1(:,ii) = pci1(:);
%     SIGMA(ii) = pd.sigma;
% [NP(ii,:),EP] = histcounts(log(r(:,ii)),linspace(-1,4.75,20),'normalization','probability');
% X2d(ii) = sum((N1-NP(ii,:)).^2./NP(ii,:));
end
% [NP,EP] = histcounts(log(r),linspace(-1,5,20),'normalization','probability');
% observed cdf
CDFo = cumsum(N1)/sum(N1);
% [H,P] = kstest2(log(r),log(b),'alpha',0.05)

    % position vector [left,bottom,width,height] between 0-1
p1 = [0.7125 0.6 0.225 0.3125];
p2 = [0.7125 0.125 0.225 0.3125];
p3 = [0.1 0.6 0.8 0.3125];
p4 = [0.1 0.125 0.8 0.3125];
% p2_5 = [0.835 0.125 0.125 0.3125];
% display best fit parameters MU & SIGMA
figure(501),hold on;%,clf;

% subplot('position',p2)
% histogram(A,15),hold on;
% plot([median(A) median(A)],[0 250],'-g')
% subplot('position',p1)
% histogram(B,15),hold on;
% plot([median(B) median(B)],[0 250],'-g')

sigma = num2str(median(A));
title(['Monte Carlo realizations'])
% [phat2,pci2] = expfit(log(cll));

% assess misfit between observed data and optimal parameter distribution
pd = makedist('exp','mu',median(MU));
pd2 = makedist('gamma','a',median(A),'b',median(B));
% pd1 = makedist('Poisson','lambda',median(LAMBDA));
pd1 = makedist('wbl','a',median(AA),'b',median(BB));
clear NP r;
parfor ii = 1:10000
    r = random(pd2,[750,1]);
    [NP(ii,:),EP] = histcounts((r),logspace(-2,6,250),'normalization','probability');
    NPnan = NP(ii,:);
    NPnan(NPnan == 0) = 0.0001;
    NP(ii,:) = NPnan;
    X2d(ii) = sum((N1-NP(ii,:)).^2./NP(ii,:));
    [H(ii),Pval(ii),kstat(ii)] = kstest2(log(r),log(b),'alpha',0.001);
%     [Ad,pval,adstat,cvad] = adtest(log(b),'distribution',pd2,'Asymptotic',0,'alpha',0.005)
end

np = median(NP,1);
chisquare = sum((N1-np).^2./np);

% expected cdf for gamma dist
% CDFe = cumsum(np)/sum(np);
CDFe = 1 - gamcdf(logspace(-2,6,249),median(A),median(B));

% expected cdf for exponential
clear r NP;
[N1,E1] = histcounts(round(log(b)),logspace(-2,6,250),'normalization','probability');% binned data
parfor ii = 1:10000
    r = random(pd1,[750,1]);
    [NP(ii,:),EP] = histcounts((r),logspace(-2,6,250),'normalization','probability');
    NPnan = NP(ii,:);
    NPnan(NPnan == 0) = 0.0001;
    NP(ii,:) = NPnan;
    X2d(ii) = sum((N1-NP(ii,:)).^2./NP(ii,:));
    [H(ii),Pval(ii),kstat(ii)] = kstest2(log(r),log(b),'alpha',0.01);
%     [Ad,pval,adstat,cvad] = adtest(log(b),'distribution',pd1,'Asymptotic',1,'alpha',0.01)
end
np = median(NP,1);
% expected cdf for exp dist
% CDFe2 = cumsum(np)/sum(np);

% [phat,pci] = gamfit(log(cll),.01);
CDFe3 = 1 - expcdf(logspace(-2,6,249),median(MU));
% CDFe2 = 1 - poisscdf(linspace(-1,6,249),median(LAMBDA));
CDFe2 = 1 - wblcdf(logspace(-2,6,249),median(AA),median(BB));

cdfmax1 = 1 - gamcdf(logspace(-2,6,249),median(CI2(2,:)),median(CI2(4,:)));
cdfmin1 = 1 - gamcdf(logspace(-2,6,249),median(CI2(1,:)),median(CI2(3,:)));

cdfmax3 = 1 - expcdf(logspace(-2,6,249),median(CI3(1,:)));
cdfmin3 = 1 - expcdf(logspace(-2,6,249),median(CI3(2,:)));
% cdfmax2 = 1 - poisscdf(linspace(-1,6,249),median(CI1(1,:))-std(CI1(1,:)));
% cdfmin2 = 1 - poisscdf(linspace(-1,6,249),median(CI1(2,:))+std(CI1(2,:)));
cdfmax2 = 1 - wblcdf(logspace(-2,6,249),median(CI1(2,:)),median(CI1(4,:)));
cdfmin2 = 1 - wblcdf(logspace(-2,6,249),median(CI1(1,:)),median(CI1(3,:)));

% subplot('position',p3)
% errorbars
% h1 = bar(linspace(-1,6,249),1 - CDFo,1,'facecolor','k');hold on;
ha = shadedErrorBarT8(logspace(-2,6,249),1 - cumsum(fx)/sum(fx),[cdfmax1; cdfmin1],1,{'color',[0.0 0.85 0.85]});
freezeColors(gca);
hold on;
h2 = plot(logspace(-2,6,249),CDFe,':','color',[0.0 0.25 0.85],'linewidth',3);
% freezeColors;
hb = shadedErrorBarT8(logspace(-2,6,249),CDFe2,[cdfmax2; cdfmin2],1,{'color',[0.85 0.65 0.0]});
freezeColors(gca);
hc = shadedErrorBarT8(logspace(-2,6,249),CDFe3,[cdfmax3; cdfmin3],1,{'color',[0.25 0.65 0.15]});
h4 = plot((logspace(-2,6,249)),CDFe2,'-.','color',[0.85 0.25 0.0],'linewidth',3);
h11 = plot((logspace(-2,6,249)),CDFe3,'--','color',[0.25 0.85 0.15],'linewidth',3);
h1 = plot(logspace(-2,6,249),1 - CDFo,'-k','linewidth',2);
% set(gca, 'Layer', 'bottom')
uistack(h1,'top')
uistack(h4,'up')
set(gca,'xscale','log')
xlim(exp([-4 3.5]))
% title({['Repose CDF'], ['\chi^2 = ' num2str(chisquare)]})
% text(2,0.5,{'Kolmogorov','Smirnov','p-value',num2str(P)})
title(['ICDF of ', num2str(numel(b)), ' Determined Repose Interval Lengths'])
ylabel('1 - Cumulative probability')
xlabel('Time')
legend([h1 h4 h2 h11],'Repose Interval Data','Weibull Fit \pm \sigma','Gamma Fit \pm \sigma','Exponential Fit \pm \sigma','Location','northeast')

numex = histcounts(b(find(idx_expl(1:numel(b))==1)),logspace(-2,6,15),'normalization','pdf')*(numel(idx_expl)/numel(b));
% X = 0.01:0.005:100;
% X = 0:0.0282:7;
X = logspace(-2,6,249);
X2 = 0:6;
figure(696),clf;
% subplot('position',p4)
% histogram([],[0 7])
hh = histogram(log(b),logspace(-2,2,15),'normalization','pdf','facecolor',[0 0 0],'facealpha',0.8);
hold on;
fx = 1./(median(B).^median(A).*gamma(median(A))).*X.^(median(A)-1).*exp(-X./median(B));
fx3 = 1/median(MU).*exp(-X./median(MU));% exponential
% fx2 = median(LAMBDA).^X2./(factorial(X2)).*exp(-median(LAMBDA));% poisson
% fx2 = random(pd1,[249,1]);
CC = 0;
fx2 = median(BB)/median(AA).*((X-median(CC))./median(AA)).^(median(BB)-1).*exp(-((X-median(CC))./median(AA)).^median(BB));
% fx = 1./(X.*pd.sigma*sqrt(2*pi)).*exp(-0.5*(log(X./median(b))./pd.sigma).^2);
bb = bar(hh.BinEdges(1:end-1)+(hh.BinEdges(2)-hh.BinEdges(1))/2,numex,1,'facecolor',[0.65 0.65 0.65]);
% histogram(log(cll(idx_expl)),linspace(1,7,15),'normalization','pdf')
h5 = plot(X,fx,'-.','color',[0.0 0.45 0.95],'linewidth',3);
h6 = plot(X,fx2,'-.','color',[0.85 0.25 0.0],'linewidth',3);
h7 = plot(X,fx3,'--','color',[0.25 0.65 0.15],'linewidth',4);
% h6 = histogram(fx2,X,'normalization','probability');
xlim(exp([-2 2]))
ylim([0 0.45])
title(['PDF for ', num2str(numel(b)), ' Determined Repose Interval Lengths'])
ylabel('Probability')
xlabel('Time')
legend([hh, bb, h6, h5, h7],'All intervals','Pre-explosion intervals','Weibull Fit','Gamma Fit','Exponential Fit','Location','northeast')
uistack(h5,'top')
set(gca,'xscale','log')
% figure(69),clf;
% histogram(log(b),linspace(-1,4.75,20),'normalization','probability')
% hold on;
% histogram(log(r),linspace(-1,4.75,20),'normalization','probability')

% test the null hypothesis that the data are from the same distribution
% plot the cdf's for visual comparison
% Y = pdf(pd,logspace(-2,5,100));
% figure();
% plot(linspace(-2,5,100),Y)
% figure(101);
% histogram(log(b),15,'normalization','probability')
% hold on;
% end
toc
%%

figure(),histogram(X2d,20,'normalization','probability')
title('\chi^2 sample distribution')




%% a model for eruption probability

% 7 day observation period
t = 7*24*60;  % [days]?
a = .010;% Omega/k --> some constant ratio
x = (a*t)/(1+a*t);
y = log(1/(1-x));
omega = 4/t;

X = 0.01:0.01:100;
Y = -log(1-X);

n = 0.1:0.1:100;
Pn = x.^n./(n.*y);
Pn = .99;

% figure(),
plot(log(n),Pn)




%% Variogram?
% figure(),clf;
for ii = 1:290
    clear p;
    for jj = 1:291-ii
    p(jj) = var(b(jj:ii:ii+jj));
%     v(ii,jj) = var(p);
    end
    v(ii) = mean(p);
    n(ii) = numel(p);
%     plot(ii,v,'.')
%     hold on;
end
figure();
scatter(1:ii,v/2,n)


%% computing expected statistics from models w/ error bounds

% computing likelihood of expl eruption
figure(987654),hh = histogram(log(b),linspace(-1,7,25),'normalization','cdf','facecolor',[0 0 0],'facealpha',0.8);
cll2 = cll;
cll2(idx_expl) = [];
ncll2 = histcounts(log(cll2),linspace(-1,7,25),'normalization','cdf');%*(numel(idx_expl)/numel(cll));
% hold on;
% histogram(log(cll2),linspace(0,7,15),'normalization','cdf','facecolor',[0.7 0.7 0],'facealpha',0.8);hold on;
figure(23),clf;
bar(hh.BinEdges(1:end-1)+(hh.BinEdges(2)-hh.BinEdges(1))/2,1 - (ncll2),1,'facecolor',[0.7 0.7 0],'facealpha',0.98);hold on;
numex = histcounts(log(cll(idx_expl)),linspace(0,7,25),'normalization','cdf');%*(numel(idx_expl)/numel(cll));
bar(hh.BinEdges(1:end-1)+(hh.BinEdges(2)-hh.BinEdges(1))/2,(numex),1,'facecolor',[0.25 0.78 0.15],'facealpha',0.8);
grid on;
title({'Comparison of Repose Intervals','Resulting in Explosive and Non-Explosive Activity'})
xlabel('Time')
ylabel('Observed Probability')
pr = (numex)./ncll2./[0,diff(numex)]/sum((numex))/2;
% pr = pr./sum(pr);
% pr(18:end) = 1;
plot(hh.BinEdges(1:end-1)+(hh.BinEdges(2)-hh.BinEdges(1))/2,pr,'-.k','linewidth',2)
legend('Non-Explosive','Explosive','Empirical Expl Likelihood','Location','southeast')
set(gca,'xTick',[0 1.1 2.3 3.4 4.5 5.6 6.59])
set(gca,'xTickLabel',{[num2str(1) ' min'], [num2str(3) ' min' ], [num2str(10) ' min' ], [num2str(30) ' min' ], [num2str(1.5) ' hr' ], [num2str(4.5) ' hr' ], [num2str(12) ' hr' ]})
jj = 0;
for ii = 0:0.25:7
    jj = jj+1;
    numa(jj) = numel(find(log(cll(idx_expl)) >= ii));% explosions
    numb(jj) = numel(find(log(cll2) >= ii));% not explosions
%     numag(jj) = 
end
ntot = numel(numa)+numel(numb);
% plot((0:0.25:ceil(max(log(cll)))),(numa./numb)/ntot,'-k')

% stair(hh.BinEdges(1:end-1)+(hh.BinEdges(2)-hh.BinEdges(1))/2,((25-cumsum(numex))./(numel(cll2)-cumsum(ncll2))-1)./2,'--k')
%%

figure(66);clf;
stairs(0:0.25:7,((numa)./numel(idx_expl)),'-.r','linewidth',2), hold on;
stairs(0:0.25:7,((numb)./numel(cll2)),':b','linewidth',2)
ylim([0 1.001])
% yyaxis right;
booyah = ((numa))./((numb));
idxnan = isnan(booyah);
booyah(idxnan) = 0;
stairs(0:0.25:5,booyah(1:21)*numel(cll2)/numel(cll),'pg','linewidth',1,'markersize',24,'markerfacecolor','g')
for ii = 1: 21
text(0.25*(ii-1)-0.05,booyah(ii)*numel(cll2)/numel(cll),[num2str(round(100*booyah(ii)*numel(cll2)/numel(cll))), '%'],'fontsize',7)
end
title('Survival Functions and Explosion Likelihood')
% plot([0 0.55],[0.1435 0.1435],':k','linewidth',2)
% text(0.15,0.25,'$\frac{N_{Expl}}{N_{Total}}$','interpreter','latex','fontsize',12)
% text(0.1,0.155,'$\frac{N_{Expl}}{N_{Total}} = 8\%$','interpreter','latex','fontsize',12)
legend({'Explosion Survival','Emergent Degassing Survival',['Emperical Expl. Likelihood ','  $\frac{N_{Expl}}{N_{Total}}$']},'interpreter','latex')
set(gca,'xTick',[0 1.1 2.3 3.4 4.5 5.6 6.59])
set(gca,'xTickLabel',{[num2str(1) ' min'], [num2str(3) ' min' ], [num2str(10) ' min' ], [num2str(30) ' min' ], [num2str(1.5) ' hr' ], [num2str(4.5) ' hr' ], [num2str(12) ' hr' ]})
xlim([0 6])
xlabel('Time')
ylabel('Probability')


