function [stack,temp_xc] = xc_helicorder(data1,data2,tLength,sps,hz,scale_opt,w,data3,data4,data5,time_vect,lag1,lag2,lag3,lag4)

% function SPECTOHELICORDER displays recording, spectrogram, and xc-vals for
% lengthy time series' as multi-line figure according to input parameters:
% 
% *this function calls another custom function XC_NORM.m*
% 
% INPUTS: 
% data1,data2 - 2 seperate channels/receivers
% tLength     - time length per line [Minutes]
% sps         - samples per second [Hz]
% hz          - vector of frequencies [Hz] to calculate power spectra at,
%               requires at least 2 inputs... i.e. [5 10]
% scale_opt   - factor to scale data by for plotting purposes
% w           - xc window size [seconds]
% 
% OUTPUTS: 
% stack       - entire dataset stack in tLength minute lines
% temp_xc     - max correlation lags calculated for each time interval (line)
% times       - quiescent times preceding activity onset
% pk          - max peak associated with activity onset
% (figure)    - Plots a helicorder style spectrogram influenced by xc_norm vals

% TO DO:
% 
% put arrows at explosions

tic
% nr = min(size(data)); % assumes number of channels does not exceed number of data points
nlines = floor(numel(data1)/(tLength*60*sps)); % number of lines in helicorder plot
data11 = reshape(data1(1:nlines*tLength*60*sps),tLength*60*sps,[]); % reshapes data for helicorder-style plot
t = reshape(time_vect(1:nlines*tLength*60*sps),tLength*60*sps,[]);
% data22 = reshape(data2(1:nlines*tLength*60*sps),tLength*60*sps,[]);
% figure(1234),clf;
% plot(0,-3) % initializes plot for text purposes
i_begin = 2; % DEFAULT = 1
i_end = nlines-2; %   DEFAULT = NLINES
    for i = i_begin:i_end
% use these if you do not have assumed lag vectors (w.r.t. channel 1)
% % %     if i > 165
    [lag4(i)] = xc_lag(data11(:,i),data5((i-1)*sps*60*tLength+1:i*sps*60*tLength),w*sps);
    [lag1(i)] = xc_lag(data11(:,i),data2((i-1)*sps*60*tLength+1:i*sps*60*tLength),w*sps);
    [lag3(i)] = xc_lag(data11(:,i),data4((i-1)*sps*60*tLength+1:i*sps*60*tLength),w*sps);
    [lag2(i)] = xc_lag(data11(:,i),data3((i-1)*sps*60*tLength+1:i*sps*60*tLength),w*sps);

    if abs(lag1(i)) > 1e+05 
        if i-1 == 0
            lag1(i) = 0;
        else
        lag1(i) = lag1(i-1);
        end
    elseif abs(lag2(i)) > 1e+05
        if i-1 == 0
            lag2(i) = 0;
        else
        lag2(i) = lag2(i-1);
        end
    elseif abs(lag3(i)) > 1e+05 
        if i-1 == 0
            lag3(i) = 0;
        else
        lag3(i) = lag3(i-1);
        end
    elseif abs(lag4(i)) > 1e+05 
        if i-1 == 0
            lag4(i) = 0;
        else
        lag4(i) = lag4(i-1);
        end
    end
% %     end

%     [~,lag3] = xc_norm(data22(:,i),data33(:,i),15*sps);
%     lags = [lag1 lag2 lag3];
%     if lag1 < 0 && lag2 < 0
%     stack(:,i) = (data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)) + data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)) + data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)) + data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)) + data11(:,i))./5;
% %     stack = reshape((data2(-lag1+1:nlines*tLength*60*sps-lag1)+data3(-lag2+1:nlines*tLength*60*sps-lag2)+data11(:))./3,tLength*60*sps,[]);
%     [temp_xc1,~] = xc_norm(data1((i-1)*sps*60*tLength+1:i*sps*60*tLength),stack(:,i),w*sps,0); % uses 15 second window
%     [temp_xc2,~] = xc_norm(data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)),stack(:,i),w*sps,0);
%     [temp_xc3,~] = xc_norm(data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)),stack(:,i),w*sps,0);
%     [temp_xc4,~] = xc_norm(data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)),stack(:,i),w*sps,0);
%     [temp_xc5,~] = xc_norm(data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)),stack(:,i),w*sps,0);
% %     [temp_xc2,temp_lag(i-i_begin+1)] = xc_norm(data2((i-1)*sps*60*tLength - lag1+1:i*sps*60*tLength - lag1),stack(:,i),w*sps);
%     temp_xc(:,i) = (temp_xc1 + temp_xc2 + temp_xc3 + temp_xc4 + temp_xc5)./5;

        stack1 = (data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)) + data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)) + data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)) + data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)))./4;
        stack2 = (data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)) + data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)) + data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)) + data11(:,i))./4;
        stack3 = (data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)) + data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)) + data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)) + data11(:,i))./4;
        stack4 = (data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)) + data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)) + data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)) + data11(:,i))./4;
        stack5 = (data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)) + data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)) + data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)) + data11(:,i))./4;
%     stack = reshape((data2(-lag1+1:nlines*tLength*60*sps-lag1)+data3(-lag2+1:nlines*tLength*60*sps-lag2)+data11(:))./3,tLength*60*sps,[]);
    [temp_xc1,~] = xc_norm(data1((i-1)*sps*60*tLength+1:i*sps*60*tLength),stack1,w*sps,0); % uses 15 second window
    [temp_xc2,~] = xc_norm(data2((i-1)*sps*60*tLength - lag1(i)+1:i*sps*60*tLength - lag1(i)),stack2,w*sps,0);
    [temp_xc3,~] = xc_norm(data3((i-1)*sps*60*tLength - lag2(i)+1:i*sps*60*tLength - lag2(i)),stack3,w*sps,0);
    [temp_xc4,~] = xc_norm(data4((i-1)*sps*60*tLength - lag3(i)+1:i*sps*60*tLength - lag3(i)),stack4,w*sps,0);
    [temp_xc5,~] = xc_norm(data5((i-1)*sps*60*tLength - lag4(i)+1:i*sps*60*tLength - lag4(i)),stack5,w*sps,0);
%     [temp_xc2,temp_lag(i-i_begin+1)] = xc_norm(data2((i-1)*sps*60*tLength - lag1+1:i*sps*60*tLength - lag1),stack(:,i),w*sps);
    temp_xc(:,i) = (temp_xc1 + temp_xc2 + temp_xc3 + temp_xc4 + temp_xc5)./5;
    stack(:,i) = (stack1 + data11(:,i))./2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---text displays sum of power levels per line---
    %     text(-5,-2*i,num2str(round(sum(sum(P([nf*i-nf+1:nf*i]+nf-nf*i_begin,:))))))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---For arrows annotating explosions---
    % [pk,loc] = findpeaks(data1(:,i).*temp_xc(:,i),'minpeakheight', 4,'minpeakdistance',60*sps);
    % plot(loc/60/sps,pk/scale_opt+2.5-2*i,'vm','markerfacecolor','m','markersize',4);
    % text(loc/60/sps-.4,pk/scale_opt+3-2*i,'<','color','m','fontsize',26,'fontweight','bold','rotation',90);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
% [~,qa] = findpeaks(temp_xc(:),'minpeakheight',0.98,'minpeakdistance',5);
% [~,qb] = findpeaks(stack(i_begin*sps*tLength*60:i_end*sps*tLength*60),'minpeakheight',5,'minpeakdistance',5);
% q = unique([qa;qb']);
% % % q = find(temp_xc(:) >= .99);
% q_times = q(2:end)-q(1:end-1);
% q_idx = find(q_times./sps/60 > 10);
% q_times = q_times(q_times./sps/60 > 10)./sps/60;
% q_idx = q(q_idx);

[~,~,T,~] = spectrogram(data11(:,1),w*sps/2,(w-1/5*w)*sps/2,hz(:)',sps);

figure(321),clf;% calling figure
imagesc(T/60,[-1 -i_end],(temp_xc'))                   % plotting PSD images
set(gca,'YDir','normal')
hold on;
    plot([0:1:tLength*sps*60-1]'./sps./60,stack./scale_opt + ones(numel(data11(:,1)),1)*[-1:-1:-i_end],'-k')
%     plot([0:1:tLength*sps*60-1]'./sps./60,temp_xc + nf*ones(numel(data11(:,1)),1)*[-i_begin:-1:-i_end] + 1.5,':g')
%     plot([0:1:tLength*sps*60-1]'./sps./60,data11(:,i_begin:i_end)./scale_opt + nf*ones(numel(data11(:,1)),1)*[-i_begin:-1:-i_end] + 1.5,'-g')
title('Xcorr-Helicorder  Sakurajima activity  07/18-07/25  2013')
xlabel('Minutes')
% ylabel('Date')
text(ones(ceil((i_end-i_begin)/2+0.00001),1)*-17,-1*[i_begin:2:i_end]',datestr(t(1,i_begin:2:i_end),'ddd HH PM'))
% text(-3,-i_begin+1.5,num2str(i_begin)) % annotate start line
% text(-3,-i_end+1.5,num2str(i_end))     % annotate end line
ylim([-i_end-1 1-i_begin]);
set(gca,'YTick',[]);

% temp_stack = stack(i_begin*sps*tLength*60:i_end*sps*tLength*60);
% for i = 1:numel(q_idx)-1
%     text(mod(q_idx(i+1),tLength*sps*60)/sps/60-.5,-2*ceil(q_idx(i+1)/tLength/sps/60)+2-2*i_begin,num2str(round(q_times(i))),'fontsize',6);
%     % to plot blue lines denoting closed periods ***IN PROGRESS***
%     %     plot(mod(q_idx(i+1),tLength*sps*60)/sps/60.*[1 1]-[1 mod(round(q_times(i)),tLength*sps*60)],[1 1].*-2*ceil(q_idx(i+1)/tLength/sps/60)+3.8-2*i_begin,'-b','linewidth',1)
%     pk(i) = max(temp_stack(q_idx(i)-12000:q_idx(i)+12000));
% end
% colormap(fliplr(hsv))
% colormap(hsv)
% colormap([1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 .7 .1; 1 0.3 0; 1 0 0])
% colormap([1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 0; 1 .7 .1; 1 0.3 0; 1 0 0])
load('yet_white.mat');
colormap(yet_white);
h = colorbar;
ylabel(h, '\fontsize{18}Mean normalized running xc value')

% figure(),clf;
% plot(q_times(1:end-1),pk,'ob','markerfacecolor','b')
% xlabel('quiescence length [minutes]')
% ylabel('peak pressure amplitude')

end
    

