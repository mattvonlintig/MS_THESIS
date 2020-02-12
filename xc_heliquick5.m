function [stack,temp_xc] = xc_heliquick5(data1,data2,data3,data4,data5,tLength,sps,w)

% function SPECTOHELICORDER displays recording, spectrogram, and xc-vals for
% lengthy time series' as multi-line figure according to input parameters:
% 
% *function calls another custom function XC_NORM.m
% *function loads custom colormap YET_WHITE.MAT
% 
% INPUTS: 
% data1,data2..-seperate channels/receivers (time-aligned)
% tLength     - time length per line [Minutes]
% sps         - samples per second [Hz]
% scale_opt   - factor to scale data by for plotting purposes
% w           - xc window size [seconds]
% time_vect   - Matlab time vector to annotate lines with datetime
% 
% OUTPUTS: 
% stack       - entire dataset stack in tLength minute lines
% temp_xc     - max correlation lags calculated for each time interval (line)

% (figure)    - Plots a helicorder style spectrogram influenced by xc_norm vals


tic
nlines = floor(numel(data1)/(tLength*60*sps)); % number of lines in helicorder plot
data11 = reshape(data1(1:nlines*tLength*60*sps),tLength*60*sps,[]); % reshapes data for helicorder-style plot
% t = reshape(time_vect(1:nlines*tLength*60*sps),tLength*60*sps,[]);
data22 = reshape(data2(1:nlines*tLength*60*sps),tLength*60*sps,[]);
data33 = reshape(data3(1:nlines*tLength*60*sps),tLength*60*sps,[]); % reshapes data for helicorder-style plot
data44 = reshape(data4(1:nlines*tLength*60*sps),tLength*60*sps,[]);
data55 = reshape(data5(1:nlines*tLength*60*sps),tLength*60*sps,[]);
clear data1 data2 data3 data4 data5;

i_begin = 1; % DEFAULT = 1
i_end = nlines-2; %   DEFAULT = NLINES
    parfor i = i_begin:i_end

        stack1 = (data22(:,i) + data33(:,i) + data44(:,i) + data55(:,i))./4;
        stack2 = (data33(:,i) + data11(:,i) + data44(:,i) + data55(:,i))./4;
        stack3 = (data22(:,i) + data11(:,i) + data44(:,i) + data55(:,i))./4;
        stack4 = (data22(:,i) + data11(:,i) + data33(:,i) + data55(:,i))./4;
        stack5 = (data22(:,i) + data11(:,i) + data44(:,i) + data33(:,i))./4;

    [temp_xc1,~] = xc_norm(data11(:,i),stack1,w*sps,0); 
    [temp_xc2,~] = xc_norm(data22(:,i),stack2,w*sps,0);
    [temp_xc3,~] = xc_norm(data33(:,i),stack3,w*sps,0);
    [temp_xc4,~] = xc_norm(data44(:,i),stack4,w*sps,0);
    [temp_xc5,~] = xc_norm(data55(:,i),stack5,w*sps,0);

    temp_xc(:,i) = (temp_xc1 + temp_xc2 + temp_xc3 + temp_xc4 + temp_xc5)./5;
    
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
stack = (data11(:,i_begin:i_end) + data22(:,i_begin:i_end) + data33(:,i_begin:i_end) + data44(:,i_begin:i_end) + data55(:,i_begin:i_end))./5;
% [~,~,T,~] = spectrogram(data11(:,1),w*sps/2,(w-1/5*w)*sps/2,1,sps);



% figure(3210),clf;% calling figure
% imagesc(T/60,[-1 -i_end],(temp_xc'))% plotting xcorr image
% set(gca,'YDir','normal')
% hold on;
% plot([0:1:tLength*sps*60-1]'./sps./60,abs(stack.^(0.5)).*sign(stack)/4 + ones(numel(data11(:,1)),1)*[-1:-1:-i_end],'-k')
% % title('Xcorr-Helicorder  Sakurajima activity  07/18-07/25  2013')
% title({['Recorded Infrasound'], ['Sakurajima 07/18-07/25  2013']},'fontsize',24)
% % xlabel('\fontsize{16}Minutes')
%  set(gca,'FontSize',16,'FontWeight','bold')
% % ylabel('Date')
% text(ones(ceil((i_end-i_begin)/2+0.00001),1)*-20,-1*[i_begin:2:i_end]',datestr(t(1,i_begin:2:i_end),'HH PM'),'fontsize',16,'fontweight','bold')
% ylim([-i_end-1 1-i_begin]);
% set(gca,'YTick',[]);
% h = colorbar;
% colormap([1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 .7 .1; 1 0.3 0; 1 0 0])
% colormap([1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 0; 1 .7 .1; 1 0.3 0; 1 0 0])
% load('yet_white.mat');
% colormap(yet_white);

% ylabel(h, '\fontsize{18}Mean normalized running xc value')

end