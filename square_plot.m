function square_plot_for_katha (uM2)
ColIx = {'PYRPYR','PYRPVB','PYRSOM','PYRVIP',...
    'SOMSOM','VIPVIP','PVBSOM',...
    'SOMVIP','PVBVIP','PVBPVB'};

% uM2  = cell(length(ColIx),2);% 2 Columns are NC pre-post learning for each combination above


% SQUARE PLOTS
lb={'PYR','PVB','SOM','VIP'};
point2inch = 1/72; % 1 point = 1/72 inch
inch2cm    = 2.54; % 1 inch = 2.54 centimeters
point2cm   = point2inch*inch2cm; % 1 point = 1/72*2.54 centimeters
%data2cm    = ps1(1)./diff(figax(1:2));
%crvlw_point = (crvlw_pix*data2cm)./point2cm
tf =0.5/point2cm;



uM3 = (uM2(:,2)-uM2(:,1));% abs change
mi = min(uM3);
ma = max(uM3);
clim = [-max(abs([mi,ma])), 0 , max(abs([mi,ma]))];
clim = [-0.5,0,0.5];
% clim = [0,0.1,0.2];
% clim = [-.1 0 .08];
posneg=1; % if bpth positive and negative values

brd=1;
bw=3;bh=3;
%bw=5;bh=5;

fa=1;
isw=3;ish=2;
bh2=(bh-2*ish)/3;
fw = 1*bw+2*brd+0*isw;
fh = 1*bh+2*brd+0*ish;

createfig=0;
f=figure('Units','centimeters','Position',[1,1,fw+2*fa,fh+2*fa],'Visible','on');
c=colormap(parula(64));
cm=32;
cf=64/(clim(3)-clim(1)); % ma*cf

ax0 = axes();
set(ax0,'Units','centimeters','Position',[fa fa fw fh],'visible','off');
axis([0 fw 0 fh]); % scale so that it has the same coordinate system as figure;
set(gca,'XTick',[1:fw],'YTick',[1:fh]);
if createfig,
    set(ax0,'Visible','on','Color','none');grid on;
end
ax = axes('Units','centimeters','Position',[fa+brd,fa+brd,bw,bh]);

set(ax,'Color','none');
hold on;



xlim([-.5 1.5]);ylim([-.5 1.5])

CC=[0,1; 1 1; 1 0; 0 0]; % coordinates of square: clockwise from top left: pyr pv som vip
list={'PYRPVB','PVBSOM','SOMVIP','PYRVIP','PYRSOM','PVBVIP'}
for li=1:length(list),
    tmp=list{li};
    cx=strmatch(tmp,ColIx); cx=cx(1); cx1=strmatch(tmp(1:3),lb); cx2=strmatch(tmp(4:6),lb);
    val=uM3(cx);
    if posneg,
        if 1% val<0,
            %cval=cm-round(val*cf);
            cval=round((val-clim(1))*cf);
        else,
            cval=cm+round(val*cf);
        end;
        %lw = (abs(cval-cm)/32)*tf+eps;
        lw = (abs(val)/max(abs(clim)))*tf+eps;
    else,
        cval=round(uM3(cx)*cf);
        lw = (cval/64)*tf;
    end
    if cval<1, cval=1; elseif cval>64, cval=64; end
    if lw<=1;lw = 1;end
    line(CC([cx1,cx2],1),CC([cx1,cx2],2),'linewidth',lw, 'color',c(cval,:));
end % for li=1:length(list),

% autocorr
[x,y]=pol2cart(45/180*pi,1);
[x,y]=pol2cart(45/180*pi,0.5);

list={'PYRPYR','PVBPVB','SOMSOM','VIPVIP'}
for li=1:length(list),
    tmp=list{li};
    
    if strcmp(tmp(1:3),'PYR'),
        xn=-x; yn=y;
    elseif strcmp(tmp(1:3),'PVB'),
        xn=x; yn=y;
    elseif strcmp(tmp(1:3),'SOM'),
        xn=x; yn=y;
    elseif strcmp(tmp(1:3),'VIP'),
        xn=-x; yn=y;
    end
    
    cx=strmatch(tmp,ColIx);cx1=strmatch(tmp(1:3),lb);
    val=uM3(cx);
    if posneg,
        if 1%val<0,
            cval=cm-round(uM3(cx)*cf);
            cval=round((val-clim(1))*cf);
        else,
            cval=cm+round(uM3(cx)*cf);
        end;
        lw = (abs(cval-cm)/32)*tf;
        lw = (abs(val)/max(abs(clim)))*tf+eps;
    else,
        cval=round(uM3(cx)*cf);
        lw = (cval/64)*tf;
    end
    if cval<1, cval=1; elseif cval>64, cval=64; end
    if lw<=1;lw = 1;end
    line(CC([cx1,cx1],1)+[0,xn]',CC([cx1,cx1],2)+[0,yn]','linewidth',lw, 'color',c(cval,:));
end

plot(CC(:,1),CC(:,2),'o','markersize',1/point2cm,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'LineWidth',2)

text(CC(1,1),CC(1,2),'PYR','HorizontalAlignment','center')
text(CC(2,1),CC(2,2),'PVB','HorizontalAlignment','center')
text(CC(3,1),CC(3,2),'SOM','HorizontalAlignment','center')
text(CC(4,1),CC(4,2),'VIP','HorizontalAlignment','center')

axis off;
set(gcf,'PaperPositionMode','auto');

dum=get(gca,'Position');
colormap(c);
hc=colorbar('Units','centimeters');
dum2=get(hc,'Position');

dum2(3)=dum2(3)/2;
dum2(4)=dum2(4)/4;
dum2(1)=dum2(1)+1.5; %dum2(1)=dum2(1)+2;
dum2(2)=dum2(2)+0;
set(gca,'Position',dum);
set(hc,'Position',dum2);%set(hc,'YLim',[-0.08 0.08],'YTick',[-0.08 0 0.08]);
axes(ax0);
text(dum2(1)+dum2(3)/2-fa,dum2(2)+dum2(4)-fa,'NC','VerticalAlignment','bottom','HorizontalAlignment','center');
%set(gca,'Color','none');
if posneg,
    set(hc,'YTick',[1,cm,64]./64,'YTickLabel',clim);
else
    set(hc,'YTick',[1,cm,64]./64,'YTickLabel',clim);
end

%xlim([0 2]);ylim([0 2])

axes(ax(1));
title('Change');

