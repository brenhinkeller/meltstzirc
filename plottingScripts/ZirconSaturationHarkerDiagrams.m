%% Load variables
cd ~/Desktop/meltstzirc/output/
n = '1F6kb3H2O';
if ~exist('mcigncn1','var') || ~exist('igncn1','var') || ~exist('suffix','var') || ~isequal(suffix,n)
    suffix = n;
    load(['mcigncn1-' suffix]);
    load(['igncn1-' suffix]);
end

%% Check out mass saturated as a function of SiO2

yElem = 'MZrn';
xElem = 'SiO2';

if ~isfield(igncn1,'MZrn') || ~isfield(mcigncn1,'MZrn')
    igncn1.MZrn = igncn1.MZr * 2.009;
    mcigncn1.MZrn = mcigncn1.MZr * 2.009;
end


xmin=43;
xmax=72;
nbins=30;
binoverlap=2;
rt=[0,1000,2000,3000,4000];

figure; hold on;
a = zeros(1,length(rt)-1);
l = cell(1,length(rt)-1);
for i=1:length(rt)-1
    test=mcigncn1.Age>rt(i)&mcigncn1.Age<rt(i+1)&mcigncn1.Elevation>-100&mcigncn1.Ff<35;

    [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);

    % Just errorbars
    a(i)=errorbar(c,m,2*e,'.','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)]);
    
%     a(i)=plot(c,m,'Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)]); % center line
%     plot(c,m+2*e,'Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)]); % upper line
%     plot(c,m-2*e,'Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)]); % lower line
    
    % Representative errorbars 
%     t = mod((1:length(c))-60,133)==0;
%     errorbar(c(t),m(t),2*e(t),'.k')
    
    l{i} = [num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma'];
end
xlabel(xElem); ylabel(yElem)
legend(a,l)
formatfigure
xlim([xmin xmax])

% ylim([0 540])
title(suffix)

%% Equivalent SiO2 for Archean samples

yElem = 'MZrn';
xElem = 'SiO2';

if ~isfield(igncn1,'MZrn') || ~isfield(mcigncn1,'MZrn')
    igncn1.MZrn = igncn1.MZr * 2.009;
    mcigncn1.MZrn = mcigncn1.MZr * 2.009;
end

xmin=47;
xmax=57;
nbins=20;
binoverlap=1;
rt=[0,541,2500,4000];

test=mcigncn1.Age>rt(1)&mcigncn1.Age<rt(2)&mcigncn1.Elevation>-100 & mcigncn1.MZr>0&mcigncn1.Ff<35;
[c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);

meq = m;
eeq = m;
for i=1:length(c)
    meq(i) = nanmean(mcigncn1.SiO2(mcigncn1.MZrn>(m(i)-2*e(i)) & mcigncn1.MZrn<(m(i)+2*e(i)) & mcigncn1.Age>rt(3) & mcigncn1.Age<rt(4)));
    eeq(i) = nansem(mcigncn1.SiO2(mcigncn1.MZrn>(m(i)-2*e(i)) & mcigncn1.MZrn<(m(i)+2*e(i)) & mcigncn1.Age>rt(3) & mcigncn1.Age<rt(4))).*sqrt(length(mcigncn1.SiO2)./length(igncn1.SiO2));
end
figure; errorbar(c,meq,2*eeq,'.')

set(gca,'ytick',[60 62 64 66 68 70])
formatfigure
xlabel('SiO2 of Phanerozoic sample')
ylabel('Average SiO2 of Archean sample saturating equivalent zircon mass')



%% Equivalent SiO2 for Phanerozoic samples

yElem = 'MZrn';
xElem = 'SiO2';

if ~isfield(igncn1,'MZrn') || ~isfield(mcigncn1,'MZrn')
    igncn1.MZrn = igncn1.MZr * 2.009;
    mcigncn1.MZrn = mcigncn1.MZr * 2.009;
end

xmin=55;
xmax=70;
nbins=20;
binoverlap=1;
rt=[2500,4000,0,541];

test=mcigncn1.Age>rt(1)&mcigncn1.Age<rt(2)&mcigncn1.Elevation>-100 & mcigncn1.MZr>0&mcigncn1.Ff<35;
[c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);

meq = m;
eeq = m;
for i=1:length(c)
    meq(i) = nanmean(mcigncn1.SiO2(mcigncn1.MZrn>(m(i)-2*e(i)) & mcigncn1.MZrn<(m(i)+2*e(i)) & mcigncn1.Age>rt(3) & mcigncn1.Age<rt(4)));
    eeq(i) = nansem(mcigncn1.SiO2(mcigncn1.MZrn>(m(i)-2*e(i)) & mcigncn1.MZrn<(m(i)+2*e(i)) & mcigncn1.Age>rt(3) & mcigncn1.Age<rt(4))).*sqrt(length(mcigncn1.SiO2)./length(igncn1.SiO2));
end
figure; errorbar(c,meq,2*eeq,'.')

% set(gca,'ytick',[60 62 64 66 68 70])
formatfigure
xlabel('SiO2 of Archean sample')
ylabel('Average SiO2 of Phanerozic sample saturating equivalent zircon mass')


%% Saturation temperature as a function of SiO2, for samples that actually saturate zircon 

yElem = 'Tsat';
xElem = 'SiO2';

xmin=40;
xmax=80;
nbins=40;
binoverlap=1;
rt=[0,1000,2000,3000,4000];
figure;
l={};
for i=1:length(rt)-1
    test=mcigncn1.Age>rt(i)&mcigncn1.Age<rt(i+1)&mcigncn1.Elevation>-100 & mcigncn1.MZr>0;

    [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);
    hold on; errorbar(c,m,2*e,'.','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)])
    xlabel(xElem); ylabel(yElem)
    formatfigure
    xlim([xmin xmax])
    l = [l,{[num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma']}];
end
legend(l)
title(suffix)

%% Crystallization temperature as a function of SiO2, for samples that actually saturate zircon 

yElem = 'Tcryst';
xElem = 'SiO2';

xmin=40;
xmax=80;
nbins=40;
binoverlap=1;
rt=[0,1000,2000,3000,4000];
figure;
l={};
for i=1:length(rt)-1
    test=mcigncn1.Age>rt(i)&mcigncn1.Age<rt(i+1)&mcigncn1.Elevation>-100 & mcigncn1.MZr>0;

    [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);
    hold on; errorbar(c,m,2*e,'.','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)])
    xlabel(xElem); ylabel(yElem)
    formatfigure
    xlim([xmin xmax])
    l = [l,{[num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma']}];
end
legend(l)
title(suffix)

%% Bulk saturation temperature as a function of SiO2

yElem = 'Tsatb';
xElem = 'SiO2';

xmin=40;
xmax=80;
nbins=40;
binoverlap=1;
rt=[0,1000,2000,3000,4000];
figure;
l={};
for i=1:length(rt)-1
    test=mcigncn1.Age>rt(i)&mcigncn1.Age<rt(i+1)&mcigncn1.Elevation>-100;

    [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);
    hold on; errorbar(c,m,2*e,'.','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)])
    xlabel(xElem); ylabel(yElem)
    formatfigure
    xlim([xmin xmax])
    l = [l,{[num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma']}];
end
legend(l)
title(suffix)

%% Final temperature as a function of SiO2

yElem = 'Tf';
xElem = 'SiO2';

xmin=40;
xmax=80;
nbins=40;
binoverlap=1;
rt=[0,1000,2000,3000,4000];
figure;
l={};
for i=1:length(rt)-1
    test=mcigncn1.Age>rt(i)&mcigncn1.Age<rt(i+1)&mcigncn1.Elevation>-100;

    [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);
    hold on; errorbar(c,m,2*e,'.','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)])
    xlabel(xElem); ylabel(yElem)
    formatfigure
    xlim([xmin xmax])
    l = [l,{[num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma']}];
end
legend(l)
title(suffix)
