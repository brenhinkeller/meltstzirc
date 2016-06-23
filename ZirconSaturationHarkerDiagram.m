%% Load variables
cd ~/Desktop/meltstzirc/output/
n = '5F6kb3H2O';
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

xmin=40;
xmax=70;
nbins=30;
binoverlap=1;
rt=[0,1000,2000,3000,4000];
figure;
l={};
for i=1:length(rt)-1
    test=mcigncn1.Age>rt(i)&mcigncn1.Age<rt(i+1)&mcigncn1.Elevation>-100&mcigncn1.Ff<35;

    [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,binoverlap);
    hold on; errorbar(c,m,e,'.','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)])
    
%     [c,m,e]=bin(mcigncn1.(xElem)(test),mcigncn1.(yElem)(test),xmin,xmax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins*10,binoverlap*10);
%     hold on; plot(c,m,'-','Color',[i/(length(rt)-1), 0, 1-i/(length(rt)-1)])

    l = [l,{[num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma']}];
end
xlabel(xElem); ylabel(yElem)
formatfigure
xlim([xmin xmax])
legend(l)

% ylim([0 540])
title(suffix)

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
