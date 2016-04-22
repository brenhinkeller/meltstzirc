%% Check out mass saturated as a function of SiO2

cd ~/Desktop/meltstzirc/output/
if ~exist('mcigncn1','var'); load mcigncn1; end
if ~exist('igncn1','var'); load igncn1; end

yElem = 'MZr';
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
