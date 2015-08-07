if ~exist('mcigncn1','var'); load mcigncn1; end
if ~exist('igncn1','var'); load igncn1; end

%% Mass of zircon saturated over time

Elem='MZr';
agemin=0;
agemax=4000;
rsi=[43,51,62,74,80];
figure;
l=[];
for i=1:length(rsi)-1
    test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,2*e,'.','Color',[1-i/(length(rsi)-1),0,i/(length(rsi)-1)])
    l=[l {[num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']}];
end
xlabel('Age (Ma)'); ylabel(Elem); legend(l);
title('All samples')
formatfigure

%% Mass of zircon saturated over time

Elem='MZr';
agemin=0;
agemax=4000;
rsi=[40 80];
figure;
l=[];
for i=1:length(rsi)-1
    test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,2*e,'.','Color',[1-i/(length(rsi)-1),0,i/(length(rsi)-1)])
    l=[l {[num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']}];
end
xlabel('Age (Ma)'); ylabel(Elem); legend(l);
title('All samples')
formatfigure


%% Zircon saturation temperature over time

Elem='Tsat';
agemin=0;
agemax=4000;
rsi=[43,51,62,74,80];
figure;
l=[];
for i=1:length(rsi)-1
    test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,2*e,'.','Color',[1-i/(length(rsi)-1),0,i/(length(rsi)-1)])
    l=[l {[num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']}];
end
xlabel('Age (Ma)'); ylabel(Elem); legend(l);
title('All samples')
formatfigure

%% Zircon saturation temperature over time for samples that actually saturate zircon

Elem='Tsat';
agemin=0;
agemax=4000;
rsi=[43,51,62,74,80];
figure;
l=[];
for i=1:length(rsi)-1
    test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100&mcigncn1.MZr>0;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,2*e,'.','Color',[1-i/(length(rsi)-1),0,i/(length(rsi)-1)])
    l=[l {[num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']}];
end
xlabel('Age (Ma)'); ylabel(Elem); legend(l);
title('Samples that actually saturate zircon')
formatfigure

%% Bulk Zircon saturation temperature over time

Elem='Tsatb';
agemin=0;
agemax=4000;
rsi=[43,51,62,74,80];
figure;
l=[];
for i=1:length(rsi)-1
    test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,2*e,'.','Color',[1-i/(length(rsi)-1),0,i/(length(rsi)-1)])
    l=[l {[num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']}];
end
xlabel('Age (Ma)'); ylabel(Elem); legend(l);
title('All samples')
formatfigure

%% Bulk zircon saturation temperature over time for samples that actually saturate zircon

Elem='Tsatb';
agemin=0;
agemax=4000;
rsi=[43,51,62,74,80];
figure;
l=[];
for i=1:length(rsi)-1
    test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100&mcigncn1.MZr>0;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,2*e,'.','Color',[1-i/(length(rsi)-1),0,i/(length(rsi)-1)])
    l=[l {[num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']}];
end
xlabel('Age (Ma)'); ylabel(Elem); legend(l);
title('Samples that actually saturate zircon')
formatfigure

