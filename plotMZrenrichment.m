%% Load variables
cd ~/Desktop/meltstzirc/output/
n = '5F6kb3H2O';
if ~exist('mcigncn1','var') || ~exist('igncn1','var') || ~exist('suffix','var') || ~isequal(suffix,n)
    suffix = n;
    load(['mcigncn1-' suffix]);
    load(['igncn1-' suffix]);
end

%% Probability of saturating zircon in basalt relative to granite

Elem='MZrn';
nbins = 20;

rsi1=[40 50];
test=mcigncn1.SiO2>rsi1(1)&mcigncn1.SiO2<rsi1(2)&mcigncn1.Elevation>-100&mcigncn1.Ff<35;
[c,m1,e1]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),0,4000,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins);
rsi2=[60 70];
test=mcigncn1.SiO2>rsi2(1)&mcigncn1.SiO2<rsi2(2)&mcigncn1.Elevation>-100&mcigncn1.Ff<35;
[~,m2,e2]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),0,4000,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins);

figure; errorbar(c, m1./m2, 1.*((e2./m2).^2+(e1./m1).^2).^0.5.*m1./m2,'.r','LineWidth',1,'MarkerSize',15)

xlabel('igncn1 Age (Ma)'); ylabel(['Mafic / Felsic ' Elem])
title(suffix)



%% Probability of saturating in a sample with < 50% SiO2 

Elem='MZr';
nbins = 20;

rsi1=[40 50];
test=mcigncn1.SiO2>rsi1(1)&mcigncn1.SiO2<rsi1(2)&mcigncn1.Elevation>-100&mcigncn1.Ff<35;
[c,m1,e1]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),0,4000,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins);
m1 = m1 * diff(rsi1);
e1 = e1 * diff(rsi1);
rsi2=[50 70];
test=mcigncn1.SiO2>rsi2(1)&mcigncn1.SiO2<rsi2(2)&mcigncn1.Elevation>-100&mcigncn1.Ff<35;
[~,m2,e2]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),0,4000,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins);
m2 = m2 * diff(rsi2);
e2 = e2 * diff(rsi2);

figure; errorbar(c, m1./m2, 2.*((e2./m2).^2+(e1./m1).^2).^0.5.*m1./m2,'.r','LineWidth',1,'MarkerSize',15)

xlabel('igncn1 Age (Ma)'); ylabel(['Mafic / Nonmafic ' Elem])
title(suffix)


