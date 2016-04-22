% plot enrichment between two different silica ranges

if ~exist('mcigncn1','var')
    load mcigncn1
end
if ~exist('igncn1','var')
    load igncn1
end

Elem='MZr';
nbins = 10;

rsi=[45 50];
test=mcigncn1.SiO2>rsi(1)&mcigncn1.SiO2<rsi(2)&mcigncn1.Elevation>-100;
[c,m1,e1]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),0,4000,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins);
rsi2=[65 70];
test=mcigncn1.SiO2>rsi2(1)&mcigncn1.SiO2<rsi2(2)&mcigncn1.Elevation>-100;
[~,m2,e2]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),0,4000,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins);



figure; errorbar(c, m1./m2, 2.*((e2./m2).^2+(e1./m1).^2).^0.5.*m1./m2,'.r')

title('MZr in basalt relative to granite')

xlabel('igncn1 Age (Ma)'); ylabel(Elem)






