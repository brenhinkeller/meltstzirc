%% Plot a variable

if ~exist('mcign','var')
    load mcign
end
if ~exist('ign','var')
    load ign
end

%% All together
test=mcign.SiO2>40&mcign.SiO2<80;
figure;

[c,m,e] = bin(mcign.Age(test),mcign.Zr(test),0,4000,length(mcign.SiO2)./length(ign.SiO2),40);
yyaxis left; errorbar(c,m,e,'.b');
ylabel('Zr')

[~, M]=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
[c,m,e] = bin(mcign.Age(test),M,0,4000,length(mcign.SiO2)./length(ign.SiO2),40);
yyaxis right; errorbar(c,m,e,'.r');
ylabel('M')

xlabel('Age (Ma)'); 
formatfigure;

%% T
n=~any(isnan([ign.CaO ign.Na2O ign.K2O ign.Al2O3 ign.SiO2 ign.FeOT, ign.MgO, ign.Zr]),2);

test=mcign.SiO2>43&mcign.SiO2<51;
T=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),T,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('43-51% SiO2'); xlabel('Age (Ma)'); ylabel('Zircon saturation temperature (C)')

test=mcign.SiO2>51&mcign.SiO2<62;
T=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),T,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('51-62% SiO2'); xlabel('Age (Ma)'); ylabel('Zircon saturation temperature (C)')

test=mcign.SiO2>62&mcign.SiO2<74;
T=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),T,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('62-74% SiO2'); xlabel('Age (Ma)'); ylabel('Zircon saturation temperature (C)')

test=mcign.SiO2>74&mcign.SiO2<80;
T=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),T,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('74-80% SiO2'); xlabel('Age (Ma)'); ylabel('Zircon saturation temperature (C)')

    
%% M

test=mcign.SiO2>43&mcign.SiO2<51;
[~, M]=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),M,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('43-51% SiO2'); xlabel('Age (Ma)'); ylabel('M')

test=mcign.SiO2>51&mcign.SiO2<62;
[~, M]=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),M,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('51-62% SiO2'); xlabel('Age (Ma)'); ylabel('M')

test=mcign.SiO2>62&mcign.SiO2<74;
[~, M]=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),M,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('62-74% SiO2'); xlabel('Age (Ma)'); ylabel('M')

test=mcign.SiO2>74&mcign.SiO2<80;
[~, M]=tzirc(mcign.CaO(test),mcign.Na2O(test),mcign.K2O(test),mcign.Al2O3(test),mcign.SiO2(test),mcign.TiO2(test),mcign.FeOT(test),mcign.MgO(test),mcign.MnO(test),mcign.P2O5(test), mcign.Zr(test));
figure; binplot(mcign.Age(test),M,0,4000,length(mcign.SiO2)./length(ign.SiO2),40,'.r')
title('74-80% SiO2'); xlabel('Age (Ma)'); ylabel('M')
