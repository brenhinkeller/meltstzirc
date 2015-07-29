
in=ign;

elems={'SiO2';'TiO2';'Al2O3';'Fe2O3';'Cr2O3';'FeO';'MnO';'MgO';'NiO';'CoO';'CaO';'Na2O';'K2O';'P2O5';'CO2';'H2O_Total';'Zr';'Kv'};

% int SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, MnO, MgO, NiO, CoO, CaO, Na2O, K2O, P2O5, CO2, H2O;


% Optional minor elements
in.Cr2O3=in.Cr*(51.99+24)/51.99/10000;
in.Cr2O3(isnan(in.Cr2O3))=0;
in.NiO=in.Ni*(58.69+24)/58.69/10000;
in.NiO(isnan(in.NiO))=0;
in.CoO=in.Co*(58.93+24)/58.93/10000;
in.CoO(isnan(in.CoO))=0;
in.MnO(isnan(in.MnO))=0;


% Start with all Fe as FeO (will use MELTS to equilibrate fO2)
in=feconversion(in);
in.FeO=in.FeOT;
in.Fe2O3=zeros(size(in.Fe2O3));

% % Set all H2O 4%
% in.H2O_Total=ones(size(in.H2O_Total))*4;

% Set undefined H2O to the average
in.H2O_Total(isnan(in.H2O_Total))=nanmean(in.H2O_Total);

% Set undefined H2O to the average
in.CO2(isnan(in.CO2))=nanmean(in.CO2);

data=zeros(length(in.SiO2),length(elems));
for i=1:length(elems)
    data(:,i)=in.(elems{i});
end

% Reject samples with missing data
test=~any(isnan(data),2);

% Reject samples with suspicious anhydrous normalizations
% test=test &~ ((sum(data,2)-in.H2O_Total)>102) &~ ((sum(data,2)-in.H2O_Total)<92);
test=test &~ (sum(data(:,1:14),2)>100) &~ (sum(data(:,1:14),2)<92);

data=data(test,:);

size(data)

exportcsv('startingcompositions.csv',data,',')


% Try running at 0.3 kbar, 1 kbar etc.