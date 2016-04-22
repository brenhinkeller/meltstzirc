% mctest.m
% Produce the the (~10^6 row) monte carlo table mcign for quick examination
% with plotmcvariable(s).m
% Uses static per-variable errors
%%

%igncn1
simitems={'Kv';'Latitude';'Longitude';'Elevation';'SiO2';'TiO2';'Al2O3';'Fe2O3';'Fe2O3T';'FeO';'FeOT';'MgO';'CaO';'Na2O';'K2O';'P2O5';'MnO';'H2O_Total';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Li';'Be';'B';'C';'CO2';'F';'Cl';'Sc';'Ti';'V';'Cr';'Co';'Ni';'Cu';'Zn';'Ga';'Zr';'Os';'Rb';'Bi';'Hg';'Ba';'Y';'Pb';'Te';'Nb';'Sr87_Sr86';'Tl';'Pt';'Sn';'Cd';'As';'Pd';'Sr';'Se';'S';'Au';'Ta';'Mo';'U';'Cs';'Sb';'Ag';'W';'Th';'Re';'Hf';'Ir';'tc1Lith';'tc1Crust';'Crust';'Vp';'Vs';'Rho';'Upper_Crust';'Upper_Vp';'Upper_Vs';'Upper_Rho';'Middle_Crust';'Middle_Vp';'Middle_Vs';'Middle_Rho';'Lower_Crust';'Lower_Vp';'Lower_Vs';'Lower_Rho';'Freeair';'Bouger';'Eustar';'Cestar';'Mbulk';'Tliq';'Tsat';'Tsatb';'Zrsat';'Zrf';'Ff';'MZr'}; % igncn1
%ign simitems={'Latitude';'Longitude';'Elevation';'SiO2';'Rb';'K2O';'Na2O';'Li';'Th';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Y';'Lu';'Zr';'Hf';'Nb';'Ta';'Mo';'W';'Cs';'U';'Pb';'Ba';'Ga';'MnO';'Sr';'P2O5';'Al2O3';'Ni';'Cr';'Co';'TiO2';'V';'Sc';'Zn';'Cu';'MgO';'Fe2O3T';'CaO';'Eustar';'tc1Crust';'tc1Lith';'Crust';'Vp';'Vs';'Rho';'Upper_Crust';'Upper_Vp';'Upper_Vs';'Upper_Rho';'Middle_Crust';'Middle_Vp';'Middle_Vs';'Middle_Rho';'Lower_Crust';'Lower_Vp';'Lower_Vs';'Lower_Rho';'Latitude';'Longitude';'Freeair';'Bouger';};
%ignf simitems={'Latitude';'Longitude';'Elevation';'SiO2';'TiO2';'Al2O3';'Fe2O3';'Fe2O3T';'FeO';'FeOT';'MgO';'CaO';'Na2O';'K2O';'P2O5';'MnO';'Loi';'H2O_Plus';'Cr2O3';'La';'NiO';'CaCO3';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Li';'Be';'B';'C';'CO2';'F';'Cl';'K';'Ca';'Mg';'Sc';'Ti';'V';'Fe';'Cr';'Mn';'Co';'Ni';'Cu';'Zn';'Ga';'Zr';'Os';'Rb';'Pb206_Pb208';'Al';'Bi';'I';'Hg';'Nd143_Nd144';'Ba';'Y';'Pb206_Pb207';'Pb';'D18O';'Te';'Hf176_Hf177';'Nb';'Pb207_Pb204';'Lu176_Hf177';'Pb206_Pb204';'Sr87_Sr86';'Os187_Os188';'Tl';'Pt';'Sn';'Cd';'Pb208_Pb204';'As';'Pd';'Sr';'Se';'S';'Au';'Os187_Os186';'Ta';'P';'Mo';'U';'Cs';'Sb';'Ag';'W';'Th';'Re';'Hf';'Ir';'Crust';'Vp';'Vs';'Rho';'Upper_Crust';'Upper_Vp';'Upper_Vs';'Upper_Rho';'Middle_Crust';'Middle_Vp';'Middle_Vs';'Middle_Rho';'Lower_Crust';'Lower_Vp';'Lower_Vs';'Lower_Rho';'Freeair';'Bouger';};
%igne? simitems={'Latitude';'Longitude';'Loc_Prec';'SiO2';'TiO2';'Al2O3';'Fe2O3';'Fe2O3T';'FeO';'FeOT';'MgO';'CaO';'Na2O';'K2O';'P2O5';'MnO';'Loi';'H2O_Plus';'Cr2O3';'La';'NiO';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Li';'Be';'B';'C';'CO2';'F';'Cl';'K';'Ca';'Mg';'Sc';'Ti';'V';'Fe';'Cr';'Mn';'Co';'Ni';'Cu';'Zn';'Ga';'Zr';'Os';'Rb';'Al';'Bi';'I';'Hg';'Nd143_Nd144';'Ba';'Y';'Pb';'D18O';'Te';'Nb';'H';'Sr87_Sr86';'Tl';'Pt';'Sn';'Cd';'As';'Pd';'Sr';'Se';'S';'Au';'Ta';'P';'Mo';'U';'Cs';'Sb';'Ag';'W';'Th';'Re';'Hf';'Ir';};


datain=zeros(length(igncn1.Age),length(simitems)+2);
uncertainty=zeros(1,length(simitems)+2);
for i=1:length(simitems)
    eval(sprintf('datain(:,%i)=igncn1.%s;', i+2, simitems{i}))
    eval(sprintf('uncertainty(%i)=igncn1.err.%s;',i+2,simitems{i}))
end
clear i

%% Produce sample weights

datain(:,1)=igncn1.Age;
agecert=igncn1.Age_Max-igncn1.Age_Min;
agecert(agecert<50)=50;
agecert(isnan(agecert))=50;
datain(:,2)=agecert;

test=~isnan(igncn1.Age)&~isnan(igncn1.Latitude)&~isnan(igncn1.Longitude)&igncn1.Elevation>-100;%&~igncn1.oibs;
data=datain(test,:);

tic;
% if isfield(igncn1,'k')
%     k=igncn1.k;
% else
    k=invweight(igncn1.Latitude(test),igncn1.Longitude(test),igncn1.Age(test));
    igncn1.k=k;
% end

prob=1./((k.*median(5./k))+1);
fprintf('Calculating sample weights: ')
toc


%% Run the monte carlo

% Number of rows to simulate
samplerows=1000000;

tic;
mcigncn1.data=NaN(samplerows,size(data,2));
toc


tic;
i=1;
while i<samplerows
    % select weighted sample of data
    r=rand(length(prob),1);
    sdata=data(prob>r,:);
    
    % Randomize ages over uncertainty interval
    r=randn(size(sdata(:,1)));
    sdata(:,1)=sdata(:,1)+r.*sdata(:,2)/2;
    
    if i+size(sdata,1)-1<=samplerows
        mcigncn1.data(i:i+size(sdata,1)-1,:)=sdata;
    else
        mcigncn1.data(i:end,:)=sdata(1:samplerows-i+1,:);
    end
    i=i+size(sdata,1);
end

% Randomize geochemical variables over uncertainty interval
mcigncn1.data=mcigncn1.data+mcigncn1.data.*repmat(uncertainty,samplerows,1).*randn(samplerows,length(simitems)+2);
mcigncn1.elements=['Age';'Age_Uncert';simitems;];


toc
mcigncn1=elementify(mcigncn1);

% save igncn1 igncn1
% save mcigncn1 mcigncn1

