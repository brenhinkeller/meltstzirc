%% Import tzirc data
cd ~/Desktop/meltstzirc/output/
 
if ~exist('igncn1','var'); load igncn1; end

suffix = '5F12kb05H2O';
name=['tzirc' suffix];
n = suffix;

% system(['grep -e ''^[0-9\.][0-9\.]*\(\t[0-9\.][0-9\.]*\)\{11\}$'' ' name '.log > ' name '.tsv']);
system(['grep -e ''^[0-9\.na][0-9\.na]*\(\t[0-9\.na][0-9\.na]*\)\{12\}$'' ' name '.log > ' name '.tsv']); % For new version with Tcryst

load([name '.tsv']);

% Make struct from input file 
% variables={'Kv','Mbulk','Tliq','Tsatb','Tf','Tsat','Zrsat','Zrf','Ff','SiO2','Zr','MZr'};
variables={'Kv','Mbulk','Tliq','Tsatb','Tf','Tsat','Zrsat','Zrf','Ff','SiO2','Zr','MZr','Tcryst'}; % New version
tzirclog=struct;
for i=1:length(variables)
   eval(['tzirclog.(variables{i})=' name '(:,i);'])
end

% Create new struct fields for imported data
% variables={'Mbulk','Tliq','Tsat','Tsatb','Tf','Zrsat','Zrf','Ff','MZr'};
variables={'Mbulk','Tliq','Tsat','Tsatb','Tf','Zrsat','Zrf','Ff','MZr','Tcryst'}; % New version
for var=variables;
    igncn1.(var{:})=NaN(size(igncn1.Kv));
end
for var=variables;
    igncn1.err.(var{:})=0.02;
end

% Parse melts struct row by row, inserting each value in the right place
for i=1:length(igncn1.Kv)
    j=find(tzirclog.Kv==igncn1.Kv(i));
    if length(j)==1
        for var=variables;
            igncn1.(var{:})(i)=tzirclog.(var{:})(j);
        end
    elseif length(j)>1
        printf('Warning: Duplicate sample number')
    end
end


%%  Produce the the (~10^6 row) monte carlo table mcign for quick examination 
% From  mctest.m
% Uses static per-variable errors

%igncn1
simitems={'Kv';'Latitude';'Longitude';'Elevation';'SiO2';'TiO2';'Al2O3';'Fe2O3';'Fe2O3T';'FeO';'FeOT';'MgO';'CaO';'Na2O';'K2O';'P2O5';'MnO';'H2O_Total';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Li';'Be';'B';'C';'CO2';'F';'Cl';'Sc';'Ti';'V';'Cr';'Co';'Ni';'Cu';'Zn';'Ga';'Zr';'Os';'Rb';'Bi';'Hg';'Ba';'Y';'Pb';'Te';'Nb';'Sr87_Sr86';'Tl';'Pt';'Sn';'Cd';'As';'Pd';'Sr';'Se';'S';'Au';'Ta';'Mo';'U';'Cs';'Sb';'Ag';'W';'Th';'Re';'Hf';'Ir';'tc1Lith';'tc1Crust';'Crust';'Vp';'Vs';'Rho';'Upper_Crust';'Upper_Vp';'Upper_Vs';'Upper_Rho';'Middle_Crust';'Middle_Vp';'Middle_Vs';'Middle_Rho';'Lower_Crust';'Lower_Vp';'Lower_Vs';'Lower_Rho';'Freeair';'Bouger';'Eustar';'Cestar';'Mbulk';'Tliq';'Tsat';'Tsatb';'Tf';'Zrsat';'Zrf';'Ff';'MZr';'Tcryst'}; % igncn1


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
k=invweight(igncn1.Latitude(test),igncn1.Longitude(test),igncn1.Age(test));
igncn1.k=k;

prob=1./((k.*median(5./k))+1);
fprintf('Calculating sample weights: ')
toc


%% Run Weighted bootstrap resampling

% Number of rows to simulate
samplerows=1000000;

tic;
fprintf('Assigning output variables: ');
mcigncn1.data=NaN(samplerows,size(data,2));
toc


tic;
fprintf('Resampling: ');
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

save(['igncn1-' suffix], 'igncn1')
save(['mcigncn1-' suffix], 'mcigncn1')


