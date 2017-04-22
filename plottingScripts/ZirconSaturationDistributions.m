% Load data
cd ~/Desktop/meltstzirc/output/

if ~exist('igncn1','var'); load igncn1; end

% name='tzircFull5F6kb';
% name='tzirc5F6kb3H2OFull';
name='tzirc1F6kb3H2OFull';
% name='tzirc1F4kb2H2OFull';
% name = 'tzirc6kb3H2OFull';

if ~exist([name '.tsv'],'file')
    % Parse log file to remove any malformed lines (without 9 tab-delimited fields)
    system(['grep -e ''^[0-9\.][0-9\.]*\(\t[0-9\.][0-9\.]*\)\{8\}$'' ' name '.log > ' name '.tsv']);
end

load([name '.tsv']);
eval(['data=' name ';']);

% Fill variables
MZr=data(:,8);
index=data(:,1);
include=zeros(size(MZr));

% Zero-out negative masses
data(MZr<0,8)=0; 
MZr(MZr<0)=0;

% Include all non-zero zircon masses and all adjacent zeros
MZr=[0;MZr;0];
index=[0;index;0];
include = MZr(2:end-1)>0 | (MZr(1:end-2)>0 & index(1:end-2)==index(2:end-1)) | (MZr(3:end)>0 & index(3:end)==index(2:end-1));

% Create the new dataset
fulltzirc.data=data(include,:);
fulltzirc.elements={'Kv';'T';'F';'M';'SiO2';'Zr';'Zrsat';'MZr';'TSat'};
fulltzirc=elementify(fulltzirc);


% Calculate Tscaled
maxkv=max(igncn1.Kv);

fulltzirc.Tscaled=NaN(size(fulltzirc.T));
fulltzirc.MZrscaled=NaN(size(fulltzirc.MZr));
fprintf('\n');

%% Tscaled
for kv=1:maxkv
    test=fulltzirc.Kv==kv;
    if sum(test)>3
        % Normalize temperature to run from 0 at first point of zircon
        % crystallization to 1 at last
        fulltzirc.Tscaled(test)=fulltzirc.T(test)-min(fulltzirc.T(test));
        fulltzirc.Tscaled(test)=1-(fulltzirc.Tscaled(test)./max(fulltzirc.Tscaled(test)));
        % Normalize zircon mass
        fulltzirc.MZrscaled(test)=fulltzirc.MZr(test)./trapz(fulltzirc.Tscaled(test),fulltzirc.MZr(test));
        
    end
end


% figure; plot(fulltzirc.Tscaled, fulltzirc.MZrscaled,'.')
% xlabel('"time"'); ylabel('Zircon amount');
% formatfigure;


%% Check out MZr distributions by silica range

maxpaths=400;
rsi=[43,51,62,74,80];
minpathlength=7;

figure;
for i=1:length(rsi)-1
    subplot(2,2,i); hold on;
    
    sitest=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1);
    path=0;
    for kv=igncn1.Kv(sitest)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(sitest))
            % Plot results for a random subset
            plot(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test));
            path=path+1;
            if path>maxpaths
                break;
            end
        end
        

    end
    ylim([0,5])
    title([num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']); 
    formatfigure
    % xlabel('"time"'); ylabel('Zircon amount');
    % formatfigure;
end

%% Check out MZr distributions by age range

maxpaths=400;
rt=[0,100,1000,2500,4000];
minpathlength=7;

figure;
for i=1:length(rt)-1
    subplot(2,2,i); hold on;
    
    agetest=igncn1.Age>rt(i)&igncn1.Age<rt(i+1);
    path=0;
    for kv=igncn1.Kv(agetest)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(agetest))
            % Plot results for a random subset
            plot(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test));
            path=path+1;
            if path>maxpaths
                break;
            end
        end
        

    end
    ylim([0,5])
    title([num2str(rt(i)) '-' num2str(rt(i+1)) ' Ma']); 
    formatfigure
    % xlabel('"time"'); ylabel('Zircon amount');
    % formatfigure;
end


%% Calculate average normalized MZr distribution as a function of scaled temperature for different silica ranges

maxpaths=70000;
minpathlength=14;
% rsi=[40,50,60,70,80];
rsi=[45,55,65,75];

l=cell(length(rsi)-1,1);
figure;
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(0,1,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test)))) && min(fulltzirc.F(test))<35
            % Plot results for a random subset
            dist = dist + interp1(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test),xi);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end

xlabel('"time"'); ylabel('Zircon amount');
set(gca,'XDir','normal');
legend(l)
formatfigure;


%% Calculate average volcanic MZr distribution as a function of scaled temperature for different silica ranges

% Calculate Tscaled
maxkv=max(igncn1.Kv);

fulltzirc.Tscaled=NaN(size(fulltzirc.T));
fulltzirc.MZrscaled=NaN(size(fulltzirc.MZr));
fprintf('\n');


for kv=1:maxkv
    test=fulltzirc.Kv==kv&fulltzirc.F>35;
    if sum(test)>1
        % Normalize temperature to run from 0 at first point of zircon
        % crystallization to 1 at last
        fulltzirc.Tscaled(test)=fulltzirc.T(test)-min(fulltzirc.T(test));
        fulltzirc.Tscaled(test)=1-(fulltzirc.Tscaled(test)./max(fulltzirc.Tscaled(test)));
        % Normalize zircon mass
        fulltzirc.MZrscaled(test)=fulltzirc.MZr(test)./trapz(fulltzirc.Tscaled(test),fulltzirc.MZr(test));
        
    end
end


minpathlength=2;
% rsi=[40,50,60,70,80];
% rsi=[45,55,65,75];
rsi=[40,80];

l=cell(length(rsi)-1,1);
figure;
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(0,1,100);
    for kv=igncn1.Kv(test)'
        t=fulltzirc.Kv==kv&fulltzirc.F>35;
        if sum(t&fulltzirc.MZr>0)>=minpathlength && all(size(fulltzirc.Tscaled(t)) == size(unique(fulltzirc.Tscaled(t)))) && min(fulltzirc.F(t))<40
            % Plot results for a random subset
            dist = dist + interp1(fulltzirc.Tscaled(t),fulltzirc.MZrscaled(t),xi);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end

xlabel('"time"'); ylabel('Zircon amount');
set(gca,'XDir','normal');
legend(l)
formatfigure;

%% Calculate average normalized MZr distribution as a function of absolute temperature for different silica ranges

maxpaths=70000;
minpathlength=14;
% rsi=[40,50,60,70,80];
rsi=[45,55,65,75];

l=cell(length(rsi)-1,1);
figure;
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(650,1000,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.T(test)) == size(unique(fulltzirc.T(test)))) && min(fulltzirc.F(test))<35
            % Plot results for a random subset
            dist = nansum([dist; interp1(fulltzirc.T(test),fulltzirc.MZr(test),xi)],1);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end

xlabel('Temperature (C)'); ylabel('Zircon amount');
set(gca,'XDir','reverse');
legend(l)
formatfigure;

%% Calculate average normalized MZr path as a function of scaled temperature for different age ranges

maxpaths=70000;
minpathlength=7;
rt=[0 541 2500 4000];
% rt=[0 4000];


l=cell(length(rt)-1,1);
figure;
for i=1:length(rt)-1
    test=igncn1.Age>rt(i)&igncn1.Age<rt(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(0,1,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test))))
            % Plot results for a random subset
            dist = dist + interp1(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test),xi);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rt(i)) ' - ' num2str(rt(i+1)) ' Ma'];
end

xlabel('"time"'); ylabel('Zircon amount');
set(gca,'XDir','normal');
legend(l)
formatfigure;



%% Calculate average MZr paths as a function of time/temperature together


maxpaths=70000;
minpathlength=3;

test=igncn1.SiO2>40&igncn1.SiO2<80;
path=0;
dist = zeros(1,200);
xi=linspace(0,1,200);
for kv=igncn1.Kv(test)'
    test=fulltzirc.Kv==kv;
    if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test)))) && min(fulltzirc.F(test))<35
        % Plot results for a random subset
        dist = dist + interp1(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test),xi);
        path=path+1;
        if path>maxpaths
            break;
        end
    end
end
dist = dist./trapz(xi,dist);
figure; plot(xi,dist);


% load watsongrowthrate
% x=watsongrowthrate(:,1);
% y=watsongrowthrate(:,2);
% y = y./trapz(x,y);
% hold on; plot(x,y);


rsi=[45,55,65,75];

l=cell(length(rsi)-1,1);
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,200);
    xi=linspace(0,1,200);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test)))) && min(fulltzirc.F(test))<35
            % Plot results for a random subset
            dist = dist + interp1(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test),xi);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end


rt=[0 541 2500 4000];

m=cell(length(rt)-1,1);
for i=1:length(rt)-1
    test=igncn1.Age>rt(i)&igncn1.Age<rt(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,200);
    xi=linspace(0,1,200);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test)))) && min(fulltzirc.F(test))<35
            % Plot results for a random subset
            dist = dist + interp1(fulltzirc.Tscaled(test),fulltzirc.MZrscaled(test),xi);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    m{i} = [num2str(rt(i)) ' - ' num2str(rt(i+1)) ' Ma'];
end

xlabel('"time"'); ylabel('Zircon amount');
set(gca,'XDir','normal');
legend([{'Average'};l;m])
formatfigure;


%% Check out F paths by silica range

maxpaths=400;
rsi=[43,51,62,74,80];
minpathlength=5;

figure;
for i=1:length(rsi)-1
    subplot(2,2,i); hold on;
    
    sitest=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1);
    path=0;
    for kv=igncn1.Kv(sitest)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(sitest))
            % Plot results for a random subset
            plot(fulltzirc.Tscaled(test),fulltzirc.F(test));
            path=path+1;
            if path>maxpaths
                break;
            end
        end
        

    end
    title([num2str(rsi(i)) '-' num2str(rsi(i+1)) '% SiO2']); 
    formatfigure
    % xlabel('"time"'); ylabel('Zircon amount');
    % formatfigure;
end

%% Calculate average F path

maxpaths=70000;
minpathlength=5;

    
test=igncn1.SiO2>40&igncn1.SiO2<80;
path=0;
Fpath = zeros(1,200);
xi=linspace(0,1,200);
for kv=igncn1.Kv(test)'
    test=fulltzirc.Kv==kv;
    if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test)))) && min(fulltzirc.F(test))<35
        % Plot results for a random subset
        Fpath = Fpath + interp1(fulltzirc.Tscaled(test),fulltzirc.F(test),xi);
        path=path+1;
        if path>maxpaths
            break;
        end
    end
end
Fpath = Fpath./path;
figure; plot(xi,Fpath);

xlabel('"time"'); ylabel('F');
formatfigure;



%% Calculate average MZr path as a function of F

maxpaths=70000;
minpathlength=5;

    
path=0;
dist = zeros(1,100);
xi=linspace(1,100,100);
for kv=igncn1.Kv'
    test=fulltzirc.Kv==kv;
    if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.F(test)) == size(unique(fulltzirc.F(test))))
        
        % Plot results for a random subset
        dist = nansum([dist; interp1(fulltzirc.F(test),fulltzirc.MZr(test),xi)],1);
        path=path+1;
        if path>maxpaths
            break;
        end
    end
end
dist = dist./trapz(xi,dist);
figure; plot(xi,dist);

xlabel('Percent melt'); ylabel('Zircon amount');
formatfigure;


%% Calculate average normalized MZr path as a function of F for different silica ranges
maxpaths=70000;
minpathlength=5;
rsi=[45,55,65,75];

l=cell(length(rsi)-1,1);
figure; 
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(1,100,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.F(test)) == size(unique(fulltzirc.F(test)))) && min(fulltzirc.F(test))<35
            
            % Plot results for a random subset
            dist = nansum([dist; interp1(fulltzirc.F(test),fulltzirc.MZr(test),xi)],1);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end

xlabel('Percent melt'); ylabel('Zircon amount');
set(gca,'XDir','reverse');
legend(l)
formatfigure;


%% Calculate average MZr path as a function of F for different silica ranges

maxpaths=70000;
minpathlength=3;
% rsi=[43,51,62,74];
rsi=[45,55,65,75];

l=cell(length(rsi)-1,1);
figure; set(gca,'ColorOrder',lines(length(rsi)-1));
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100&igncn1.Age<2500;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(1,100,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.F(test)) == size(unique(fulltzirc.F(test)))) && min(fulltzirc.F(test))<35
            F=fulltzirc.F(test);       
            MZrn=fulltzirc.MZr(test)./([100; F(1:end-1)]-F)*2.009;
            % Plot results for a random subset
            dist = nansum([dist; interp1(F,MZrn,xi)],1);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./path;
%     dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end

for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100&igncn1.Age>2500;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(1,100,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.F(test)) == size(unique(fulltzirc.F(test)))) && min(fulltzirc.F(test))<35
            F=fulltzirc.F(test);
            MZrn=fulltzirc.MZr(test)./([100; F(1:end-1)]-F)*2.009;
            % Plot results for a random subset
            dist = nansum([dist; interp1(F,MZrn,xi)],1);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./path;
%     dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end

xlabel('Percent melt'); ylabel('Zircon saturated (ug/g/%)');
set(gca,'XDir','reverse');
legend(l)
formatfigure;

%% Calculate average MZr path as a function of T for different silica ranges

maxpaths=70000;
minpathlength=3;
% rsi=[43,51,62,74];
rsi=[45,55,65,75];

l=cell(length(rsi)-1,1);
figure; set(gca,'ColorOrder',lines(length(rsi)-1));
for i=1:length(rsi)-1
    test=igncn1.SiO2>rsi(i)&igncn1.SiO2<rsi(i+1)&igncn1.Elevation>-100;
    
    path=0;
    dist = zeros(1,100);
    xi=linspace(650,1000,100);
    for kv=igncn1.Kv(test)'
        test=fulltzirc.Kv==kv;
        if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.T(test)) == size(unique(fulltzirc.T(test)))) && min(fulltzirc.F(test))<35
            % Plot results for a random subset
            dist = nansum([dist; interp1(fulltzirc.T(test),fulltzirc.MZr(test)*2.009,xi)],1);
            path=path+1;
            if path>maxpaths
                break;
            end
        end
    end
    dist = dist./path;
%     dist = dist./trapz(xi,dist);
    hold on; plot(xi,dist);
    l{i} = [num2str(rsi(i)) ' - ' num2str(rsi(i+1)) ' % SiO2'];
end


xlabel('Temperature (C)'); ylabel('Zircon saturated (ug/g/10C)');
set(gca,'XDir','reverse');
legend(l)
formatfigure;