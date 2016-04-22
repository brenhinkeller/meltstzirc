% Load data
cd ~/Desktop/meltstzirc/output/

if ~exist('igncn1','var'); load igncn1; end

name='tzircFull5F6kb';

% system(['grep -e ''^[0-9\.][0-9\.]*\(\t[0-9\.][0-9\.]*\)\{8\}$'' ' name '.log > ' name '.tsv']);

load([name '.tsv']);




% Fill variables
MZr=tzircFull5F6kb(:,8);
index=tzircFull5F6kb(:,1);
include=zeros(size(MZr));

% Zero-out negative masses
tzircFull5F6kb(MZr<0,8)=0; 
MZr(MZr<0)=0;

% Include all non-zero zircon masses and all adjacent zeros
MZr=[0;MZr;0];
index=[0;index;0];
include = MZr(2:end-1)>0 | (MZr(1:end-2)>0 & index(1:end-2)==index(2:end-1)) | (MZr(3:end)>0 & index(3:end)==index(2:end-1));

% Create the new dataset
fulltzirc.data=tzircFull5F6kb(include,:);
fulltzirc.elements={'Kv';'T';'F';'M';'SiO2';'Zr';'Zrsat';'MZr';'TSat'};
fulltzirc=elementify(fulltzirc);

figure; plot(fulltzirc.T, fulltzirc.MZr,'.')
set(gca, 'xdir','reverse')
xlabel('Temperature'); ylabel('Zircon mass saturated');
formatfigure

%%
maxkv=max(fulltzirc.Kv);

fulltzirc.Tscaled=NaN(size(fulltzirc.T));
fulltzirc.MZrscaled=NaN(size(fulltzirc.MZr));
fprintf('\n');

figure; hold on;

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
    
    % Display progress
    if mod(kv,100)==0
        bspstr=repmat('\b',1,floor(log10(kv-100))+1);
        fprintf(bspstr)
        fprintf('%i',kv)
    end
end

figure; plot(fulltzirc.Tscaled, fulltzirc.MZrscaled,'.')
xlabel('"time"'); ylabel('Zircon amount');
formatfigure;


%% Check out MZr paths by silica range

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

%% Check out MZr paths by age range

maxpaths=400;
rt=[0,100,1000,2500,4000];
minpathlength=5;

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

%% Calculate average MZr path

maxpaths=70000;
minpathlength=5;

    
test=igncn1.SiO2>40&igncn1.SiO2<80;
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
figure; plot(xi,dist);

xlabel('"time"'); ylabel('Zircon amount');
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
Fpath = zeros(1,100);
xi=linspace(0,1,100);
for kv=igncn1.Kv(test)'
    test=fulltzirc.Kv==kv;
    if sum(test)>minpathlength && rand < (10 * maxpaths / sum(test)) && all(size(fulltzirc.Tscaled(test)) == size(unique(fulltzirc.Tscaled(test))))
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

