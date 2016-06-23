%% Mass of zircon saturated over time, combined

cd ~/Desktop/meltstzirc/output/
% n = {'5F1kb3H2O','5F1kb1H2O','5F6kb3H2O','5F12kb3H2O','5F12kb1H2O','5F12kb01H2O'};
% n = {'5F1kb3H2O','5F4kb2H2O','5F12kb2H2O','5F12kb01H2O'};
% n = {'5F1kb3H2O','5F8kb3H2O','5F12kb3H2O','5F12kb01H2O'};

% n = {'5F4kb2H2O','5F12kb4H2O','5F12kb2H2O','5F12kb1H2O'};
n = {'5F2kb4H2O','5F4kb2H2O','5F12kb2H2O','5F12kb1H2O'};


Elem='MZrn';
rsi=[40 80];
agemin=0;
agemax=4000;

figure;
for j = 1:length(n)
    if ~exist('mcigncn1','var') || ~exist('igncn1','var') || ~exist('suffix','var') || ~isequal(suffix,n{j})
        suffix = n{j};
        load(['mcigncn1-' suffix]);
        load(['igncn1-' suffix]);
    end 
    if ~isfield(igncn1,'MZrn') || ~isfield(mcigncn1,'MZrn')
        igncn1.MZrn = igncn1.MZr * 2.009;
        mcigncn1.MZrn = mcigncn1.MZr * 2.009;
    end

    test=mcigncn1.SiO2>rsi(1)&mcigncn1.SiO2<rsi(2)&mcigncn1.Elevation>-100;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),40);
    hold on; errorbar(c,m,e,'.','Color',[1-j/(length(n)),0,j/(length(n))])
    
%     [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),400,10);
%     hold on; plot(c,m,'-','Color',[1-j/(length(n)),0,j/(length(n))]) 
end

xlabel('Age (Ma)'); ylabel(Elem); legend(n);
formatfigure


%% Calculate zircon-to-crust mapping factor / corrected zircon age spectra

% n = {'5F2kb4H2O','5F4kb2H2O','5F12kb4H2O','5F12kb2H2O','5F12kb1H2O'};
n = {'5F2kb4H2O','5F4kb2H2O','5F12kb2H2O','5F12kb1H2O'};


Elem='MZr';
rsi=[40 80];
agemin=0;
agemax=4000;

figure; hold on;
if ~exist('BelousovaSpectrum','var'); load BelousovaSpectrum; end
plot(BelousovaSpectrum(:,1), BelousovaSpectrum(:,2)./max(BelousovaSpectrum(1:50,2)));

for j = 1:length(n)
    if ~exist('mcigncn1','var') || ~exist('igncn1','var') || ~exist('suffix','var') || ~isequal(suffix,n{j})
        suffix = n{j};
        load(['mcigncn1-' suffix]);
        load(['igncn1-' suffix]);
    end
    
    for i=1:length(rsi)-1
        test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100&mcigncn1.Ff<70;
        [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),400,10);
    end
    
    CorrectedSpectrum=BelousovaSpectrum;
    for k=1:length(CorrectedSpectrum)
        [~,i]=min(abs(c-CorrectedSpectrum(k,1))); % For the closest sample in age
        CorrectedSpectrum(k,2)=CorrectedSpectrum(k,2)/m(i); % Apply the correction
    end
    
    plot(CorrectedSpectrum(:,1), CorrectedSpectrum(:,2)./max(CorrectedSpectrum(1:50,2)))
end


xlabel('Age (Ma)'); ylabel('Sample density'); legend([{'Uncorrected'}, n]);
set(gca,'YTick',[]);
formatfigure
% xlim([0,4500])
% ylim([0,1.15])



