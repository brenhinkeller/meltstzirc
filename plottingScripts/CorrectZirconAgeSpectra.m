%% Calculate Voice age spectrum
%     [N,Age] = ksdensity(voice.Best_Age,0:1:4500,'bandwidth',10);
%     VoiceSpectrum.Age = Age;
%     VoiceSpectrum.N = N;

    if ~exist('VoiceSpectrum','var'); load VoiceSpectrum;  end
    RawSpectrum = VoiceSpectrum;
    

%% Apply correction to zircon age spectrum
    
%     ns = {'5F2kb4H2O','5F4kb2H2O','5F12kb2H2O','5F12kb1H2O'};
    ns = {'5F2kb4H2O','5F12kb1H2O'};

	figure;
    hold on;
    plot(RawSpectrum.Age, RawSpectrum.N./max(RawSpectrum.N(1:120)))
    
    for n = ns

	% Load data
        if ~exist('mcigncn1','var') || ~exist('igncn1','var') || ~exist('suffix','var') || ~isequal(suffix,n{:})
            suffix = n{:};
            load(['mcigncn1-' suffix]);
            load(['igncn1-' suffix]);
        end
        
	% Mass of zircon saturated over time, combined
        Elem='MZr';
        agemin=0;
        agemax=4000;
        rsi=[40 80];
        for i=1:length(rsi)-1
            test=mcigncn1.SiO2>rsi(i)&mcigncn1.SiO2<rsi(i+1)&mcigncn1.Elevation>-100&mcigncn1.Ff<70;
            [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(Elem)(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),400,10);
        end
        
	% Start with plain uncorrected spectrum
        CorrectedSpectrum=RawSpectrum;
        
	% Linear interpolation
        ClosestM = interp1(c,m,CorrectedSpectrum.Age);
        CorrectedSpectrum.N = CorrectedSpectrum.N./ClosestM;
        
	% Present-day-normalized comparison
        plot(CorrectedSpectrum.Age, CorrectedSpectrum.N./max(CorrectedSpectrum.N(1:120)))
        
    end

    xlabel('Age (Ma)'); ylabel('Adjusted Zircon Abundance');
    set(gca,'YTick',[]);
    legend([{'Original'} ns])
    formatfigure
    xlim([0,4000])
    ylim([0,1.4])

