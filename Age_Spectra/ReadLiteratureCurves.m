%% Define range
tmax=4500;
nsteps=450;
step=tmax/nsteps;
bincenters=((step:step:tmax)-step/2)';


%% Digitize Belousova curve
% 
% %RGB values of colors used in the figure
% colors=[165  32  35 % Red (curve we want)
%         138 188  65 % Green
%          59  59  59 % Grey 
%         255 255 255 % Black
%           0   0   0 % White
% ];
% colors=int16(colors);
% 
% % Corresponding output values
% types=[1;0;0;0;0]; % We only want the red
% 
% % Read the plot
% [x,y]=digitizePlotCurveLineBounds('BelousovaCurveImage.jpg', colors, types, 0, 4500, 0, 1);
% 
% % Smooth the results and rescale
% y=smooth(y,10);
% y=y/max(y);
% 
% % interpolate to make output arrtmax=4500;
% nBelInt=interp1(x,y,bincenters);
% figure; plot(bincenters,nBelInt);
% 
% nBel=nBelInt-[nBelInt(2:end); 0];


%% Digitize Dhuime curve

%RGB values of colors used in the figure
colors=[209  24  31 % Red (curve we want)
          0  88 168 % Blue
        226 117  88 % Grey 
        255 255 255 % Black
          0   0   0 % White
];
colors=int16(colors);

% Corresponding output values
types=[1;0;0;0;0]; % We only want the red

% Read the plot
[x,y]=digitizePlotCurveLineBounds('DhuimeCurveImage.png', colors, types, 0, 4500, 0, 1, 'HorizontalFirst');

% Smooth the results
x=smooth(x,20);

% interpolate to make output array
nDhuimeInt=interp1(x,y,bincenters);
figure; plot(bincenters,nDhuimeInt,'k');

nDhuime=nDhuimeInt-[nDhuimeInt(2:end); 0];

% %% Apply model age correction to Dhuime
% 
% nDhuimeCorrected = nDhuime.*Tdmtocrust;
% nDhuimeCorrected = nDhuimeCorrected./sum(nDhuimeCorrected);
% nDhuimeIntCorrected = cumsum(nDhuimeCorrected,'reverse');
% hold on; plot(bincenters,nDhuimeIntCorrected)

%% Dhuime with 50% of modern crust at 4.5 Ga

nDhuime50i = nDhuime./sum(nDhuime)./2;
nDhuime50i(end)=0.5;
nDhuime50iInt = cumsum(nDhuime50i,'reverse');
hold on; plot(bincenters,nDhuime50iInt)

% However, we are hesitant to reccomend these corrections as a simple
% adjustment to existing crustal growth curves, as we believe there are
% more fundamental issues at stake regarding/concerning the unavoidable choice of
% assumptions/boundary conditions in such models (supplementary discussion).

% Particularly, while such models typically account for magmatic reworking
% within the crust, they typically do not 

% Particularly, zirconium content and M value only uniquely constrain
% zircon crystallization mass in the case of closed-system solidification.
% While typically acceptible for most crustal magmatic rocks today, this
% assumption may not be justified in early Earth history.

% No initial crust
% Include reworking but no recycling of crust to the mantle 



%% Dhume with 30% coverage of flotation crust at 4.5 Ga and 1 continental volume net crust destruction / recycling to mantle


nDhuime100i = nDhuime./sum(nDhuime);
nDhuime100i(end)=1;
nDestruction = ones(size(nDhuime))./length(nDhuime);

nDhuime100irInt = cumsum(nDhuime100i-nDestruction,'reverse');
hold on; plot(bincenters,nDhuime100irInt)

% Ditto and exponential recycling to mantle with 1 Gyr e-folding time

nDestruction = flip(exp(-bincenters/700));
nDestruction = nDestruction./sum(nDestruction); % Normalize
hold on; plot(bincenters,cumsum(nDhuime100i-nDestruction,'reverse'))

% Ditto and exponential recycling to mantle with 4.5 Gyr e-folding time

nDestruction = flip(exp(-bincenters/1400));
nDestruction = nDestruction./sum(nDestruction); % Normalize
hold on; plot(bincenters,cumsum(nDhuime100i-nDestruction,'reverse'))


%% Dhume with 50% coverage of flotation crust at 4.5 Ga and linear recycling


nDhuime160i = nDhuime./sum(nDhuime);
nDhuime160i(end)=6/3;
nDestruction = ones(size(nDhuime))./length(nDhuime)*6/3;

% hold on; plot(bincenters, cumsum(nDhuime160i-nDestruction,'reverse'))

% Ditto and exponential recycling to mantle with 1 Gyr e-folding time
nDestruction = flip(exp(-bincenters/700));
nDestruction = nDestruction./sum(nDestruction).*6/3; % Normalize
hold on; plot(bincenters,cumsum(nDhuime160i-nDestruction,'reverse'))

% Ditto and exponential recycling to mantle with 4.5 Gyr e-folding time
nDestruction = flip(exp(-bincenters/1400));
nDestruction = nDestruction./sum(nDestruction).*6/3; % Normalize
hold on; plot(bincenters,cumsum(nDhuime160i-nDestruction,'reverse'))

%% Dhume with 75% coverage of flotation crust at 4.5 Ga and linear recycling


nDhuime160i = nDhuime./sum(nDhuime);
nDhuime160i(end)=2.5;
nDestruction = ones(size(nDhuime))./length(nDhuime)*2.5;

% hold on; plot(bincenters, cumsum(nDhuime160i-nDestruction,'reverse'))

% % Ditto and exponential recycling to mantle with 1 Gyr e-folding time
% nDestruction = flip(exp(-bincenters/700));
% nDestruction = nDestruction./sum(nDestruction).*2.5; % Normalize
% hold on; plot(bincenters,cumsum(nDhuime160i-nDestruction,'reverse'))

% % Ditto and exponential recycling to mantle with 4.5 Gyr e-folding time
% nDestruction = flip(exp(-bincenters/4500));
% nDestruction = nDestruction./sum(nDestruction).*2.5; % Normalize
% hold on; plot(bincenters,cumsum(nDhuime160i-nDestruction,'reverse'))

%%
legend('Dhuime, 2012','with 5km anorthosite, no recycling','with 10km anorthosite, linear recycling','with 10km anorthosite, exponential recycling l=0.7Gyr','with 10km anorthosite, exponential recycling l=1.4Gyr','with 20km anorthosite, exponential recycling l=0.7Gyr','with 20km anorthosite, exponential recycling l=1.4Gyr');
hold on; plot([0 4500],[1,1],'--k')
xlabel('Age (Ma)'); ylabel('Crust fraction relative to present day');
ylim([0 2.05])

