%% Define range and load data
tmax=4500;
nsteps=450;
step=tmax/nsteps;
bincenters=((step:step:tmax)-step/2)';

% % Digitize Dhuime curve
% 
% %RGB values of colors used in the figure
% colors=[209  24  31 % Red (curve we want)
%           0  88 168 % Blue
%         226 117  88 % Grey 
%         255 255 255 % Black
%           0   0   0 % White
% ];
% colors=int16(colors);
% 
% % Corresponding output values
% types=[1;0;0;0;0]; % We only want the red
% 
% % Read the plot
% [x,y]=digitizePlotCurveLineBounds('DhuimeCurveImage.png', colors, types, 0, 4500, 0, 1, 'HorizontalFirst');
% 
% % Smooth the results
% x=smooth(x,20);
% 
% % interpolate to make output array
% nDhuimeInt=interp1(x,y,bincenters);
% figure; plot(bincenters,nDhuimeInt,'k');
% 
% nDhuime=nDhuimeInt-[nDhuimeInt(2:end); 0];

if ~exist('nDhuime','var'); load nDhuime; end;
figure;

%% Dhuime crustal growth curve modified assuming 0.5 modern crustal volumes of flotation crust at 4.5 Ga

nDhuime50i = nDhuime./sum(nDhuime)./2;
nDhuime50i(end)=0.5;
nDhuime50iInt = cumsum(nDhuime50i,'reverse');
hold on; plot(bincenters,nDhuime50iInt)


%% With 10km global flotation crust at 4.5 Ga

nDhuime10i = nDhuime./sum(nDhuime);
nDhuime10i(end)=1;

% With linear recycling
nDestruction = ones(size(nDhuime))./length(nDhuime);
hold on; plot(bincenters,cumsum(nDhuime10i-nDestruction,'reverse'))

% As above but with exponential recycling to mantle with 0.7 Gyr e-folding time
nDestruction = flip(exp(-bincenters/700));
nDestruction = nDestruction./sum(nDestruction); % Normalize
hold on; plot(bincenters,cumsum(nDhuime10i-nDestruction,'reverse'))

% As above but with exponential recycling to mantle with 1.4 Gyr e-folding time
nDestruction = flip(exp(-bincenters/1400));
nDestruction = nDestruction./sum(nDestruction); % Normalize
hold on; plot(bincenters,cumsum(nDhuime10i-nDestruction,'reverse'))


%% With 20km global flotation crust at 4.5 Ga

nDhuime20i = nDhuime./sum(nDhuime);
nDhuime20i(end)=6/3;

% % With linear recycling
% nDestruction = ones(size(nDhuime))./length(nDhuime)*6/3;
% hold on; plot(bincenters, cumsum(nDhuime20i-nDestruction,'reverse'))

% With exponential recycling to mantle with 0.7 Gyr e-folding time
nDestruction = flip(exp(-bincenters/700));
nDestruction = nDestruction./sum(nDestruction).*6/3; % Normalize
hold on; plot(bincenters,cumsum(nDhuime20i-nDestruction,'reverse'))

% With exponential recycling to mantle with 1.4 Gyr e-folding time
nDestruction = flip(exp(-bincenters/1400));
nDestruction = nDestruction./sum(nDestruction).*6/3; % Normalize
hold on; plot(bincenters,cumsum(nDhuime20i-nDestruction,'reverse'))

%%
legend('Dhuime, 2012','with 5km anorthosite, no recycling','with 10km anorthosite, linear recycling','with 10km anorthosite, exponential recycling l=0.7Gyr','with 10km anorthosite, exponential recycling l=1.4Gyr','with 20km anorthosite, exponential recycling l=0.7Gyr','with 20km anorthosite, exponential recycling l=1.4Gyr');
hold on; plot([0 4500],[1,1],'--k')
xlabel('Age (Ma)'); ylabel('Crust fraction relative to present day');
ylim([0 2.05])

