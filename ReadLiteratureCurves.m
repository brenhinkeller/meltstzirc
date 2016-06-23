%% Define range
tmax=4500;
nsteps=450;
step=tmax/nsteps;
bincenters=((step:step:tmax)-step/2)';


%% Digitize Belousova curve

%RGB values of colors used in the figure
colors=[165  32  35 % Red (curve we want)
        138 188  65 % Green
         59  59  59 % Grey 
        255 255 255 % Black
          0   0   0 % White
];
colors=int16(colors);

% Corresponding output values
types=[1;0;0;0;0]; % We only want the red

% Read the plot
[x,y]=digitizePlotCurveLineBounds('BelousovaCurveImage.jpg', colors, types, 0, 4500, 0, 1);

% Smooth the results and rescale
y=smooth(y,10);
y=y/max(y);

% interpolate to make output arrtmax=4500;
nBelInt=interp1(x,y,bincenters);
figure; plot(bincenters,nBelInt);

nBel=nBelInt-[nBelInt(2:end); 0];


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
hold on; plot(bincenters,nDhuimeInt,'r');

nDhuime=nDhuimeInt-[nDhuimeInt(2:end); 0];

