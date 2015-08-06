 %% Import tzirc data
if ~exist('igncn1','var'); load igncn1; end

% Make struct from input file 
for var={'Kv','Mbulk','Tliq','Tsat','Tsatb','Zrsat','Zrf','Ff','SiO2','Zrbulk','MZr'};
    eval(['tzirclog.' var{:} '=' var{:}])
end

% Create new struct fields for imported data
variables={'Mbulk','Tliq','Tsat','Tsatb','Zrsat','Zrf','Ff','MZr'};
for var=variables;
    igncn1.(var{:})=NaN(size(igncn1.Kv));
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
    
    
    
%%
t=igncn1.MZr>0;
figure; plot(igncn1.SiO2(~t), igncn1.MZr(~t),'.r');
hold on; plot(igncn1.SiO2(t), igncn1.MZr(t),'.b');
xlabel('Bulk SiO2 (%)'); ylabel('Mass of zircon saturated (ug/g)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.Zr(~t), igncn1.MZr(~t),'.r');
hold on; plot(igncn1.Zr(t), igncn1.MZr(t),'.b');
xlabel('Bulk Zr (ppm)'); ylabel('Mass of zircon saturated (ug/g)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.Tsatb(~t), igncn1.MZr(~t),'.r');
hold on; plot(igncn1.Tsatb(t), igncn1.MZr(t),'.b');
xlabel('Bulk zircon saturation temp (C)'); ylabel('Mass of zircon saturated (ug/g)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.Tsatb(~t), igncn1.Tsat(~t),'.r');
hold on; plot(igncn1.Tsatb(t), igncn1.Tsat(t),'.b');
xlabel('Bulk zircon saturation temp (C)'); ylabel('Adjusted zircon saturation temp (C)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.Tsat(~t), igncn1.MZr(~t),'.r');
hold on; plot(igncn1.Tsat(t), igncn1.MZr(t),'.b');
xlabel('Adjusted zircon saturation temp (C)'); ylabel('Mass of zircon saturated (ug/g)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.Tliq(~t), igncn1.Tsat(~t),'.r');
hold on; plot(igncn1.Tliq(t), igncn1.Tsat(t),'.b');
xlabel('Liquidus temperature (C)'); ylabel('Adjusted zircon saturation temp (C)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.Tliq(~t), igncn1.Tsatb(~t),'.r');
hold on; plot(igncn1.Tliq(t), igncn1.Tsatb(t),'.b');
xlabel('Liquidus temperature (C)'); ylabel('Bulk zircon saturation temp (C)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.SiO2(~t), igncn1.Tsatb(~t),'.r');
hold on; plot(igncn1.SiO2(t), igncn1.Tsatb(t),'.b');
xlabel('Bulk SiO2 (%)'); ylabel('Bulk zircon saturation temp (C)');
formatfigure

t=igncn1.MZr>0;
figure; plot(igncn1.SiO2(~t), igncn1.Tsat(~t),'.r');
hold on; plot(igncn1.SiO2(t), igncn1.Tsat(t),'.b');
xlabel('Bulk SiO2 (%)'); ylabel('Adjusted zircon saturation temp (C)');
formatfigure
