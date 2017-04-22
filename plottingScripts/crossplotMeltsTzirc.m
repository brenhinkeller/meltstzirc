
    
%% Zircon saturation systematics cross-plots
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