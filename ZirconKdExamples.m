% Percent crystals
X = 0:0.01:100;
% Percent melt
F = 100-X;

figure;
l = {};
for kd = [0.05 0.1 0.2 1];
    hold on; plot(F,100.*kd.*X./(kd.*X+F))
    l = [l; {['Bulk Kd: ' num2str(kd)]}];
end
set(gca,'xdir','reverse')
xlabel('Percent melt'); ylabel('Percent of total Zr partitioned into solid')
legend(l)
formatfigure;