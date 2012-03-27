figure(3)
clf;
set(3,'Name','Relative Phase Diagrams');
subplot(2,2,1)
phaseSpace = phaseSpace_spgl1;
phaseSpace(~isfinite(phaseSpace)) = 10;
imagesc(delta,rho,log10(phaseSpace)')
caxis([-7,0]);
set(gca,'YDir','normal')
xlabel('Undersampling, \delta = n / N');
ylabel('Sparsity, \rho = k / n');
colorbar;
title('SPGl1');
subplot(2,2,3)
phaseSpace = phaseSpace_ist;
phaseSpace(~isfinite(phaseSpace)) = 10;
imagesc(delta,rho,log10(phaseSpace)')
caxis([-7,0]);
set(gca,'YDir','normal')
xlabel('Undersampling, \delta = n / N');
ylabel('Sparsity, \rho = k / n');
colorbar;
title('IST');
subplot(2,2,4)
phaseSpace = phaseSpace_amp;
phaseSpace(~isfinite(phaseSpace)) = 10;
imagesc(delta,rho,log10(phaseSpace)')
caxis([-7,0]);
set(gca,'YDir','normal')
colorbar
xlabel('Undersampling, \delta = n / N');
ylabel('Sparsity, \rho = k / n');
title('AMP');
