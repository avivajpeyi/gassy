load spiral_out.mat out
plot(out.X, out.Y)
title('1M object sprilaing into r 15M object')
xlabel('X [units]')
ylabel('Y [units]')
saveas(gcf,'orbit','png')


load waveform_data.mat out
plot(out.h_plus)
hold on
plot(out.h_cross)
title('GW Strain from inspiral')
xlabel('Time [s]')
ylabel('Strain')
saveas(gcf,'strain','png')
hold off