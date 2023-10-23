clear

pH4_B1 = [99.86712633 100 90.94448856];
err_4B = [8.232547361 2.191165779 2.938348373];
pH9_B1 = [100 49.34476956 33.47918137];
err_9B = [10.96842443 2.182214705 3.673989034];
pH4_G2 = [95.84154922 100 96.22647317];
err_4G = [2.505998945 1.252665961 1.442760582];
pH9_G2 = [100 1.93764831 1.979107002];
err_9G = [0.281777495 0.255867609 0.176645583];

TP = [0 24 48];

figure
    hold on
    %plot(TP,pH4_B1)
    errorbar(TP,pH4_B1,err_4B)
    %plot(TP,pH9_B1)
    errorbar(TP,pH9_B1,err_9B)
    xlim([0 48])
    xticks([0 12 24 36 48])
    ylim([0 100])
    xlabel('Time (hrs)')
    ylabel('Normalized RFU')
    legend('pH 4','pH 9')
    title('AFB_1 degradation by buffered medium')
    
figure
    hold on
    %plot(TP,pH4_G2)
    errorbar(TP,pH4_G2,err_4G)
    %plot(TP,pH9_G2)
    errorbar(TP,pH9_G2,err_9G)
    xlim([0 48])
    xticks([0 12 24 36 48])
    ylim([0 100])
    xlabel('Time (hrs)')
    ylabel('Normalized RFU')
    legend('pH 4','pH 9')
    title('AFG_2 degradation by buffered medium')
    