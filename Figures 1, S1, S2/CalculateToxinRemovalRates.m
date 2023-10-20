clear
load('pH_ERY_RFUtoToxin')

%Time of analysis (how long the plate was run) in hours
tf = 48;

%Dimensions of the plate; i=rows, j=columns
i = 8;
j = 12;

dt = 5/60; %Time between read points in hours
T_BGg = mean(FLT1(7:8,10,1:Nr),1); %pH 7.0 G2
T_BGg4 = mean(FLT1(5:6,10,1:Nr),1); %pH 4.0 G2 
T_BGg5 = mean(FLT1(5:6,11,1:Nr),1); %pH 5.0 G2
T_BGg6 = mean(FLT1(5:6,12,1:Nr),1); %pH 6.0 G2
T_BGg8 = mean(FLT1(7:8,11,1:Nr),1); %pH 8.0 G2
T_BGg8t = mean(FLT1(7:8,12,1:Nr),1); %pH 8.0 in Tri-HCl G2
T_BGg9 = mean(FLT1(8,5:7,1:Nr),2); %pH 9.0 G2
T_BGb = mean(FLT2(4,10,1:Nr),1); %pH 7.0 B1
T_BGb4 = mean(FLT2(2:3,10,1:Nr),1); %pH 4.0 B1
T_BGb5 = mean(FLT2(2:3,11,1:Nr),1); %pH 5.0 B1
T_BGb6 = mean(FLT2(2:3,12,1:Nr),1); %pH 6.0 B1
T_BGb8 = mean(FLT2(4,11,1:Nr),1); %pH 8.0 B1
T_BGb8t = mean(FLT2(4,12,1:Nr),1); %pH 8.0 in Tri-HCl B1
T_BGb9 = mean(FLT2(8,2:4,1:Nr),2); %pH 9.0 B1

%Calculate percent degradation
pdegn = zeros(8,12);
pdeg_minn = zeros(8,12);
FLTn = zeros(8,12,Nr);
for q = 2:4
    for v = 2:9
        FLTn(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb(1,1,1)*T_BGb(1,1,1:Nr));
        start = FLTn(q,v,1);
        final = FLTn(q,v,Nr);
        min_toxin = min(FLTn(q,v,1:Nr));
        pdegn(q,v) = (start - final)/start;
        pdeg_minn(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 2:9
        FLTn(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg(1,1,1)*T_BGg(1,1,1:Nr));
        start = FLTn(q,v,1);
        final = FLTn(q,v,Nr);
        min_toxin = min(FLTn(q,v,1:Nr));
        pdegn(q,v) = (start - final)/start;
        pdeg_minn(q,v) = (start - min_toxin)/start;
    end
end
%disp(pdegn)
disp(pdeg_minn)

pdegpH = zeros(8,12);
pdeg_minpH = zeros(8,12);
FLTpH = zeros(8,12,Nr);
for q = 2:4
    for v = 3
        FLTpH(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb4(1,1,1)*T_BGb4(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 2:4
    for v = 4
        FLTpH(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb5(1,1,1)*T_BGb5(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 2:4
    for v = 5
        FLTpH(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb6(1,1,1)*T_BGb6(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 2:4
    for v = 7
        FLTpH(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb8(1,1,1)*T_BGb8(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 2:4
    for v = 8
        FLTpH(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb8t(1,1,1)*T_BGb8t(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 2:4
    for v = 9
        FLTpH(q,v,:) = FLT2(q,v,1:Nr)./(1/T_BGb9(1,1,1)*T_BGb9(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 3
        FLTpH(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg4(1,1,1)*T_BGg4(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 4
        FLTpH(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg5(1,1,1)*T_BGg5(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 5
        FLTpH(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg6(1,1,1)*T_BGg6(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 7
        FLTpH(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg8(1,1,1)*T_BGg8(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 8
        FLTpH(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg8t(1,1,1)*T_BGg8t(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
for q = 5:7
    for v = 9
        FLTpH(q,v,:) = FLT1(q,v,1:Nr)./(1/T_BGg9(1,1,1)*T_BGg9(1,1,1:Nr));
        start = FLTpH(q,v,1);
        final = FLTpH(q,v,Nr);
        min_toxin = min(FLTpH(q,v,1:Nr));
        pdegpH(q,v) = (start - final)/start;
        pdeg_minpH(q,v) = (start - min_toxin)/start;
    end
end
%disp(pdegn)
disp(pdeg_minpH)

%FLT_subcell = zeros(8,12);
%FLT_subcell(2,4,1:Nr) = FLT(2,4,1:Nr) - FLT(2,9,1:Nr);

for q = 2:7
    for v = 3:9
         percent_tox(q,v,:) = mean(FLTpH(q,v,:),1)./max(mean(FLTpH(q,v,:),1))*100;
         per_deg (q,v) = min(percent_tox(q,v,:));
    end
end
for q = 2:7
    for v = 2
         percent_tox(q,v,:) = mean(FLTn(q,v,:),1)./max(mean(FLTn(q,v,:),1))*100;
         per_deg (q,v) = min(percent_tox(q,v,:));
    end
end
for q = 2:7
    for v = 6
         percent_tox(q,v,:) = mean(FLTn(q,v,:),1)./max(mean(FLTn(q,v,:),1))*100;
         per_deg (q,v) = min(percent_tox(q,v,:));
    end
end

%% Figures

figure %Plotting degradation as a function of loss of toxin conc. over time
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,1,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:3,10,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:3,11,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:3,12,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(4,10,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(4,11,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(4,12,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(8,2:3,1:Nr),2),2),15))
    xlabel('Time (hrs)')
    ylabel('AFB_1 Toxin Conc. (ug/ml)')
    xlim([0 tf]) 
    xticks([0 12 24 36 48])
    legend('Z medium','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFB_1 at different pH')
    
figure %Plotting degradation as a function of loss of toxin conc. over time
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,1,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:6,10,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:6,11,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:6,12,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(7:8,10,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(7:8,11,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(7:8,12,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(8,4:5,1:Nr),2),2),15))
    xlabel('Time (hrs)')
    ylabel('AFG_2 Toxin Conc. (ug/ml)')
    xlim([0 tf]) 
    xticks([0 12 24 36 48])
    legend('Z medium','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFG_2 at different pH')

% figure %Plotting degradation as a function of loss of toxin conc. over time
%     hold on
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,2,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(2:4,3,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(2:4,4,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(2:4,5,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,6,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(2:4,7,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(2:4,8,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(2:4,9,1:Nr),1),2),15))
%     xlabel('Time (hrs)')
%     ylabel('AFB_1 Toxin Conc. (ug/ml)')
%     xlim([0 tf]) 
%     xticks([0 12 24 36 48 60 72])
%     legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
%     title('AFB_1 degradation by R. ery filtrate at different pH')
    
figure %Plotting degradation as a function of loss of toxin conc. over time
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,2,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,3,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,4,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,5,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,6,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,7,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,8,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(2:4,9,1:Nr),1),2),15))
    xlabel('Time (hrs)')
    ylabel('AFB_1 Toxin Conc. (ug/ml)')
    xlim([0 tf]) 
    xticks([0 12 24 36 48])
    legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFB_1 degradation by R. ery filtrate at different pH')
    
figure %Plotting degradation as a function of loss of toxin conc. over time
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,2,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,3,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,4,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,5,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,6,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,7,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,8,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,9,1:Nr),1),2),15))
    xlabel('Time (hrs)')
    ylabel('AFB_1 Toxin Conc. (ug/ml)')
    xlim([0 tf]) 
    xticks([0 12 24 36 48])
    legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFB_1 degradation by R. ery filtrate at different pH')
%     
% figure %Plotting degradation as a function of loss of toxin conc. over time
%     hold on
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,2,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(5:7,3,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(5:7,4,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(5:7,5,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,6,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(5:7,7,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(5:7,8,1:Nr),1),2),15))
%     plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTpH(5:7,9,1:Nr),1),2),15))
%     xlabel('Time (hrs)')
%     ylabel('AFG_2 Toxin Conc. (ug/ml)')
%     xlim([0 tf]) 
%     xticks([0 12 24 36 48 60 72])
%     legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
%     title('AFG_2 degradation by R. ery filtrate at different pH')
    
figure %Plotting degradation as a function of loss of toxin conc. over time
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,2,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,3,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,4,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,5,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,6,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,7,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,8,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLTn(5:7,9,1:Nr),1),2),15))
    xlabel('Time (hrs)')
    ylabel('AFG_2 Toxin Conc. (ug/ml)')
    xlim([0 tf]) 
    xticks([0 12 24 36 48])
    legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFG_2 degradation by R. ery filtrate at different pH')   

figure %Plotting degradation as a function of loss of toxin conc. over time
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,2,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,3,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,4,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,5,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,6,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,7,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,8,1:Nr),1),2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,9,1:Nr),1),2),15))
    xlabel('Time (hrs)')
    ylabel('AFG_2 Toxin Conc. (ug/ml)')
    xlim([0 tf]) 
    xticks([0 12 24 36 48])
    legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFG_2 degradation by R. ery filtrate at different pH')       
    
figure
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,2,1:Nr),1)./max(mean(FLT1(5:7,2,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,3,1:Nr),1)./max(mean(FLT1(5:7,3,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,4,1:Nr),1)./max(mean(FLT1(5:7,4,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,5,1:Nr),1)./max(mean(FLT1(5:7,5,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,6,1:Nr),1)./max(mean(FLT1(5:7,6,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,7,1:Nr),1)./max(mean(FLT1(5:7,7,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,8,1:Nr),1)./max(mean(FLT1(5:7,8,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:7,9,1:Nr),1)./max(mean(FLT1(5:7,9,1:Nr),1))*100,2),15))
    xlabel('Time (hrs)')
    ylabel('Toxin remaining (%)')
    ylim([0 100])
    xlim([0 tf]) 
    xticks([0 12 24 36 48 60 72])
    legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFG_2 degradation by R. ery filtrate at different pH')
    
figure
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,2,1:Nr),1)./max(mean(FLT2(2:4,2,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,3,1:Nr),1)./max(mean(FLT2(2:4,3,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,4,1:Nr),1)./max(mean(FLT2(2:4,4,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,5,1:Nr),1)./max(mean(FLT2(2:4,5,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,6,1:Nr),1)./max(mean(FLT2(2:4,6,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,7,1:Nr),1)./max(mean(FLT2(2:4,7,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,8,1:Nr),1)./max(mean(FLT2(2:4,8,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:4,9,1:Nr),1)./max(mean(FLT2(2:4,9,1:Nr),1))*100,2),15))
    xlabel('Time (hrs)')
    ylabel('Toxin remaining (%)')
    ylim([0 100])
    xlim([0 tf]) 
    xticks([0 12 24 36 48 60 72])
    legend('Normal filtrate','pH 4.0','pH 5.0','pH 6.0','pH 7.0','pH 8.0','pH 8.0','pH 9.0')
    title('AFB_1 degradation by R. ery filtrate at different pH')
    
figure
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(5:6,10,1:Nr),1)./max(mean(FLT1(5:6,10,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(7:8,10,1:Nr),1)./max(mean(FLT1(7:8,10,1:Nr),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT1(8,4:5,1:Nr),2)./max(mean(FLT1(8,4:5,1:Nr),2))*100,2),15))
    xlabel('Time (hrs)')
    ylabel('Normalized RFU (%)')
    ylim([0 100])
    xlim([0 tf]) 
    xticks([0 12 24 36 48 60 72])
    legend('pH 4.0','pH 7.0','pH 9.0')
    title('AFG_2 at different pH')
    
figure
    hold on
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(2:3,10,1:Nr),1)./(mean(FLT2(2:3,10,1),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(4,10,1:Nr),1)./(mean(FLT2(4,10,1),1))*100,2),15))
    plot(dt*(0:Nr-1),smooth(shiftdim(mean(FLT2(8,2:3,1:Nr),2)./(mean(FLT2(8,2:3,1),2))*100,2),15))
    xlabel('Time (hrs)')
    ylabel('Normalized RFU (%)')
    ylim([0 400])
    xlim([0 tf]) 
    xticks([0 12 24 36 48 60 72])
    legend('pH 4.0','pH 7.0','pH 9.0')
    title('AFB_1 at different pH')

save('pH_ERY_degeff')