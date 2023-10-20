clear

%This is for the conversion of RFU data into toxin concentration for AFG2

load('pH_ERY_07282022')

%Size of well matrix
q = 8;
v = 12;
%Matrix of zeros the size of well matrix and number of time reads
FLT1 = zeros(q,v,Nr);
FLT2 = zeros(q,v,Nr);

%AFG2
FLT1 = 54*(FL1- 27*ones(8,12,577))./(60027*ones(8,12,577) - FL1);

%AFB1
for q = 1:8
    for v = 1:12
        for Nr = 1:577
            FLT2(q,v,Nr) = 300*(FL2(q,v,Nr)- 275)/(51975 - FL2(q,v,Nr));
        end
    end
end
%Simplier way to write the same code as above
% FLT2 = 54*(FL1- 27*ones(8,12,865))./(60027*ones(8,12,865) - FL1);

%Check of random read/time point
disp(FLT1(2,2,200))

%Saves the data in this file name, must change it each time
save('pH_ERY_RFUtoToxin')