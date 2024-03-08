% Efficient Portfolio using quadprog function

clear
clc

tic

%% Question 6
%%5 Factor Proxies

FactorsFivetable = readtable('Fama-French 5 Factor Proxies.xlsx','Range', 'B1:K127');
FactorsFivedata=table2array(FactorsFivetable)
FactorsFiveRF= 0.0036
FactorsFive= FactorsFivedata(4:end, 3:2:end) / 100;
FactorsFive=FactorsFive 

FactorsFive = rmmissing(FactorsFive)

up=max(mean(FactorsFive))
[FivemeanGMVP,FivestddevGMVP,FivemeanTP,FivestdDevTP,Fivevarpor_r,Fiveretport_r,FiveretMeansTarget,EFretFive,EFvarFive]=Factors (FactorsFive,0.35,0.04)

figure(1);
plot(sqrt(Fivevarpor_r),FiveretMeansTarget,'b-');
    hold on;
    plot(sqrt(EFvarFive),EFretFive,'g-')
    cal_xFive = [0,FivestdDevTP];
    cal_yFive = [0,FivemeanTP];
    stdDevExtendedFive = FivestdDevTP * 1.8; % Example: extend the line 50% further out
    meanTPExtendedFive = (FivemeanTP / FivestdDevTP) * stdDevExtendedFive; % Use the slope to calculate the new y-value
    plot([0, stdDevExtendedFive], [0, meanTPExtendedFive], 'r--');% Extended CAL with dashed line for differentiation
    plot([FivestdDevTP, stdDevExtendedFive], [FivemeanTP, meanTPExtendedFive], 'k-');
    plot(FivestdDevTP,FivemeanTP, 'g.', 'MarkerSize', 20)
    plot(FivestddevGMVP,FivemeanGMVP, 'm.', 'MarkerSize', 20)

    xlabel('Standard Deviation (Risk)');
    ylabel('Excess Returns');
    title('FF 5-Factor Proxy Efficient Frontier with Tangent Portfolio and CML');
    legend('MVF','5-Factor Proxy Efficient Frontier','CML','Portfolio Efficient Frontier');
    text(FivestdDevTP,FivemeanTP , 'Tangency Portfolio', 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    text(FivestddevGMVP,FivemeanGMVP,'GMVP','VerticalAlignment','bottom', 'HorizontalAlignment','left' )
%% Calling my function can return GMVP value, don't need this part


%% 3 Factor Proxy
FactorsThreetable = readtable('Fama-French 3 Factor Proxies.xlsx','Range', 'B3:G436');
FactorsThreedata=table2array(FactorsThreetable)
FactorsThreeRF= 0.0036
FactorsThree= FactorsThreedata(4:end, 2:2:end);
FactorsThree=FactorsThree

FactorsThree = rmmissing(FactorsThree); 

[ThreemeanGMVP,ThreestddevGMVP,ThreemeanTP,ThreestdDevTP,Threevarpor_r,Threeretport_r,ThreeretMeansTarget,EFretThree,EFvarThree]=Factors (FactorsThree,0.008,-0.001)

figure(2);
plot(sqrt(Threevarpor_r),Threeretport_r,'b-');
    hold on;
    plot(sqrt(EFvarThree),EFretThree,'g-')
    cal_xThree = [0,ThreestdDevTP];
    cal_yThree = [0,ThreemeanTP];
    stdDevExtendedThree = ThreestdDevTP * 2; % Example: extend the line 50% further out
    meanTPExtendedThree = (ThreemeanTP / ThreestdDevTP) * stdDevExtendedThree; % Use the slope to calculate the new y-value
    plot([0, stdDevExtendedThree], [0, meanTPExtendedThree], 'r--');% Extended CAL with dashed line for differentiation
    plot([ThreestdDevTP, stdDevExtendedThree], [ThreemeanTP, meanTPExtendedThree], 'k-');
    plot(ThreestdDevTP,ThreemeanTP, 'g.', 'MarkerSize', 20)
    plot(ThreestddevGMVP,ThreemeanGMVP,'m.', 'MarkerSize', 20)

    xlabel('Standard Deviation (Risk)');
    ylabel('Excess Returns');
    title('FF 3-Factor Proxy Efficient Frontier with Tangent Portfolio and CML');
    legend('MVF','3-Factor Proxy Efficient Frontier','CML','Portfolio Efficient Frontier');
    text(ThreestdDevTP,ThreemeanTP , 'Tangency Portfolio', 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    text(ThreestddevGMVP,ThreemeanGMVP,'GMVP','VerticalAlignment','bottom', 'HorizontalAlignment','left');

