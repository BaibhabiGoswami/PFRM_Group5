% Efficient Portfolio using quadprog function

clear
clc

tic
%% Question Four
data = readtable('17_Industry_Portfolios.CSV')
% Slicing for value-weighted monthly returns
data=data(1:1170,:)
data.Var1=string(data.Var1)
data.Var1 = datetime(data.Var1, 'InputFormat', 'yyyyMM');
startIndex = find(data.Var1 >= datetime(1964, 1, 1));
pdata=data(startIndex:end, :);
matrixData = table2array(pdata(:,2:end));
matrixData = matrixData / 100;
FactorsFivetable = readtable('F-F_Research_Data_5_Factors_2x3.csv')
FactorsFivedata=table2array(FactorsFivetable)/100
FactorsFiveRF=FactorsFivedata(3:722,7)
FactorsFiveRF=mean(FactorsFiveRF)


[meanGMVP,stddevGMVP,meanTP,stdDevTP,varpor_r,retpor_r,retMeansTarget,EFret,EFvar]=RiskyRiskFree (matrixData,FactorsFiveRF,0.03,-0.015)
figure(1);
plot(sqrt(varpor_r),retpor_r,'b-');
    hold on;
    plot(sqrt(EFvar),EFret,'g-')
    cal_x = [0,stdDevTP];
    cal_y = [0,meanTP-FactorsFiveRF];
    stdDevExtended = stdDevTP * 3; % Example: extend the line 50% further out
    meanTPExtended = ((meanTP-FactorsFiveRF) / stdDevTP) * stdDevExtended; % Use the slope to calculate the new y-value
    plot([0, stdDevExtended], [0, meanTPExtended], 'r--');% Extended CAL with dashed line for differentiation
    plot([stdDevTP, stdDevExtended], [meanTP-FactorsFiveRF, meanTPExtended], 'k-');
    plot(stdDevTP,meanTP-FactorsFiveRF, 'g.', 'MarkerSize', 20)

    xlabel('Standard Deviation (Risk)');
    ylabel('Excess Return');
    title('Risky Assets with Risk-Free Asset');
    legend('MVF','Capital Market Line','Risky-Riskfree Assets Efficient Frontier','Risky-Assets Efficient Frontier');
    text(stdDevTP, meanTP-FactorsFiveRF, 'Tangency Portfolio', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');

%% Question Five
% 5 Factors

FactorsFivetable = readtable('F-F_Research_Data_5_Factors_2x3.csv')
FactorsFivedata=table2array(FactorsFivetable)/100
FactorsFiveRF=FactorsFivedata(3:722,7)
FactorsFiveRF=mean(FactorsFiveRF)
FactorsFive=FactorsFivedata(3:722,2:6)
FactorsFive=FactorsFive + FactorsFiveRF
%mean1=mean(FactorsFive)%%% mean factors even underperform the risk-free?
%% Five Factors
up=max(mean(FactorsFive))
[FivemeanGMVP,FivestddevGMVP,FivemeanTP,FivestdDevTP,Fivevarpor_r,Fiveretport_r,FiveretMeansTarget,EFretFive,EFvarFive]= Factors (FactorsFive,0.014,-0.001)
figure(2);
plot(sqrt(Fivevarpor_r),FiveretMeansTarget,'b-');
    hold on;
    plot(sqrt(EFvarFive),EFretFive,'g-')
    cal_xFive = [0,FivestdDevTP];
    cal_yFive = [0,FivemeanTP];
    stdDevExtendedFive = FivestdDevTP * 3; % Example: extend the line 50% further out
    meanTPExtendedFive = (FivemeanTP / FivestdDevTP) * stdDevExtendedFive; % Use the slope to calculate the new y-value
    plot([0, stdDevExtendedFive], [0, meanTPExtendedFive], 'r--');% Extended CAL with dashed line for differentiation
    plot([FivestdDevTP, stdDevExtendedFive], [FivemeanTP, meanTPExtendedFive], 'k-');
    plot(FivestdDevTP,FivemeanTP, 'g.', 'MarkerSize', 20)

    xlabel('Standard Deviation (Risk)');
    ylabel('Excess Return');
    title('F-F Five Factors');
    legend('MVF','Risky-Assets Efficient Frontier','Capital Market Line','5 Factors Risky-Riskfree Assets Efficient Frontier');
    text(FivestdDevTP,FivemeanTP , 'Tangency Portfolio', 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    
%% 3 Factors
FactorsThreetable = readtable('F-F_Research_Data_Factors.CSV')
FactorsThreedata=table2array(FactorsFivetable)/100
FactorsThreeRF=FactorsThreedata(3:722,5)
FactorsThreeRF=mean(FactorsThreeRF)
FactorsThree=FactorsThreedata(3:722,2:4)+FactorsThreeRF

[ThreemeanGMVP,ThreestddevGMVP,ThreemeanTP,ThreestdDevTP,Threevarpor_r,Threeretport_r,ThreeretMeansTarget,EFretThree,EFvarThree]=Factors (FactorsThree,0.014,-0.001)

figure(3);
plot(sqrt(Threevarpor_r),ThreeretMeansTarget,'b-');
    hold on;
    plot(sqrt(EFvarThree),EFretThree,'g-')
    cal_xThree = [0,ThreestdDevTP];
    cal_yThree = [0,ThreemeanTP];
    stdDevExtendedThree = ThreestdDevTP * 2; % Example: extend the line 50% further out
    meanTPExtendedThree = (ThreemeanTP / ThreestdDevTP) * stdDevExtendedThree; % Use the slope to calculate the new y-value
    plot([0, stdDevExtendedThree], [0, meanTPExtendedThree], 'r--');% Extended CAL with dashed line for differentiation
    plot([ThreestdDevTP, stdDevExtendedThree], [ThreemeanTP, meanTPExtendedThree], 'k-');
    plot(ThreestdDevTP,ThreemeanTP, 'g.', 'MarkerSize', 20)

    xlabel('Standard Deviation (Risk)');
    ylabel('Excess Return');
    title('F-F Three Factors');
    legend('MVF','Risky-Assets Efficient Frontier','Capital Market Line','3 Factors Risky-Riskfree Assets Efficient Frontier');
    text(FivestdDevTP,FivemeanTP , 'Tangency Portfolio', 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
    %% Comparing three portfolios
    figure(4);
    p1_1=plot(sqrt(varpor_r),retpor_r,'b-');
    hold on;
    p1_2=plot(sqrt(EFvar),EFret,'g-')
    p1_3=plot([0, stdDevExtended], [0, meanTPExtended], 'r--');
    p1_4=plot([stdDevTP, stdDevExtended], [meanTP-FactorsFiveRF, meanTPExtended], 'k-');
    p1_5=plot(stdDevTP,meanTP-FactorsFiveRF, 'g.', 'MarkerSize', 20);
    p1_6=plot(stddevGMVP,meanGMVP,'r.','MarkerSize', 20);

    p2_1=plot(sqrt(Fivevarpor_r),FiveretMeansTarget,'b-');
    p2_2=plot(sqrt(EFvarFive),EFretFive,'g-')
    p2_3=plot([0, stdDevExtendedFive], [0, meanTPExtendedFive], 'r--');% Extended CAL with dashed line for differentiation
    p2_4=plot([FivestdDevTP, stdDevExtendedFive], [FivemeanTP, meanTPExtendedFive], 'k-');
    p2_5=plot(FivestdDevTP,FivemeanTP, 'g.', 'MarkerSize', 20)
    p2_6=plot(FivestddevGMVP,FivemeanGMVP,'r.','MarkerSize', 20);

    p3_1=plot(sqrt(Threevarpor_r),ThreeretMeansTarget,'b-');
    p3_2=plot(sqrt(EFvarThree),EFretThree,'g-')
    p3_3=plot([0, stdDevExtendedThree], [0, meanTPExtendedThree], 'r--');% Extended CAL with dashed line for differentiation
    p3_4=plot([ThreestdDevTP, stdDevExtendedThree], [ThreemeanTP, meanTPExtendedThree], 'k-');
    p3_5=plot(ThreestdDevTP,ThreemeanTP, 'g.', 'MarkerSize', 20)
    p3_6=plot(ThreestddevGMVP,ThreemeanGMVP,'r.','MarkerSize', 20);

    xlabel('Standard Deviation (Risk)');
    ylabel('Excess Return');
    title('Comparation For Industries Returns and Factors');
    legend([p1_1,p1_2,p1_3,p1_4,p2_2,p2_4,p3_2,p3_4,p1_5,p1_6], {'MVF','Risky-assets Efficient Frontier','CML','Portfolio Efficient Frontier','5 Factors Risky Efficient Frontier','5 Factors Efficient Frontier','3 Factors Risky Efficient Frontier','3 Factors Efficient Frontier','TP','GMVP'});
    
