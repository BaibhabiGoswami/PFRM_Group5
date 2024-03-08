clear
clc

% Risk-free Rate
Factorsdata = readtable('F-F_Research_Data_Factors 3.CSV','ReadVariableNames',true);
Factorsdata=Factorsdata(1:1171,:)
% Select Risk-free Rate and divide the perid to In-sample and Out-of-smaple
%19630101-20231201
Factorsdata.Properties.VariableNames{1} = 'Date';
Factorsdata.Date=string(Factorsdata.Date)
Factorsdata.Date = datetime(Factorsdata.Date, 'InputFormat', 'yyyyMM');
InsampleRfStart = find(Factorsdata.Date >= datetime(1963, 1, 1),1,'first');
InsampleRfEnd = find(Factorsdata.Date < datetime(1993, 1, 1),1,'last');
InsampleRf=Factorsdata(InsampleRfStart:InsampleRfEnd, 5);
OutsampleRfStart=find(Factorsdata.Date >= datetime(1993, 1, 1),1,"first");
OutsampleRfEnd=find(Factorsdata.Date <= datetime(2022, 12, 1),1,"last");
OutsampleRf=Factorsdata(OutsampleRfStart:OutsampleRfEnd, 5);
InsampleRf=table2array(InsampleRf)/100
OutsampleRf=table2array(OutsampleRf)/100

data = readtable('17_Industry_Portfolios.CSV')
% Slicing for value-weighted monthly returns
% Normal return for 17 industry
data=data(1:1170,:)
data.Var1=string(data.Var1)
data.Var1 = datetime(data.Var1, 'InputFormat', 'yyyyMM');
InsampleStart = find(data.Var1 >= datetime(1963, 1, 1));
InsampleEnd = find(data.Var1 >= datetime(1992, 12, 1));
Insample=data(InsampleStart:InsampleEnd, :);
OutsampleStart=find(data.Var1 >= datetime(1993, 1, 1));
OutsampleEnd=find(data.Var1 <= datetime(2022, 12, 1),1,"last");
Outsample=data(OutsampleStart:OutsampleEnd, :);
IS_returns = table2array(Insample(:,2:end))/100;
OS_returns= table2array(Outsample(:,2:end))/100;

IS_retMeans = mean( IS_returns ) ;
IS_mvarcov  = cov( IS_returns )  ; 

OS_retMeans = mean( OS_returns ) ;
OS_mvarcov  = cov( OS_returns )  ; 

% Excess Return
IS_xreturns = IS_returns - InsampleRf;
IS_xretMeans = mean( IS_xreturns ) ; 
IS_xmvarcov  = cov( IS_xreturns) ;

OS_xreturns = OS_returns - OutsampleRf ;
OS_xretMeans = mean( OS_xreturns ) ; 
OS_xmvarcov  = cov( OS_xreturns) ;

% 5 factors data:
FactorsFivetable = readtable('FF Factor 5 new series.csv')

%%
FactorsFivetable.Properties.VariableNames{2} = 'Date';
FactorsFivetable.Date=string(FactorsFivetable.Date)
FactorsFivetable.Date = datetime(FactorsFivetable.Date, 'InputFormat', 'yyyyMM');
InsampleRfStart = find(FactorsFivetable.Date >= datetime(1963, 1, 1),1,'first');
InsampleRfEnd = find(FactorsFivetable.Date <= datetime(1992, 12, 1),1,'last');
InsampleFactors=FactorsFivetable(InsampleRfStart:InsampleRfEnd, 3:7);
OutsampleRfStart=find(FactorsFivetable.Date >= datetime(1993, 1, 1),1,"first");
OutsampleRfEnd=find(FactorsFivetable.Date <= datetime(2022, 12, 1),1,"last");
OutsampleFactors=FactorsFivetable(OutsampleRfStart:OutsampleRfEnd,3:7 );
InsampleFactors=table2array(InsampleFactors)/100
OutsampleFactors=table2array(OutsampleFactors)/100

retMeans_ISFactors = mean(InsampleFactors) ;
mvarcov_ISFactors  = cov( InsampleFactors )  ; 
retMeans_OSFactors = mean( OutsampleFactors ) ;
mvarcov_OSFactors  = cov( OutsampleFactors )  ; 

%% 17 Industries portfolios in sample and out of sample
% the mean TP using this function corresponds to the normal return;
% for excess return, risk free should be substracted from normal return;
% But for BS, the TP return corresponds to the excess return
% in sample:
[meanGMVP_IS1701,stddevGMVP_IS1701,meanTP_IS1701,stdDevTP_IS1701,varpor_r_IS1701,retport_r_IS1701,retMeansTarget_IS1701,filteredRet_IS1701,filteredVar_IS1701] = RiskyRiskFree (IS_returns,InsampleRf,0.005,0.022)
[BSretMeans_IS17,BSvcov_IS17]=BS(IS_xreturns,IS_xmvarcov,IS_xretMeans);
[meanGMVP_IS1702,stddevGMVP_IS1702,meanTP_IS1702,stdDevTP_IS1702,varpor_r_IS1702,retport_r_IS1702,retMeansTarget_IS1702,filteredRet_IS1702,filteredVar_IS1702] = BSEF (BSretMeans_IS17,IS_mvarcov,0.001,0.010)
[meanGMVP_IS1703,stddevGMVP_IS1703,meanTP_IS1703,stdDevTP_IS1703,varpor_r_IS1703,retport_r_IS1703,retMeansTarget_IS1703,filteredRet_IS1703,filteredVar_IS1703] = BSEF (BSretMeans_IS17,BSvcov_IS17,0.001,0.008)

% out of sample:
[meanGMVP_OS1701,stddevGMVP_OS1701,meanTP_OS1701,stdDevTP_OS1701,varpor_r_OS1701,retport_r_OS1701,retMeansTarget_OS1701,filteredRet_OS1701,filteredVar_OS1701] = RiskyRiskFree (OS_returns,OutsampleRf,0.005,0.016)
[BSretMeans_OS17,BSvcov_OS17]=BS(OS_xreturns,OS_xmvarcov,OS_xretMeans);
[meanGMVP_OS1702,stddevGMVP_OS1702,meanTP_OS1702,stdDevTP_OS1702,varpor_r_OS1702,retport_r_OS1702,retMeansTarget_OS1702,filteredRet_OS1702,filteredVar_OS1702] = BSEF (BSretMeans_OS17,OS_mvarcov,0.001,0.012)
[meanGMVP_OS1703,stddevGMVP_OS1703,meanTP_OS1703,stdDevTP_OS1703,varpor_r_OS1703,retport_r_OS1703,retMeansTarget_OS1703,filteredRet_OS1703,filteredVar_OS1703] = BSEF (BSretMeans_OS17,BSvcov_OS17,0.0013,0.010)
%%
% 5 Factors:
% in sample:
[meanGMVP_ISF01,stddevGMVP_ISF01,meanTP_ISF01,stdDevTP_ISF01,varpor_r_ISF01,retport_r_ISF01,retMeansTarget_ISF01,filteredRet_ISF01,filteredVar_ISF01] = Factors (InsampleFactors,0.0025,0.0045)
[BSretMeans_ISF,BSvcov_ISF]=BS(InsampleFactors,mvarcov_ISFactors,retMeans_ISFactors)
[meanGMVP_ISF02,stddevGMVP_ISF02,meanTP_ISF02,stdDevTP_ISF02,varpor_r_ISF02,retport_r_ISF02,retMeansTarget_ISF02,filteredRet_ISF02,filteredVar_ISF02] = BSEF (BSretMeans_ISF,mvarcov_ISFactors,0.0025,0.0035)
[meanGMVP_ISF03,stddevGMVP_ISF03,meanTP_ISF03,stdDevTP_ISF03,varpor_r_ISF03,retport_r_ISF03,retMeansTarget_ISF03,filteredRet_ISF03,filteredVar_ISF03] = BSEF (BSretMeans_ISF,BSvcov_ISF,0.0025,0.0033)

% out of sample:
[meanGMVP_OSF01,stddevGMVP_OSF01,meanTP_OSF01,stdDevTP_OSF01,varpor_r_OSF01,retport_r_OSF01,retMeansTarget_OSF01,filteredRet_OSF01,filteredVar_OSF01] = Factors (OutsampleFactors,0.0025,0.005)
[BSretMeans_OSF,BSvcov_OSF]=BS(OutsampleFactors,mvarcov_OSFactors,retMeans_OSFactors)
[meanGMVP_OSF02,stddevGMVP_OSF02,meanTP_OSF02,stdDevTP_OSF02,varpor_r_OSF02,retport_r_OSF02,retMeansTarget_OSF02,filteredRet_OSF02,filteredVar_OSF02] = BSEF (BSretMeans_OSF,mvarcov_OSFactors,0.0025,0.0045)
[meanGMVP_OSF03,stddevGMVP_OSF03,meanTP_OSF03,stdDevTP_OSF03,varpor_r_OSF03,retport_r_OSF03,retMeansTarget_OSF03,filteredRet_OSF03,filteredVar_OSF03] = BSEF (BSretMeans_OSF,BSvcov_OSF,0.0025,0.0043)

%% Industries Sharpe ratios of TP and GMVP (Times-series of out-of-sample returns)
% In-sample GMVP weights:
N=17
IS_17weightsGMVP = (IS_mvarcov^(-1)*ones(N,1));
IS_17weightsGMVP = IS_17weightsGMVP./sum(IS_17weightsGMVP);
% Times-series of out-of-sample excess returns on the GMVP
Out_17returnsGMVP = OS_xreturns*IS_17weightsGMVP ;  % notice I use excess returns to calculate excess returns on GMVP
Out_17rbarGMVP= IS_17weightsGMVP' * OS_xretMeans'
Out_17varGMVP=IS_17weightsGMVP'*OS_xmvarcov*IS_17weightsGMVP
Out_17sigGMVP=sqrt(Out_17varGMVP)

% Calculate the weights and time-series of excess returns on the Tangency portfolio (TP)
IS_17weightsTP=IS_xmvarcov^(-1) * IS_xretMeans';
IS_17weightsTP = IS_17weightsTP./sum(IS_17weightsTP);

% Times-series of out-of-sample excess returns on the TP
Out_17returnsTP = OS_xreturns*IS_17weightsTP ; % notice I use excess returns to calculate excess returns TP strategy
Out_17rbarTP= IS_17weightsTP' * OS_xretMeans'
Out_17varTP=IS_17weightsTP'*OS_xmvarcov*IS_17weightsTP
Out_17sigTP=sqrt(Out_17varTP)

% Industry Sharpe ratios of TP and GMVP (Times-series of out-of-sample returns)
Out_17SRTP  = mean(Out_17returnsTP)/std(Out_17returnsTP);
Out_17SRGMVP = mean(Out_17returnsGMVP)/std(Out_17returnsTP);

%%% Sharpe ratios of TP and GMVP in smaple--17 industry
IS_17SRTP = (meanTP_IS1701)/stdDevTP_IS1701
IS_17SRGMVP = (meanGMVP_IS1701)/stddevGMVP_IS1701

% GMVP and TP In sample returns
IS_17returnsGMVP=IS_xreturns*IS_17weightsGMVP;
IS_17rbarGMVP= IS_17weightsGMVP' * IS_xretMeans';
IS_17varGMVP=IS_17weightsGMVP'*IS_xmvarcov*IS_17weightsGMVP;
IS_17sigGMVP=sqrt(IS_17varGMVP);

IS_17returnsTP=IS_xreturns*IS_17weightsTP;
IS_17rbarTP= IS_17weightsTP' * IS_xretMeans';
IS_17varTP=IS_17weightsTP'*IS_xmvarcov*IS_17weightsTP;
IS_17sigTP=sqrt(IS_17varTP)

fprintf('Out-of-sample Sharpe Ratio for TP: %f\n', Out_17SRTP);
fprintf('Out-of-sample Sharpe Ratio for GMVP: %f\n', Out_17SRGMVP);
fprintf('Out-sample average return for GMVP: %f\n', Out_17rbarGMVP);
fprintf('Out-sample standard deviation for GMVP: %f\n', Out_17sigGMVP);
fprintf('Out-sample average return for TP: %f\n', Out_17rbarTP);
fprintf('Out-sample standard deviation for TP: %f\n', Out_17sigTP);

fprintf('In-sample Sharpe Ratio for TP: %f\n', IS_17SRTP);
fprintf('In-sample Sharpe Ratio for GMVP: %f\n', IS_17SRGMVP);
fprintf('In-sample average return for GMVP: %f\n', IS_17rbarGMVP);
fprintf('In-sample standard deviation for GMVP: %f\n', IS_17sigGMVP);
fprintf('In-sample average return for TP: %f\n', IS_17rbarTP);
fprintf('In-sample standard deviation for TP: %f\n', IS_17sigTP);

%% Constructed Out-of-sample portfolio using BS TP weights

% Calculate the weights and time-series of excess returns on the Tangency
% portfolio (TP) in BS
IS_BS17weightsTP=BSvcov_IS17 * BSretMeans_IS17';
IS_BS17weightsTP = IS_BS17weightsTP./sum(IS_BS17weightsTP);

% Times-series of out-of-sample excess returns on the TP
Out_BS17returnsTP = OS_xreturns*IS_BS17weightsTP ; % notice I use excess returns to calculate excess returns TP strategy
Out_BS17rbarTP= IS_BS17weightsTP' * OS_xretMeans'
Out_BS17varTP=IS_BS17weightsTP'*OS_xmvarcov*IS_BS17weightsTP
Out_BS17sigTP=sqrt(Out_BS17varTP)

% Industry Sharpe ratios of TP and GMVP (Times-series of out-of-sample returns)
Out_BS17SRTP  = mean(Out_BS17returnsTP)/std(Out_BS17returnsTP);
%Out_17SRGMVP = mean(Out_17returnsGMVP)/std(Out_17returnsTP);

%%% Sharpe ratios of TP and GMVP in smaple--17 industry
IS_BS17SRTP = (meanTP_IS1703)/stdDevTP_IS1703

IS_BS17returnsTP=IS_xreturns*IS_BS17weightsTP;
IS_BS17rbarTP= IS_BS17weightsTP' * BSretMeans_IS17';
IS_BS17varTP=IS_BS17weightsTP'* BSvcov_IS17 *IS_BS17weightsTP;
IS_BS17sigTP=sqrt(IS_BS17varTP)

fprintf('Out-of-sample Sharpe Ratio for TP: %f\n', Out_BS17SRTP);
fprintf('Out-sample average return for TP: %f\n', Out_BS17rbarTP);
fprintf('Out-sample standard deviation for TP: %f\n', Out_BS17sigTP);

fprintf('In-sample Sharpe Ratio for TP: %f\n', IS_BS17SRTP);
fprintf('In-sample average return for TP: %f\n', IS_BS17rbarTP);
fprintf('In-sample standard deviation for TP: %f\n', IS_BS17sigTP);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the weights and time-series of excess returns on the global minimum variance portfolio
N=5
% For fctors
IS_FweightsGMVP = (mvarcov_ISFactors^(-1)*ones(N,1));
IS_FweightsGMVP = IS_FweightsGMVP./sum(IS_FweightsGMVP);
% Times-series of out-of-sample excess returns on the GMVP
Out_FreturnsGMVP = OutsampleFactors*IS_FweightsGMVP ;  % notice I use excess returns to calculate excess returns on GMVP
Out_FrbarGMVP= IS_FweightsGMVP' * retMeans_OSFactors'
Out_FvarGMVP=IS_FweightsGMVP'*mvarcov_OSFactors*IS_FweightsGMVP
Out_FsigGMVP=sqrt(Out_FvarGMVP)

% Calculate the weights and time-series of excess returns on the Tangency portfolio (TP)
IS_FweightsTP=mvarcov_ISFactors^(-1) * retMeans_ISFactors';
IS_FweightsTP = IS_FweightsTP./sum(IS_FweightsTP);

% Times-series of out-of-sample excess returns on the TP
Out_FreturnsTP = OutsampleFactors*IS_FweightsTP ; % notice I use excess returns to calculate excess returns TP strategy
Out_FrbarTP= IS_FweightsTP' * retMeans_OSFactors'
Out_FvarTP=IS_FweightsTP'*mvarcov_OSFactors*IS_FweightsTP
Out_FsigTP=sqrt(Out_FvarTP)

% Factors Sharpe ratios of TP and GMVP (Times-series of out-of-sample returns)
Out_FSRTP  = mean(Out_FreturnsTP )/std(Out_FreturnsTP );
Out_FSRGMVP = mean(Out_FreturnsGMVP)/std(Out_FreturnsGMVP);

%%% Sharpe ratios of TP and GMVP in smaple--17 industry
IS_FSRTP = (meanTP_ISF01)/stdDevTP_ISF01
IS_FSRGMVP = (meanGMVP_ISF01)/stddevGMVP_ISF01

% GMVP and TP In sample returns
IS_FreturnsGMVP=InsampleFactors*IS_FweightsGMVP;
IS_FrbarGMVP= IS_FweightsGMVP' * retMeans_ISFactors';
IS_FvarGMVP=IS_FweightsGMVP'*mvarcov_ISFactors*IS_FweightsGMVP;
IS_FsigGMVP=sqrt(IS_FvarGMVP);

IS_FreturnsTP=InsampleFactors*IS_FweightsTP;
IS_FrbarTP= IS_FweightsTP' * retMeans_ISFactors';
IS_FvarTP=IS_FweightsTP'*mvarcov_ISFactors*IS_FweightsTP;
IS_FsigTP=sqrt(IS_FvarTP);

% Print results for GMVP
fprintf('--- GMVP ---\n');
fprintf('In-sample Mean Return: %f\n', IS_FrbarGMVP);
fprintf('In-sample Standard Deviation: %f\n', IS_FsigGMVP);
fprintf('In-sample Sharpe Ratio: %f\n', IS_FSRGMVP);
fprintf('Out-of-sample Mean Return: %f\n', Out_FrbarGMVP);
fprintf('Out-of-sample Standard Deviation: %f\n', Out_FsigGMVP);
fprintf('Out-of-sample Sharpe Ratio: %f\n', Out_FSRGMVP);

% Print results for TP
fprintf('--- TP ---\n');
fprintf('In-sample Mean Return: %f\n', IS_FrbarTP);
fprintf('In-sample Standard Deviation: %f\n', IS_FsigTP);
fprintf('In-sample Sharpe Ratio: %f\n', IS_FSRTP);
fprintf('Out-of-sample Mean Return: %f\n', Out_FrbarTP);
fprintf('Out-of-sample Standard Deviation: %f\n', Out_FsigTP);
fprintf('Out-of-sample Sharpe Ratio: %f\n', Out_FSRTP);


%% Test of differences in Sharpe ratio
startTime=tic;  % Start timer
IndustryTP = [IS_17returnsTP, Out_17returnsTP] ;
[rejected, pval, tstat] = robustSharpe(IndustryTP) ;
elapsedTime =toc;  % End timer and display elapsed time
fprintf('Elapsed Time: %.2f seconds\n', elapsedTime);
%%
IndustryGMVP = [IS_17returnsGMVP, Out_17returnsGMVP];
[rejected_InGMVP, pval_InGMVP, tstat_InGMVP] = robustSharpe(IndustryGMVP) ;
%%
FactorsTP = [IS_FreturnsTP, Out_FreturnsTP] ;
[rejected_FTP, pval_FTP, tstat_FTP] = robustSharpe(FactorsTP) ;
%%
FactorsGMVP = [IS_FreturnsGMVP, Out_FreturnsGMVP] ;
[rejected_FGMVP, pval_FGMVP, tstat_FGMVP] = robustSharpe(FactorsGMVP) ;
%%
BSIndustryTP=[IS_BS17returnsTP,Out_BS17returnsTP];
[rejected_BSIndustryTP, pval_BSIndustryTP, tstat_BSIndustryTP] = robustSharpe(BSIndustryTP) ;
%% JobsonKorkie test
% T1=360,[IS_17returnsTP, Out_17returnsTP]
T1=360;
covm1=cov(IS_17returnsTP, Out_17returnsTP);
covar1=covm1(1,2);
tstat1=(Out_17sigTP*IS_17rbarTP-IS_17sigTP*Out_BS17rbarTP)/(2/(T1*(IS_17sigTP^2*Out_17sigTP^2-IS_17sigTP*Out_17sigTP*covar1)))^0.5
% T1=360,[IS_17returnsGMVP, Out_17returnsGMVP]
covm2=cov(IS_17returnsGMVP, Out_17returnsGMVP);
covar2=covm2(1,2);
tstat2=(Out_17sigGMVP*IS_17rbarGMVP-IS_17sigGMVP*Out_17rbarGMVP)/(2/(T1*(IS_17sigGMVP^2*Out_17sigGMVP^2-IS_17sigGMVP*Out_17sigGMVP*covar2)))^0.5
%T1=360,[IS_FreturnsTP, Out_FreturnsTP]
covm3=cov(IS_FreturnsTP, Out_FreturnsTP);
covar3=covm3(1,2);
tstat3=(Out_FsigTP*IS_FrbarTP-IS_FsigTP*Out_FrbarTP)/(2/(T1*(IS_FsigTP^2*Out_FsigTP^2-IS_FsigTP*Out_FsigTP*covar3)))^0.5
%[IS_FreturnsGMVP, Out_FreturnsGMVP]
covm4=cov(IS_FreturnsGMVP, Out_FreturnsGMVP);
covar4=covm4(1,2);
tstat4=(Out_FsigGMVP*IS_FrbarGMVP-IS_FsigGMVP*Out_FrbarGMVP)/(2/(T1*(IS_FsigGMVP^2*Out_FsigGMVP^2-IS_FsigGMVP*Out_FsigGMVP*covar4)))^0.5
%[IS_BS17returnsTP,Out_BS17returnsTP]
covm5=cov(IS_BS17returnsTP, Out_BS17returnsTP);
covar5=covm5(1,2);
tstat5=(Out_BS17sigTP*IS_BS17rbarTP-IS_BS17sigTP*Out_BS17rbarTP)/(2/(T1*(IS_BS17sigTP^2*Out_BS17sigTP^2-IS_BS17sigTP*Out_BS17sigTP*covar5)))^0.5

%%
fprintf('--- IS_17GMVP ---\n');
% Assume excess_returns is your vector of excess returns data
[IS17mean_ret, IS17std_dev, IS17skew, IS17kurt, IS17jb_stat, IS17jb_pvalue, IS17autocorr_p] = calculateExcessReturnStats(IS_17returnsGMVP);

% Display the results
fprintf('Mean Return: %.4f%%\n', IS17mean_ret);
fprintf('Standard Deviation: %.4f%%\n', IS17std_dev);
fprintf('Skewness: %.4f\n', IS17skew);
fprintf('Kurtosis: %.4f\n', IS17kurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', IS17jb_stat, IS17jb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', IS17autocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', IS17autocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', IS17autocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', IS17autocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', IS17autocorr_p(5));

fprintf('--- IS_17TP ---\n');
% Assume excess_returns is your vector of excess returns data
[IS17TPmean_ret, IS17TPstd_dev, IS17TPskew, IS17TPkurt, IS17TPjb_stat, IS17TPjb_pvalue, IS17TPautocorr_p] = calculateExcessReturnStats(IS_17returnsTP);

% Display the results
fprintf('Mean Return: %.4f%%\n', IS17TPmean_ret);
fprintf('Standard Deviation: %.4f%%\n', IS17TPstd_dev);
fprintf('Skewness: %.4f\n', IS17TPskew);
fprintf('Kurtosis: %.4f\n', IS17TPkurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', IS17TPjb_stat, IS17TPjb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', IS17TPautocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', IS17TPautocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', IS17TPautocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', IS17TPautocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', IS17TPautocorr_p(5));

fprintf('--- IS_FGMVP ---\n');
% Assume excess_returns is your vector of excess returns data
[ISFmean_ret, ISFstd_dev, ISFskew, ISFkurt, ISFjb_stat, ISFjb_pvalue, ISFautocorr_p] = calculateExcessReturnStats(IS_FreturnsGMVP);

% Display the results
fprintf('Mean Return: %.4f%%\n', ISFmean_ret);
fprintf('Standard Deviation: %.4f%%\n', ISFstd_dev);
fprintf('Skewness: %.4f\n', ISFskew);
fprintf('Kurtosis: %.4f\n', ISFkurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', ISFjb_stat, ISFjb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', ISFautocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', ISFautocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', ISFautocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', ISFautocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', ISFautocorr_p(5));

fprintf('--- In_FTP ---\n');
% Assume excess_returns is your vector of excess returns data
[ISFTPmean_ret, ISFTPstd_dev, ISFTPskew, ISFTPkurt, ISFTPjb_stat, ISFTPjb_pvalue, ISFTPautocorr_p] = calculateExcessReturnStats(IS_FreturnsTP);

% Display the results
fprintf('Mean Return: %.4f%%\n', ISFTPmean_ret);
fprintf('Standard Deviation: %.4f%%\n', ISFTPstd_dev);
fprintf('Skewness: %.4f\n', ISFTPskew);
fprintf('Kurtosis: %.4f\n', ISFTPkurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', ISFTPjb_stat, ISFTPjb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', ISFTPautocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', ISFTPautocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', ISFTPautocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', ISFTPautocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', ISFTPautocorr_p(5));

%%
fprintf('--- OS_17GMVP ---\n');
% Assume excess_returns is your vector of excess returns data
[OS17mean_ret, OS17std_dev, OS17skew, OS17kurt, OS17jb_stat, OS17jb_pvalue, OS17autocorr_p] = calculateExcessReturnStats(Out_17returnsGMVP);

% Display the results
fprintf('Mean Return: %.4f%%\n', OS17mean_ret);
fprintf('Standard Deviation: %.4f%%\n', OS17std_dev);
fprintf('Skewness: %.4f\n', OS17skew);
fprintf('Kurtosis: %.4f\n', OS17kurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', OS17jb_stat, OS17jb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', OS17autocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', OS17autocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', OS17autocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', OS17autocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', OS17autocorr_p(5));

fprintf('--- OS_17TP ---\n');
% Assume excess_returns is your vector of excess returns data
[OS17TPmean_ret, OS17TPstd_dev, OS17TPskew, OS17TPkurt, OS17TPjb_stat, OS17TPjb_pvalue, OS17TPautocorr_p] = calculateExcessReturnStats(Out_17returnsTP);

% Display the results
fprintf('Mean Return: %.4f%%\n', OS17TPmean_ret);
fprintf('Standard Deviation: %.4f%%\n', OS17TPstd_dev);
fprintf('Skewness: %.4f\n', OS17TPskew);
fprintf('Kurtosis: %.4f\n', OS17TPkurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', OS17TPjb_stat, OS17TPjb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', OS17TPautocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', OS17TPautocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', OS17TPautocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', OS17TPautocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', OS17TPautocorr_p(5));

fprintf('--- OS_FGMVP ---\n');
% Assume excess_returns is your vector of excess returns data
[OSFmean_ret, OSFstd_dev, OSFskew, OSFkurt, OSFjb_stat, OSFjb_pvalue, OSFautocorr_p] = calculateExcessReturnStats(Out_FreturnsGMVP);

% Display the results
fprintf('Mean Return: %.4f%%\n', OSFmean_ret);
fprintf('Standard Deviation: %.4f%%\n', OSFstd_dev);
fprintf('Skewness: %.4f\n', OSFskew);
fprintf('Kurtosis: %.4f\n', OSFkurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', OSFjb_stat, OSFjb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', OSFautocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', OSFautocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', OSFautocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', OSFautocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', OSFautocorr_p(5));

fprintf('--- Out_FTP ---\n');
% Assume excess_returns is your vector of excess returns data
[OSFTPmean_ret, OSFTPstd_dev, OSFTPskew, OSFTPkurt, OSFTPjb_stat, OSFTPjb_pvalue, OSFTPautocorr_p] = calculateExcessReturnStats(Out_FreturnsTP);
% Display the results
fprintf('Mean Return: %.4f%%\n', OSFTPmean_ret);
fprintf('Standard Deviation: %.4f%%\n', OSFTPstd_dev);
fprintf('Skewness: %.4f\n', OSFTPskew);
fprintf('Kurtosis: %.4f\n', OSFTPkurt);
fprintf('Jarque-Bera Test Statistic: %.4f, p-value: %.4f\n', OSFTPjb_stat, OSFTPjb_pvalue);
fprintf('Autocorrelations at lag 1: %.4f\n', OSFTPautocorr_p(1));
fprintf('Autocorrelations at lag 2: %.4f\n', OSFTPautocorr_p(2));
fprintf('Autocorrelations at lag 3: %.4f\n', OSFTPautocorr_p(3));
fprintf('Autocorrelations at lag 6: %.4f\n', OSFTPautocorr_p(4));
fprintf('Autocorrelations at lag 12: %.4f\n', OSFTPautocorr_p(5));

%%
% plot
% Highlight the EF:
figure(1);
p1=plot(sqrt(varpor_r_IS1701),retport_r_IS1701,'b-', 'LineWidth', 2);
hold on;
p2=plot(sqrt(varpor_r_IS1702),retport_r_IS1702,'k--', 'LineWidth', 2);
p3=plot(sqrt(varpor_r_IS1703),retport_r_IS1703,'r--', 'LineWidth', 2);

p7=plot(sqrt(varpor_r_OS1701),retport_r_OS1701,'g-', 'LineWidth', 2);
p8=plot(sqrt(varpor_r_OS1702),retport_r_OS1702,'c:', 'LineWidth', 2);
p9=plot(sqrt(varpor_r_OS1703),retport_r_OS1703,'m:', 'LineWidth', 2);

scatter(stddevGMVP_IS1701,meanGMVP_IS1701 , 50, 'b', 'filled')
scatter(stdDevTP_IS1701,meanTP_IS1701, 50, 'r', 'filled')
scatter(stddevGMVP_IS1702,meanGMVP_IS1702 , 50, 'b', 'filled')
scatter(stdDevTP_IS1702,meanTP_IS1702, 50, 'r', 'filled')
scatter(stddevGMVP_IS1703,meanGMVP_IS1703 , 50, 'b', 'filled')
scatter(stdDevTP_IS1703,meanTP_IS1703, 50, 'r', 'filled')

scatter(stddevGMVP_OS1701,meanGMVP_OS1701 , 50, 'b', 'filled')
scatter(stdDevTP_OS1701,meanTP_OS1701, 50, 'r', 'filled')
scatter(stddevGMVP_OS1702,meanGMVP_OS1702 , 50, 'b', 'filled')
scatter(stdDevTP_OS1702,meanTP_OS1702, 50, 'r', 'filled')
scatter(stddevGMVP_OS1703,meanGMVP_OS1703 , 50, 'b', 'filled')
scatter(stdDevTP_OS1703,meanTP_OS1703, 50, 'r', 'filled')

scatter(Out_17sigGMVP,Out_17rbarGMVP , 50, 'g', 'filled')
scatter(Out_17sigTP,Out_17rbarTP, 50, 'g', 'filled')
scatter(Out_BS17sigTP,Out_BS17rbarTP, 50, 'm', 'filled')

cal_yOS17 = [0,stdDevTP_OS1701];
cal_xOS17 = [0,meanTP_OS1701];
stdDevExtendedOS17 = stdDevTP_OS1701 * 1.6; % Example: extend the line 50% further out
meanTPExtendedOS17 = (meanTP_OS1701 / stdDevTP_OS1701) * stdDevExtendedOS17; % Use the slope to calculate the new y-value
plot([0, stdDevExtendedOS17], [0, meanTPExtendedOS17], 'r--');% Extended CAL with dashed line for differentiation
p10=plot([stdDevTP_OS1701, stdDevExtendedOS17], [meanTP_OS1701, meanTPExtendedOS17], 'k-');

text(stddevGMVP_IS1701,meanGMVP_IS1701, ' MVP', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_IS1701,meanTP_IS1701,'TP','VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_IS1702,meanTP_IS1702,'TP','VerticalAlignment','bottom', 'HorizontalAlignment','left');
%text(stdDevTP_IS1703,meanTP_IS1703,'TP','VerticalAlignment','top', 'HorizontalAlignment','left');
text(stddevGMVP_OS1701,meanGMVP_OS1701, ' MVP', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_OS1701,meanTP_OS1701,'TP','VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_OS1702,meanTP_OS1702,'TP','VerticalAlignment','bottom', 'HorizontalAlignment','left');
text(Out_17sigGMVP,Out_17rbarGMVP,'Con.GMVP','VerticalAlignment','top', 'HorizontalAlignment','left');
text(Out_17sigTP,Out_BS17rbarTP,'Con.TP','VerticalAlignment','bottom', 'HorizontalAlignment','left');
text(Out_BS17sigTP,Out_BS17rbarTP,'Con.BS TP','VerticalAlignment','top', 'HorizontalAlignment','left')

xlabel('Standard Deviation (Risk)');
ylabel('Excess Return');
title('17 Industry Portfolio In-Sample VS Out-Sample Efficient Frontier');
legend([p1,p2,p3,p7,p8,p9,p10], {'In Sample MVF','In Sample BS Shrunk Mean','In Sample BS Shrunk Mean Variance','Out Sample MVF','Out Sample BS Shrunk Mean','Out Sample BS Shrunk Mean Variance','Portfolio Efficient Frontier'});

%% Plot:
% plot
% Highlight the EF:
figure(2);
p1=plot(sqrt(varpor_r_ISF01),retport_r_ISF01,'b-', 'LineWidth', 2);
hold on;
p2=plot(sqrt(varpor_r_ISF02),retport_r_ISF02,'k--', 'LineWidth', 2);
p3=plot(sqrt(varpor_r_ISF03),retport_r_ISF03,'r--', 'LineWidth', 2);

p7=plot(sqrt(varpor_r_OSF01),retport_r_OSF01,'g-', 'LineWidth', 2);
p8=plot(sqrt(varpor_r_OSF02),retport_r_OSF02,'c:', 'LineWidth', 2);
p9=plot(sqrt(varpor_r_OSF03),retport_r_OSF03,'m:', 'LineWidth', 2);

scatter(stddevGMVP_ISF01,meanGMVP_ISF01 , 50, 'b', 'filled')
scatter(stdDevTP_ISF01,meanTP_ISF01, 50, 'r', 'filled')
scatter(stddevGMVP_ISF02,meanGMVP_ISF02 , 50, 'b', 'filled')
scatter(stdDevTP_ISF02,meanTP_ISF02, 50, 'r', 'filled')
scatter(stddevGMVP_ISF03,meanGMVP_ISF03 , 50, 'b', 'filled')
scatter(stdDevTP_ISF03,meanTP_ISF03, 50, 'r', 'filled')

scatter(stddevGMVP_OSF01,meanGMVP_OSF01 , 50, 'b', 'filled')
scatter(stdDevTP_OSF01,meanTP_OSF01, 50, 'r', 'filled')
scatter(stddevGMVP_OSF02,meanGMVP_OSF02 , 50, 'b', 'filled')
scatter(stdDevTP_OSF02,meanTP_OSF02, 50, 'r', 'filled')
scatter(stddevGMVP_OSF03,meanGMVP_OSF03 , 50, 'b', 'filled')
scatter(stdDevTP_OSF03,meanTP_OSF03, 50, 'r', 'filled')

p11=scatter(Out_FsigGMVP,Out_FrbarGMVP , 50, 'g', 'filled')
p12=scatter(Out_FsigTP,Out_FrbarTP, 50, 'm', 'filled')

cal_OSF = [0,stdDevTP_OSF01];
cal_OSF = [0,meanTP_OSF01];
stdDevExtendedOSF = stdDevTP_OSF01 * 1.2; % Example: extend the line 50% further out
meanTPExtendedOSF = (meanTP_OSF01 / stdDevTP_OSF01) * stdDevExtendedOSF; % Use the slope to calculate the new y-value
plot([0, stdDevExtendedOSF], [0, meanTPExtendedOSF], 'r--');% Extended CAL with dashed line for differentiation
p10=plot([stdDevTP_OSF01, stdDevExtendedOSF], [meanTP_OSF01, meanTPExtendedOSF], 'k-');

text(stddevGMVP_ISF01,meanGMVP_ISF01, ' MVP', 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_ISF01,meanTP_ISF01,'TP','VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_ISF02,meanTP_ISF02,'TP','VerticalAlignment','top', 'HorizontalAlignment','right');
%text(stdDevTP_ISF03,meanTP_ISF03,'BS Returns and BS Var TP','VerticalAlignment','top', 'HorizontalAlignment','left');
text(stddevGMVP_OSF01,meanGMVP_OSF01, ' MVP', 'VerticalAlignment','top', 'HorizontalAlignment','right');
text(stdDevTP_OSF01,meanTP_OSF01,'TP','VerticalAlignment','bottom', 'HorizontalAlignment','right');
text(stdDevTP_OSF02,meanTP_OSF02,'TP','VerticalAlignment','top', 'HorizontalAlignment','left');
text(Out_FsigGMVP,Out_FrbarGMVP,'Constructed GMVP','VerticalAlignment','bottom', 'HorizontalAlignment','left');
text(Out_FsigTP,Out_FrbarTP,'Constructed TP','VerticalAlignment','top', 'HorizontalAlignment','left');

xlabel('Standard Deviation (Risk)');
ylabel('Excess Return');
title('5 Factors Portfolio In-Sample In-Sample VS Out-Sample Efficient Frontiers');
legend([p1,p2,p3,p7,p8,p9,p10], {'In Sample MVF','In Sample BS Shrunk Mean','In Sample BS Shrunk Mean Variance','Out Sample MVF','Out Sample BS Shrunk Mean','Out Sample BS Shrunk Mean Variance','Portfolio Efficient Frontier'});