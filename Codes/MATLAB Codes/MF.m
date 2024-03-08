 mutualFundData = readtable('MFA.xlsx', 'ReadVariableNames', false); 
 returns = mutualFundData{:, 1};    
 riskFreeRate = 0.0036;  
  
%market index 
marketIndexData = readtable('MARKET1.xlsx', 'ReadVariableNames', false); 
marketReturn = mean(marketIndexData{:,:}, 2);  
meanReturn = mean(returns) * 12; 
meanMarketReturn = mean(marketReturn) * 12;  
stdReturn = std(returns) * sqrt(12); 
   
X = [ones(length(marketReturn), 1), marketReturn - riskFreeRate];  
Y = returns - riskFreeRate;  
b = regress(Y, X);  
alphaSample = b(1);  
betaSample = b(2);  

%Performance Measures
sharpeRatio = (meanReturn - riskFreeRate) / stdReturn;  
treynorsMeasure = (meanReturn - riskFreeRate) / betaSample;  
jensenAlpha = alphaSample; % Jensen's Alpha is the intercept from the regression  
informationRatio = (meanReturn - meanMarketReturn) / std(returns - marketReturn);  
  
% Display  
fprintf('Sharpe Ratio: %f\n', sharpeRatio);  
fprintf('Treynors Measure: %f\n', treynorsMeasure);  
fprintf('Jensens Alpha: %f\n', jensenAlpha);  
fprintf('Information Ratio: %f\n', informationRatio);  
fprintf('Beta: %f\n', betaSample)
fprintf('Expected Annual Return: %f%%\n', meanReturn)
fprintf('Annual Standard Deviation: %f%%\n', stdReturn)