function [mvpReturn,mvpVolatility,mvpSharpeRatio,tpReturn,tpVolatility,tpSharpeRatio]= MVPTP(retpor,varpor,xretMeans,xmvarcov,retMeans,mvarcov,riskFreeRate,legendLabel)
    % Find MVP (minimum volatility)
    [~, mvpIndex] = min(varpor);
    mvpReturn = retpor(mvpIndex);
    mvpVolatility = sqrt(varpor(mvpIndex));
    mvpSharpeRatio = (mvpReturn-riskFreeRate)/mvpVolatility;

    % Find TP:analytical solution for weights of tangency portfolio 
    % weights unscaled
    weightsTP = xmvarcov^(-1)*xretMeans';
    %weights scaled to sum to 1
    weightsTP = weightsTP./sum(weightsTP);
    % Tangency portfolio statistics
    tpReturn= retMeans*weightsTP ;
    varTP  = weightsTP'*mvarcov*weightsTP ;
    tpVolatility  = sqrt(varTP) ;
    tpSharpeRatio=tpReturn/tpVolatility


% Print details for MVP
    fprintf('%s - Minimum Volatility Portfolio (MVP):\n',legendLabel);
    fprintf('Return: %f\n', mvpReturn);
    fprintf('Volatility: %f\n', mvpVolatility);
    fprintf('Sharpe Ratio: %f\n\n', mvpSharpeRatio);

    % Print details for TP
    fprintf('%s - Tangency Portfolio (TP):\n', legendLabel);
    fprintf('Return: %f\n', tpReturn);
    fprintf('Volatility: %f\n', tpVolatility);
    fprintf('Sharpe Ratio: %f\n\n', tpSharpeRatio);

end