function [filteredVar,filteredRet,weights] = EfficientFrontier(expectedReturns,covMatrix,up)
    [n,col]=size(expectedReturns)
    low = min(expectedReturns);
    %up = max(expectedReturns)
    step = (up-low)/100
    retMeansTarget = (low:step:up)
    mm=length(retMeansTarget)
    varpor=zeros(1,mm);
    weights = zeros(col,mm);
    c = zeros(col,1);
    Aeq = [expectedReturns;ones(1,col)];
    filteredRet = [];
    filteredStdDev = [];
    filteredVar=[];
    
    for i=1:mm
        beq = [retMeansTarget(i);1];
        vlb = -ones(col,1)*10000
        vub = ones(col,1)*10000;
        x0 = repmat(1/col, col, 1);
        options = []; 
        options = optimset('display', 'on','largescale', 'off');
        x = quadprog(covMatrix,c,[],[],Aeq,beq,vlb,vub,x0,options);
        varpor(i) = x'*covMatrix*x;
        weights(:,i) = x;
    end
    retport=expectedReturns*weights
    stdDevpor=sqrt(varpor)

    % Identifying the minimum standard deviation and filtering

    [minStdDev, minIndex] = min(stdDevpor);
    minStdDevReturn = retport(minIndex);

    % Filter for portfolios with expected returns >= minStdDevReturn
    for i = 1:mm
        if retport(i) >= minStdDevReturn
            filteredRet = [filteredRet, retport(i)];
            filteredStdDev = [filteredStdDev, stdDevpor(i)];
        end
    end
    filteredVar=filteredStdDev.^2
end