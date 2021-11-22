function tform = registerRigid(moving, fixed, maxIterations, ...
    w, tolerance, verbose)

if(isSimMode())
    printer = vision.internal.MessagePrinter.configure(verbose);
end
X = fixed;
Y = moving;

D = size(X, 2);
M = size(Y, 1);
N = size(X, 1);

% normalize to zero mean
xmean = mean(X);
ymean = mean(Y);

if(isSimMode())
    X = X - xmean;
    Y = Y - ymean;
else
    for i = 1 : D
        X(:,i) = X(:,i) - xmean(i);
        Y(:,i) = Y(:,i) - ymean(i);
    end
end
% Initialize rotation matrix to identity, translation vector to 0 and
% scaling parameter to 1, i.e. R = I, t = 0
R = cast(eye(D), 'like', Y);
t = cast([0 0 0]', 'like', Y);

sigma2    = cast(0, 'like', Y);

for col = 1 : D
    if(isSimMode())
        sigma2 = sigma2 + sum(sum((X(:,col) - Y(:,col)').^2 ));
    else
        for j = 1 : N
            sigma2 = sigma2 + sum((X(j, col) - Y(:,col)).^2);
        end
    end
end
sigma2 = sigma2 / (D*N*M);

nIter  = 0;
negativeLogLikelihood      = cast(0, 'like', X);

transformedY      = Y ;
% EM optimization, repeat until convergence
while (nIter < maxIterations)
    negativeLogLikelihoodPrev = negativeLogLikelihood;
    
    % E-step: Compute P
    [ P1, Pt1, Px, negativeLogLikelihood ] = computeEStep(X, transformedY, sigma2, w);
    
    ntol            = abs((negativeLogLikelihood-negativeLogLikelihoodPrev)/negativeLogLikelihood);
    
    % M-step : Solve R, t, sigma2
    Np  = sum(P1);
    mux = X'*Pt1 / Np;
    muy = Y'*P1 / Np;
    
    A   = Px'*Y - Np*(mux*muy');
    rcondVal = rcond(A);
    if(isnan(rcondVal) || rcondVal<eps)
        if(isSimMode())
            printer.printMessage('vision:pointcloud:cpdStopCondIllMatrix');
        end
        break;
    end
    [U,~,V]    = svd(A);
    C          = eye(size(U,2));
    C(end,end) = sign(det(U*V'));
    R = U*C*V';
    
    t = mux - R*muy;
    if(isSimMode())
        X1  = X - mux';
        Y1  = Y - muy';
        sigma2 = abs(( sum(sum((X1.^2).*Pt1 )) - (trace(A'*R)^2)/sum(sum((Y1.^2).*P1)) )/(Np*D));
        transformedY      = Y*R' + t';
    else
        sum1 = cast(0, 'like', X);
        sum2 = cast(0, 'like', X);
        for c = 1 : D
            sum1 = sum1 + sum(((X(:, c)-mux(c)).^2).*Pt1);
            sum2 = sum2 + sum(((Y(:, c)-muy(c)).^2).*P1);
        end
        sigma2 = abs(( sum1 - (trace(A'*R)^2)/sum2 )/(Np*D));
        transformedY      = Y*R';
        
        for col = 1 : D
            transformedY(:, col) = transformedY(:, col) + t(col);
        end
    end
    nIter  = nIter + 1;
    if(isSimMode())
        printer.linebreak;
        printer.print('--------------------------------------------\n');
        printer.printMessage('vision:pointcloud:cpdIteration','Rigid', nIter);
        printer.printMessage('vision:pointcloud:cpdCurrentCovariance', mat2str(sigma2));
        printer.printMessage('vision:pointcloud:cpdCurrentVar', mat2str(t, 6), mat2str(R', 6));
        printer.printMessage('vision:pointcloud:cpdCurrentFcn', mat2str(ntol));
    end
    if(ntol <= tolerance)
        break;
    end
    
end
if(isSimMode())
    if(nIter>0)
        printer.print('--------------------------------------------\n');
        printer.printMessage('vision:pointcloud:cpdTotalIterations', nIter);
    end
end
t = t+xmean' - (R*ymean');


tform = rigid3d(R', t');
end