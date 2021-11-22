function [ tform, movingReg, rmse ] = pcregistercpd(moving, fixed, varargin)
%PCREGISTERCPD Register two point clouds using CPD algorithm.
%   tform = PCREGISTERCPD(moving, fixed) returns a transformation that
%   registers the moving point cloud with the fixed point cloud using the
%   Coherent Point Drift (CPD) algorithm. moving and fixed are pointCloud
%   objects. Returned tform is a rigid3d object, affine3d object or a
%   displacement field. By default pcregistercpd returns non-rigid
%   transformation.
%
%   Consider downsampling point clouds using pcdownsample before using
%   PCREGISTERCPD to improve accuracy and efficiency of registration.
%
%   [tform, movingReg] = PCREGISTERCPD(moving, fixed) additionally
%   returns the transformed  point cloud, movingReg, that is aligned with
%   the fixed point cloud.
%
%   [..., rmse] = PCREGISTERCPD(...) returns the root mean squared error
%   of the Euclidean distance between the aligned point clouds.
%
%   [...] = PCREGISTERCPD(..., Name, Value) specifies additional
%   name-value pairs described below:
%
%   'Transform'              Specify transform type estimated during
%                            registration as one of: 'Rigid', 'Affine',
%                            'Nonrigid'.
%
%                            Default: 'Nonrigid'
%
%   'OutlierRatio'           A positive scalar less than 1, specifying the
%                            expected ratio of noise and outlier points
%                            present in the point clouds. Higher outlier
%                            ratio reduces the influence of noise and
%                            outliers.
%
%                            Default: 0.1
%
%   'MaxIterations'          A positive integer specifying the maximum
%                            number of iterations.
%
%                            Default: 20
%
%   'Tolerance'              Tolerance criteria for CPD convergence. The
%                            algorithm stops when absolute percentage
%                            change in log likelihood function values
%                            measured between consecutive iterations
%                            reaches or falls below this value. Typical
%                            value is small to increase likelihood of a
%                            better alignment.
%
%                            Default: 1e-5
%
%   'InteractionSigma'       Positive scalar defined as a standard
%                            deviation of a Gaussian filter. This value
%                            controls the amount of interaction between
%                            points. The greater the interaction between
%                            points, the more a coherent motion is observed
%                            in the point cloud, approaching the same
%                            translation for each point. Use smaller values
%                            to produce a more localized deformation in the
%                            output displacement field. Use larger values
%                            to increase the interaction between points.
%                            Typical values range from 1.5 to 3.
%
%                            This parameter is used only with the
%                            'Nonrigid' transform.
%
%                            Default: 2.0
%
%   'SmoothingWeight'        Positive scalar that defines the weight used
%                            to smooth the motion of the point cloud. Use a
%                            larger weight to produce a more coherent
%                            motion in the output displacement field.
%                            Typical values range from 0.1 to 10.
%
%                            This parameter is used only with the
%                            'Nonrigid' transform.
%
%                            Default: 3.0
%
%   'Verbose'                Set true to display progress information.
%
%                            Default: false
%
%   Note
%   ----
%   Both MaxIterations and Tolerance are used as stopping criteria.
%   The algorithm stops when stopping condition of either of the two values
%   is met i.e, when algorithm iterations reaches MaxIterations or when
%   absolute percentage change in log likelihood function reaches or falls
%   below Tolerance value.
%
%   Class Support
%   -------------
%   moving and fixed must be pointCloud objects.  tform is a rigid3d object
%   when used for rigid registration, an affine3d object when used for
%   affine registration and an array of size same as 'moving.Location' for
%   non-rigid registration. movingReg is a pointCloud object.
%
%   Example : Align two point clouds using non-rigid registration.
%   ---------------------------------------------------------------
%   handData = load('hand3d.mat');
%
%   % Downsample point clouds to accelerate the registration process
%   fixedDownsampled = pcdownsample(handData.fixed, 'gridAverage', 0.03);
%   movingDownsampled = pcdownsample(handData.moving, 'gridAverage', 0.03);
%
%   % Visualize input point clouds
%   figure
%   pcshowpair(movingDownsampled, fixedDownsampled, 'MarkerSize', 50)
%   title('Before Registration')
%
%   % Perform non-rigid registration
%   tform = pcregistercpd(movingDownsampled, fixedDownsampled);
%   movingReg = pctransform(movingDownsampled, tform);
%
%   % Visualize results
%   figure
%   pcshowpair(movingReg, fixedDownsampled, 'MarkerSize', 50)
%   title('After Registration')
%
%   See also pointCloud, pctransform, pcdownsample, rigid3d, affine3d, 
%            pcshow, pcshowpair, pcregistericp, pcregisterndt, pcdenoise,
%            pcmerge

% Copyright 2018-2019 The MathWorks, Inc.
%
% References
% ----------
% [1] A. Myronenko, X. Song, "Point set registration: Coherent point drift",
% IEEE Trans. Pattern Anal. Mach. Intell., vol. 32, no. 12, pp. 2262-2275,
% Dec. 2010.

%#codegen

narginchk(2,16);

% Validate arguments
[movingValid, fixedValid, movingValidPointIndices, tformType, paramsStruct] = validateAndParseInputArgs(moving, fixed, varargin{:});

% Invoke appropriate registration method
if(strcmpi('Rigid', tformType))
    tform = registerRigid(movingValid.Location, fixedValid.Location, paramsStruct.MaxIterations,...
        paramsStruct.OutlierRatio, paramsStruct.Tolerance, paramsStruct.Verbose);
elseif(strcmpi('Affine', tformType))
    tform = registerAffine(movingValid.Location, fixedValid.Location, paramsStruct.MaxIterations,...
        paramsStruct.OutlierRatio, paramsStruct.Tolerance, paramsStruct.Verbose);
else
    tform = registerNonrigid(movingValid.Location, fixedValid.Location, paramsStruct.MaxIterations,...
        paramsStruct.OutlierRatio, paramsStruct.Tolerance, paramsStruct.InteractionSigma, ...
        paramsStruct.SmoothingWeight, paramsStruct.Verbose);
    
    % Return displacement fields in same size as moving.Location
    tempTform = nan(size(moving.Location));
    ndims = size(tempTform,3);
    isOrganized = ndims > 1;
    if(isOrganized)
        offset = size(tempTform,1) * size(tempTform,2);
    else
        offset = size(tempTform,1);
    end
    nchannels = size(tform,2);
    for col = 1 : nchannels
        tempTform( movingValidPointIndices + (col-1)*offset ) = tform(:, col);
    end
    tform = cast(tempTform, 'like', moving.Location);
end

if nargout >= 2
    movingReg = pctransform(moving, tform);
end

if nargout >= 3
    % rmse is based on point-to-point correspondence
    [~, dists] = multiQueryKNNSearchImpl(removeInvalidPoints(movingReg), fixedValid.Location, 1);
    rmse = sqrt(mean(dists(:)));
    rmse = cast(rmse, 'like', moving.Location);
end
end

%======================================================================
% function to perform CPD Rigid registration
%======================================================================
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

%======================================================================
% function to perform CPD Affine registration
%======================================================================
function tform = registerAffine(moving, fixed, maxIterations, ...
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

% Initialize affine transformation matrix to identity and translation
% vector to 0, i.e. B = I, t = 0
B = cast(eye(D), 'like', Y);
t = cast([0 0 0]', 'like', Y);

sigma2     = cast(0, 'like', Y);

for col = 1 : D
    if(isSimMode())
        sigma2 = sigma2 + sum(sum( (X(:,col) - Y(:,col)').^2 ));
    else
        for j = 1 : N
            sigma2 = sigma2 + sum((X(j, col) - Y(:, col)).^2);
        end
    end
end
sigma2 = sigma2 / (D*N*M);

nIter  = 0;

negativeLogLikelihood      = cast(0, 'like', X);

transformedY      = Y;

if(~isSimMode())
    Y1 = zeros(size(Y), 'like', Y);
end
% EM optimization, repeat until convergence
while (nIter < maxIterations)
    
    negativeLogLikelihoodPrev  = negativeLogLikelihood;
    
    % E-step: Compute P
    [ P1, Pt1, Px, negativeLogLikelihood ] = computeEStep(X, transformedY, sigma2, w);
    ntol   = abs((negativeLogLikelihood-negativeLogLikelihoodPrev)/negativeLogLikelihood);
    
    % M-step : Solve B, t, sigma2
    Np     = sum(P1);
    mux    = X'*Pt1 / Np;
    muy    = Y'*P1 / Np;
    if(isSimMode())
        X1     = X - mux';
        Y1     = Y - muy';
        
        C     = ((Y1'.*(P1'))*Y1);
    else
        for col = 1 : D
            Y1(:,col) = Y(:,col) - muy(col);
        end
        
        C     = Y1;
        for i = 1 : D
            C(:, i) = C(:, i).*P1;
        end
        C = C'*Y1;
    end
    rcondVal = rcond(C);
    if(isnan(rcondVal) || rcondVal<eps)
        if(isSimMode())
            printer.printMessage('vision:pointcloud:cpdStopCondIllMatrix');
        end
        break;
    end
    A      = Px'*Y - Np*(mux*muy');
    B      = A/C;
    t      = mux - B*muy;
    if(isSimMode())
        sigma2 = abs(( sum(sum( (X1.^2).*Pt1 )) - trace(A*B'))/(Np*D));
        transformedY      = Y*B' + t';
    else
        sum1 = cast(0, 'like', X);
        for c = 1 : D
            sum1 = sum1 + sum(((X(:, c)-mux(c)).^2).*Pt1);
        end
        
        sigma2 = abs(( sum1 - trace(A*B'))/(Np*D));
        
        transformedY      = Y*B';
        for col = 1 : D
            transformedY(:, col) = transformedY(:, col) + t(col);
        end
    end
    
    
    nIter  = nIter + 1;
    if(isSimMode())
        printer.linebreak;
        printer.print('--------------------------------------------\n');
        printer.printMessage('vision:pointcloud:cpdIteration','Affine', nIter);
        printer.printMessage('vision:pointcloud:cpdCurrentCovariance', mat2str(sigma2));
        printer.printMessage('vision:pointcloud:cpdCurrentVar', mat2str(t, 6), mat2str(B', 6));
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
t = t+xmean' - (B*ymean');

Bt = [ B, t(:); 0 0 0 1];
tform = affine3d(Bt');
end

%======================================================================
% function to perform CPD Nonrigid registration
%======================================================================
function tform = registerNonrigid(moving, fixed, maxIterations, ...
    w, tolerance, beta, lambda, verbose)

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

% Initialization
W = zeros(size(Y), 'like', Y);

sigma2     = cast(0, 'like', Y);

for col = 1 : D
    if(isSimMode())
        sigma2 = sigma2 + sum(sum((X(:,col) - Y(:,col)').^2 ));
    else
        for j = 1 : N
            sigma2 = sigma2 + sum((X(j, col) - Y(:, col)).^2);
        end
    end
end
sigma2 = sigma2 / (D*N*M);

nIter  = 0;
negativeLogLikelihood      = cast(0, 'like', X);
G      = zeros(M, M, 'like', Y);
if(~isSimMode())
    A      = zeros(M, M, 'like', Y);
end

for col = 1 : D
    if(isSimMode())
        G = G + (Y(:,col) - Y(:,col)').^2;
    else
        for j = 1 : M
            G(:, j) = G(:, j) + (Y(j, col) - Y(:, col)).^2;
        end
    end
end
G = exp((-1/(2*beta*beta))*G);

transformedY = Y;

if(~isSimMode())
    P1Y = zeros(size(Y), 'like', Y);
end

% EM optimization, repeat until convergence
while (nIter < maxIterations)
    negativeLogLikelihoodPrev = negativeLogLikelihood;
    
    % E-step: Compute P
    [ P1, Pt1, Px, negativeLogLikelihood ] = computeEStep(X, transformedY, sigma2, w);
    
    negativeLogLikelihood      = negativeLogLikelihood+lambda/2*trace(W'*G*W);
    ntol   = abs((negativeLogLikelihood-negativeLogLikelihoodPrev)/negativeLogLikelihood);
    
    % M-step: Solve G, W, sigma2
    if(isSimMode())
        A = (P1.*G + lambda*sigma2*eye(M));
    else
        for i = 1 : M
            A(:, i) = G(:, i).*P1;
            A(i, i) = A(i, i) + lambda*sigma2;
        end
    end
    rcondVal = rcond(A);
    if(isnan(rcondVal) || rcondVal<eps)
        if(isSimMode())
            printer.printMessage('vision:pointcloud:cpdStopCondIllMatrix');
        end
        break;
    end
    if(isSimMode())
        W      = A\(Px - P1.*Y);
    else
        for i = 1 : D
            P1Y(:, i) = Px(:, i) - P1.*Y(:, i);
        end
        W      = A\(P1Y);
    end
    
    Np     = sum(P1);
    transformedY      = Y + G*W;
    if(isSimMode())
        sigma2 = (sum(sum((X.^2).*Pt1)) - 2*trace(Px'*transformedY) + sum(sum((transformedY.^2).*P1)))/(Np*D) ;
    else
        sum1 = cast(0, 'like', X);
        sum3 = cast(0, 'like', X);
        sum2 = cast(0, 'like', X);
        for c = 1 : D
            sum1 = sum1 + sum((X(:, c).^2).*Pt1);
            sum2 = sum2 + sum(Px(:, c).*transformedY(:, c));
            sum3 = sum3 + sum((transformedY(:, c).^2).*P1);
        end
        sigma2 = (sum1 - 2*sum2 + sum3)/(Np*D) ;
    end
    nIter  = nIter + 1;
    if(isSimMode())
        printer.linebreak;
        printer.print('--------------------------------------------\n');
        printer.printMessage('vision:pointcloud:cpdIteration','Nonrigid', nIter);
        printer.printMessage('vision:pointcloud:cpdCurrentCovariance', mat2str(sigma2));
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
if(isSimMode())
    tform = G*W + (xmean-ymean);
else
    tform = G*W;
    for i = 1 : D
        tform(:, i) = tform(:, i) + (xmean(i)-ymean(i));
    end
end
end

%======================================================================
% function to compute P, i.e, E-step of CPD registration
%======================================================================
function [ P1, Pt1, Px, negativeLogLikelihood ] = computeEStep( X, Y, sigma2, w )
% The elements m, n of the matrix P are defined as:
%
%                      exp(C1*(Xn-Ym)^2)
% P(m,n) = ----------------------------------------,
%            sum-over-k(exp(C1*(Xn-Yk)^2)) + C2
%
% where C1 = -1/(2*sigma2), M is the number of points in Y, N is the number
% of points in X and D is the dimensionality of points (ex D=3 in case of 3
% dimensional points),
% and   C2 = ( (M*w)/(N*(1-w)) )*((2*pi*sigma2)^(D/2))
%
% The outputs are: P1 = P*1, Pt1 = P'*1, Px = P*X,
%  where 1 is a column vector of all ones.


D = size(X, 2);
M = size(Y, 1);
N = size(X, 1);

% Compute C2 as given in above equation
c   = power(2*pi*sigma2, D/2) *( w /(1-w) *M/N );
if(isSimMode())
    % Compute elements of P matrix
    pMatrix = zeros(M, N, 'like', X);
    for col = 1 : D
        pMatrix = pMatrix + (X(:,col)' - Y(:,col)).^2;
    end
    pMatrix     = exp(pMatrix*(-1/(2*sigma2)));
    
    % Compute Pt1, P1, Px
    pMatrixColSum  = sum(pMatrix);
    Pt1         = (pMatrixColSum./(pMatrixColSum +c))';
    pMatrix     = pMatrix./(pMatrixColSum +c);
    P1          = sum(pMatrix, 2);
    
    Px  = zeros(M, D, 'like', X);
    for col = 1 : D
        Px(:, col)  = pMatrix*X(:, col);
    end
    negativeLogLikelihood  = -sum(log(pMatrixColSum +c)) + D*N*log(sigma2)/2;
else
    % Codegen
    P1  = zeros(M,1, 'like', X);
    Px  = zeros(M, D, 'like', X);
    Pt1 = zeros(N, 1, 'like', X);
    
    for j = 1 : N
        pMatrixCol = (X(j,1) - Y(:,1)).^2;
        for col = 2 : D
            pMatrixCol = pMatrixCol + (X(j,col) - Y(:,col)).^2;
        end
        
        pMatrixCol = exp(pMatrixCol.*(-1/(2*sigma2)));
        pMatrixColSum =  sum(pMatrixCol);
        
        pMatrixCol = pMatrixCol./(c+pMatrixColSum);
        P1 = P1 + pMatrixCol;
        Px = Px + pMatrixCol*X(j, :);
        Pt1(j) = pMatrixColSum;
    end
    
    negativeLogLikelihood  = cast(-sum(log(Pt1 + c)) + D*N*log(sigma2)/2, 'like', X);
    Pt1 = Pt1./(Pt1+c);
end
end

%==========================================================================
% Parameter validation
%==========================================================================
function [outMoving, outFixed, movingValidPointIndices, tformType, paramsStruct] = validateAndParseInputArgs(inMoving, inFixed, varargin)

validateattributes(inMoving, {'pointCloud'}, {'scalar'});
validateattributes(inFixed, {'pointCloud'}, {'scalar'});

% Remove invalid points
[outMoving, movingValidPointIndices] = removeInvalidPoints(inMoving);
outFixed = removeInvalidPoints(inFixed);

% Check for minimum points
if outFixed.Count < 2 || outMoving.Count < 2
    coder.internal.error('vision:pointcloud:cpdNotEnoughPoints');
end

% Convert to double internally if both are not of same type
if (~isa(outFixed.Location, class(outMoving.Location)))
    if(isa(outFixed.Location, 'single'))
        outFixed = pointCloud(double(outFixed.Location));
    end
    if(isa(outMoving.Location, 'single'))
        outMoving = pointCloud(double(outMoving.Location));
    end
end

% Set input parser
if(isSimMode())
    % Simulation
    defaults = struct(...
        'Transform',             'Nonrigid', ...
        'OutlierRatio',          0.1, ...
        'Tolerance',             1e-5, ...
        'MaxIterations',         20, ...
        'InteractionSigma',      2, ...
        'SmoothingWeight',       3, ...
        'Verbose',               false);
    
    parser               = inputParser;
    parser.CaseSensitive = false;
    
    parser.addParameter('Transform', defaults.Transform, ...
        @(x)validateattributes(x,{'char', 'string'}, {'nonempty'}));
    parser.addParameter('OutlierRatio', defaults.OutlierRatio, ...
        @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative' , '>=', 0, '<', 1}));
    parser.addParameter('Tolerance', defaults.Tolerance, ...
        @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
    parser.addParameter('MaxIterations', defaults.MaxIterations, ...
        @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
    
    parser.addParameter('InteractionSigma', defaults.InteractionSigma, ...
        @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'positive'}));
    parser.addParameter('SmoothingWeight', defaults.SmoothingWeight, ...
        @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'positive'}));
    
    parser.addParameter('Verbose', defaults.Verbose, ...
        @(x)validateattributes(x, {'logical'}, {'nonsparse', 'scalar', 'nonempty'}));
    
    parser.parse(varargin{:});
    paramsStruct           = parser.Results;
    
    paramsStruct.Transform = validatestring(paramsStruct.Transform, {'Rigid', 'Affine', 'Nonrigid'}, mfilename, 'Transform');
    
    if(~strcmpi(paramsStruct.Transform, 'Nonrigid'))
        % Throw warning if 'InteractionSigma' or 'SmoothingWeight' NV pair
        % is passed when 'Transform' is not 'Nonrigid'.
        defaultInteractionSigma = strfind(parser.UsingDefaults, 'InteractionSigma');
        if(sum([defaultInteractionSigma{:}]) == 0)
            warning(message('vision:pointcloud:cpdNVParameterAvailableForNonrigidTransformOnly', 'InteractionSigma'));
        end
        defaultSmoothingWeight = strfind(parser.UsingDefaults, 'SmoothingWeight');
        if(sum([defaultSmoothingWeight{:}]) == 0)
            warning(message('vision:pointcloud:cpdNVParameterAvailableForNonrigidTransformOnly', 'SmoothingWeight'));
        end
    end
    tformType = paramsStruct.Transform;
    paramsStruct = rmfield(paramsStruct, 'Transform');
else
    % Codegen
    defaults = struct(...
        'Transform',             'Nonrigid', ...
        'OutlierRatio',          0.1, ...
        'Tolerance',             1e-5, ...
        'MaxIterations',         20, ...
        'InteractionSigma',      2, ...
        'SmoothingWeight',       3);
    
    pvPairs = struct(...
        'Transform',        uint32(0),...
        'OutlierRatio',     uint32(0),...
        'Tolerance',        uint32(0),...
        'MaxIterations',    uint32(0),...
        'InteractionSigma', uint32(0),...
        'SmoothingWeight',  uint32(0));
    
    properties =  struct( ...
        'CaseSensitivity', false, ...
        'StructExpand',    true, ...
        'PartialMatching', false);
    
    optarg = eml_parse_parameter_inputs(pvPairs, properties, varargin{:});
    
    tformType = eml_get_parameter_value(optarg.Transform, defaults.Transform, varargin{:});
    
    paramsStruct.OutlierRatio = eml_get_parameter_value(optarg.OutlierRatio, defaults.OutlierRatio, varargin{:});
    paramsStruct.Tolerance = eml_get_parameter_value(optarg.Tolerance, defaults.Tolerance, varargin{:});
    paramsStruct.MaxIterations = eml_get_parameter_value(optarg.MaxIterations, defaults.MaxIterations, varargin{:});
    paramsStruct.InteractionSigma = eml_get_parameter_value(optarg.InteractionSigma, defaults.InteractionSigma, varargin{:});
    paramsStruct.SmoothingWeight = eml_get_parameter_value(optarg.SmoothingWeight, defaults.SmoothingWeight, varargin{:});
    
    % Since Codegen doesnot support vision.internal.MessagePrinter, setting
    % verobse to false.
    paramsStruct.Verbose = false;
    
    validateattributes(tformType, {'char', 'string'}, {'nonempty'})
    validateattributes(paramsStruct.OutlierRatio, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative' , '>=', 0, '<', 1})
    validateattributes(paramsStruct.Tolerance, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'})
    validateattributes(paramsStruct.MaxIterations, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'})
    
    validateattributes(paramsStruct.InteractionSigma, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'positive'})
    validateattributes(paramsStruct.SmoothingWeight, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'positive'})
    
    tformType = validatestring(tformType, {'Rigid', 'Affine', 'Nonrigid'}, mfilename, 'Transform');
    
    if(~strcmpi(tformType,'Nonrigid'))
        if( 0 ~= optarg.InteractionSigma)
            coder.internal.warning('vision:pointcloud:cpdNVParameterAvailableForNonrigidTransformOnly', 'InteractionSigma');
        end
        if( 0 ~= optarg.SmoothingWeight)
            coder.internal.warning('vision:pointcloud:cpdNVParameterAvailableForNonrigidTransformOnly', 'SmoothingWeight');
        end
    end
end
end

function flag = isSimMode()
flag = isempty(coder.target);
end