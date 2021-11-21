% LSG-CPD: CPD with Local Surface Geometry
% Author: Weixiao Liu, Hongtao Wu 
% Johns Hopkins University
% Nov 6th, 2020

function [tform] = LSGCPD(source, target, varargin)

% --------------------------- varargin ------------------------------------
% 判断varargin变量的规模
flag = cellfun(@isequal, varargin, repmat({'maxIter'}, size(varargin)));
%any函数判断，flag数组存在非0则为真
if any(flag)
    %idx矩阵为flag上移一行
    idx = circshift(flag,1);
    parm.maxIter = varargin{idx};
else
    parm.maxIter = 50;
end
% B = REPMAT(A,M,N)，the size of B is [size(A,1)*M, size(A,2)*N].
flag = cellfun(@isequal, varargin, repmat({'tolerance'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.tolerance = varargin{idx};
else
    parm.tolerance = 1e-3;
end

flag = cellfun(@isequal, varargin, repmat({'outlierRatio'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.w = varargin{idx};
else
    parm.w = 0.1;
end

flag = cellfun(@isequal, varargin, repmat({'truncationThreshold'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.truncate_threshold = varargin{idx};
else
    parm.truncate_threshold = 0.19;
end

flag = cellfun(@isequal, varargin, repmat({'optimizationIter'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.opti_maxIter = varargin{idx};
else
    parm.opti_maxIter = 2;
end

flag = cellfun(@isequal, varargin, repmat({'optimizationTolerance'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.opti_tolerance = varargin{idx};
else
    parm.opti_tolerance = 1e-3;
end

flag = cellfun(@isequal, varargin, repmat({'neighbours'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.neighbours = varargin{idx};
else
    parm.neighbours = 10;
end

flag = cellfun(@isequal, varargin, repmat({'maxPlaneRatio'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    parm.alimit = varargin{idx};
else
    parm.alimit = 30;
end
%前面这部分都是判断实验的各种参数
% --------------------------- Load Data -----------------------------------
flag = cellfun(@isequal, varargin, repmat({'dataType'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    %参数形式是array则output的值为1
    if strcmp(varargin{idx}, 'array')
        parm.output = 1;
        X = source;
        Y = target;
    else
%参数形式为pointcloud，output值为0
        if strcmp(varargin{idx}, 'pointCloud')
            parm.output = 0;
            X = source.Location; % source
            Y = target.Location; % target (GMM)
        else
            error('Input data type unsupported.')
        end
    end
else
    parm.output = 0;
    X=source.Location;
    Y=target.Location;
%    X = gpuArray(source.Location); % source
%    Y = gpuArray(target.Location); % target (GMM)
end

% ----------------------volume of bounding cube----------------------------
V = (max(Y(:, 1)) - min(Y(:, 1))) *  ...
    (max(Y(:, 2)) - min(Y(:, 2))) * ...
    (max(Y(:, 3)) - min(Y(:, 3))); 

% ------------------------ Confidence Filtering ---------------------------
% 判断参数是否有confidenceFilter
flag = cellfun(@isequal, varargin, repmat({'confidenceFilter'},...
    size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    if strcmp(varargin{idx}, 'true')
        % Depth of points
        X_depth = X(:, 3);
        Y_depth = Y(:, 3);
        
        % Confidence of points
        confidence_X = confidence_filter(X_depth);
        confidence_Y = confidence_filter(Y_depth);
        
        % Truncate the points if confidence < threshold
        X = X(confidence_X > parm.truncate_threshold, :);
        Y = Y(confidence_Y > parm.truncate_threshold, :);
        % Truncate the corresponding confidence
        confidence_X = confidence_X(confidence_X > parm.truncate_threshold, :);
        confidence_Y = confidence_Y(confidence_Y > parm.truncate_threshold, :);
        
        % Size of points
%         size第二个参数返回该维度的值
        N = size(X, 1);
        M = size(Y, 1);
        
        % Compute surface normal and variation
        [Normal_cpu, Curvature] = findPointNormals(gather(Y), parm.neighbours);
        
        % pi(m)
        f_Y = confidence_Y ./ sum(confidence_Y);
    else
        
        % Size of points
        N = size(X, 1);
        M = size(Y, 1);
        
        % Confidence is simply one if not weighted
        confidence_X = ones(N, 1);
        
        % Compute surface normal and variation
        [Normal_cpu, Curvature] = findPointNormals(target.Location, ...
                                                   parm.neighbours);
        
        % pi(m)
        f_Y = single(ones(size(Y, 1), 1)) ./ size(Y, 1);
    end
else
    % Size of points
    N = size(X, 1);
    M = size(Y, 1);
    
    % Confidence is simply one if not weighted
    confidence_X = ones(N, 1);
    
    % Compute surface normal and variation
    [Normal_cpu, Curvature] = findPointNormals(target.Location, ...
                                               parm.neighbours);
    
    % pi(m)
    f_Y = single(ones(size(Y, 1), 1)) ./ size(Y, 1);
end
Normal =Normal_cpu;
% Normal = gpuArray(Normal_cpu);

% ------------------ Compute centroid and xform to centroid ---------------
flag = cellfun(@isequal, varargin, repmat({'xform2center'}, size(varargin)));
parm.mean_xform = 0;
if any(flag)
    idx = circshift(flag,1);
    if strcmp(varargin{idx}, 'true')
        parm.mean_xform = 1;
        xmean = mean(X);
        ymean = mean(Y);
        X = X - xmean;
        Y = Y - ymean;
%         xmean = gather(xmean);
%         ymean = gather(ymean);
    end
end


% --------------- Pre-calculations and transpose of locations ---------------
Y_Normal = arrayfun(@(a1, a2, a3, b1, b2, b3) a1 * b1 + a2 * b2 + a3 * b3, ...
    Y(:, 1), Y(:, 2), Y(:, 3), Normal(:, 1), Normal(:, 2), Normal(:, 3));

X = X';
Y = Y';
X_X2 = vecnorm(X).^2;
Y_Y2 = (vecnorm(Y).^2)';
X_Y = X_X2 + Y_Y2;
% 
X_cpu = X;
% X_cpu = gather(X);
Y_cpu = Y;
% Y_cpu = gather(Y);

E1 = [0 0 0; 0 0 -1; 0 1 0];
E2 = [0 0 1; 0 0 0; -1 0 0];
E3 = [0 -1 0; 1 0 0; 0 0 0];
E1X_cpu = E1 * X;
E2X_cpu = E2 * X;
E3X_cpu = E3 * X;


% 
% E1X_cpu = gather(E1 * X);
% E2X_cpu = gather(E2 * X);
% E3X_cpu = gather(E3 * X);

% --------------- Initialize sigma2 ---------------
flag = cellfun(@isequal, varargin, repmat({'sigma2'}, size(varargin)));
if any(flag)
    idx = circshift(flag,1);
    sigma2 = varargin{idx};
else
    sigma2 = 0;
    for i = 1:3
        sigma2 = sigma2 + sum(sum((X(i, :)' - Y(i, :)).^2));
    end
    sigma2 = sigma2 / (3 * M * N);
end

% Pre-calculations---------------------------------------------------------
parm.lambda = 0.2; % lambda
a = max(2 ./ (1 + exp(parm.lambda .* (3 - 1 ./ Curvature))) - 1, 0) .* parm.alimit;
vol = (a + 1) .^ (1 / 2);
f_Y = f_Y .* vol;
% estimate outlier weight--------------------------------------------------
w0 = V * parm.w * reshape(f_Y, 1, []) * ...
    single((2 * pi * sigma2) ^ (- 3 / 2) .* (a + 1) .^ (1 / 2)); 
w0 = w0 / (1 - parm.w + w0);
wn = reshape(1 - (1 - w0) .* confidence_X, 1, []);
f_X = (1 - wn) ./ wn;
F_matrix = f_X .* f_Y;

invSigma_flatten_const = zeros(M, 9); % M x 9
y_invSigma_const = zeros(M, 3); % M x 3
y_invSigma_y_const = zeros(M, 1);

for m = 1:M
    invSigma_const_m = a(m) .* Normal_cpu(m, :)' * Normal_cpu(m, :) + eye(3);
    invSigma_flatten_const(m, :) = reshape(invSigma_const_m, 1, []);
    y_invSigma_const(m, :) = Y_cpu(:, m)' * invSigma_const_m;
    y_invSigma_y_const(m) = Y_cpu(:, m)' * invSigma_const_m * Y_cpu(:, m);
end
invSigma_flatten_const = single(invSigma_flatten_const);
% invSigma_flatten_const = gpuArray(single(invSigma_flatten_const));
y_invSigma_const = single(y_invSigma_const);
% y_invSigma_const = gpuArray(single(y_invSigma_const));
y_invSigma_y_const = single(y_invSigma_y_const);
% y_invSigma_y_const = gpuArray(single(y_invSigma_y_const));

% -------------------------------EM-process--------------------------------
% initialize EM parameters
iter = 0;
loglikelihood = 0;
R = eye(3);
t = [0; 0; 0];

% start EM iterations
while iter <= parm.maxIter
    loglikelihood_prev = loglikelihood;
    
    % -----------------------------E-step----------------------------------
    C = single((2 * pi * sigma2) ^ (3 / 2) * (1 / V));
    c = -1 / (2 * sigma2);
      RX =R * X_cpu;
%     RX = gpuArray(R * X_cpu);
    [P, M_0, M_1, M_2] = E_step(RX, Y, Normal, t, a, X_Y, ...
                                F_matrix, Y_Normal, N, C, c,...
                                invSigma_flatten_const, ...
                                y_invSigma_const, y_invSigma_y_const);
    
    % -----------------------------M-step----------------------------------
    [R, t] = NewtonSE3(R, t, M_0, M_1, X_cpu, N, E1X_cpu, E2X_cpu, E3X_cpu, ...
                       parm.opti_maxIter, parm.opti_tolerance);
    iter = iter + 1;
    
    % ----------------Shrinking and Convergence Checking-------------------
    [loglikelihood, sigma2] = Shrink_step(R, t, X_cpu, P, M_0, M_1, M_2, N);  
    % Check convergency
    if abs(loglikelihood - loglikelihood_prev) / loglikelihood < ...
            parm.tolerance || loglikelihood < 1e-5
        break
    end
    
end
if parm.mean_xform == 1
    t = t + ymean' - (R * xmean');
end

% return transform

if parm.output == 1
    tform = [R, t; 0 0 0 1];
else
    tform = rigid3d(R', t');
end
end

%---------------------------Utility-Functions------------------------------
function [P, M_0, M_1, M_2] = E_step(RX, Y, Normal, t, a, X_Y, ...
    F_matrix, Y_Normal, N, C, c,...
    invSigma_flatten_const, y_invSigma_const, y_invSigma_y_const)

P = F_matrix .* exp(c * (a .* ((Normal * RX + Normal * t - Y_Normal) .^ 2) + ...
    (X_Y + t' * t + 2 * (t' * RX - Y' * RX - Y' * t))));
denominator = sum(P) + C;
P = P ./ denominator;

M_0_flatten = P' * invSigma_flatten_const; % N x 9
M_0 = gather(reshape(M_0_flatten', 3, 3, N));
M_1 = gather(P' * y_invSigma_const); % N x 3
M_2 = gather(P' * y_invSigma_y_const); % N

end

%--------------------------------------------------------------------------
function [g_gradient, H] = GradientSE3(R, t, M_0, M_1, X, N, E1X, E2X, E3X)

M_1_flatten = reshape(M_1', 1, 3 * N);
gX_flatten = reshape(R * X + t, 3 * N, 1);

g_E1_X = reshape(R * E1X + t, 1, 3, N); % 1 * 3 * N
g_E2_X = reshape(R * E2X + t, 1, 3, N); % 1 * 3 * N
g_E3_X = reshape(R * E3X + t, 1, 3, N); % 1 * 3 * N
g_E4_X = reshape(R(:, 1), 1, 3);
g_E5_X = reshape(R(:, 2), 1, 3);
g_E6_X = reshape(R(:, 3), 1, 3);

g_E1_X_flatten = reshape(g_E1_X, 3 * N, 1); % 3N * 1
g_E2_X_flatten = reshape(g_E2_X, 3 * N, 1); % 3N * 1
g_E3_X_flatten = reshape(g_E3_X, 3 * N, 1); % 3N * 1
g_E4_X_flatten = reshape(repmat(R(:, 1), 1, N), 3 * N, 1); % 3N * 1
g_E5_X_flatten = reshape(repmat(R(:, 2), 1, N), 3 * N, 1); % 3N * 1
g_E6_X_flatten = reshape(repmat(R(:, 3), 1, N), 3 * N, 1); % 3N * 1

gE1X_M_0 = reshape(pagemtimes(g_E1_X, M_0), 1, 3 * N);
gE2X_M_0 = reshape(pagemtimes(g_E2_X, M_0), 1, 3 * N);
gE3X_M_0 = reshape(pagemtimes(g_E3_X, M_0), 1, 3 * N);
gE4X_M_0 = reshape(pagemtimes(g_E4_X, M_0), 1, 3 * N);
gE5X_M_0 = reshape(pagemtimes(g_E5_X, M_0), 1, 3 * N);
gE6X_M_0 = reshape(pagemtimes(g_E6_X, M_0), 1, 3 * N);

g_gradient = 2 .* ([gE1X_M_0; gE2X_M_0; gE3X_M_0; gE4X_M_0; gE5X_M_0; ...
                    gE6X_M_0] * gX_flatten  - (M_1_flatten * ...
                    [g_E1_X_flatten, g_E2_X_flatten, g_E3_X_flatten,...
                     g_E4_X_flatten, g_E5_X_flatten, g_E6_X_flatten])');

H = 2 .* ([gE1X_M_0; gE2X_M_0; gE3X_M_0; gE4X_M_0; gE5X_M_0; gE6X_M_0] * ...
          [g_E1_X_flatten, g_E2_X_flatten, g_E3_X_flatten, ...
           g_E4_X_flatten, g_E5_X_flatten, g_E6_X_flatten]);
end

%--------------------------------------------------------------------------
function [R, t] = NewtonSE3(R, t, M_0, M_1, X, N, E1X, E2X, E3X, maxIter, tolerance)
iter = 1;
while iter <= maxIter
    % Calculate gradient and Hessian Matrix
    [g_gradient, H] = GradientSE3(R, t, M_0, M_1, X, N, E1X, E2X, E3X);
    % Check optimality
%     norm默认为二范数
    if norm(g_gradient) <= tolerance
        break
    else
    % Newton update
        x_opti = - (pinv(((1 / 2) .* (H + H'))) * eye(6)) * g_gradient;
        X_opti = [0 -x_opti(3) x_opti(2) x_opti(4); ...
            x_opti(3) 0 -x_opti(1) x_opti(5); ...
            -x_opti(2) x_opti(1) 0 x_opti(6);
            0 0 0 0];
        g = [R t; 0 0 0 1] * expm(X_opti);
        R = g(1 : 3, 1 : 3);
        t = g(1 : 3, 4);
    end
    iter = iter + 1;
end

end

%--------------------------------------------------------------------------
function [loglikelihood, sigma2] = Shrink_step(R, t, X, P, M_0, M_1, M_2, N)

gX = R * X + t;
gX_T = reshape(gX, 1, 3, N);
gX   = reshape(gX, 3, 1, N);

% Calculate Log-likelihood
loglikelihood = sum(pagemtimes(gX_T, pagemtimes(M_0, gX)));
loglikelihood = loglikelihood - 2 * reshape(gX, 3*N, 1)' * reshape(M_1', 3*N, 1);
loglikelihood = loglikelihood + sum(M_2);

sum_P = sum(sum(P));

% Shrink local covariances
sigma2 = gather(loglikelihood / (3 * sum_P));
end

%--------------------------------------------------------------------------

function [confidence] = confidence_filter(depth)
% 将depth的三个维度赋三个权重
p1 = 0.002203;
p2 = -0.001028;
p3 = 0.0005351;
min_depth = 0.4;

error = p1 * depth.^2 + p2 * depth + p3;
confidence = (p1 * min_depth + p2 * min_depth + p3) ./ error;
end

%--------------------------------end---------------------------------------


