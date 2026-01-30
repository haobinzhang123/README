function [gbest, gbestval, con] = DE12(fhd, dim, N, Max_Gen, VRmin, VRmax, varargin)
% --- 1. Boundary Mapping ---
if length(VRmin) == 1
    lb = VRmin * ones(dim, 1);
    ub = VRmax * ones(dim, 1);
else
    lb = VRmin(:);
    ub = VRmax(:);
end
FEsmax = Max_Gen * N;

% --- 2. Parameters ---
F = 0.5;
CR = 0.5;
rarc = 2.6;
p = 0.11;
H = 6;
NPinit = N;
NPmin = 4;
count_eval = 0;
count_iter = 1;

% --- 3. Population Initialization ---
X = lb + (ub - lb) .* rand(dim, N);

fitness = feval(fhd, X, varargin{:});
count_eval = count_eval + N;

[fitness, fidx] = sort(fitness);
X = X(:, fidx);

improved_flag = false(1, N); % Initialize improvement flag

Asize = round(rarc * N);
A = [];
nA = 0;

MF = F * ones(H, 1);
MCR = CR * ones(H, 1);
MCR(H) = 0.9;
MF(H) = 0.9;
iM = 1;

A(:, nA + 1) = X(:, 1);
Afitness(nA + 1) = fitness(1);
nA = nA + 1;

V = X;
U = X;
S_CR = zeros(1, N);
S_F = zeros(1, N);
S_df = zeros(1, N);

Chy = cauchyrnd(0, 0.1, N + 200);
iChy = 1;

con = [];
con(1) = fitness(1);

% --- 4. Main Loop ---
while count_eval < FEsmax
    
    % --- Linear Decreasing Weights ---
    WL = ceil(nA / 2);
    if WL > 0
        weights = (WL : -1 : 1)'; 
        weights = weights / sum(weights);
        Xsel = A(:, 1:WL);
        xmean = Xsel * weights;
    else
        xmean = mean(X, 2);
    end
    % -------------------
    
    pbest = 1 + floor(max(2, round(p * N)) * rand(1, N));
    r = floor(1 + H * rand(1, N));
    
    CR = MCR(r)' + 0.1 * randn(1, N);
    CR((CR < 0) | (MCR(r)' == -1)) = 0;
    CR(CR > 1) = 1;
    
    F = zeros(1, N);
    for i = 1 : N
        while F(i) <= 0
            F(i) = MF(r(i)) + Chy(iChy);
            iChy = mod(iChy, numel(Chy)) + 1;
        end
    end
    F(F > 1) = 1;
    
    PA = [X, A];
    sizePA = size(PA, 2);
    
    % --- Dynamic Mutation Strategy Allocation Logic ---
    mutation_idx = zeros(1, N); 
    
    if count_iter == 1
        mutation_idx = randi([1, 2], 1, N);
    else
        mutation_idx = randi([1, 2], 1, N);
        not_imp_indices = find(~improved_flag);
        
        if ~isempty(not_imp_indices)
            sorted_not_imp = sort(not_imp_indices, 'descend'); 
            
            n_select = round(N / 10);
            if n_select < 1, n_select = 1; end
            
            num_bad = min(length(sorted_not_imp), n_select);
            bad_indices = sorted_not_imp(1:num_bad);
            
            diffs = X(:, bad_indices) - X(:, 1);
            dists = sum(diffs.^2, 1);
            
            [~, dist_sort_idx] = sort(dists, 'descend');
            sorted_bad_indices = bad_indices(dist_sort_idx);
            
            mid_pt = ceil(num_bad / 2);
            idx_far = sorted_bad_indices(1 : mid_pt);       
            idx_near = sorted_bad_indices(mid_pt+1 : end);
            
            mutation_idx(idx_far) = 3;
            mutation_idx(idx_near) = 4;
        end
    end
    % ---------------------------
    
    for i = 1 : N
        r1 = floor(1 + N * rand);
        while i == r1
            r1 = floor(1 + N * rand);
        end
        
        r2 = floor(1 + sizePA * rand);
        while i == r1 || r1 == r2
            r2 = floor(1 + sizePA * rand);
        end
        
        strat = mutation_idx(i);
        
        if strat == 1
            V(:, i) = X(:, i) + F(i) .* (X(:, pbest(i)) - X(:, i)) + F(i) .* (X(:, r1) - PA(:, r2));
        elseif strat == 2
            V(:, i) = X(:, i) + F(i) .* (xmean - X(:, i)) + F(i) .* (X(:, r1) - PA(:, r2));
        elseif strat == 3
            vec_A = 3 * rand .* (0.5*rand) .* (X(:, 1) - X(:, i));
            vec_B = 3 * randn(dim, 1) .* cos((pi * i) / N) .* (X(:, 1) - X(:, i));
            vec_C = 3 * X(:, i);
            WB = [vec_C, vec_A, vec_B];
            kernel = 1/size(WB, 2) * ones(1, size(WB, 2)); 
            P_res = conv2(WB, kernel, 'valid');
            V(:, i) = P_res;
        elseif strat == 4
            r_rand = rand * (1 - count_iter/Max_Gen + 0.05); 
            levy_vec = Levy(dim)'; 
            step_size = r_rand .* rand(dim, 1) .* levy_vec .* X(:, i);
            newpos_temp = X(:, i) + step_size;
            rotation_angle = 0.5* rand * pi * (1 - count_iter/Max_Gen + 0.05);
            rotation_matrix = eye(dim);
            for ri = 1 : dim - 1
                rotation_matrix(ri, ri+1) = -sin(rotation_angle);
                rotation_matrix(ri+1, ri) = sin(rotation_angle);
                rotation_matrix(ri, ri) = cos(rotation_angle);
                rotation_matrix(ri+1, ri+1) = cos(rotation_angle);
            end
            rotated_diff = rotation_matrix * (newpos_temp - X(:, i));
            V(:, i) = X(:, i) + rotated_diff;
        end
        
        % Boundary Correction
        for j = 1 : dim
            if V(j, i) < lb(j)
                V(j, i) = 0.5 * (lb(j) + X(j, i));
            end
            if V(j, i) > ub(j)
                V(j, i) = 0.5 * (ub(j) + X(j, i));
            end
        end
        
        jrand = floor(1 + dim * rand);
        for j = 1 : dim
            if rand < CR(i) || j == jrand
                U(j, i) = V(j, i);
            else
                U(j, i) = X(j, i);
            end
        end
    end
    
    % Evaluation
    fu = feval(fhd, U, varargin{:});
    count_eval = count_eval + N;
    
    nS = 0;
    for i = 1 : N
        if fu(i) < fitness(i)
            nS = nS + 1;
            S_CR(nS) = CR(i);
            S_F(nS) = F(i);
            S_df(nS) = abs(fu(i) - fitness(i));
            X(:, i) = U(:, i);
            fitness(i) = fu(i);
            
            improved_flag(i) = true;
            
            if nA < Asize
                A(:, nA + 1) = X(:, i);
                Afitness(nA + 1) = fu(i);
                nA = nA + 1;
            else
                ri = floor(1 + Asize * rand);
                A(:, ri) = X(:, i);
                Afitness(ri) = fu(i);
            end
        else
            improved_flag(i) = false; 
            
            if fu(i) == fitness(i)
                X(:, i) = U(:, i);
            end
        end
    end
    
    if nS > 0
        w = S_df(1 : nS) ./ sum(S_df(1 : nS));
        if all(S_CR(1 : nS) == 0)
            MCR(iM) = -1;
        elseif MCR(iM) ~= -1
            MCR(iM) = sum(w .* S_CR(1 : nS) .* S_CR(1 : nS)) / sum(w .* S_CR(1 : nS));
        end
        MF(iM) = sum(w .* S_F(1 : nS) .* S_F(1 : nS)) / sum(w .* S_F(1 : nS));
        iM = mod(iM, H - 1) + 1;
    end
    
    % Sort Population
    [fitness, fidx] = sort(fitness);
    X = X(:, fidx);
    
    improved_flag = improved_flag(fidx);
    
    % Update NP (LPSR)
    NPinit_val = NPinit; 
    N = round(NPinit_val - (NPinit_val - NPmin) * count_eval / FEsmax);
    % Must ensure N is not smaller than NPmin
    N = max(N, NPmin); 
    
    fitness = fitness(1 : N);
    X = X(:, 1 : N);
    U = X;
    V = X;
    
    improved_flag = improved_flag(1 : N);
    
    % Resize Archive
    [Afitness, Ax] = sort(Afitness);
    A = A(:, Ax);
    Asize = round(rarc * N);
    if nA > Asize
        nA = Asize;
        A = A(:, 1 : Asize);
        Afitness = Afitness(1 : Asize);
    end
    disp(fitness(1))
    
    count_iter = count_iter + 1;
    con(count_iter) = fitness(1);
end

% --- 5. Finalize Output ---
Global_score = fitness(1);
Global_pos = X(:, 1);

gbest = Global_pos';
gbestval = Global_score;

% If iteration count is used as the criterion
actual_iters = length(con);
if actual_iters > Max_Gen
    indices = round(linspace(1, actual_iters, Max_Gen));
    con = con(indices);
elseif actual_iters < Max_Gen
    con = [con, repmat(con(end), 1, Max_Gen - actual_iters)];
end

end

% --- Helper Functions ---

function r = cauchyrnd(varargin)
a = 0.0;
b = 1.0;
n = 1;
if(nargin >= 1)
    a = varargin{1};
    if(nargin >= 2)
        b = varargin{2};
        b(b <= 0) = NaN;
        if(nargin >= 3)
            n = [varargin{3:end}];
        end
    end
end
r = cauchyinv(rand(n), a, b);
end

function x = cauchyinv(p, varargin)
a = 0.0;
b = 1.0;
if(nargin >= 2)
    a = varargin{1};
    if(nargin == 3)
        b = varargin{2};
        b(b <= 0) = NaN;
    end
end
p(p < 0 | 1 < p) = NaN;
x = a + b .* tan(pi * (p - 0.5));
if(numel(p) == 1), p = repmat(p, size(x)); end
x(p == 0) = -Inf;
x(p == 1) = Inf;
end

function o = Levy(d)
beta = 1.5;
sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
    (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
u = randn(1, d) * sigma;
v = randn(1, d);
step = u ./ abs(v).^(1 / beta);
o = step;
end