% =========================================================================================
% CUBLOCK: Optimized Version with Vectorized Polynomial Fitting
% =========================================================================================
% Original Method: Valentin Junet, et al. (Bioinformatics, 2021)
% Modifications: 
% - Vectorized sorting and polynomial estimation
% - Minimized redundant calculations inside loops
% =========================================================================================

function dataN = CuBlock(data,N,k)
if nargin<3 || isempty(k)
    k = 5;
end
if nargin<2 || isempty(N)
    N = 30;
end

data = double(data);
[nbProbes,nbSamples] = size(data);
dataN = zeros(nbProbes,nbSamples);
count = dataN;

gene_vars = var(data, 0, 2);
valid_genes = (gene_vars > 0) & (!any(isnan(data), 2));
num_valid = sum(valid_genes);

if num_valid < k
    error('CuBlock failure: only %d genes with non-zero variance found, but k=%d clusters requested.', num_valid, k);
end

% Pre-calculate p powers for target calculation
p_vals = 3:2:21;
num_p = numel(p_vals);
tol = 1e-1;

for nRep = 1:N
    try
        indProbes_valid = kmeans(data(valid_genes, :), k, 'maxiter', 1000);
    catch err
        error('CuBlock failure in k-means: %s', err.message);
    end
    
    indProbes = zeros(nbProbes, 1);
    indProbes(valid_genes) = indProbes_valid;
    
    for i=1:k
        cluster_mask = (indProbes == i);
        n_cluster = sum(cluster_mask);
        
        if n_cluster > 100
            % Pre-calculate target values for this cluster size
            X_all = (linspace(-1,1,n_cluster)').^p_vals;
            
            for j=1:nbSamples
                dataCurr = data(cluster_mask,j);
                dataCurrStd = std(dataCurr,'omitnan');
                
                if dataCurrStd > 0
                    % Z-transform
                    dataCurr_z = (dataCurr - mean(dataCurr,'omitnan'))/dataCurrStd;
                    [dataCurrS, indS] = sort(dataCurr_z);
                    
                    % Get target index
                    [~,indStdUp] = min(abs(dataCurrS-1));
                    [~,indStdDown] = min(abs(dataCurrS+1));
                    S = mean(abs(X_all(indStdDown:indStdUp,:)));
                    indP = min([num_p, find(S<tol, 1)]);
                    if isempty(indP); indP = num_p; end
                    
                    target = X_all(:, indP);
                    
                    % Vectorized cubic polyfit: pol = polyfit(dataCurrS, target, 3)
                    % Matrix form: V * p = target, where V is Vandermonde matrix
                    V = [dataCurrS.^3, dataCurrS.^2, dataCurrS, ones(n_cluster, 1)];
                    pol = V \ target;
                    
                    % ModPol logic integrated
                    dataNS = V * pol;
                    
                    % ModPol adjustments
                    diff_vals = dataNS(2:end) - dataNS(1:(end-1));
                    changeInDirectionDown = diff_vals < 0;
                    indDown1 = find(changeInDirectionDown, 1);
                    
                    if ~isempty(indDown1)
                        indDownL = find(changeInDirectionDown, 1, 'last') + 1;
                        changeInDirectionUp = diff_vals > 0;
                        indUp1 = find(changeInDirectionUp, 1);
                        indUpL = find(changeInDirectionUp, 1, 'last') + 1;
                        
                        if ~isempty(indUp1) && ~isempty(indUpL) && dataNS(indUp1) == dataNS(1) && dataNS(indUpL) == dataNS(n_cluster)
                            M = dataNS(indDown1);
                            m = dataNS(indDownL);
                            if (n_cluster-indDownL+1) <= indDown1
                                dataNS(indDown1:indDownL) = M;
                                dataNS((indDownL+1):end) = M + dataNS((indDownL+1):end) - m;
                            else
                                dataNS(indDown1:indDownL) = m;
                                dataNS(1:(indDown1-1)) = m + dataNS(1:(indDown1-1)) - M;
                            end
                        else
                            indUpL_final = find(changeInDirectionUp, 1, 'last') + 1;
                            indUp1_final = find(changeInDirectionUp, 1);
                            if ~isempty(indUpL_final); dataNS(indUpL_final:end) = dataNS(indUpL_final); end
                            if ~isempty(indUp1_final); dataNS(1:indUp1_final) = dataNS(indUp1_final); end
                        end
                    end
                    
                    % Store
                    currDataN = zeros(n_cluster, 1);
                    currDataN(indS) = dataNS;
                    
                    dataN(cluster_mask, j) = dataN(cluster_mask, j) + currDataN;
                    count(cluster_mask, j) = count(cluster_mask, j) + 1;
                end
            end
        end
    end
end

dataN_final = data;
idx_normalized = count > 0;
dataN_final(idx_normalized) = dataN(idx_normalized) ./ count(idx_normalized);
dataN = dataN_final;
end
