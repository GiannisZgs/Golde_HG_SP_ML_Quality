function [M,varargout] = batch_intra_sensor_metrics(beats,metric,varargin)
    % beat: [L x N], where L = length of each beat, N = number of beats
    [~, N] = size(beats);

    M = zeros(N, N);  % Pairwise metric matrix
    if (strcmp(metric,'xcorr') | strcmp(metric,'jsd'))
        O = zeros(N, N); % Optimal lags matrix or sqrt(JSD) matrix
    end
    % Normalize each beat beforehand
    for i = 1:N
        for j = i+1:N  % Only upper triangle, avoid i == j
            b1 = beats(:, i);
            b2 = beats(:, j);
            
            switch metric
                case 'xcorr'
                    b1 = b1 - mean(b1);
                    b1 = b1 / norm(b1);

                    b2 = b2 - mean(b2);
                    b2 = b2 / norm(b2);

                    % Cross-correlation
                    [xcorr_vals,lags] = xcorr(b1, b2, 'coeff');
                    [max_corr, max_idx] = max(xcorr_vals);

                    % Symmetric assignment
                    M(i, j) = max_corr;
                    O(i,j) = lags(max_idx);
                    
                case 'mse'
                    diff = b1 - b2;
                    M(i,j) = mean(diff.^2);
                case 'nrmse'
                    diff = b1 - b2;
                    rmse = sqrt(mean(diff.^2));
                    range_val = max([b1; b2]) - min([b1; b2]);
                    M(i,j) = rmse / range_val;
                case 'cosine'
                    b1n = b1 - mean(b1);
                    b2n = b2 - mean(b2);
                    cos_sim = dot(b1n, b2n) / (norm(b1n) * norm(b2n));
                    M(i,j) = 1 - cos_sim;
                case 'jsd'
                    nbins = varargin{1};
                    [p1, edges] = histcounts(b1, nbins, 'Normalization', 'probability');
                    p2 = histcounts(b2, edges, 'Normalization', 'probability');
                    eps_val = 1e-10;
                    p1 = p1 + eps_val;
                    p2 = p2 + eps_val;
                    m = 0.5 * (p1 + p2);
                    kl1 = sum(p1 .* log(p1 ./ m));
                    kl2 = sum(p2 .* log(p2 ./ m));
                    jsd = 0.5 * (kl1 + kl2);
                    M(i,j) = jsd;
                    O(i,j) = sqrt(jsd);
                case 'dtw'
                    M(i,j) = dtw(b1, b2);
                case 'emd'
                    %Earth mover's distance
                    p1 = b1 / sum(b1);
                    p2 = b2 / sum(b2);
                    emd_val = sum(abs(cumsum(p1) - cumsum(p2)));
                    max_emd_1 = sum(abs(cumsum(p1) - cumsum(sort(p1, 'descend'))));
                    max_emd_2 = sum(abs(cumsum(p2) - cumsum(sort(p2, 'descend'))));
                    max_emd = max(max_emd_1, max_emd_2);
                    M(i,j) = emd_val / max_emd;
                otherwise
                    error('Unknown metric');
            end
            M(j,i) = M(i,j);
            if (strcmp(metric,'xcorr') | strcmp(metric,'jsd'))
                O(i,j) = O(j,i);
            end
        end
        M(i,i) = 0;
        if (strcmp(metric,'xcorr') | strcmp(metric,'jsd'))
            O(i,i) = O(i,i);
        end
    end
    if (strcmp(metric,'xcorr') | strcmp(metric,'jsd'))
        varargout{1} = O;
        if strcmp(metric,'jsd')
            varargout{2} = p1;
            varargout{3} = p2;
            varargout{4} = edges;
        end
    end
    if strcmp(metric,'emd')
        varargout{1} = p1;
        varargout{2} = p2;
    end
end