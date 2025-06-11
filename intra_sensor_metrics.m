function D = intra_sensor_metrics(beats, metric, varargin)
    [~, N] = size(beats);
    D = zeros(N);
    for i = 1:N
        for j = i+1:N
            b1 = beats(:,i);
            b2 = beats(:,j);
            switch metric
                case 'mse'
                    diff = b1 - b2;
                    D(i,j) = mean(diff.^2);
                case 'nrmse'
                    diff = b1 - b2;
                    rmse = sqrt(mean(diff.^2));
                    range_val = max([b1; b2]) - min([b1; b2]);
                    D(i,j) = rmse / range_val;
                case 'cosine'
                    b1n = b1 - mean(b1);
                    b2n = b2 - mean(b2);
                    cos_sim = dot(b1n, b2n) / (norm(b1n) * norm(b2n));
                    D(i,j) = 1 - cos_sim;
                case 'js'
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
                    D(i,j) = sqrt(jsd);
                case 'dtw'
                    D(i,j) = dtw(b1, b2);
                case 'emd'
                    p1 = b1 / sum(b1);
                    p2 = b2 / sum(b2);
                    D(i,j) = sum(abs(cumsum(p1) - cumsum(p2)));
                case 'corr'
                    D(i,j) = corr(b1, b2);
                otherwise
                    error('Unknown metric');
            end
            D(j,i) = D(i,j);
        end
        D(i,i) = 0;
    end
end