function [xcorr_vals, lags, max_corr, lag_at_max] = cross_corr_beats(beat1, beat2)
    % Ensure the inputs are column vectors
    beat1 = beat1(:);
    beat2 = beat2(:);
    
    % Check that both beats are the same length
    if length(beat1) ~= length(beat2)
        error('Heartbeat segments must be of the same length');
    end

    % Compute normalized cross-correlation
    [xcorr_vals, lags] = xcorr(beat1, beat2, 'coeff');

    % Find maximum correlation value and its lag
    [max_corr, max_idx] = max(xcorr_vals);
    lag_at_max = lags(max_idx);

    % Optional: Plot the cross-correlation
    figure;
    plot(lags, xcorr_vals, 'LineWidth', 2);
    xlabel('Lag (samples)');
    ylabel('Normalized Cross-Correlation');
    title(sprintf('Max Correlation = %.3f at Lag = %d samples', max_corr, lag_at_max));
    grid on;
end