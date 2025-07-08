function p_significance_on_plot(pvals, y_offset, group_positions)
    % pvals: upper triangle of p-value matrix
    % y_offset: scalar offset above the max value of each box
    % group_positions: x positions of each group (usually 1:nGroups)
    
    n = size(pvals, 1);
    for i = 1:n
        for j = i+1:n
            p = pvals(i,j);
            if ~isnan(p) && p < 0.05
                % Decide how many stars
                if p < 0.001
                    stars = '***';
                elseif p < 0.01
                    stars = '**';
                else
                    stars = '*';
                end

                % Get y position (just above the boxes)
                y = max(get(gca,'YLim')) - y_offset * (j - i);
                x1 = group_positions(i);
                x2 = group_positions(j);

                % Plot horizontal line and star annotation
                plot([x1 x2], [y y], 'k', 'LineWidth', 1.5)
                text(mean([x1 x2]), y + 0.01*y, stars, 'HorizontalAlignment', 'center', 'FontSize', 12)
            end
        end
    end
end