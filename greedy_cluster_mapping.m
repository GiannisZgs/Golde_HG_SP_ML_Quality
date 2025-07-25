function [acc,macro_f1_score,micro_f1,assignments] = greedy_cluster_mapping(y_true, y_pred)
%Greedy cluster mapping using the Hungarian (Munkres) algorithm 

    shape_cond = any((size(y_pred) == size(y_true))==0); 
    if shape_cond
        y_true = y_true';
        assert(shape_cond)
    end
    
    n = max(max(y_true), max(y_pred));
    cost = zeros(n);
    for i = 1:n
        for j = 1:n
            cost(i,j) = -sum((y_true == i) & (y_pred == j));
        end
    end
    % Solve assignment problem using assignmunkres
    [assignments, ~] = assignmunkres(cost,1e-6);
    new_pred = zeros(size(y_pred));
    for j = 1:n
        if assignments(j) > 0
            new_pred(y_pred == j) = assignments(j);
        end
    end
    acc = sum(new_pred == y_true) / length(y_true);
    unique_vals = unique(y_true);
    f1_scores = [];
    tp = sum(new_pred == y_true);
    fp_fn = length(y_true) - tp;
    micro_f1 = 2*tp / (2*tp + fp_fn);
    for i = 1:length(unique_vals)
        true_class = unique_vals(i);
        where_class_true = (y_true == true_class);
        where_class_pred = (new_pred == true_class);
        
        per_class_tp = 0;
        per_class_fn = 0;
        per_class_fp = 0;
        for j = 1:length(where_class_true)
            if where_class_true(j) && where_class_pred(j)
                per_class_tp = per_class_tp + 1;
            elseif where_class_true(j) && ~where_class_pred(j)
                per_class_fn = per_class_fn + 1;
            elseif ~where_class_true(j) && where_class_pred(j)
                per_class_fp = per_class_fp + 1;
            end
        end
        f1_scores = [f1_scores 2*per_class_tp/(2*per_class_tp + per_class_fp + per_class_fn)]; 
    end
    macro_f1_score =  mean(f1_scores);
end
