%syndrome based puncture pattern search
function [pattern, G_best] = find_puncture_pattern(H_base, len, tests)
    tic
    [~, n] = size(H_base);
    G_best = {};
    pattern = [];
%     tests = 100000;
    for test = 1:tests
        % test erasure pattern
        temp = randperm(n);
        era_pos= temp(1:len);
        G = {};
        it = 0;
        curr_pattern = [];
        
        while ~isempty(era_pos)
            it = it+1;
            
            input = zeros(1, n);
            input(era_pos) = 1;

            syndrome = H_base*input';
            corr_pos = [];

            for pos = era_pos
                codes = H_base(:, pos);
                era_numbers = syndrome(find(codes));
                temp = find(era_numbers==1);
                if ~isempty(temp)
                    corr_pos = [corr_pos [pos; length(temp)]];
                end
            end
            
            if isempty(corr_pos)
                break;
            end
            [~, perm] = sort(corr_pos(2,:), 'descend');
            corr_pos = corr_pos(:, perm);
            G = {G{:}, corr_pos};
            curr_pattern = [curr_pattern corr_pos(1,:)];
            era_pos = setdiff(era_pos, corr_pos(1,:));
        end
        if isempty(era_pos)
            % compare with currently best pattern
            if isempty(pattern) || pattern_worse(G_best, G)
                G_best = G;
                pattern = curr_pattern;
            end 
        end
    end
    toc 
end

function result = pattern_worse(G1, G2)
    result = 0;
    
    if length(G1) > length(G2)
        result = 1;
    elseif length(G1) == length(G2)
        equal_length = 1;
        for it = 1:length(G1)
            if size(G1{it},2) < size(G2{it},2)
                result = 1;
                equal_length = 0;
                break;
            elseif size(G1{it},2) > size(G2{it},2)
                result = 0;
                equal_length = 0;
                break;
            end
        end
        if equal_length == 1
            for it = 1:length(G1)
                temp1 = G1{it};
                temp2 = G2{it};
                if sum(temp1(2,:)>1) < sum(temp2(2,:)>1)
                    result = 1;
                    break;
                end
            end
        end
    end
   
end

