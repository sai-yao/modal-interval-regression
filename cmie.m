function [p_lower,p_upper,p_lower_value,p_upper_value] = cmie(dataset,x,alpha,bandwidth)

SORT = dataset;
L = length(SORT);
number_alpha = length(alpha);
Q = zeros(L,2);
p_lower = zeros(1,number_alpha);
p_upper = zeros(1,number_alpha);
p_lower_value = zeros(number_alpha,1);
p_upper_value = zeros(number_alpha,1);

for i = 1:L
    [Q(i,1),Q(i,2)] = cdfe(dataset,x,SORT(i,2),bandwidth);
end

for k = 1:number_alpha
    Modal_Interval_Size = Inf;
    for i = 1:L
        if Q(i,1)>1-alpha(k)
            break;
        end
        for j = i:L
            if Q(j,1)-Q(i,1)>=alpha(k)
                Interval_Size = SORT(j,2) - SORT(i,2);
                if Modal_Interval_Size > Interval_Size
                    Modal_Interval_Size = Interval_Size;
                    p_lower(1,k) = Q(i,1);
                    p_upper(1,k) = Q(j,1);
                    p_lower_value(k,1) = Q(i,2);
                    p_upper_value(k,1) = Q(j,2);
                end
                break;
            end
        end
    end   
end

end