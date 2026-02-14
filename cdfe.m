function [quantile,value] = cdfe(dataset, x, y, bandwidth)
sigma = 1;
SET = dataset;
L=length(SET);
Weight = zeros(L,2);
for i = 1 : L
    Weight(i,1) = exp( -abs((SET(i,1) - x)/bandwidth)^2 / (2 * sigma^2));
    Weight(i,2) = SET(i,2);
end
total_weight = 0;
for i = 1 : L
    total_weight = total_weight + Weight(i,1);
end
accumulate_weight = 0;
Weight = sortrows(Weight,2);
for i = 1 : L
    if Weight(i,2) <= y
        accumulate_weight = accumulate_weight + Weight(i,1);
    else
        break;
    end
end
quantile = accumulate_weight / total_weight;
value = y;

end