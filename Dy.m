function [D]=Dy(m)
diaga_values = -ones(m, 1);
diagb_values = ones(m-1, 1);
diaga_indices = 1:m;
diagb_indices = 1:(m-1);
values = [diaga_values; diagb_values];
row_indices = [diaga_indices, diagb_indices];
col_indices = [diaga_indices, diagb_indices+1];
D= sparse(row_indices, col_indices, values, m, m);
D(end, 1) = 1;
end