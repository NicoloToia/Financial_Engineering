function [ratesSet_shift] = shift_rates(ratesSet)
% shift_rates: parallel shift of market data contained in the ratesSet by 1bp
%
% INPUT
% ratesSet      : market data for bid&ask rates stored as follow
%                -> ratesSet.depos      : Matix bid&ask Depos
%                -> ratesSet.futures    : Matix bid&ask Futures
%                -> ratesSet.swaps      : Matix bid&ask Swaps
%
% OUTPUT
% ratesSet_shift: new rates data shifted by 1bp

% Initialize the output struct with the same fields as the input struct
ratesSet_shift = ratesSet;

% Get the number of fields in the struct
num_fields = numel(fieldnames(ratesSet));

% Iterate over each field of the struct
for field_index = 1:num_fields
    % Extract the name of the field
    field_name = fieldnames(ratesSet);
    
    % Extract the matrix associated with that field
    matrix = ratesSet.(field_name{field_index});
    
    % Shift each value of the matrix by adding +0.0001 (1bp)
    shifted_matrix = matrix + 0.0001;
    
    % Assign the shifted matrix to the corresponding field in the output struct
    ratesSet_shift.(field_name{field_index}) = shifted_matrix;
end

end