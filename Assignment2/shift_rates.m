function [ratesSet_shift] = shift_rates(ratesSet)
    % Inizializza la struct di output con gli stessi campi della struct di input
    ratesSet_shift = ratesSet;

    % Ottieni il numero di campi nella struct
    num_fields = numel(fieldnames(ratesSet));

    % Itera su ciascun campo della struct
    for field_index = 1:num_fields
        % Estrai il nome del campo
        field_name = fieldnames(ratesSet);
        
        % Estrai la matrice associata a quel campo
        matrix = ratesSet.(field_name{field_index});
        
        % Shifta ogni valore della matrice aggiungendo +0.0001
        shifted_matrix = matrix + 0.0001;
        
        % Assegna la matrice shiftata al campo corrispondente nella struct di output
        ratesSet_shift.(field_name{field_index}) = shifted_matrix;
    end
end
