function [rates_Set_shift] = shift_rates(ratesSet, datesSet, bucketDate,bp)

    % shift rates by 1bp only in the bucket date
    
    DeposDate = 3; 
    FutureDate = 7;
    SwapDate = 2; 
    dates = [
    datesSet.settlement;
    datesSet.depos(1:DeposDate);
    datesSet.futures(1:FutureDate,2);
    datesSet.swaps(SwapDate:end)];
    
    %dates = datetime(dates, 'ConvertFrom', 'datenum');
    rates_Set_shift = ratesSet;
    
    bucketDate = datenum(bucketDate);
    
    %if bucket date in datesset is present, extract the index -> DA METTERE A POSTO CON I BUSDAYS
    if ismember(bucketDate, dates)
        position = ismember(bucketDate, dates);
        idx = find(position);
        %if idx is less than DeposDate, sum 1bp to the deposit rate corrispondent to the bucket date
        if idx <= DeposDate
        rates_Set_shift.depos(idx) = ratesSet.depos(idx) + bp;
        end 
        if idx > DeposDate+FutureDate
        rates_Set_shift.swaps(idx) = ratesSet.swaps(idx) + bp;
        end 
        %if idx is between DeposDate and DeposDate+FutureDate, sum 1bp to the future rate corrispondent to the bucket date
        if idx > DeposDate && idx <= DeposDate+FutureDate
        rates_Set_shift.futures(idx) = ratesSet.futures(idx) + bp;
        end
    else
        %find the previous and following date in dateset
        

    end
    
    end

end