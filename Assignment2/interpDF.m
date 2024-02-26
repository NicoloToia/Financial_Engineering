function [discount]=interpDF(targetDate, prevDate, nextDate, prevdisc ,nextdisc)

%   interpolation betwean two dates given a certain target date
         
if targetDate == prevDate
     discount=prevdisc;
else 
    discount=interp1([prevDate,nextDate],[prevdisc,nextdisc],targetDate);
end

end
 

