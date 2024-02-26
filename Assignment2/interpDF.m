function [discount]=interpDF(targetDate, prevDate, nextDate, prevdisc ,nextdisc)

%   interpolation betwean two dates given a certain target date

discount=interp1([prevDate,nextDate],[prevdisc,nextdisc],targetDate);


end