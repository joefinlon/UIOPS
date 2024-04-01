function area=single_area(diameter, habit, iNewHabit)
% 
% Used to calculate the particle area using A-D relations
% This is in cgs units
%  
% Created on Feb 14, 2014 by Will Wu
%
area = 0;
if 0 == iNewHabit
    switch char(habit)
        case {'l','o'}  % Linear and Oriented
            if diameter<0.03
                a=0.0696;
                b=1.5;
            else
                a=0.0512; %0.000907;
                b=1.414;
            end
            area = a *diameter^b;
        case 'h'  % Plate
            a=0.65;
            b=2.0;
            area = a *diameter^b;
        case {'i','a','s','t'}       % Inregular 
            a=0.2285;
            b=1.88;
            area = a *diameter^b;
        case {'g'}       % Graupel 
            a=0.5;
            b=2.0;
            area = a *diameter^b;
        case 'd'       % Dendrite  
            a=0.21;
            b=1.76;
            area = a *diameter^b;
    end
else
    switch char(habit)
        case {'l'}  % Linear
            if diameter<0.03
                a=0.0696;
                b=1.5;
            else
                a=0.0512; %0.000907;
                b=1.414;
            end
            area = a *diameter^b;
        case 'h'  % Plate
            a=0.65;
            b=2.0;
            area = a *diameter^b;
        case {'i','a','s','t','o'} % Irregular, aggregate, sphere, tiny, bad 
            a=0.2285;
            b=1.88;
            area = a *diameter^b;
        case {'g'}       % Graupel 
            a=0.5;
            b=2.0;
            area = a *diameter^b;
        case 'd'       % Dendrite  
            a=0.21;
            b=1.76;
            area = a *diameter^b;
    end
end

area = area*100; % Change unit into mm^2

end