global g;
global opt;
load("data.mat");
g = table2array(data);
B0 = rand(4,1);
parms = zeros(4,4);

for opt = 1:4
    for i = 1:2000
    [B,fval,exitflag,output] = fminsearch(@model,B0);
    if sum(B0 - B) == 0
        break;
    else
        B0 = B;
        k = i;
    end
    end
parms(opt, :) = B;
end
params = array2table(parms, 'VariableNames',{'a','b','c','d'});


