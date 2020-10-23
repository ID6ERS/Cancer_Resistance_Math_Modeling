global g; 
%global declaration is needed because functions supplied to fminsearch can only have 1 argument.

load("data.mat");
g = table2array(Red);
B0 = [0.03; 50000]; % Intial guess
for i = 1:2000
    [B,fval,exitflag,output] = fminsearch(@model,B0);
    if sum(B0 - B) == 0 
        break;
    else
        B0 = B;
        k = i;
    end
end
Bred = B;
Ered = fval;
g = table2array(Green);
B0 = [0.03; 50000];
for i = 1:2000
    [B,fval,exitflag,output] = fminsearch(@model,B0);
    if sum(B0 - B) == 0
        break;
    else
        B0 = B;
        k = i;
    end
end
Bgreen = B;
Egreen = fval;

