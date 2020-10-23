function S = model(B) % returns the euclidean norm between data and fit
global g;
t = g(:, 1);
g1 = g(:, 2:end);
x_in_list = g1(1,:);
res = t;
for x_init = x_in_list
    [T,Sv] = ode23s(@DifEq, t, x_init);
    res = [res Sv];
end

    function dS = DifEq(t, x)
        xdot = x*B(1)*(1-(x/B(2)));
        dS = xdot;
    end
S = norm(res-g);
end