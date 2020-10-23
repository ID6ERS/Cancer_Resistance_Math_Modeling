function S = model(B)
global g;
global opt;
t = g(:, 1);
g1 = g(:, 2:end);
n = size(g1,2);
x_in_list = zeros(2,1);
for i = 1:2:n
    x_in_list = [x_in_list g1(1, i:(i+1))'];
end
x_in_list = x_in_list(:, 2:end);
res = t;
for i = 1:size(x_in_list,2)
    x_init = x_in_list(:,i);
    [T,Sv] = ode23s(@DifEq, t, x_init);
    if (size(Sv,1) == 61)
        res = [res Sv];
    else
        res = [res g(:, (2*i-1):(2*i))];
    end
end

      function dS = DifEq(t, x)
          xdot = zeros(2,1);
          R = x(1);
          G = x(2);
          rf = R/(R+G);
          gf = G/(R+G);
          if (opt == 1)
              Kr = 119590/(1+gf);
              Kg = 100520/(1+rf);
          elseif (opt == 2)
              Kr = 119590*gf/(1+gf);
              Kg = 100520/(1+rf);
          elseif (opt == 3)
              Kr = 119590/(1+gf);
              Kg = 100520*rf/(1+rf);
          else
              Kr = 119590*gf/(1+gf);
              Kg = 100520*rf/(1+rf);
          end
   
          xdot(1) = R*0.0262.*(1-(B(1)*R/Kr + B(2)*G/Kg));
          xdot(2) = G*0.0377.*(1-(B(3)*R/Kg + B(4)*G/Kg));
          dS = xdot;
      end
S = norm(res-g);
end