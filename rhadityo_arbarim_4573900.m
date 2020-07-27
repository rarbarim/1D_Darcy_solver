% single phase flow solver across a sandpack
% R B Arbarim : 4573900

%case sand A over sand B

clear all 
%input parameters
opt = str2double(input('press 1 for A over B, press 2 for B over A = ','s'));

if opt==1
    Da = 2e-4; %m
    Db = 4e-4; %m
    por_a = 0.1;
    por_b =0.26;
    za = 0.2; %m
    zb = 0.3; %m
    
elseif opt==2
    Db = 2e-4; %m
    Da = 4e-4; %m
    por_b = 0.1;
    por_a =0.26;
    zb = 0.2; %m
    za = 0.3; %m
end
    
N = 100; % number of grid cells
a = zeros (N,1);
b = zeros (N,1);
Da = 2e-4; %m
Db = 4e-4; %m
rho = 1000; % [kg/m3]
g = 9.81; %m/s^2
mu = 1; %cp
zw = 0.1; %m
L = za + zb; %m
Lt = L + zw; %m

dz = L/N;
na = za/dz;
nb = zb/dz;
pressure = 0; %Pa
phi_top = pressure + rho*g*Lt; %Pa
phi_bot = 0; %Pa
Q = zeros(N,1);
A = zeros(N,N);
k = zeros(N,1);
ka = por_a^3*Da^2/(180*(1-por_a)^2); %employ 180 to constitute tortuosity
kb = por_b^3*Db^2/(180*(1-por_b)^2); %employ 180 to constitute tortuosity

%Question - Part 1 Analytical
kh = (za/2 + zb/2)/(za/2/ka + zb/2/kb);

for i = 1 : 2
    if i == 1
        alpha(i) = kh/(mu*za*(za/2 + zb/2));
        beta(i) = ka/(mu*za^2/2);
    else
        alpha(i) = kb/(mu*zb^2/2);
        beta(i) = kh/(mu*zb*(za/2 + zb/2));
    end
end

for i = 1 : 2
    X(i,i) = alpha(i) + beta(i);
    if i == 1
        X(i,i+1) = -alpha(i) ;
        Y(i) = beta(i)*phi_top;
    else
        X(i,i-1) = -beta(i);
        Y(i) = alpha(i)*phi_bot;
    end
end
potential = X\Y';
pot_A = potential(1); %potential at centre of sand A
pot_B = potential(2); %potential at centre of sand B
gradA = (phi_top - pot_A)/(za/2); %Pa/m
gradB = (pot_B - phi_bot)/(zb/2); %Pa/m
Ca = phi_top - gradA*(za+zb);
Cb = phi_bot - gradB*0;
h = (dz*(N:-1:1))';
z = [0.5:-0.1:0];

 for i = 1 : length(h)
     if h(i) >= zb
         pot_anl(i) = gradA*h(i) + Ca;
     else
         pot_anl(i) = gradB*h(i) + Cb;
     end
     P(i) = pot_anl(i) - rho*g*h(i);
 end
 
%Question - Part 2 Numerical
for i = 1 : N
    if i*dz <= za
        k(i) = ka;
    else
        k(i) = kb;
    end    
end
    
for i = 1 : N
    if i == 1
    b(i) = k(i)/(mu*(dz^2)/2);
    else
    b(i) = 2*k(i).*k(i-1)/(k(i) + k(i-1))/(mu*dz^2);
    end
    if i == N
    a(i) = k(i)/(mu*(dz^2)/2);
    else
    a(i) = 2*k(i).*k(i+1)/(k(i) + k(i+1))/(mu*dz^2);
    end
end

for i = 1 : N
    if i == 1
        Q(i) = b(i)*phi_top;
    elseif i == N
        Q(i) = a(i)*phi_bot;
    end
    A(i,i) = a(i) + b(i);
    if i == 1
        A(i,i+1)= -a(i);
    elseif i == N
        A(i,i-1) = -b(i);
    else
        A(i,i+1)= -a(i);
        A(i,i-1)= -b(i);
    end
end

pot = A\Q;

for i = 1 : N
    p(i,1) = pot(i) - rho*g*h(i);
end

hold on
plot(p,h); %plot pressure vs height for numerical approach
plot(pot,h); %plot potential vs height for numerical approach
plot(P,h,'--');
plot(pot_anl,h,'--');
xlabel('potential/pressure (Pa)');
ylabel('height (m)');
legend('pressure(numerical)','potential (numerical)','pressure (analytical)','potential (analytical)');
hold off

