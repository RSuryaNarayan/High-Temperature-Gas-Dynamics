# -*- coding: utf-8 -*-
# +
import Pkg
Pkg.add("NLsolve")
Pkg.add("Plots")
using NLsolve
using Plots
##Oxygen
#R_A2 = 8.314*1000/32; #the molar gas constant for Oxygen
#Θ_d = 59500; #dissociation temperature in Kelvin
#ρ_d = 150000; #dissociation density
#Nitrogen
R_A2 = 8.314*1000/28; #the molar gas constant for Nitrogen
Θ_d = 113000; #dissociation temperature in Kelvin
ρ_d = 130000; #dissociation density
γ = 1.33; #  for calculation of initial guess
c_p = 3*R_A2;; #  for calculation og initial guess

function state(h,P,t,αo)
    t_c = t/Θ_d;
    P_c = P./(ρ_d*R_A2*Θ_d)#reduced Pressure
    h_c = h./(R_A2*Θ_d)#reduced enthalpy
    function f!(F, x)
        F[1] =  h_c-((4+x[1])*x[2]+x[1])
        F[2] =  x[1]-sqrt(1/(1+(P_c/(x[2]*exp(-1/x[2])))))
    end
    results = nlsolve(f!, [α0;t_c])
    x = results.zero[1]
    y = results.zero[2]*Θ_d
    r = P/((1+x)*R_A2*y)
    a = sqrt(R_A2*(y)*(x*(1-x^2)*(1+2*(y/Θ_d)) + (8+3*x-x^2)*(y/Θ_d)^2)/(x*(1-x) + 3*(2-x)*(y/Θ_d)^2))
    if y<0 || r<0 ||  x>1
        println(y,r,x)
        error("EOS solver crashed")
    end
    vec = [x,y,r,h,P,a]
    return vec
end
#MAIN PROGRAM
#upstream properties
u_1 =collect(700:100:4000);#re-entry speed
T_a = 500;#upstream temperature
h_a = 1004*T_a;#upstream enthalpy
p_a = 0.10132500;#upstream Pressure
ρ_a = p_a/(R_A2*T_a);#upstream density
#variable initialization for initial guess
u_a = u_1[1];
m_a =  u_a/sqrt(1.33*R_A2*T_a);#Mach number assuming ideal gas
ρ_r = (2+(γ-1)*m_a^2)/(γ+1)/m_a^2; #ideal gas density ratio as initial guess
t0 = u_1[1]^2/2/c_p;#ideal gas temperature ratio
α0 = 0.05;#initial alpha
#declare result arrays
l = length(u_1);
p_r = zeros(l);
t_r = zeros(l);
r_r = zeros(l);
h_r = zeros(l);
a = zeros(l);
u_b = zeros(l);
M_b = zeros(l);
#high temperature gas dynamics solver
for i = 1:l
    global α0,t0,ρ_r
    u_a = u_1[i];
    results =zeros(0);
    err = 0;
    h_b = 0; p_b = 0;
    if i>1
        ρ_r = 1/r_r[i-1];
    end
    for k = 1:6
        #ρ_r = R_A2*T_a/u_1[1]^2;
        p_b=p_a + (ρ_a*(u_a^2))*(1 - (ρ_r));
        h_b=h_a + (0.5*(u_a^2))*(1-(ρ_r^2));
        #println(h_b);
        results = state(h_b,p_b,t0,α0);
        ρ_rn = ρ_a/results[3];
        err = abs((ρ_rn- ρ_r)/ρ_r);
        ρ_r =  ρ_rn;
        #println(p_b,h_b);
        if err < 1e-3
            #println(k)
            break;
        elseif k==6
            error("err=",err)
            error("Not converged. Try increasing k")
        end
    end
    out = [err,h_b,p_b];
    if results[2]>7000
        error("Temperature outside IDG bounds")
    end
    p_r[i] = results[5]/p_a;
    t_r[i] = results[2]/T_a;
    r_r[i] = results[3]/ρ_a;
    h_r[i] = results[4]/h_a;
    a[i] = results[1];
    M_b[i]=results[6];
    u_b[i] = ρ_a*u_a/results[3];
    α0 = results[1];
    t0 = results[2];
end
println(p_a);
display(plot(u_1 ./ 1000,t_r, seriestype = :scatter,grid = true, title = "Temperature Ratio ",label = p_a));
#savefig("<Enter location where you want to save>\\Temperature Ratio");
display(plot(u_1 ./ 1000,p_r,seriestype = :scatter,grid = true, title = "Pressure Ratio ",label = p_a));
#savefig("<Enter location where you want to save>\\Pressure Ratio");
display(plot(u_1 ./ 1000,a,seriestype = :scatter,grid = true, title = "Dissociation ",label = p_a));
#savefig("<Enter location where you want to save>\\Dissociaton Plot");
display(plot(u_1 ./ 1000,r_r ,seriestype = :scatter, grid = true, title = "Density Ratio",label = p_a));
#savefig("<Enter location where you want to save>\\Density Ratio");
display(plot(u_1 ./ 1000,h_r,seriestype = :scatter,grid = true, title = "Enthaply Ratio ",label = p_a));
#savefig("<Enter location where you want to save>\\Enthalpy ratio");
display(plot(u_1 ./ 1000,u_b,seriestype = :scatter,grid = true, title = "Downstream Velocity  ",label = p_a));
#savefig("<Enter location where you want to save>\\Post Shock Velcoity");
display(plot(u_1 ./ 1000,M_b,seriestype = :scatter,grid = true, title = "Downstream Mach number  ",label = p_a));
#savefig("<Enter location where you want to save>\\Downstream Mach Number");
# -
