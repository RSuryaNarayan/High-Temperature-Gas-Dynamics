# -*- coding: utf-8 -*-
# # +
# import Pkg
# import DataFrames
# Pkg.add("NLsolve")
# Pkg.add("Plots")
# Pkg.add("DataFrames")
using Crayons.Box
using NLsolve
using Plots
##GAS PROPERTIES
#Oxygen
#R_A2 = 8.314*1000/32; #the molar gas constant for Oxygen
#Θ_d = 59500; #dissociation temperature in Kelvin
#ρ_d = 150000; #dissociation density
#Nitrogen
R_A2 = 8.314*1000/28; #the molar gas constant for Nitrogen
Θ_d = 113000; #dissociation temperature in Kelvin
ρ_d = 130000; #dissociation density
γ = 1.33; #  for calculation of initial guess
c_p = 3*R_A2;; #  for calculation og initial guess
##state function
function state(P,T)
    t_c = T/Θ_d
    P_c = P/(ρ_d*R_A2*Θ_d)#reduced Pressure
    α = sqrt(1/(1+((P_c)/(t_c*exp(-1/t_c)))))
    h = R_A2*Θ_d*((4+α)*t_c + α)
    r = P/((1+α)*R_A2*T)
    a = sqrt(R_A2*T*(α*(1-α^2)*(1+2*t_c) + (8+3*α-α^3)*t_c^2)/(α*(1-α) + 3*(2-α)*t_c^2))
    s = R_A2*(3*log(t_c) + α*(1-2log(α))-(1-α)log(1-α)-(1+α)log(r/ρ_d))
    if α<0 || r<0 ||  α>1
        error(t_c,r,α)
        error("EOS solver crashed")
    end
    vec = [h,s,r,α,T,P,a]
    return vec
end
##function to initialize guess
function initialize_guess(h,s)
    s_c = s/R_A2#reduced entropy
    h_c = h/(R_A2*Θ_d);#reduced enthalpy
    #set the  default guess
    a0 =0.55;
    t0 =0.07
    r0 = 6;
    #conditional guesses
    if h_c<0.5
        a0 = 0.075;
        t0 =0.05;
        r0 =5;
    elseif h_c>1.125
        a0 = 0.8;
        t0 =0.1;
        r0 =6;
    elseif h_c >=0.5
        a0 =0.1;
        t0=0.07;
        r0 = 4;
    elseif h_c>=0.6 && s_c<7.35
        a0 =0.16;
        t0 =0.09;
        r0 =3;
    end
    guess = [a0,t0,r0];
    return guess
end
##Mach function calculator for High-Temperature Gas Dynamics
#NON-LINEAR EQUATION SYSTEM FOR STATE
#x[1]-->α_*, x[2]-->T/Θ_d, x[3]-->log10(ρ_d/ρ)
function mach(h,s,flag,a0,t0,r0)
    s_c = s/R_A2#reduced entropy
    h_c = h/(R_A2*Θ_d)#reduced enthalpy
    if a0>1e-5
        function f!(F, x)
            F[1] =  h_c-((4+x[1])*x[2]+x[1])
            F[2] =  exp(x[1]-s_c) - (x[1]^(2*x[1]))*((1-x[1])^(1-x[1]))*((1/10^x[3])^(1+x[1]))/(x[2]^3)
            k = ((10^x[3]))*exp(-1/x[2])
            F[3] = x[1] - 0.5*(sqrt(k^2+4k)-k)
        end
        results = nlsolve(f!, [a0,t0,r0])
        x = results.zero[1]#dissociation
        y = results.zero[2]#temperature ratio
        z = ρ_d/(10^results.zero[3])#density
        p = (1+x)*R_A2*z*y*Θ_d
        a = sqrt(R_A2*y*Θ_d*(x*(1-x^2)*(1+2*y) + (8+3*x-x^3)*y^2)/(x*(1-x) + 3*(2-x)*y^2))
        if flag==1
            return a
        else
            return vec= [p,z,y*Θ_d,a,x,h,s]
        end
    else
        T = h/(4*R_A2)
        ρ = ρ_d*exp(-(s_c-3*log(T/Θ_d)))
        p = R_A2*T*ρ
        a = sqrt(4*R_A2*T/3)
        if flag==1
            return a
        else
            return vec = [p,ρ,T,a,0,h,s]
        end
    end
end
##throat property solver
function throat(h0,s0)
    h = h0*0.9;
    guess = initialize_guess(h,s0);
    a0 = guess[1];
    t0 = guess[2];
    r0 = guess[3];
    for i = 1:1000
        function f!(F,x)
            F[1] = sqrt(2*(h0-x[1]))-mach(x[1],s0,1,a0,t0,r0)
        end
        h = nlsolve(f!,[h],iterations = 10,show_trace=false)
        h=h.zero[1];
        res = mach(h,s0,2,a0,t0,r0);
        a0 = res[5];
        t0 = res[3]/Θ_d;
        r0 = log10(ρ_d/res[2]);
    end
    res = mach(h,s0,2,a0,t0,r0);
    err = abs(sqrt(2*(h0-h))-res[4]);
    println("")
    println(LIGHT_GREEN_FG,"Error in throat_property calculation = ", err);
    return res
end
##Property solver at a given area
function area_find(A_sup)
    ref_i=0;
    h = throat_properties[6];
    guess = initialize_guess(h,s_0);
    a0 = guess[1];
    t0 = guess[2];
    r0 = guess[3];
    h = throat_properties[6];
    stp = 0.01; erro = 0;
    for i = 1:200
        #h = throat_properties[6]*(100-i-1)*0.01;
        h = h-throat_properties[6]*stp;
        u = sqrt(2*(h_0-h));
        res = mach(h,s_0,2,a0,t0,r0);
        ρ = res[2];
        a0 = res[5];
        t0 = res[3]/Θ_d;
        r0 = log10(ρ_d/ρ);
        A = mass_flow_rate/(ρ*u);
        err=(A-A_sup)
        if err*erro<0
            stp = -stp*0.1;
        end
        #println("al=",A);
        if abs(err)<1e-3
            res= [mach(h,s_0,2,a0,t0,r0); u];
            return res;
            println("");
            println(RED_FG,"Supersonic Calculation error-stage-1 = ", err)
            #ref_i = (100-i-1);
            break;
        elseif i ==200
            println("");
            println(RED_FG,"Not converged.Increase Max(i)",err)
            #ref_i = (100-i-1);
            break;           
        end
        erro = err;
    end
end
## Main program for Nozzle flow
##Reservoir properties
p_0 = 121590;#reservoir pressure in pascal
T_0 = 5280.27;#reservoir temperature in kelvin
results = state(p_0,T_0);
s_0 = results[2];#constant entropy
h_0 = results[1];#reservoir enthalpy
##Calculating throat parameters
A_t = 1;#throat cross-sectional area
throat_properties = throat(h_0,s_0);
mass_flow_rate = throat_properties[2]*throat_properties[4]*A_t;
A_sup = 20 ;
other_properties = area_find(A_sup);
#println(other_properties);
println("");
println(BOLD,WHITE_FG,RED_BG("ALL PROPERTIES EXPRESSED IN SI UNITS"))
println(YELLOW_FG("PROPERTIES AT RESERVOIR"));
df1 = DataFrames.DataFrame(PRESSURE = results[6], DENSITY = results[3],TEMPERATURE = results[5], SONIC_VELOCITY = results[7], DISSOCIATION = results[4], ENTHALPY = results[1],ENTROPY = results[2]);
println(df1);
println(YELLOW_FG("PROPERTIES AT THROAT"));
df2 = DataFrames.DataFrame(PRESSURE = throat_properties[1], DENSITY = throat_properties[2],TEMPERATURE = throat_properties[3], SONIC_VELOCITY = throat_properties[4], DISSOCIATION = throat_properties[5], ENTHALPY = throat_properties[6],ENTROPY = throat_properties[7]);
println(df2);
println(YELLOW_FG("PROPERTIES AT EXIT: SUPERSONIC AREA RATIO = $A_sup"))
df3 = DataFrames.DataFrame(PRESSURE = other_properties[1], DENSITY = other_properties[2],TEMPERATURE = other_properties[3], SONIC_VELOCITY = other_properties[4], DISSOCIATION = other_properties[5], ENTHALPY = other_properties[6],MACH = other_properties[8]/other_properties[4]);
println(df3);
# -
