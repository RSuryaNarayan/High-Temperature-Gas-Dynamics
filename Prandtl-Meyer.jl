# -*- coding: utf-8 -*-
# +
#PREAMBLE
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.1
#     language: julia
#     name: julia-1.4
# ---
import Pkg
Pkg.add("NLsolve")
Pkg.add("Plots")
Pkg.add("DifferentialEquations")
Pkg.add("CSV")
Pkg.add("DataFrames")
using NLsolve
using DifferentialEquations
using Plots
using CSV, DataFrames
# About this code module
# Julia Program to solve Prandtl-Meyer flow considering high temperature effects
# GAS PROPERTIES
#comment out the appropriate section/lines to select
#Oxygen
#R_A2 = 8.314*1000/32; #the molar gas constant for Oxygen
#Θ_d = 59500; #dissociation temperature in Kelvin
#ρ_d = 150000; #dissociation density
#flag = 1; #to display on the plot
#Nitrogen
R_A2 = 8.314*1000/28; #the molar gas constant for Oxygen
global Θ_d = 113000; #dissociation temperature in Kelvin
global ρ_d = 130000; #dissociation density
#flag =2; #to display on the plot
#SOLUTION
#STATE EQUATION for high temperature gas dynamics in-terms of pressure and temperature input
#Intended to calculate inlet properties
function state(P,T)
#x[1]-->α_*, x[2]-->log10(ρ_d/ρ)
    function g!(G,x)
        G[1] = P - (1+x[1])*R_A2*T*ρ_d/10^x[2]
        k = (10^(x[2]))*exp(-Θ_d/T)
        G[2] = x[1] - 0.5*(sqrt(k^2 + 4k) -k)
    end
    res = nlsolve(g!,[0.5,6])
    x = res.zero[1];#dissociation
    y = res.zero[2];#density ratio
    s = R_A2*((3*log(T/Θ_d) + (x*(1-2*log(x)))-((1-x)*log(1-x)) - (1+x)*log(1/10^y)))
    a = sqrt(R_A2*(T/Θ_d)*Θ_d*(x*(1-x^2)*(1+2*(T/Θ_d)) + (8+3*x-x^2)*(T/Θ_d)^2)/(x*(1-x) + 3*(2-x)*(T/Θ_d)^2))
    h = R_A2*((4+x)*T + x*Θ_d)
    vec = [s,a,h,T/Θ_d,log10(R_A2*ρ_d*Θ_d/P),ρ_d/10^y,x]
    return vec
end
#Mach function calculator for High-Temperature Gas Dynamics
#NON-LINEAR EQUATION SYSTEM FOR STATE
#x[1]-->α_*, x[2]-->T/Θ_d, x[3]-->log10(ρ_d/ρ)
function mach(h,s,flag,a0,t0,r0)
    s_c = s/R_A2
    h_c = h/(R_A2*Θ_d)
    if a0>=1e-5
        function f!(F, x)
            F[1] =  h_c-((4+x[1])*x[2]+x[1])
            F[2] =  exp(x[1]-s_c) - (x[1]^(2*x[1]))*((1-x[1])^(1-x[1]))*((1/10^x[3])^(1+x[1]))/(x[2]^3)
            k = ((10^x[3]))*exp(-1/x[2])
            F[3] = x[1] - 0.5*(sqrt(k^2+4k)-k)
        end
        results = nlsolve(f!,[a0,t0,r0])
        x = results.zero[1]
        y = results.zero[2]
        z = ρ_d/(10^results.zero[3])
        p = (1+x)*R_A2*z*y*Θ_d
        a = sqrt(R_A2*y*Θ_d*(x*(1-x^2)*(1+2*y) + (8+3*x-x^3)*y^2)/(x*(1-x) + 3*(2-x)*y^2))
        if flag==1
            return a
        else
            vec= [p,z,y*Θ_d,a,x]
        end
    else
        T = h/(4*R_A2)
        ρ = ρ_d*exp(-(s_c-3*log(T/Θ_d)))
        p = R_A2*T*ρ
        a = sqrt(4*R_A2*T/3)
        if flag==1
            return a
        else
            vec = [p,ρ,T,a,0]
        end
    end
end
#Main program for Prandtl-Meyer Expansion
#Upstream properties
p_a = 1215.90;#upstream pressure in pascal
T_a = 4140;#upstream temperature in kelvin
M_a = 1.05;#upstream Mach number
state_info = state(p_a,T_a);
s_a = state_info[1];#inlet entropy
u_a = M_a*state_info[2];#inlet velocity
ρ_a = state_info[6];#inlet density
alpha_a = state_info[7];#inlet dissociation
h_a = state_info[3];#inlet enthalpy
h_t_a = h_a + 0.5*u_a^2;#inlet stagnaion enthalpy
##Prandtl-Meyer Solver-Iteration 1
f(u,p,t) = 2(u-h_t_a)/sqrt((2*(h_t_a - u)/(mach(u,s_a,1,p[1], p[2],p[3]))^2)-1)
u0 = h_a
tspan = (0.0,178.0*pi/180)#theta in degrees
p = [alpha_a,T_a/Θ_d,log10(ρ_d/ρ_a)];
function affect!(integrator)
    guess = mach(integrator.u,s_a,2,p[1],p[2],p[3]);
    p[1] = guess[5];
    p[2] = guess[3]/113000;
    p[3] = log10(ρ_d/guess[2]);
end
function condition_warn(u,t,integrator)
    res = mach(u,s_a,2,p[1],p[2],p[3]);
    println(res[5])
    res[1]<1;
end
update_guess = PeriodicCallback(affect!,0.01,initial_affect = true);
pressure_drop = DiscreteCallback(condition_warn, terminate!);
cb = CallbackSet(update_guess,pressure_drop);
prob = ODEProblem(f,u0,tspan,p,callback =cb,saveat =0.01);
sol = solve(prob,Tsit5())#fifth order Tsitouras integration
#calculate other parameters
θ = (sol.t).*180 ./pi;
println("Max deflection angle is ",last(θ))
h_b = sol.u;
s_b = s_a;
#extract properties
p_b = zeros(length(h_b));
ρ_b = zeros(length(h_b));
u_b = zeros(length(h_b));
T_b = zeros(length(h_b));
a_b = zeros(length(h_b));
alpha_b = zeros(length(h_b));
#initialize the first element of the array
ans_vec = mach(h_b[1],s_b,2,0.05, 0.01,4.5);#call the function to extract the property vector
p_b[1] = ans_vec[1];
ρ_b[1] = ans_vec[2];
T_b[1] = ans_vec[3];
a_b[1] = ans_vec[4];
alpha_b[1] = ans_vec[5];
u_b[1] = sqrt(2*(h_t_a-h_b[1]));
for i = 2:length(h_b)
    ans_vec = mach(h_b[i],s_b,2,alpha_b[i-1], T_b[i-1]/Θ_d,log10(ρ_d/ρ_b[i-1]));#call the function to extract the property vector
    p_b[i] = ans_vec[1];
    ρ_b[i] = ans_vec[2];
    T_b[i] = ans_vec[3];
    a_b[i] = ans_vec[4];
    alpha_b[i] = ans_vec[5];
    u_b[i] = sqrt(2*(h_t_a-h_b[i]));
end
#Post-processing and plotting
display(plot(θ,T_b./T_a, title="Temperature ratio"));#temperature ratio
display(plot(θ,ρ_b./ρ_a,title="Density ratio"));#density ratio
display(plot(θ,p_b,title="Pressure ratio"));#pressure ratio
display(plot(θ,u_b./a_b,title="Mach Number"));#Mach number
display(plot(θ,alpha_b,title="Dissociation"));#Mach number
# -
