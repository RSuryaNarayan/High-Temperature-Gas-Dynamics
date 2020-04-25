##PREAMBLE
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
## About this code module
#Uses Ideal Dissociating gas model (LightHill) for computation
#Generates Equilibrium composition and properties
#Plots Mollier Diagram for specific gases namely Oxygen and Nitrogen
##Some preliminary packages to enable plotting and display them in the notebook
import Pkg
Pkg.add("GR")
Pkg.add("Plots")
using Plots
gr()
##GAS PROPERTIES
#comment out the appropriate section/lines to select
##Oxygen
R_A2 = 8.314*1000/16; #the molar gas constant for Oxygen
Θ_d = 59500; #dissociation temperature in Kelvin
ρ_d = 150000; #dissociation density
flag = 1; #to display on the plot
##Nitrogen
#R_A2 = 8.314*1000/28; #the molar gas constant for Oxygen
#Θ_d = 113000; #dissociation temperature in Kelvin
#ρ_d = 130000; #dissociation density
#flag =2; #to display on the plot
##PLOT 1: Degree of dissociation with Temperature
t = 0:0.005:0.16; #here t = T/Θ_d, Θ_d being the dissociation temperature
plot();
xlabel!("Temperature(K)");
ylabel!("Degree of dissociation");
for i = 4:7
    rho_r = 10^i; #here rho_r = ρ_d/ρ;
    k = rho_r.*exp.(-1 ./ t);
    alpha_star = (-k .+ sqrt.(k.^2 .+ 4 .* k)) ./ 2;
    T = t.*Θ_d; #compute the absolute temperature range from the temperature fraction
    ρ = ρ_d ./ rho_r #compute the density from the density fraction
    plot!(T,alpha_star, title = "Degree of dissociation vs. Temperature", label = "Density=$ρ");
end
#Generating the constant enthalpy lines on the alpha plot
t = 0.04:0.005:0.12; #here t = T/Θ_d, Θ_d being the dissociation temperature
for i = 1:5
    h_c = 0.5 +(i-1)*0.25; #h_c = h/(R_A2*Θ_d)
    alpha_star = (h_c .- 4 .*t) ./ (1 .+ t);
    T = t.*Θ_d; #compute the absolute temperature range from the temperature fraction
    H = h_c*R_A2*Θ_d; #compute the absolute enthalpy from the enthalpy fraction
    plot!(T,alpha_star, linestyle = :dot, label = "Enthalpy=$H");
end
display(plot!(legend = :outertopright));
#Dotted lines indicate lines of constant enthalpy
##PLOT 2: Mollier diagram
plot();
if flag ==1
    title!("Mollier Diagram for Oxygen");
else
    title!("Mollier Diagram for Nitrogen");
end
xlabel!("Entropy(kJ/K)");
ylabel!("Enthalpy(GJ/kg)");
#fix rho
t = 0.05:0.01:0.1; #here t = T/Θ_d, Θ_d being the dissociation temperature
for i = 4:8
    rho_r = 10^i; #here rho_r = ρ_d/ρ;
    k = rho_r.*exp.(-1 ./ t);
    alpha_star = (-k .+ sqrt.(k.^2 .+ 4 .* k)) ./ 2;
    h_c = (4 .+ alpha_star).*t + alpha_star; #h_c = h/(R_A2*Θ_d)
    s = 3 .*log.(t) + alpha_star.*(1 .- 2 .*log.(alpha_star)) .- (1 .- alpha_star).*log.(1 .- alpha_star) .+ (1 .+ alpha_star).*log(rho_r);#s = S/RA2
    H = h_c.*(R_A2*Θ_d); #compute the absolute enthalpy from the enthalpy fraction
    S = R_A2.*s; #compute the absolute entropy from the entropy fraction
    ρ = ρ_d ./ rho_r #compute the density from the density fraction
    plot!(S,H, linestyle = :dot, label = "Density=$ρ",legend = :none);
end
#fix t
rho_r= 4:0.005:8.5; #here rho_r = ρ_d/ρ;
for i = 1:6
    t = (i-1)*0.01 + 0.05; #here t = T/Θ_d, Θ_d being the dissociation temperature
    k = (10 .^rho_r).*exp(-1/t);
    alpha_star = (-k .+ sqrt.(k.^2 .+ 4 .* k)) ./ 2;
    h_c = (4 .+ alpha_star).*t + alpha_star; #h_c = h/(R_A2*Θ_d)
    s = 3 .*log(t) .+ alpha_star.*(1 .- 2 .*log.(alpha_star)) .- (1 .- alpha_star).*log.(1 .- alpha_star) .+ (1 .+ alpha_star).*log(10).*(rho_r);#s = S/RA2
    H = h_c.*(R_A2*Θ_d); #compute the absolute enthalpy from the enthalpy fraction
    S = R_A2.*s; #compute the absolute entropy from the entropy fraction
    ρ = ρ_d ./ rho_r #compute the density from the density fraction
    T = t.*Θ_d; #compute the absolute temperature range from the temperature fraction
    if i==2;
            temp = 0.06*Θ_d;
            plot!(S,H, label = "Temperature=$temp", legend = :none);
    else
    plot!(S,H, label = "Temperature=$T", legend = :none);
    end
end
#fix pressure
t = 0.05:0.005:0.1; #here t = T/Θ_d, Θ_d being the dissociation temperature
for i = 1:8
    k = 5+(i-1)*0.5;
    alpha_star = sqrt.(1 ./ (1 .+ ((10^-k) ./ (t.*exp.(-1 ./ t)))));
    h_c = (4 .+ alpha_star).*t .+ alpha_star; #h_c = h/(R_A2*Θ_d)
    s = 3 .*log.(t) .+ alpha_star.*(1 .- 2 .*log.(alpha_star)) .- (1 .- alpha_star).*log.(1 .- alpha_star) .+ (1 .+ alpha_star).*log.((alpha_star.^2 ./ (1 .- alpha_star)) .* (exp.(1 ./ t)));#s = S/RA2
    H = h_c.*(R_A2*Θ_d); #compute the absolute enthalpy from the enthalpy fraction
    S = R_A2.*s; #compute the absolute entropy from the entropy fraction
    T = t.*Θ_d; #compute the absolute temperature range from the temperature fraction
    P = ρ_d*R_A2*Θ_d/(10^k); #compute the pressure from the value of k
    plot!(S,H, legend = :none,label = "Pressure=$P",legendfontvalign = :bottom);
end
#fix dissociation
t = 0.05:0.005:0.1;#here t = T/Θ_d, Θ_d being the dissociation temperature
for i = 1:10
    alpha_star = 0.1+(i-1)*0.1;
    h_c = (4 .+ alpha_star).*t .+ alpha_star;#h_c = h/(R_A2*Θ_d)
    s = 3 .*log.(t) .+ alpha_star.*(1 .- 2 .*log.(alpha_star)) .- (1 .- alpha_star).*log.(1 .- alpha_star) .+ (1 .+ alpha_star).*log.((alpha_star.^2 ./ (1 .- alpha_star)) .* (exp.(1 ./ t)));#s = S/RA2
    H = h_c.*(R_A2*Θ_d); #compute the absolute enthalpy from the enthalpy fraction
    S = R_A2.*s; #compute the absolute entropy from the entropy fraction
    plot!(S,H, linestyle = :dot, label = "Dissociation = $alpha_star", legend = :none, legendfontvalign = :bottom);
end
display(plot!());
