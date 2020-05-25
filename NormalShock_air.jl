# -*- coding: utf-8 -*-
import Pkg
Pkg.add("NLsolve")
Pkg.add("DataFrames")
Pkg.add("Plots")
using NLsolve
using Crayons.Box
using DataFrames
using Plots
##equation of state
# T = T(P,ρ)
function TGAS3n(p,Rho)
#
#     INPUTS FOR SUBROUTINE :
#
#     P= PRESSURE, IN NEWTONS/M**2.
#     RHO= DENSITY, IN KG/M**3.
#    OUTPUT :
#
#     T= TEMPERATURE, IN KELVIN.
#
#p = 1.0133e05*100;Rho = 1.292*2;
r0 =1.292e00;  p0 =1.0133e05;  t0 =273.15e00;    gascon=287.06e00;
y = log10(Rho./r0);
x = log10(p./p0);
if( abs(y+4.5e00)<4.5e-02 )
    iflag = 0;
    rsave = Rho;
    ym = y;
    y = -4.5e00 + 4.5e-02;
    yhigh = y;
    Rho =(10 .^y).*r0;
    jflag = -1;
elseif( abs(y+0.5e00)<5.0e-03 )
    iflag = 1;
    rsave = Rho;
    ym = y;
    y = -0.5e00 + 0.5e-02;
    yhigh = y;
    Rho =(10 .^y).*r0;
    jflag = -1;
else
    iflag = -1;
end
opt = 3;
while (opt == 3)
    # LABEL three
    z1 = x - y;
    if( y>-0.5e00 )
        if( z1<=0.25e00 )
            t = p./(Rho.*gascon);
        elseif( z1>1.00e00 )
            if( z1>1.45e00 )
                if( z1>2.3e00 )
                    error("OUTSIDE OF VALIDITY RANGE OF CURVE FIT")
                end
                gas1 = -1.66249e00 - 8.91113e-02.*y;
                gas2 =(4.11648e00+8.78093e-02.*y).*z1;
                gas3 =(-3.09742e-03+1.99879e-03.*z1+6.85472e-05.*y).*y.*y;
                gas4 =(-1.84445e00-7.50324e-03.*y+3.05784e-01.*z1).*z1.*z1;
                gas5 = 1.11555e01 + 1.32100e00.*y;
                gas6 =(-1.71236e01-1.29190e00.*y).*z1;
                gas7 =(6.28124e-02-3.07949e-02.*z1+1.57743e-03.*y).*y.*y;
                gas8 =(8.63804e00+3.07809e-01.*y-1.42634e00.*z1).*z1.*z1;
                gas9 = exp(1.330611e02+8.979635e00.*y-7.265298e01.*z1-2.449009e00.*y.*z1);
            else
                gas1 = 8.06492e-01 + 9.91293e-02.*y;
                gas2 =(-1.70742e00-2.28264e-01.*y).*z1;
                gas3 =(5.03500e-03-6.13927e-03.*z1+1.69824e-04.*y).*y.*y;
                gas4 =(3.02351e00+1.31574e-01.*y-1.12755e00.*z1).*z1.*z1;
                gas5 = -1.17930e-01 - 2.12207e-01.*y;
                gas6 =(1.36524e00+4.05886e-01.*y).*z1;
                gas7 =(-1.88260e-02+1.65486e-02.*z1-5.11400e-04.*y).*y.*y;
                gas8 =(-2.10926e00-1.89881e-01.*y+8.79806e-01.*z1).*z1.*z1;
                gas9 = exp(1.959604e02-4.269391e01.*y-1.734931e02.*z1+3.762898e01.*y.*z1);
            end
            tnon = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8) ./ (1 .+ gas9);
            t =(10 .^tnon).*t0;
        else
            gas1 = -1.54141e-03 + 6.58337e-04.*y;
            gas2 =(9.82201e-01-3.85028e-03.*y).*z1;
            gas3 =(1.23111e-04-4.08210e-04.*z1+2.13592e-05.*y).*y.*y;
            gas4 =(3.77441e-02+4.56963e-03.*y-2.35172e-02.*z1).*z1.*z1;
            tnon = gas1 + gas2 + gas3 + gas4;
            t =(10 .^tnon).*t0;
        end
    elseif( y>-4.5e00 )
        if( z1<=0.25e00 )
            t = p./(Rho.*gascon);
        elseif( z1>0.95e00 )
            if( z1<=1.45e00 )
                gas1 = -5.12404e00 - 2.84740e-01.*y;
                gas2 =(1.54532e01+4.52475e-01.*y).*z1;
                gas3 =(-1.22881e-02+8.56845e-03.*z1-3.25256e-04.*y).*y.*y;
                gas4 =(-1.35181e01-1.68725e-01.*y+4.18451e00.*z1).*z1.*z1;
                gas5 = 7.52564e00 + 8.35238e-01.*y;
                gas6 =(-1.95558e01-1.23393e00.*y).*z1;
                gas7 =(3.34510e-02-2.34269e-02.*z1+4.81788e-04.*y).*y.*y;
                gas8 =(1.71779e01+4.54628e-01.*y-5.09936e00.*z1).*z1.*z1;
                gas9 = exp(6.148442e01-1.828123e01.*y-5.468755e01.*z1+1.562500e01.*y.*z1);
            elseif( z1>2.05e00 )
                if( z1>2.50e00 )
                    error("OUTSIDE OF VALIDITY RANGE OF CURVE FIT")
                end
                gas1 = -1.27244e01 - 1.66684e00.*y;
                gas2 =(1.72708e01+1.45307e00.*y).*z1;
                gas3 =(-3.64515e-02+1.90463e-02.*z1+4.80787e-04.*y).*y.*y;
                gas4 =(-6.97208e00-3.04323e-01.*y+9.67524e-01.*z1).*z1.*z1;
                gas5 = 7.71330e00 + 5.08340e-01.*y;
                gas6 =(-9.82110e00-4.49138e-01.*y).*z1;
                gas7 =(-9.41787e-04-2.40293e-03.*z1-8.28450e-04.*y).*y.*y;
                gas8 =(4.16530e00+9.63923e-02.*y-5.88807e-01.*z1).*z1.*z1;
                gas9 = exp(-1.092654e03-3.05312e02.*y+4.656243e02.*z1+1.312498e02.*y.*z1);
            else
                gas1 = -1.23779e01 - 1.14728e00.*y;
                gas2 =(2.41382e01+1.38957e00.*y).*z1;
                gas3 =(-3.63693e-02+2.24265e-02.*z1-3.23888e-04.*y).*y.*y;
                gas4 =(-1.42844e01-4.06553e-01.*y+2.87620e00.*z1).*z1.*z1;
                gas5 = 4.40782e00 + 1.33046e00.*y;
                gas6 =(-1.15405e01-1.59892e00.*y).*z1;
                gas7 =(5.30580e-02-3.10376e-02.*z1+4.77650e-04.*y).*y.*y;
                gas8 =(8.57309e00+4.71274e-01.*y-1.96233e00.*z1).*z1.*z1;
                gas9 = exp(1.4075e02-6.499992e00.*y-7.75e01.*z1+5.0e00.*y.*z1);
            end
            tnon = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8) ./ (1 .+ gas9);
            t =(10 .^tnon).*t0;
        else
            gas1 = 2.03910e-02 + 7.67310e-03.*y;
            gas2 =(8.48581e-01-2.93086e-02.*y).*z1;
            gas3 =(8.40269e-04-1.47701e-03.*z1+3.13687e-05.*y).*y.*y;
            gas4 =(2.67251e-01+2.37262e-02.*y-1.41973e-01.*z1).*z1.*z1;
            tnon = gas1 + gas2 + gas3 + gas4;
            t =(10 .^tnon).*t0;
        end
    elseif( z1>0.25e00 )
        if( z1>0.95e00 )
            if( z1<=1.4e00 )
                gas1 = -8.12952e00 - 8.28637e-01.*y;
                gas2 =(2.26904e01+1.41132e00.*y).*z1;
                gas3 =(-2.98633e-02+2.70066e-02.*z1-2.28103e-04.*y).*y.*y;
                gas4 =(-1.91806e01-5.78875e-01.*y+5.62580e00.*z1).*z1.*z1;
                gas5 = -3.99845e00 + 2.26369e-01.*y;
                gas6 =(2.52876e00-7.28448e-01.*y).*z1;
                gas7 =(1.09769e-02-1.83819e-02.*z1-1.51380e-04.*y).*y.*y;
                gas8 =(2.99238e00+3.91440e-01.*y-2.04463e00.*z1).*z1.*z1;
                gas9 = exp(-3.887015e01-2.908228e01.*y+4.070557e01.*z1+2.682347e01.*y.*z1);
            elseif( z1>1.95e00 )
                if( z1>2.60e00 )
                    error("OUTSIDE OF VALIDITY RANGE OF CURVE FIT")
                end
                gas1 = -2.33271e01 - 1.89958e00.*y;
                gas2 =(3.21440e01+1.68622e00.*y).*z1;
                gas3 =(-4.42123e-02+2.82629e-02.*z1+6.63272e-04.*y).*y.*y;
                gas4 =(-1.38645e01-3.40976e-01.*y+2.04466e00.*z1).*z1.*z1;
                gas5 = 8.35474e00 + 1.71347e00.*y;
                gas6 =(-1.60715e01-1.63139e00.*y).*z1;
                gas7 =(4.14641e-02-2.30068e-02.*z1+1.53246e-05.*y).*y.*y;
                gas8 =(8.70275e00+3.60966e-01.*y-1.46166e00.*z1).*z1.*z1;
                gas9 = exp(1.115884e02-6.452606e00.*y-5.337863e01.*z1+2.026986e00.*y.*z1);
            else
                gas1 = -1.98573e01 - 1.67225e00.*y;
                gas2 =(3.76159e01+2.10964e00.*y).*z1;
                gas3 =(-3.40174e-02+2.31712e-02.*z1-9.80275e-05.*y).*y.*y;
                gas4 =(-2.22215e01-6.44596e-01.*y+4.40486e00.*z1).*z1.*z1;
                gas5 = -5.36809e00 + 2.41201e-01.*y;
                gas6 =(-1.25881e00-8.62744e-01.*y).*z1;
                gas7 =(-3.79774e-03-7.81335e-03.*z1-3.80005e-04.*y).*y.*y;
                gas8 =(5.58609e00+3.78963e-01.*y-1.81566e00.*z1).*z1.*z1;
                gas9 = exp(2.08e01-2.56e01.*y+1.0e00.*z1+1.80e01.*y.*z1);
            end;
            tnon = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8) ./ (1 .+ gas9);
            t =(10. ^tnon).*t0;
        else
            gas1 = 1.23718e-01 + 1.08623e-02.*y;
            gas2 =(2.24239e-01-8.24608e-02.*y).*z1;
            gas3 =(-1.17615e-03-1.87566e-03.*z1-1.19155e-04.*y).*y.*y;
            gas4 =(1.18397e00+6.48520e-02.*y-5.52634e-01.*z1).*z1.*z1;
            tnon = gas1 + gas2 + gas3 + gas4;
            t =(10. ^tnon).*t0;
        end
    else
        t = p./(Rho.*gascon);
    end
    if( iflag<0 )
            opt = 1;
            return [p,Rho,t]
    end
    if( iflag!=0 )
        if( jflag<0 )
            thigh = t;
            y = -0.5e00 - 0.5e-02;
            ylow = y;
            Rho =(10 .^y).*r0;
            jflag = 0;
            opt = 3;
        elseif( jflag==0 )
            tlow = t;
            opt = 2;
            t = tlow +(thigh-tlow)./(yhigh-ylow).*(ym-ylow);
            Rho = rsave;
            return [p,Rho,t]
        else
            opt = 1;
            return [p,Rho,t]
        end
    elseif( jflag<0 )
        thigh = t;
        y = -4.5e00 - 4.5e-02;
        ylow = y;
        Rho =(10 .^y).*r0;
        jflag = 0;
        opt = 3;
    elseif( jflag==0 )
        tlow = t;
        opt = 2;
        t = tlow +(thigh-tlow)./(yhigh-ylow).*(ym-ylow);
        Rho = rsave;
        return [p,Rho,t]
    else
        opt = 1;
        return [p,Rho,t]
    end
end
end
##equation of state
#h=h(P,ρ)
function TGAS4n(p,Rho,h)
#     SUBROUTINE TGAS4 (P,RHO)
#
#
#     INPUTS FOR SUBROUTINE :
#
#     P = PRESSURE, IN NEWTONS/M.**2.
#     RHO = DENSITY, IN KG/M**3.
#
#     OUTPUT :
#
#     H = SPECIFIC ENTHALPY, IN (M/SEC)**2.
#p = 10147000;Rho = 4.0000;
r0 =1.292e00;  p0=1.0133e05;
y = log10(Rho./r0);
x = log10(p./p0);
if( abs(y+4.5e00)<4.5e-02 )
    iflag = 0;
    rsave = Rho;
    ym = y;
    y = -4.5e00 + 4.5e-02;
    yhigh = y;
    Rho =(10.0 .^y).*r0;
    jflag = -1;
elseif( abs(y+0.5e00)<5.0e-03 )
    iflag = 1;
    rsave = Rho;
    ym = y;
    y = -0.5e00 + 0.5e-02;
    yhigh = y;
    Rho =(10.0 .^y).*r0;
    jflag = -1;
else
    iflag = -1;
end
opt = 3;
while (opt == 3)
    #opt = 3
    z1 = x - y;
    if( y>-0.5e00 && y<=3.0e00)
        if( z1<=0.1e00 )
            gamm = 1.4017e00;
        elseif( z1>1.05e00 )
            if( z1>1.60e00 )
                if( z1>2.30e00 )
                    error("OUTSIDE OF VALIDITY RANGE OF CURVE FIT")
                end
                gas1 = 9.21537e-01 - 2.39670e-01.*y;
                gas2 =(1.30714e00+3.42990e-01.*y).*z1;
                gas3 =(-2.18847e-02+1.36691e-02.*z1-4.90274e-04.*y).*y.*y;
                gas4 =(-1.20916e00-1.10206e-01.*y+3.087920e-01.*z1).*z1.*z1;
                gas5 = -6.77089e00 - 6.90476e-02.*y;
                gas6 =(8.18168e00-9.52708e-02.*y).*z1;
                gas7 =(2.98487e-02-1.78706e-02.*z1+6.28419e-04.*y).*y.*y;
                gas8 =(-3.07662e00+6.60408e-02.*y+3.38590e-01.*z1).*z1.*z1;
                gas9 = exp(1.5916669e02+3.976192e01.*y-7.966199e01.*z1-1.66667e01.*y.*z1);
            else
                gas1 = -2.67593e-01 - 1.87457e-01.*y;
                gas2 =(5.07693e00+2.72286e-01.*y).*z1;
                gas3 =(1.04541e-02-1.42211e-02.*z1+6.38962e-04.*y).*y.*y;
                gas4 =(-5.08520e00-7.81935e-02.*y+1.58711e00.*z1).*z1.*z1;
                gas5 = 2.87969e00 + 3.9009e-01.*y;
                gas6 =(-8.06179e00-5.51250e-01.*y).*z1;
                gas7 =(-1.01903e-02+1.35906e-02.*z1-8.97772e-04.*y).*y.*y;
                gas8 =(7.29592e00+1.83861e-01.*y-2.15153e00.*z1).*z1.*z1;
                gas9 = exp(1.828573e02-3.428596e01.*y-1.51786e02.*z1+2.976212e01.*y.*z1);
            end
            gamm = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8)./(1 .+ gas9);
        else
            gas1 = -9.67488e01 + 2.05296e-01.*y;
            gas2 =(2.69927e02-1.92887e00.*y).*z1;
            gas3 =(3.78392e-01-3.24965e-01.*z1-3.61036e-03.*y).*y.*y;
            gas4 =(-2.46711e02+1.54416e00.*y+7.48760e01.*z1).*z1.*z1;
            gas5 = 9.81502e01 - 2.05448e-01.*y;
            gas6 =(-2.69913e02+1.93052e00.*y).*z1;
            gas7 =(-3.78527e-01+3.24832e-01.*z1+3.66182e-03.*y).*y.*y;
            gas8 =(2.46630e02-1.54646e00.*y-7.48980e01.*z1).*z1.*z1;
            gas9 = exp(-2.659865e01+1.564631e00.*y+2.312926e01.*z1-1.360543e00.*y.*z1);
            gamm = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8)./(1 .- gas9);
        end
    elseif( y>-4.5e00 && y<=-0.5e00 )
        if( z1<=0.1e00 )
            gamm = 1.399e00;
        elseif( z1>0.95e00 )
            if( z1<=1.50e00 )
                gas1 = -7.36684e00 - 1.13247e00.*y;
                gas2 =(2.47879e01+1.99625e00.*y).*z1;
                gas3 =(-4.91630e-02+4.16673e-02.*z1-6.58149e-04.*y).*y.*y;
                gas4 =(-2.32990e01-8.59418e-01.*y+7.19016e00.*z1).*z1.*z1;
                gas5 = -2.42647e00 + 5.57912e-01.*y;
                gas6 =(-2.03055e00-1.22031e00.*y).*z1;
                gas7 =(3.74866e-02-3.39278e-02.*z1+5.21042e-04.*y).*y.*y;
                gas8 =(7.75414e00+6.08488e-01.*y-3.68326e00.*z1).*z1.*z1;
                gas9 = exp(8.077385e01-1.273807e01.*y-6.547623e01.*z1+1.190475e01.*y.*z1);
            elseif( z1>2.00e00 )
                if( z1>2.5e00 )
                    error("OUTSIDE OF VALIDITY RANGE OF CURVE FIT")
                end
                gas1 = -3.77766e00 - 5.53738e-01.*y;
                gas2 =(6.60834e00+4.87181e-01.*y).*z1;
                gas3 =(-2.11045e-02+9.67277e-03.*z1-2.19420e-04.*y).*y.*y;
                gas4 =(-2.94754e00-1.02365e-01.*y+4.39620e-01.*z1).*z1.*z1;
                gas5 = 4.05813e01 + 3.25692e00.*y;
                gas6 =(-4.79583e01-2.53660e00.*y).*z1;
                gas7 =(9.06436e-02-3.47578e-02.*z1+1.00077e-03.*y).*y.*y;
                gas8 =(1.89040e01+4.94114e-01.*y-2.48554e00.*z1).*z1.*z1;
                gas9 = exp(5.34718e02+7.495657e01.*y-2.219822e02.*z1-3.017229e01.*y.*z1);
            else
                gas1 = 4.31520e-01 - 2.83857e-01.*y;
                gas2 =(2.27791e00+3.99159e-01.*y).*z1;
                gas3 =(-1.29444e-02+8.78724e-03.*z1-1.60583e-04.*y).*y.*y;
                gas4 =(-1.84314e00-1.28136e-01.*y+4.45362e-01.*z1).*z1.*z1;
                gas5 = -1.03883e01 - 3.58718e-01.*y;
                gas6 =(1.35068e01+1.87268e-01.*y).*z1;
                gas7 =(-4.28184e-03-9.52016e-04.*z1-4.10506e-05.*y).*y.*y;
                gas8 =(-5.63894e00-1.45626e-03.*y+7.39915e-01.*z1).*z1.*z1;
                gas9 = exp(2.949221e02+1.368660e01.*y-1.559335e02.*z1-3.787766e00.*y.*z1);
            end
            gamm = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8)./(1 .+ gas9);
        else
            gas1 = -1.33083e02 - 9.98707e00.*y;
            gas2 =(3.94734e02+2.35810e01.*y).*z1;
            gas3 =(1.43957e00-1.43175e00.*z1+1.77068e-05.*y).*y.*y;
            gas4 =(-3.84712e02-1.36367e01.*y+1.24325e02.*z1).*z1.*z1;
            gas5 = 1.34486e02 + 9.99122e00.*y;
            gas6 =(-3.94719e02-2.35853e01.*y).*z1;
            gas7 =(-1.43799e00+1.43039e00.*z1+1.44367e-04.*y).*y.*y;
            gas8 =(3.84616e02+1.36318e01.*y-1.24348e02.*z1).*z1.*z1;
            gas9 = exp(-2.141444e01+1.381584e00.*y+2.039473e01.*z1-1.315789e00.*y.*z1);
            gamm = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8)./(1 .- gas9);
        end
    elseif( y>=-7.0e00 && y<=-4.5e00 )
        if( z1>0.10e00 )
            if( z1>0.85e00 )
                if( z1<=1.30e00 )
                    gas1 = -1.05745e01 - 1.93693e00.*y;
                    gas2 =(3.07202e01+3.35578e00.*y).*z1;
                    gas3 =(-7.79965e-02+6.68790e-02.*z1-9.86882e-04.*y).*y.*y;
                    gas4 =(-2.60637e01-1.42391e00.*y+7.23223e00.*z1).*z1.*z1;
                    gas5 = -1.86342e01 + 2.41997e-02.*y;
                    gas6 =(3.20880e01-7.46914e-01.*y).*z1;
                    gas7 =(3.75161e-02-4.10125e-02.*z1+5.74637e-04.*y).*y.*y;
                    gas8 =(-1.69985e01+5.39041e-01.*y+2.56253e00.*z1).*z1.*z1;
                    gas9 = exp(2.768567e02+2.152383e01.*y-2.164837e02.*z1-1.394837e01.*y.*z1);
                elseif( z1>1.95e00 )
                    if( z1>2.60e00 )
                        error("OUTSIDE OF VALIDITY RANGE OF CURVE FIT")
                    end
                    gas1 = -8.32595e00 - 3.50219e-01.*y;
                    gas2 =(1.36455e01+3.59350e-01.*y).*z1;
                    gas3 =(-3.70109e-03+3.30836e-03.*z1+1.10018e-04.*y).*y.*y;
                    gas4 =(-6.49007e00-8.38594e-02.*y+1.02443e00.*z1).*z1.*z1;
                    gas5 = -3.08441e01 - 1.49510e00.*y;
                    gas6 =(3.00585e01+9.19650e-01.*y).*z1;
                    gas7 =(-3.60024e-02+1.02522e-02.*z1-4.68760e-04.*y).*y.*y;
                    gas8 =(-9.33522e00-1.35228e-01.*y+8.92634e-01.*z1).*z1.*z1;
                    gas9 = exp(8.800047e01-1.679356e01.*y-3.333353e01.*z1+8.465574e00.*y.*z1);
                else
                    gas1 = 6.17584e-01 - 2.40690e-01.*y;
                    gas2 =(1.95904e00+3.41644e-01.*y).*z1;
                    gas3 =(-1.01073e-02+6.77631e-03.*z1-1.15922e-04.*y).*y.*y;
                    gas4 =(-1.68951e00-1.10932e-01.*y+4.26058e-01.*z1).*z1.*z1;
                    gas5 = -1.34222e01 - 5.43713e-01.*y;
                    gas6 =(1.81528e01+3.95928e-01.*y).*z1;
                    gas7 =(-7.41105e-03+1.67768e-03.*z1-3.32714e-06.*y).*y.*y;
                    gas8 =(-7.97425e00-5.80593e-02.*y+1.12448e00.*z1).*z1.*z1;
                    gas9 = exp(8.677803e01-8.370349e00.*y-4.074084e01.*z1+7.407405e00.*y.*z1);
                end
                gamm = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8)./(1 .+ gas9);
            else
                gas1 = 2.53908e02 + 1.01491e02.*y;
                gas2 =(-3.87199e02-1.54304e02.*y).*z1;
                gas3 =(7.28532e00-8.04378e00.*z1-1.82577e-03.*y).*y.*y;
                gas4 =(9.86233e01+4.63763e01.*y+2.18994e01.*z1).*z1.*z1;
                gas5 = -2.52423e02 - 1.01445e02.*y;
                gas6 =(3.87210e02+1.54298e02.*y).*z1;
                gas7 =(-7.2773e00+8.04277e00.*z1+2.28399e-03.*y).*y.*y;
                gas8 =(-9.87576e01-4.63883e01.*y-2.19438e01.*z1).*z1.*z1;
                gas9 = exp(-11.0e00+2.0e00.*y+11.0e00.*z1-2.0e00.*y.*z1);
                gamm = gas1 + gas2 + gas3 + gas4 +(gas5+gas6+gas7+gas8)./(1 .- gas9);
            end
        else
            gamm = 1.3986e00;
        end
    end
    h = gamm./(gamm-1.0e00).*p./Rho;
    if( iflag<0 )
        opt = 1;
        return [p,Rho,h]
    end
    if( iflag!=0 )
        if( jflag<0 )
            hhigh = h;
            y = -0.5e00 - 0.5e-02;
            ylow = y;
            Rho =(10 .^y).*r0;
            jflag = 0;
            opt = 3;
        elseif( jflag==0 ) ;
            hlow = h;
            opt = 2;
            h = hlow +(hhigh-hlow)./(yhigh-ylow).*(ym-ylow);
            Rho = rsave;
            return [p,Rho,h]
        else
            opt = 1;
            return [p,Rho,h]
        end
    elseif( jflag<0 )
        hhigh = h;
        y = -4.5e00 - 4.5e-02;
        ylow = y;
        Rho =(10 .^y).*r0;
        jflag = 0;
        opt = 3;
    elseif( jflag==0 )
        hlow = h;
        opt = 2;
        h = hlow +(hhigh-hlow)./(yhigh-ylow).*(ym-ylow);
        Rho = rsave;
        return [p,Rho,h]
    else
        opt = 1;
        return [p,Rho,h]
    end
end
end
##Normal Shock solver
function normal_shock(P1,ρ1,u1,Pg,ρg)
    state = TGAS4n(P1,ρ1,0)
    h1= state[3]
    h_0 = h1 + u1^2/2
    m_dot = ρ1*u1
    P_0 = P1+ρ1*u1^2
    #x[1]-->P2, x[2]-->ρ2
    function f!(F,x)
        F[1] = x[2]*(P_0-x[1])-(m_dot)^2
        s = TGAS4n(x[1],x[2],0)
        h = s[3]
        F[2] = x[2] - m_dot/sqrt(2*(h_0-h))
    end
    #initial guess
   # R = 287;
   # γ = 1.4;
   # T1 = P1/ρ1/R;
   # M = u1/sqrt(γ*R*T1)
   # ρg = ρ1*2.4*M^2/(.4*M^2+2)
   # Pg = P1*(2.8*M^2-.4)/2.4
    #solve
    res = nlsolve(f!,[Pg,ρg],show_trace = false)
    P2 = res.zero[1]
    ρ2 = res.zero[2]
    u2 = m_dot/ρ2
    state = TGAS3n(P2,ρ2)
    #println(state)
    T2 = state[3]
    h2 = h_0 - u2^2/2
    vec = [P2,ρ2,T2,u2,h2]
    return vec
end
## Program main
#upstream properties
P1 = 101325; #upstream pressure
ρ1 = 1.225;#upstream density
#u1 = 1800;#upstream velocity
u1 =collect(800:10:6000)
l = length(u1);
P2 = zeros(l)
ρ2 = zeros(l)
T2 = zeros(l)
h2 = zeros(l)
u2 = zeros(l)
s = TGAS3n(P1,ρ1);
T1 = s[3];
h1 = TGAS4n(P1,ρ1,0);
h1 = h1[3]
#downstream properties
for i = 1:l
    if i==1
        R = 287;
        γ = 1.4;
        T1 = P1/ρ1/R;
        M = u1[i]/sqrt(γ*R*T1)
        ρg = ρ1*2.4*M^2/(.4*M^2+2)
        Pg = P1*(2.8*M^2-.4)/2.4
        downstream = normal_shock(P1,ρ1,u1[i],Pg,ρg);
        P2[i] = downstream[1];
        ρ2[i] = downstream[2];
        T2[i] = downstream[3];
        u2[i] = downstream[4];
        h2[i] = downstream[5];
    else
        downstream = normal_shock(P1,ρ1,u1[i],P2[i-1],ρ2[i-1]);
        P2[i] = downstream[1];
        ρ2[i] = downstream[2];
        T2[i] = downstream[3];
        u2[i] = downstream[4];
        h2[i] = downstream[5];
    end
end
println("")
println(BOLD,WHITE_FG,RED_BG("ALL PROPERTIES EXPRESSED IN SI UNITS"))
println(YELLOW_FG(" UPSTREAM PROPERTIES"));
df1 = DataFrames.DataFrame(PRESSURE = P1, DENSITY = ρ1,TEMPERATURE = T1, ENTHALPY = h1, VELOCITY = u1);
println(df1);
println(YELLOW_FG("DOWNSTREAM PROPERTIES"));
df2 = DataFrames.DataFrame(PRESSURE = P2, DENSITY = ρ2,TEMPERATURE = T2, ENTHALPY = h2, VELOCITY = u2);
println(df2);
#Plot graphs
display(plot(u1 ./ 1000,T2./T1, seriestype = :scatter,grid = true, title = "Temperature Ratio "));
#savefig("<Enter location where you want to save>\\Temperature Ratio");
display(plot(u1 ./ 1000,P2./P1,seriestype = :scatter,grid = true, title = "Pressure Ratio "));
#savefig("<Enter location where you want to save>\\Dissociaton Plot");
display(plot(u1 ./ 1000,ρ2./ρ1 ,seriestype = :scatter, grid = true, title = "Density Ratio"));
#savefig("<Enter location where you want to save>\\Density Ratio");
display(plot(u1 ./ 1000,h2./h1,seriestype = :scatter,grid = true, title = "Enthaply Ratio "));
#savefig("<Enter location where you want to save>\\Enthalpy ratio");
display(plot(u1 ./ 1000,u2,seriestype = :scatter,grid = true, title = "Downstream Velocity  "));
#savefig("<Enter location where you want to save>\\Post Shock Velcoity");
# -
