
import Plots
import LinearAlgebra
using Plots
using LinearAlgebra

#data example
x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
y = [0.6761488864859304; 0.6345697680852508; 0.6396283580587062; 0.6132010027973919; 0.5906142598705267; 0.5718728461471725; 0.5524549902830562;
    0.538938885654085; 0.5373495476994958; 0.514904589752926; 0.49243437874655027]

function jacobnewton(gc)
    ab = gc
    global x
    global y
    iter = 0
    F = zeros(2)
    J = zeros(2,2)
    while iter <= 100
        C = 0
        D = 0
        E = 0 
        G = 0
        H = 0
        I = 0
        for i = 1:11
            C = C + 2*((ab[1]/(ab[2]+x[i])) - y[i])*(1/(ab[2]+x[i]))
            D = D + 2*((ab[1]/(ab[2]+x[i])) - y[i])*(-ab[1]/((ab[2]+x[i])^2))
            E = E + 2*(1/((ab[2]+x[i]))^2)
            G = G + 2*((y[i]*ab[2] + x[i]*y[i] - 2*ab[1]))/((ab[2]+x[i])^3)
            H = H + -2*(2*ab[1]+y[i]*(-x[i] - ab[2]))/((ab[2]+x[i])^3)
            I = I + -2*ab[1]*(2*y[i]*ab[2] + 2*x[i]*y[i]-3*ab[1])/((ab[2]+x[i])^4)
        end
        F[1] = C
        F[2] = D
        J[1,1], J[1,2], J[2,1], J[2,2] = E, G, H, I
        iter +=1
        ab = ab - LinearAlgebra.inv(J)*F
    end
    return ab
end

function set()
    gc = [1;1]
    nt = jacobnewton(gc)
    global x
    s = []
    for i = 1:11
        r = nt[1] / (nt[2]+x[i])
        append!(s, r)
    end
    return s
end

Plots.plot!(x , set(), label = "Approximation")
Plots.scatter!(x , y, label = "Data")