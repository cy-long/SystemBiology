using Distributions, Plots, FFTW

# Initialization
d1 = 0.1; d2 = 0.05; b = 0.1; p1 = 0.25; p2 = 0.05; N = 3000;
r = 2b - d2;
K = r/(2b);
p = p1 + p2 + b;

function g1(n, m)
    d1*n 
end

function f2(n, m)
    2b*(m/N)*(N-n-m)
end

function g2(n, m)
    2p2*(n*m)/N + d2*m
end

function f1(n, m)
    2p1*(m*n)/N
end


# Gillespie Simulation
n = Float64[0]; 
m = Float64[0];
t = Float64[0];
Nsample = 10^6;

i = 1; 
n[i] = 100; m[i] = 200;


for run in 1:Nsample
    u = rand(Uniform(0,1), 4)
    t1 = (1/g1(n[i], m[i])) * log(1/u[1])
    t2 = (1/f2(n[i], m[i])) * log(1/u[2])
    t3 = (1/g2(n[i], m[i])) * log(1/u[3])
    t4 = (1/f1(n[i], m[i])) * log(1/u[4])
    tt = minimum([t1, t2, t3, t4])

    if(t1 == tt)
        push!(n, n[i] - 1)
        push!(m, m[i])
    elseif (t2 == tt)
        push!(n, n[i])
        push!(m, m[i] + 1)
    elseif (t3 == tt)
        push!(n, n[i])
        push!(m, m[i] - 1)
    elseif (t4 == tt)
        push!(n, n[i] + 1)
        push!(m, m[i] - 1)
    end
    push!(t, t[i] + tt)
    i += 1
end

# Spectrum analysis
k = abs.(fft(n[1:10000]))
plot(k[1:10])
plot(n)


# Visulization
# p1 = plot(t, [n m], title = "Indi-Based Model",
    # lw =2, legend = :topright, label=["predator" "prey"],
    # xlabel = "Time", ylabel = "Abundances")

# p1 = plot(t[500000: 1000000], n[500000: 1000000], title = "Indi-Based Model")



# Amplitude - N relationship
N = [300, 3000, 30000]
A1 = sqrt(var(n[Integer(Nsample/2):Nsample]))
#similarly A2 A3
# A = [A1, A2, A3]

# plot(N, A, xaxis=:log, yaxis=:log, legend = :topleft, xlabel = "N", ylabel = "Amplitude", label = "slope = 0.476")
# plot!(N, A, seriestype=:scatter, label = "")

# lN = log.(N)
# lA = log.(A)

# m = (3*sum(lN .* lA) - sum(lN)sum(lA))/(3*sum(lN .* lN) - sum(lN)^2)
# b = (sum(lA) - m*sum(lN))/3