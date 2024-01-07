using DifferentialEquations
using Plots
using LaTeXStrings



Notes = "N25_K_red_larger_omega_on"
#  Try my stuff
ρ = 0.25
γ = 2
ξ = 0.025
K = 0.003


λ = 1 + ρ*exp(γ)

# parameters for plotting comparison with PDE
ω_on = 0.003


cd("/Users/mbcx4dh9/Documents/Woolner,S_Jensen,O/PhD/Julia_code/ODE solver")
mkpath("ODE_Data_$Notes")
cd("ODE_Data_$Notes")

touch("run_parameters.txt")
save_params = open("run_parameters.txt","w")
write(save_params, "ρ γ ξ λ K \n")
write(save_params, string(ρ, " ", γ," ", ξ," ", λ," ", K))
close(save_params)

function model!(du, u, p, t)
    # u in order bu, bl, z
    ρ, γ, ξ, K = p
    du[3]=(u[1]-u[2]-K*u[3])/(ξ+u[1]+u[2]) #dz/dt
    Fu=ρ*exp(γ*(1-du[3]))
    Fl=ρ*exp(γ*(1+du[3]))
#

    du[1] = 1/2-1/2*u[1]*(1+Fu) # db+/dt
    du[2] = 1/2-1/2*u[2]*(1+Fl) # db-/dt
end




u0 = [0.3, 0.3, 0.05] # initial conditions, b+, b-, z
p=(ρ,γ,ξ,K)


EndTime = 100000.0
# Time
tspan = (0.0, EndTime*ω_on)                # time start and end
prob=ODEProblem(model!,u0,tspan,p) # set up problem for ODE solver
sol = solve(prob)  # solutions saved in sol

# sol can be plotted


save_points = LinRange(0, EndTime*ω_on, Int(EndTime*100))
A = Array(sol(save_points)) # dense solutions saved therefore can access intermediate timepoints


Dz_Dt = (A[1,:] .- A[2,:] .- K.*A[3,:])./(ξ .+ A[1,:] .+ A[2,:])

plot(save_points/ω_on, A[3,:]/ω_on, linewidth = 2, c = :black)
plot!(legend = false, grid = false, dpi=300, xguidefontsize = 18, yguidefontsize = 18, thickness_scaling = 2)
xlabel!(L"\bar{t}")
ylabel!(L"\bar{z}")
savefig("Pole position")

# other plots can be constructed with the solutions held in A as required
