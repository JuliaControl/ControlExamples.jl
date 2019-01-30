using ControlSystems, OrdinaryDiffEq

"""
    GainSchedulingSimulator
# Fields:
    s::Simulator
    f = (x,p,t) -> x
    y = (x,t) -> y
    r = (x,t) -> reference
    e = (x,t) -> r(x,t) .- x # Error
    controllers::Vector{StateSpace}
    conditions::Vector{Function}
"""
struct GainSchedulingSimulator <: ControlSystems.AbstractSimulator
    s::Simulator
    f
    y
    r
    e
    controllers::Vector{StateSpace}
    conditions::Vector{T} where T <: Function
    feedforward::Vector{Float64}
end
"""
    GainSchedulingSimulator(P,r,controllers::AbstractVector{LTISystem},conditions::AbstractVector{Function (x,y,r)->Bool}; inputfun=(u,t)->u, feedforward=zeros())
`r` is a function (x,t)-> reference signal
`controllers` is a list of LTIsystems representing the controllers
`conditions` is a list of function that accepts state, output and reference and returns a Bool indicating which controller is active. At least one function should return true for every input. If more than one condition is true, the index of the first condition in the list is chosen.
´feedforward´ is a vector of constants that are added to the control signal calculated by the active controller, i.e., `u = C*x + D*e + feedforward[active_index]`
For a usage example, see notebook in example folder.
```
"""
function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc};
    inputfun=(u,t)->u,
    feedforward=zeros(length(controllers))) where Tu <: StateSpace where Tc <: Function
    s = Simulator(P)
    pinds = 1:P.nx # Indices of plant-state derivative
    r = (x,t) -> ri(x[pinds],t)
    y(x,t) = s.y(x[pinds],t)
    y(sol::ODESolution,t) = P.C*sol(t)[pinds,:]
    e = (x,t) -> r(x,t) .- y(x,t)
    f = function(der,x,p,t)
        xyr = (x[pinds],y(x,t),r(x,t))
        index = findfirst(c->c(xyr...), conditions)
        @assert index > 0 "No condition returned true"
        et = e(x,t)
        ind = P.nx+1 # First index of currently updated controller
        for c in controllers
            c.nx == 0 && continue # Gain only, no states
            inds = ind:(ind+c.nx-1) # Controller indices
            der[inds] = c.A*x[inds] + c.B*et # Controller dynamics, driven by error
            ind += c.nx
        end
        c = controllers[index] # Active controller
        cind = P.nx + (index == 1 ? 1 : 1+sum(i->controllers[i].nx, 1:index-1)) # Get index of first controller state
        der[pinds] = P.A*x[pinds] # System dynamics
        u = c.C*x[cind:(cind+c.nx-1)] .+ c.D*et .+ feedforward[index] # Form control signal
        size(inputfun(P.B*u,t)), size(der[pinds])
        der[pinds] .+= vec(inputfun(P.B*u,t)) # Add input from active controller to system dynamics
        der
    end
    GainSchedulingSimulator(s,f,y,r,e,controllers,conditions,feedforward)
end

function GainSchedulingSimulator(P,ri,controllers::AbstractVector{Tu},conditions::AbstractVector{Tc};
    inputfun=(u,t)->u,
    feedforward=zeros(length(controllers))) where Tu <: LTISystem where Tc <: Function
    GainSchedulingSimulator(P,ri,ss.(controllers),conditions; inputfun=inputfun, feedforward=feedforward)
end
# ============================================================================================


"""
    sol = solve(s::AbstractSimulator, x0, tspan,  args...; kwargs...)
Simulate the system represented by `s` from initial state `x0` over time span `tspan = (t0,tf)`.
`args` and `kwargs` are sent to the `solve` function from `OrdinaryDiffEq`, e.g., `solve(s, x0, tspan,  Tsit5(), reltol=1e-5)` solves the problem with solver [`Tsit5()`](http://docs.juliadiffeq.org/stable/solvers/ode_solve.html) and relative tolerance 1e-5.

See also `Simulator` `lsim`
"""

DiffEqBase.solve(s::GainSchedulingSimulator, x0, tspan, solver=Tsit5(), args...; kwargs...) = solve(ODEProblem(s.f,vcat(x0, zeros(sum(c->c.nx, s.controllers))),tspan), solver, args...; kwargs...)
