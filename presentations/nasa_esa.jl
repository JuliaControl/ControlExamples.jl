### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8fe0da54-e41f-45a2-9d79-3df73932ee40
begin
	using Revise, Pkg
	Pkg.activate()
	using ControlSystems, Plots, PlutoUI, Random, LinearAlgebra
end

# ╔═╡ 0a084289-12b1-49b8-bed7-e59b192b5c0a
let
	using Symbolics
	@variables δ k2 γ
	s  = tf("s")
	Gc = -k2 * (s + δ) / (s + γ)
	H  = tf(1, [1, 1])
	sys_cl = feedback(H, Gc)
end

# ╔═╡ 78c746e9-1348-47ae-83f7-d3fb1dfcc117
begin
	using Optim, Statistics
	using Plots.Measures
end

# ╔═╡ 4e5b9b1a-9e06-490b-932a-b311354962bc
begin
	using MonteCarloMeasurements: Particles, StaticParticles, pmean, pmaximum, unsafe_comparisons, Uniform
	unsafe_comparisons(true, verbose=false)
end;

# ╔═╡ 7eaa0e06-f1d1-4e6d-9080-c33603b4d62f
# ╠═╡ disabled = true
#=╠═╡
begin
	using ControlSystemIdentification
	w = 2pi .* exp10.(LinRange(-3, log10(0.5), 500)) # Frequency vector for plots
	G0 = tf(1, [10, 1]) # The true system, 10ẋ = -x + u
	Gtemp = c2d(G0, 1)  # discretize with a sample time of 1s
	G0
end
  ╠═╡ =#

# ╔═╡ 1569ee96-9567-11ec-2f28-419213ea3fab
html"<button onclick='present()'>present</button>"

# ╔═╡ 37d05643-4906-4df7-9905-bd6a88f82a21
md"""
# JuliaControl
## Outline

- Who am I?
- What does JuliaControl offer?
    - [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
    - [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl)
    - [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl)
    - [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)
    - [ModelingToolkit.jl](https://mtk.sciml.ai/dev/)
- Interaction with the Julia ecosystem

![logo](https://avatars.githubusercontent.com/u/10605979?s=400&u=7b2efdd404c4db3b3f067f04c305d40c025a8961&v=4)
"""

# ╔═╡ 7dbef7ff-437d-4f1f-ad3e-c03f1394107e
md"""
# Who am I?

My name is **Fredrik Bagge Carlson**, I live in Lund in southern **Sweden**.

I got my MSc and PhD from Dept. Automatic Control at **Lund University**.

My background is within **robotics, control and system identification**.

I enjoy writing **software toolboxes** 😄 [baggepinnen@github](https://github.com/baggepinnen)

I work with simulation and control at **Julia Computing** (I represent myself here, all opinions are my own).
"""

# ╔═╡ 484dc785-59b3-4e5d-918f-6b71aa2d057c
md"""
# What is JuliaControl?
[ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) is a toolbox for analyzing and designing control systems
- Similar in spirit and scope to Control systems toolbox
- Linear systems
- Transfer functions, statespace systems, frequency response, LQR, LQG, Kalman filter, Bode and Nyquist plots, root locus, PID controllers
- [Paper describing the toolbox](https://portal.research.lu.se/en/publications/controlsystemsjl-a-control-toolbox-in-julia)
> FBC., Fält, Heimerson & Troeng, (2021). ControlSystems.jl: A Control Toolbox in Julia. 60th IEEE Conference on Decision and Control (CDC)https://doi.org/10.1109/CDC45484.2021.9683403

[RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl) contains facilities for robust analysis and design
- Named statespace systems where states, inputs and outputs are accessible by names rather than indices. This also facilitates creating complicated feedback interconnections.
- An interface to [DescriptorSystems.jl](https://github.com/andreasvarga/DescriptorSystems.jl).
- Robust/optimal design methods such as $H_{\infty}$, $H_{2}$, LQG and Glover-McFarlane.
- Robustness-related metrics such as $\nu$-gap, `ncfmargin`, `diskmargin` etc.
- Uncertainty modeling with the $M\Delta$ framework (and more).
- Model augmentation.

[ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl) contains LTI identification
- Time domain
- Frequency domain
- Subspace-based methods (N4SID etc.)
- Prediction-error method (PEM)
[SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl) contains basic C-code generation.

##### JuliaControl is currently *not*
- Open-loop optimal control
- Nonlinear control (except basic simulation)

##### Main Contributors
![Mattias](https://portal.research.lu.se/files-asset/7619687/_MPP9977.jpg?w=160&f=webp)
![Albin](https://portal.research.lu.se/files-asset/110428020/albinh2_cropped_lowres.jpg?w=160&f=webp)
![Olof](https://portal.research.lu.se/files-asset/7213474/OlofTroeng.jpg?w=160&f=webp)
![Fredrik](https://portal.research.lu.se/files-asset/114575837/bagge.jpg?w=160&f=webp)

Started by Jim Christ in 2016
"""

# ╔═╡ a76b10ce-e373-4361-ba27-9995b5c3f730
md"""
## Look and feel
"""

# ╔═╡ fbc35df0-c8ae-4704-996e-9c5dd4b27bd2
md"#### Transfer functions:"

# ╔═╡ 3cd2a0cf-c2e5-494d-abe3-9474dda7fc84
tf(1, [1,2,1])

# ╔═╡ 02e3f332-c2fd-4bce-93a1-373e3b8eb125
md"#### Discrete-time systems:"

# ╔═╡ 4b4ee402-19e4-44a4-8423-dd3835aba538
c2d(tf(1, [1,2,1]), 0.1)

# ╔═╡ 2b543ce0-2acf-4c53-b9e9-fd310f96ffed
md"#### Zero-pole-gain representation:"

# ╔═╡ f488eb97-17a1-44da-8b7a-cc0f236404ac
zpk(tf(1, [1,2,1]))

# ╔═╡ be97646c-68e2-46bb-b66b-ab3a7a432196
md"#### State-space systems:"

# ╔═╡ 8fe8e18a-7975-44d5-81d1-c3ed35e27a7e
let
	A = [1 0.1; 0 1]
	B = [0; 1]
	C = [1 0]
	ss(A,B,C,0, 0.1)
end

# ╔═╡ d1fa95a9-92d6-441f-9b00-2f14cc5e96e0
md"#### Support for time delays:"

# ╔═╡ 717b3cfb-4ee9-409e-a8fe-74e12adf65c6
Pdelay = tf(1, [1,2,1]) * delay(1.5)

# ╔═╡ 468ef1b4-f0d2-41a9-81ef-1b747a52caaf
plot(step(Pdelay, 12))

# ╔═╡ 458fe3d3-671c-40ed-9a52-c2d9e7a945b1
md"#### Support for nonlinearities:"

# ╔═╡ aa01b91c-d77a-4f23-8f74-a6f543b1f7d5
md"""
saturation: = $(@bind sat Slider(0.5:0.5:5, default=1, show_value=true))
"""

# ╔═╡ 9b9a1365-59d9-490f-bd9c-78c8f4d73957
let
	using ControlSystems: saturation
	P  = tf(1, [1,2,1])
	C  = saturation(sat) * tf(10) 					# Saturated P controller
	T  = feedback(P*C)
	CS = feedback(C, P)
	stepsim = step([T; CS], 5) 						# Step-response simulation
	plot(stepsim, ylim=[(0, 1.3) (-5,5)], ylabel=["y" "u"])
end

# ╔═╡ 232de9ed-9ca9-4e4e-9344-c9306755aec5
md"#### Basic design:"

# ╔═╡ 6c90abe5-e0b7-444c-8282-7ac52592964d
md"""
``T_s`` = $(@bind Ts Slider(exp10.(-2:0.1:0), default=0.1, show_value=true))
``R`` = $(@bind r Slider(exp10.(-3:0.2:3), default=10, show_value=true))
"""

# ╔═╡ 798829dd-af76-48df-891b-fe273aa7dd6d
let
	A   = [1 Ts; 0 1]
	B   = [0; 1]
	C   = [1 0]
	sys = ss(A,B,C,0,Ts)
	Q   = I
	R   = Ts*r*I
	L   = lqr(Discrete,A,B,Q,R) # lqr(sys,Q,R) can also be used
	
	u(x,t) = -L*x .+ 1.5(t ≥ 2.5)   # Form control law (u is a function of t and x), a constant input disturbance is affecting the system from t ≧ 2.5
	t      = 0:Ts:8
	x0     = [1, 0]
	res    = lsim(sys,u,t,x0=x0)
	plot(res, lab=["Input + disturbance" "Position" "Velocity"], xlabel="Time [s]", plotx = true, plotu=true, ploty=false)
end

# ╔═╡ cca92baf-988a-4f17-a472-9e0870d4239d
md"""
## Do we offer something unique?
System types are generic w.r.t. element types
- Systems with symbolic coefficients, uncertain parameters, `BigFloat` etc.
- Delay systems
- GPU support
- Automatic differentiation
"""

# ╔═╡ 33e2e0a2-653c-46f2-a8d7-3dc98dfd0f20
md"#### Symbolics"

# ╔═╡ 9d90fae3-5cfc-4ff1-a6aa-112e25c98ef0
md""" #### GPU
```julia
using CUDA
A   = cu([-1.2 0; 0.1 -0.4]);  B = cu([1.2; 0])
C   = cu([1 0]);               D = cu([0])
sys = HeteroStateSpace(A, B, C, D)
```
"""

# ╔═╡ 3fe2ae0a-256a-468e-bd59-e80cd00f9c8e
md"#### Delay systems and High Precision"

# ╔═╡ c8321df5-a9d5-49ee-b1c5-f1ecbca20663
let
	s = tf("s")
	sys = feedback(BigFloat(1)/s, delay(1))
	impresp = impulse(sys, BigFloat.(0:0.1:10))
	@show typeof(impresp)
	plot(impresp)
end

# ╔═╡ da1cdab6-c7c0-4b20-af71-c4bb6f258e3d
md"#### Automatic differentiation"

# ╔═╡ 5a8fcedb-8c04-4ba3-b0ff-107e7e763d42
function cost_function(params)
	gain, pole = params
	P = tf(1, [1,1])
	C = tf(gain, [1, pole])
	closed_loop = feedback(P*C)
	y = step(closed_loop, 0:0.1:5).y
	sum(abs2, y .- 1)
end

# ╔═╡ caa06a90-0434-4d79-ae19-08c06254762d
let
	using ForwardDiff
	params = [1.0, 1.0]
	∇cost = ForwardDiff.gradient(cost_function, params)
end

# ╔═╡ 42b8af96-db21-4f2c-bb79-5a8570a710fa
md"""
## Example: Smith predictor for delay system

This example designs a controller for a plant with a time delay using a Smith predictor. The plant is given by
```math
\dfrac{1}{s + 1}e^{-s\tau} = P_0 e^{-s\tau}
```

and the control architecture looks like this
```
                ┌──────┐              ┌─────────────┐
r               │      │          u   │             │
───⊕──⊕────────►│  C0  ├───────────┬─►│ P0*exp(-sτ) ├─┐y
   ▲  ▲         │      │           │  │             │ │
  -│  │         └──────┘           │  └─────────────┘ │
   │  │                            │                  │
   │  │ ┌──────────┐    ┌──────┐   │                  │
   │  │ │          │    │      │   │                  │
   │  └─┤exp(-sτ)-1│◄───┤  P0  │◄──┘                  │
   │    │          │    │      │                      │
   │    └──────────┘    └──────┘                      │
   │                                                  │
   └──────────────────────────────────────────────────┘
```
The benefit of this approach is that the controller $C_0$ can be designed for the nominal plant $P_0$ without time delay, and still behave well in the presence of the delay. 

We now set up the nominal system and PI controller
"""

# ╔═╡ f4c42628-42a8-45a1-9e8c-3637a9433d98
md"τ = $(@bind τ Slider(0.2:0.2:8, default=8, show_value=true))"

# ╔═╡ 14bfe60e-fe07-4774-9055-4b76806ca4d8
let
	## Set up the nominal system and PI controller
	P0 = ss(-1, 1, 1, 0)
	
	# PI controller for nominal system P0
	# To verify the pole placement, use, e.g., dampreport(feedback(P0, C0))
	ω0 = 2
	ζ = 0.7
	_, C0 = placePI(P0, ω0, ζ)
	
	## Setup delayed plant + Smith predictor-based
	#  controller for a given delay τ
	P = delay(τ) * P0
	C_sp = feedback(C0, (1.0 - delay(τ))*P0)
	
	## Plot the closed loop response 
	# Reference step at t = 0 and load disturbance step at t = 15
	G = [feedback(P*C_sp, 1) feedback(P, C_sp)]
	fig_timeresp = plot(lsim(G, t -> [1; t >= 15], 0:0.1:40),  title="τ = $τ", ylims=(-0.1, 2.1))
	
	## Plot the frequency response of the predictor part and compare to a negative delay
	C_pred = feedback(1, C0*(ss(1.0) - delay(τ))*P0)
	fig_bode = bodeplot([C_pred, delay(-τ)], exp10.(-1:0.001:1), ls=[:solid :solid :dash :dash], title="", ylims=[(1e-2, 1e2) (0, 900)])
	plot!(yticks=[0.1, 1, 10], sp=1)
	plot!(yticks=0:180:1080, sp=2)
	
	## Check the Nyquist plot
	# Note that the Nyquist curve encircles -1 for τ > 2.99
	fig_nyquist = nyquistplot(P * C_sp, exp10.(-1.5:3e-4:2), title="τ = $τ", xlims=(-20, 10), ylims=(-10,10))
	plot(fig_timeresp, fig_bode, fig_nyquist, layout=(1,3), size=(1200, 500), margin=4mm)
end

# ╔═╡ 047be3bb-357b-4a2e-abac-c39b797ae4ab
md"""
## Example: Controller optimization with plant uncertainty
We will demonstrate an example of controller optimization with a plant model containing **uncertain parameters**.
We start by defining a nominal model, visalize it's Bode diagram and perform a time-domain simulation. The model we'll consider is on the form
$Y(s) = P(s)U(s)$ and of second order. This kind of model may represent something like a moving mass attached to a flexible element.

$P(s) = \dfrac{p\omega}{s^2 + 2\zeta\omega s + \omega^2}$

The controller will be a standard PID controller with a second-order filter on the form

$C(s) = \left(k_P + \dfrac{k_I}{s} + k_D s \right) \dfrac{1}{(0.05s + 1)^2}$

### Initial design, no uncertainty
"""

# ╔═╡ 77e162be-11a9-46f9-a7bd-c9f8c49f0313
begin
	p = 1  
	ζ = 0.3
	ω = 1  
	P = ss(tf([p*ω], [1, 2ζ*ω, ω^2]))
end

# ╔═╡ a1070752-64fb-427f-962d-721db16a44b9
begin
	Ω = exp10.(-2:0.04:3)
	kp,ki,kd =  1, 0.1, 0.1 # controller parameters
	
	C  = ss(pid(; kp, ki, kd)*tf(1, [0.05, 1])^2) # Construct a PID controller with filter
	G  = feedback(P*C) # Closed-loop system
	S  = 1/(1 + P*C)   # Sensitivity function
	Gd = c2d(G, 0.1)   # Discretize the system
	res = step(Gd,15)  # Step-response
	
	mag = bode(S, Ω)[1][:]
	plot(res, title="Time response", layout = (1,3), legend=:bottomright)
	plot!(Ω, mag, title="Sensitivity function", xscale=:log10, yscale=:log10, subplot=2, legend=:bottomright, ylims=(3e-2, Inf))
	nyquistplot!(P*C, Ω, sp=3, ylims=(-2.1,1.1), xlims=(-2.1,1.2), size=(1200,400))
end

# ╔═╡ 3ad0c65d-2431-45a8-a022-9c49cbee996c
md"
#### Design by optimization
Next, we try to tune the controller using optimization. We therefore package the creation of the systems above in a function that takes in controller parameters and outputs a cost.

In order to promote robustness, we place a constraint on the maximum magnitude $M_S$ of the sensitivity function $S$.
"

# ╔═╡ 717e664b-9d46-4729-9aa6-0e6c32471111
begin
	const Msc = 1.2 # Constraint on Ms
	
	function systems(P, params)
	    kp,ki,kd = exp.(params)
	    C   = ss(pid(; kp, ki, kd)*tf(1, [0.05, 1])^2)
	    G   = feedback(P*C) # Closed-loop system
	    S   = 1/(1 + P*C)   # Sensitivity function
	    CS  = C*S           # Noise amplification
	    Gd  = c2d(G, 0.1)   # Discretize the system
	    res = step(Gd,15)  # Step-response
	    C, G, S, CS, res
	end
	
	function cost(P, params)
	    C, G, S, CS, res = systems(P, params)
	    Ms     = maximum(bode(S, Ω)[1]) # max sensitivity
	    perf   = mean(abs, 1 .- res.y)
	    robust = (Ms > Msc ? 10000(Ms-Msc) : zero(eltype(params)))    
	    noise  = sum(bode(CS, Ω[end-30:end])[1])
	    100perf + robust + 0.002noise
	end
	
	params  = log.([1,0.1,0.1]) # Initial guess (optimize in log space to force positive parameters)
	optres = Optim.optimize(p->cost(P, p), params, Optim.Options(
	    show_trace = true,
	    show_every = 50,
	))

	function plot_optimized(P, labparams...)
	    fig = plot(layout=(1,3), size=(1200,400), bottommargin=2mm)
	    for (i,(lab, params)) = enumerate(labparams)
	        C, G, S, CS, r = systems(P, params)
	        mag = bode(S, Ω)[1][:]
	        plot!(r; title="Time response", subplot=1,
				lab, legend=:bottomright,
				fillalpha=0.05, linealpha=0.8, seriestype=:path)
	        plot!(Ω, mag; title="Sensitivity function",
				xscale=:log10, yscale=:log10, subplot=2,
				lab, legend=:bottomright, fillalpha=0.05, linealpha=0.8)
			nyquistplot!(P*C, Ω; Ms_circles=Msc, sp=3, ylims=(-2.1,1.1), xlims=(-2.1,1.2), lab, points=true, seriescolor=i)
	    end
	    hline!([Msc], l=(:black, :dash), subplot=2, lab="", ylims=(9e-2, Inf))
	    fig
	end
	
	## We can now perform the same computations as above to visualize the found controller	
	plot_optimized(P, "Optimized"=>optres.minimizer, "Initial"=>params)
end

# ╔═╡ a16cbfbc-3ea3-4891-bf3a-65304fca3fc1
md"""
### Design with uncertainty

The next step is to add uncertainty to the system. Lets say all the parameters $(p, ζ, ω)$ are associated with a Gaussian uncertainty. We can create such parameters using [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl)
"""

# ╔═╡ 74859413-5696-4ed3-8015-9794f6bbff98
begin
	±(m,s) = m + s*Particles(50)
	..(a,b) = Particles(50, Uniform(a,b))
	
	pᵤ = 1 ± 0.3
	ζᵤ = 0.1..0.35
	ωᵤ = 1 ± 0.05
	Pᵤ = ss(tf([pᵤ*ωᵤ], [1, 2ζᵤ*ωᵤ, ωᵤ^2]))
end

# ╔═╡ 557d94ce-89a8-4f92-9e17-cbccb0ecb1d6
plot_optimized(Pᵤ, "Optimized"=>optres.minimizer)

# ╔═╡ 72b702ac-7d7c-4d4d-93b4-9af4657cda5e
md"If we now visualize the found controller, we find that it might violate the constraint we placed on $M_s$"

# ╔═╡ 465b2854-5e02-4389-b917-f3d9ecb619bc
md"""
#### Robust optimization
To alleviate this, we perform the optimization again, this time with the uncertain system.
In order to do this, we must change our cost function slightly so that it outputs a scalar instead of an uncertain numer
"""

# ╔═╡ 06414cfb-b613-4f5f-a0e2-8dec15e8c6e4
cost(Pᵤ, params)

# ╔═╡ 9834d057-887e-4875-b023-bf763296fe0e
function cost_uncertain(P, params)
    C, G, S, CS, res = systems(P, params)
    Ms = maximum(bode(S, Ω)[1])     |> pmaximum# max sensitivity
    perf = mean(abs, 1 .- res.y)    |> pmean
    robust = (Ms > Msc ? 10000(Ms-Msc) : zero(eltype(params)))    
    noise = sum(bode(CS, Ω[end-30:end])[1]) |> pmean
    100perf + robust + 0.002noise
end

# ╔═╡ ef6804c9-58a2-496b-876e-253618678af4
cost_uncertain(Pᵤ, params)

# ╔═╡ c2554c31-d09a-40f5-aec1-95e6d8ab145c
md"We can now perform the optimization and visualization again"

# ╔═╡ e8e80d7d-3476-4825-9799-9a083bad8aeb
# ╠═╡ disabled = true
#=╠═╡
begin
	optresu = optimize(p->cost_uncertain(Pᵤ, p), optres.minimizer, Optim.Options(
	    show_trace        = true,
	    show_every        = 50,
	))
	plot_optimized(Pᵤ, "Optimized"=>optresu)
end
  ╠═╡ =#

# ╔═╡ 4dbad986-5975-40f2-aad0-8795e53d0a92
md"""
This time, all realizations are satisfying the constraints! The controller achieved this by being less aggressive. The controller parameters ($k_P,k_I,k_D$) are given below
"""

# ╔═╡ c7d9f4f3-e01a-425f-9037-4807cea0ab08
#=╠═╡
(nominal=exp.(optres.minimizer), robust=exp.(optresu.minimizer))
  ╠═╡ =#

# ╔═╡ 6476dc63-52e0-4304-a849-190e82f55a8f
md"""
# [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl)
An extension to ControlSystems.jl containing
- Uncertainty modeling and $\mu$-analysis
- Named systems
- Extra LQG functionality
-  $\mathcal{H}_\infty$ and $\mathcal{H}_2$ optimal design.
"""

# ╔═╡ dc4b13f1-9cc8-49e7-a30c-ba5903de2dc9
md"""
# System identification
The art and science of estimating models of dynamical systems.

[ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl)
- Black-box linear models.
- Time and frequency domain.

## Example: temperature model
A typical model for a temperature-controlled system is 
```math
\tau \dot T = -T  + Bu + c
```
where $T$ is the temperature, $u$ the control signal and $c$ a constant offset, e.g., related to the temperature surrounding the controlled system. The time constant $\tau$ captures the relation between stored energy and the resistance to heat flow and determines how fast the temperature is changing. This system can be written on transfer-function form like (omitting $c$)
```math
Y(s) = \dfrac{B}{\tau s + 1}U(s)
```
This is a simple first-order transfer function which can be estimated with, e.g., the functions [`subspaceid`](https://baggepinnen.github.io/ControlSystemIdentification.jl/dev/ss/#ControlSystemIdentification.subspaceid) or [`plr`](https://baggepinnen.github.io/ControlSystemIdentification.jl/dev/tf/#ControlSystemIdentification.plr). To illustrate this, we create such a system and simulate some data from it.
"""

# ╔═╡ 54b9a02c-d707-42af-9661-f26a57640009
# ╠═╡ disabled = true
#=╠═╡
begin
	u = sign.(sin.((0:0.01:20) .^ 2))' # sample a control input for identification
	y, t, x = lsim(ss(Gtemp), u) # Simulate the true system to get test data
	yn = y .+ 0.2 .* randn.() # add measurement noise
	data = iddata(yn, u, t[2] - t[1]) # create a data object
	plot(data)
end
  ╠═╡ =#

# ╔═╡ 738e2006-b560-4481-9a50-679f6026ccd6
md"""
We see that the data we're going to use for identification is a chirp input. Chirps are excellent for identification as they have a well defined and easily controllable interval of frequencies for identification. We start by inspecting the coherence plot to ensure that the data is suitable for identification of a linear system
"""

# ╔═╡ e6a9a177-8e4f-4e1e-b31e-d5d2ff482700
#=╠═╡
coherenceplot(data, hz=true)
  ╠═╡ =#

# ╔═╡ eb5cea54-1efb-451b-8c6f-76dc3debe99d
md"""
The coherence is high for all frequencies spanned by the chirp, after which it drops significantly. This implies that we can only ever trust the identified model to be accurate up to the highest frequency that was present in the chirp input.

Next we set the parameters for the estimation, the numerator and denominator have one parameter each, so we set $n_a = n_b = 1$ and estimate two models.
"""

# ╔═╡ dd4e6b96-fd36-41ff-8e92-be8eda68162d
#=╠═╡
begin
	na, nb = 1, 1 # system order and number of parameters in the numerator
	Gh = subspaceid(data, na, r=6, zeroD=true).sys
	Gh2, noise_model = plr(data, na, nb, 1) # try another identification method
	Gh, Gh2
end
  ╠═╡ =#

# ╔═╡ 88382bc4-53e0-4211-a8e7-aec464cae474
md"""
We can plot the results in several different ways:
"""

# ╔═╡ 759f6471-908e-4e08-9e0b-a56c4d0ece96
#=╠═╡
d2c(Gh) # Transform to continuous time
  ╠═╡ =#

# ╔═╡ 2b05b91e-abd1-44d1-8016-32d78f02b1c7
#=╠═╡
begin
	bp = bodeplot(Gtemp, w, lab = "G (true)", hz = true, l = 3)
	bodeplot!(Gh, w, lab = "subspace", hz = true)
	bodeplot!(Gh2, w, lab = "plr", hz = true, ticks = :default)
	
	sp = plot(step(Gtemp, 150), lab="G (true)")
	plot!(step(Gh, 150), lab = "subspace")
	plot!(step(Gh2, 150), lab = "plr", ticks = :default)
	hline!([1], primary = false, l = (:black, :dash))
	
	lp = plot(lsim(ss(Gtemp), u), lab="G (true)")
	plot!(lsim(ss(Gh), u), lab = "subspace")
	plot!(lsim(ss(Gh2), u), lab = "plr", ticks = :default)
	plot!(data.t, yn[:], lab = "Estimation data", alpha=0.5)
	
	plot(bp, sp, lp, layout = @layout([[a b]; c]), size=(1000,800))
end
  ╠═╡ =#

# ╔═╡ 62219eaf-6b01-4e06-ba7b-5247b081ad49
md"""
# [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)
- Tools for working with systems with symbolic coefficients.
- C-code generation of controllers and filters.
"""

# ╔═╡ d35f2cf1-8111-4bc1-9e77-fdd0db770bca
md"""
# Summary
- [JuliaControl](https://github.com/JuliaControl/) tries to offer a fully-featured environment for analysis, design and estimation of linear control systems.
- Don't hesitate to reach out!
    - [Issue](https://github.com/JuliaControl/ControlSystems.jl/issues)
    - [Discourse](https://discourse.julialang.org/tag/control)
    - [`#control`](https://julialang.slack.com/archives/CKH1UTZT9) on slack
    - PRs are welcome, including tutorials and examples.
- If you use the software, be sure to say hi and mention what you like and what can be improved.

##### Thank you 😃
"""

# ╔═╡ 22a258e7-71a0-433b-b4b9-cbfd8a3d4a84
html"""<style>
main {
    max-width: 60%;
}
"""

# ╔═╡ 61341695-7a95-4164-a519-7fe347a56738
PlutoUI.TableOfContents()

# ╔═╡ Cell order:
# ╟─1569ee96-9567-11ec-2f28-419213ea3fab
# ╠═8fe0da54-e41f-45a2-9d79-3df73932ee40
# ╟─37d05643-4906-4df7-9905-bd6a88f82a21
# ╟─7dbef7ff-437d-4f1f-ad3e-c03f1394107e
# ╟─484dc785-59b3-4e5d-918f-6b71aa2d057c
# ╟─a76b10ce-e373-4361-ba27-9995b5c3f730
# ╟─fbc35df0-c8ae-4704-996e-9c5dd4b27bd2
# ╠═3cd2a0cf-c2e5-494d-abe3-9474dda7fc84
# ╟─02e3f332-c2fd-4bce-93a1-373e3b8eb125
# ╠═4b4ee402-19e4-44a4-8423-dd3835aba538
# ╟─2b543ce0-2acf-4c53-b9e9-fd310f96ffed
# ╠═f488eb97-17a1-44da-8b7a-cc0f236404ac
# ╟─be97646c-68e2-46bb-b66b-ab3a7a432196
# ╠═8fe8e18a-7975-44d5-81d1-c3ed35e27a7e
# ╟─d1fa95a9-92d6-441f-9b00-2f14cc5e96e0
# ╠═717b3cfb-4ee9-409e-a8fe-74e12adf65c6
# ╠═468ef1b4-f0d2-41a9-81ef-1b747a52caaf
# ╟─458fe3d3-671c-40ed-9a52-c2d9e7a945b1
# ╠═9b9a1365-59d9-490f-bd9c-78c8f4d73957
# ╟─aa01b91c-d77a-4f23-8f74-a6f543b1f7d5
# ╟─232de9ed-9ca9-4e4e-9344-c9306755aec5
# ╠═798829dd-af76-48df-891b-fe273aa7dd6d
# ╟─6c90abe5-e0b7-444c-8282-7ac52592964d
# ╟─cca92baf-988a-4f17-a472-9e0870d4239d
# ╟─33e2e0a2-653c-46f2-a8d7-3dc98dfd0f20
# ╠═0a084289-12b1-49b8-bed7-e59b192b5c0a
# ╟─9d90fae3-5cfc-4ff1-a6aa-112e25c98ef0
# ╟─3fe2ae0a-256a-468e-bd59-e80cd00f9c8e
# ╠═c8321df5-a9d5-49ee-b1c5-f1ecbca20663
# ╟─da1cdab6-c7c0-4b20-af71-c4bb6f258e3d
# ╠═5a8fcedb-8c04-4ba3-b0ff-107e7e763d42
# ╠═caa06a90-0434-4d79-ae19-08c06254762d
# ╟─42b8af96-db21-4f2c-bb79-5a8570a710fa
# ╟─f4c42628-42a8-45a1-9e8c-3637a9433d98
# ╠═14bfe60e-fe07-4774-9055-4b76806ca4d8
# ╟─047be3bb-357b-4a2e-abac-c39b797ae4ab
# ╠═77e162be-11a9-46f9-a7bd-c9f8c49f0313
# ╠═a1070752-64fb-427f-962d-721db16a44b9
# ╟─3ad0c65d-2431-45a8-a022-9c49cbee996c
# ╠═78c746e9-1348-47ae-83f7-d3fb1dfcc117
# ╠═717e664b-9d46-4729-9aa6-0e6c32471111
# ╟─a16cbfbc-3ea3-4891-bf3a-65304fca3fc1
# ╠═4e5b9b1a-9e06-490b-932a-b311354962bc
# ╠═74859413-5696-4ed3-8015-9794f6bbff98
# ╠═557d94ce-89a8-4f92-9e17-cbccb0ecb1d6
# ╟─72b702ac-7d7c-4d4d-93b4-9af4657cda5e
# ╟─465b2854-5e02-4389-b917-f3d9ecb619bc
# ╠═06414cfb-b613-4f5f-a0e2-8dec15e8c6e4
# ╠═9834d057-887e-4875-b023-bf763296fe0e
# ╠═ef6804c9-58a2-496b-876e-253618678af4
# ╟─c2554c31-d09a-40f5-aec1-95e6d8ab145c
# ╠═e8e80d7d-3476-4825-9799-9a083bad8aeb
# ╟─4dbad986-5975-40f2-aad0-8795e53d0a92
# ╠═c7d9f4f3-e01a-425f-9037-4807cea0ab08
# ╠═6476dc63-52e0-4304-a849-190e82f55a8f
# ╟─dc4b13f1-9cc8-49e7-a30c-ba5903de2dc9
# ╠═7eaa0e06-f1d1-4e6d-9080-c33603b4d62f
# ╠═54b9a02c-d707-42af-9661-f26a57640009
# ╟─738e2006-b560-4481-9a50-679f6026ccd6
# ╠═e6a9a177-8e4f-4e1e-b31e-d5d2ff482700
# ╟─eb5cea54-1efb-451b-8c6f-76dc3debe99d
# ╠═dd4e6b96-fd36-41ff-8e92-be8eda68162d
# ╟─88382bc4-53e0-4211-a8e7-aec464cae474
# ╠═759f6471-908e-4e08-9e0b-a56c4d0ece96
# ╠═2b05b91e-abd1-44d1-8016-32d78f02b1c7
# ╟─62219eaf-6b01-4e06-ba7b-5247b081ad49
# ╟─d35f2cf1-8111-4bc1-9e77-fdd0db770bca
# ╠═22a258e7-71a0-433b-b4b9-cbfd8a3d4a84
# ╠═61341695-7a95-4164-a519-7fe347a56738
