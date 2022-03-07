### A Pluto.jl notebook ###
# v0.18.1

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
using Revise, ControlSystems, Plots, PlutoUI, Random

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
begin
	using ControlSystemIdentification
	w = 2pi .* exp10.(LinRange(-3, log10(0.5), 500)) # Frequency vector for plots
	G0 = tf(1, [10, 1]) # The true system, 10ẋ = -x + u
	Gtemp = c2d(G0, 1)  # discretize with a sample time of 1s
	G0
end

# ╔═╡ 1569ee96-9567-11ec-2f28-419213ea3fab
html"<button onclick='present()'>present</button>"

# ╔═╡ 37d05643-4906-4df7-9905-bd6a88f82a21
md"""
# JuliaControl
## Outline

- Who am I?
- What is control?
    - The curse of uncertainty
    - The magic of feedback
    - A little on linear dynamical systems
- What does JuliaControl offer?
    - [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
    - [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl)
    - [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl)
    - [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)

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

# ╔═╡ a26ef44f-b7d7-4811-b235-2054b190214b
md"""
# What is control?
The theory and practice of making dynamical systems behave the way we want.

Canonical examples:
- Control the speed of a car.
- Keep an airplane in the air.
- Keep the indoor temperature at a comfortable level.

Controlled systems are everywhere:
- Your body has countless of tightly regulated systems.
- Electrical amplifiers are control systems.
- Financial markets are regulated.
- Server and network queues may not grow too long.

Control theory is often a hidden technology
- Always present
- Never talked about
    - (Unless it fails)
"""

# ╔═╡ 13965d46-319a-4921-89b8-a2f27304ca69
md"""
# Why do we need control?
### Uncertainty.

Our knowledge of the world is inaccurate, models are inaccurate, sensors are noisy, the world is changing.

**Uncertainty is the fundamental motivation for feedback** .
"""

# ╔═╡ b53f9a9e-366e-422a-bf1d-124b77c501f4
md"""
# The magic of feedback
Feedback implies making and acting on **measurements**.

Feedback allows you to make accurate systems out of inaccurate components:

### Early example: Vacuum tube amplifiers
Vacuum tubes were used for early telecommunication.

A vacuum tube is *highly nonlinear* and the amplification *changes with temperature*. How can you possbly make a good sounding amplifier out of such a component?

Enter feedback...
"""

# ╔═╡ 8fa558af-0b87-4cd1-93d3-b1f5710a594b
md"""
Let's say the vacuum tube has an amplification gain $g$, which may be time varying and depend on the input voltage etc.

```
 r  ┌───┐ y
───►│ g ├──►
    └───┘
```
$y = gr$

When connected in *negative feedback*, the gain becomes

```
r     ┌───┐   y
───+─►│ g ├─┬───►
  -│  └───┘ │
   │        │feedback
   │  ┌───┐ │
   └──┤ k │◄┘
      └───┘
```
```math
\begin{aligned}
y &= g(r-ky) \\
y &= \dfrac{g}{1 + kg}r
\end{aligned}
```



- By controlling the *feedback gain* $k$, we control the overall amplification.
- As long as $g$ is very large, **the amplification remains constant and linear even if $g$ varies**!

Want to learn more about the history of control?
[Karl Johan Åström. The fascinating history and success of feedback control (Youtube seminar)](https://www.youtube.com/watch?v=R-h66PrQ808&ab_channel=NTNUCybernetics)
"""

# ╔═╡ 7a871898-e00b-4f10-ab55-58058758c3cf
md"""
## Example: Suppress disturbances
**A temperature controlled room** may be influenced by disturbances $d$
- People entering (about 100W per person)
- Open windows
- Sun shining in

The task of the control system is to *regulate* the temperature, keep it constant.
```
               d
               │
r      ┌───┐u  ▼  ┌──────┐temp
────+─►│ K ├───+─►│ Room ├─┬─►
   -│  └───┘      └──────┘ │
    │                      │
    └──────────────────────┘
      measurement feedback
```
$\text{temp} = \dfrac{RK}{1 + RK}r + \dfrac{R}{1 + RK}d$

-  $R = R(s)$: room transfer function
-  $K = K(s)$: controller transfer function

By designing the controller $K$, we can shape the response of the temperature to both changes in reference $r$ and disturbances $d$ (but not independently).

Control theory deals with designing $K$
"""

# ╔═╡ fbac90ff-299b-4512-9fe5-2294baeab12f
md"""
## Caveats
Feedback can make stable systems **unstable**.
"""

# ╔═╡ 9fa1ed83-6463-46ef-9b5d-db828b99cfae
md"Feedback gain $K$ $(@bind K Slider(0:12, show_value=true))"

# ╔═╡ 53e19605-57b7-4e3c-8285-5c3808685408
md"$H(s) = \dfrac{1}{0.1s^3 + s^2 + s + 1}$"

# ╔═╡ 9e372c3b-90e1-4c59-81d2-44a39d8165e6
begin
	H = tf(1, [0.1,1,1,1])
	closed_loop = feedback(H, K)
	plot(step(closed_loop, 25), ylims=(-0.3, 1.3))
end

# ╔═╡ 347f3a17-dd66-48e7-b97b-b6585d643a40
md"""
The field of control theory is very concerned with **stability and robustness** of dynamical systems.
"""

# ╔═╡ 601b389c-4bad-45d9-b079-1ddfdbc9b9a2
md"""
# On linear dynamical systems
Linear time-invariant systems (LTI) are common in control.
- Correspond to linear systems of ODEs.
- Can be analyzed algebraically ([Laplace transform](https://en.wikipedia.org/wiki/Laplace_transform)).
- Nonlinear systems are often *linearized*.

Control systems are often tasked with *regulation*, keeping a variable at a fixed set point.

Linear theory often works well if we stay close to the point of linearization!
"""

# ╔═╡ 484dc785-59b3-4e5d-918f-6b71aa2d057c
md"""
# What is JuliaControl?
[ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) is a toolbox for analyzing and designing control systems
- Similar in spirit and scope to Matlab control systems toolbox
- Linear systems
- Transfer functions, statespace systems, frequency response, LQR, LQG, Kalman filter, Bode and Nyquist plots, root locus, PID controllers
- [Paper describing the toolbox](https://portal.research.lu.se/en/publications/controlsystemsjl-a-control-toolbox-in-julia)
> FBC., Fält, Heimerson & Troeng, (2021). ControlSystems.jl: A Control Toolbox in Julia. 60th IEEE Conference on Decision and Control (CDC)https://doi.org/10.1109/CDC45484.2021.9683403

##### JuliaControl is currently *not*
- Open-loop optimal control
- Nonlinear control

##### Main Contributors
![Mattias](https://portal.research.lu.se/files-asset/7619687/_MPP9977.jpg?w=160&f=webp)
![Albin](https://portal.research.lu.se/files-asset/110428020/albinh2_cropped_lowres.jpg?w=160&f=webp)
![Olof](https://portal.research.lu.se/files-asset/7213474/OlofTroeng.jpg?w=160&f=webp)
![Fredrik](https://portal.research.lu.se/files-asset/114575837/bagge.jpg?w=160&f=webp)

Started by Jim Christ in 2016
"""

# ╔═╡ cca92baf-988a-4f17-a472-9e0870d4239d
md"""
## Do we offer something unique?
System types are generic w.r.t. element types
- Systems with symbolic coefficients, uncertain parameters, `BigFloat` etc.
- Delay systems
- GPU support
"""

# ╔═╡ 9d90fae3-5cfc-4ff1-a6aa-112e25c98ef0
md""" #### GPU
```julia
A   = cu([-1.2 0; 0.1 -0.4]);  B = cu([1.2; 0])
C   = cu([1 0]);               D = cu([0])
sys = HeteroStateSpace(A, B, C, D)
```
#### Symbolics
```julia
@syms δ k2 γ
Gc = -k2 * (s + δ) / (s + γ)
sys_cl = feedback(H, Gc)
```
#### Delay systems and High Precision
```julia
sys = feedback(BigFloat(1)/s, delay(1))
impulse(sys, BigFloat.(0:0.1:10))
```
"""

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
	τ = 8
	P = delay(τ) * P0
	C_sp = feedback(C0, (1.0 - delay(τ))*P0)
	
	## Plot the closed loop response 
	# Reference step at t = 0 and load disturbance step at t = 15
	G = [feedback(P*C_sp, 1) feedback(P, C_sp)]
	fig_timeresp = plot(lsim(G, t -> [1; t >= 15], 0:0.1:40),  title="τ = $τ")
	
	## Plot the frequency response of the predictor part and compare to a negative delay
	C_pred = feedback(1, C0*(ss(1.0) - delay(τ))*P0)
	fig_bode = bodeplot([C_pred, delay(-τ)], exp10.(-1:0.002:0.4), ls=[:solid :solid :dash :dash], title="")
	plot!(yticks=[0.1, 1, 10], sp=1)
	plot!(yticks=0:180:1080, sp=2)
	
	## Check the Nyquist plot
	# Note that the Nyquist curve encircles -1 for τ > 2.99
	fig_nyquist = nyquistplot(C_sp * P, exp10.(-1:1e-4:2), title="τ = $τ")
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
	optres = optimize(p->cost(P, p), params, Optim.Options(
	    show_trace        = true,
	    show_every        = 50,
	))

	function plot_optimized(P, params, res)
	    fig = plot(layout=(1,3), size=(1200,400), bottommargin=2mm)
	    for (i,params) = enumerate((params, res.minimizer))
	        C, G, S, CS, r = systems(P, params)
	        mag = bode(S, Ω)[1][:]
	        plot!(r, title="Time response", subplot=1,
				lab= i==1 ? "Initial" : "Optimized", legend=:bottomright,
				fillalpha=0.05, linealpha=0.8, seriestype=:path)
	        plot!(Ω, mag, title="Sensitivity function",
				xscale=:log10, yscale=:log10, subplot=2,
				lab= i==1 ? "Initial" : "Optimized",
				legend=:bottomright, fillalpha=0.05, linealpha=0.8)
			nyquistplot!(P*C, Ω, Ms_circles=Msc, sp=3, ylims=(-2.1,1.1), xlims=(-2.1,1.2), lab= i==1 ? "Initial" : "Optimized", points=true, seriescolor=i)
	    end
	    hline!([Msc], l=(:black, :dash), subplot=2, lab="", ylims=(9e-2, Inf))
	    fig
	end
	
	## We can now perform the same computations as above to visualize the found controller	
	plot_optimized(P, params, optres)
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
	pu = 1   ± 0.3
	ζu = 0.1..0.35
	ωu = 1   ± 0.05
	Pu = ss(tf([pu*ωu], [1, 2ζu*ωu, ωu^2]))
end

# ╔═╡ 557d94ce-89a8-4f92-9e17-cbccb0ecb1d6
plot_optimized(Pu, params, optres)

# ╔═╡ 72b702ac-7d7c-4d4d-93b4-9af4657cda5e
md"If we now visualize the found controller, we find that it might violate the constraint we placed on $M_s$"

# ╔═╡ 465b2854-5e02-4389-b917-f3d9ecb619bc
md"""
#### Robust optimization
To alleviate this, we perform the optimization again, this time with the uncertain system.
In order to do this, we must change our cost function slightly so that it outputs a scalar instead of an uncertain numer
"""

# ╔═╡ 06414cfb-b613-4f5f-a0e2-8dec15e8c6e4
cost(Pu, params)

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
cost_uncertain(Pu, params)

# ╔═╡ c2554c31-d09a-40f5-aec1-95e6d8ab145c
md"We can now perform the optimization and visualization again"

# ╔═╡ e8e80d7d-3476-4825-9799-9a083bad8aeb
begin
	optresu = optimize(p->cost_uncertain(Pu, p), optres.minimizer, Optim.Options(
	    show_trace        = true,
	    show_every        = 50,
	))
	plot_optimized(Pu, params, optresu)
end

# ╔═╡ 4dbad986-5975-40f2-aad0-8795e53d0a92
md"""
This time, all realizations are satisfying the constraints! The controller achieved this by being less aggressive. The controller parameters ($k_P,k_I,k_D$) are given below
"""

# ╔═╡ c7d9f4f3-e01a-425f-9037-4807cea0ab08
(nominal=exp.(optres.minimizer), robust=exp.(optresu.minimizer))

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
begin
	u = sign.(sin.((0:0.01:20) .^ 2))' # sample a control input for identification
	y, t, x = lsim(ss(Gtemp), u) # Simulate the true system to get test data
	yn = y .+ 0.2 .* randn.() # add measurement noise
	data = iddata(yn, u, t[2] - t[1]) # create a data object
	plot(data)
end

# ╔═╡ 738e2006-b560-4481-9a50-679f6026ccd6
md"""
We see that the data we're going to use for identification is a chirp input. Chirps are excellent for identification as they have a well defined and easily controllable interval of frequencies for identification. We start by inspecting the coherence plot to ensure that the data is suitable for identification of a linear system
"""

# ╔═╡ e6a9a177-8e4f-4e1e-b31e-d5d2ff482700
coherenceplot(data, hz=true)

# ╔═╡ eb5cea54-1efb-451b-8c6f-76dc3debe99d
md"""
The coherence is high for all frequencies spanned by the chirp, after which it drops significantly. This implies that we can only ever trust the identified model to be accurate up to the highest frequency that was present in the chirp input.

Next we set the parameters for the estimation, the numerator and denominator have one parameter each, so we set $n_a = n_b = 1$ and estimate two models.
"""

# ╔═╡ dd4e6b96-fd36-41ff-8e92-be8eda68162d
begin
	na, nb = 1, 1 # system order and number of parameters in the numerator
	Gh = subspaceid(data, na, r=6, zeroD=true).sys
	Gh2, noise_model = plr(data, na, nb, 1) # try another identification method
	Gh, Gh2
end

# ╔═╡ 88382bc4-53e0-4211-a8e7-aec464cae474
md"""
We can plot the results in several different ways:
"""

# ╔═╡ 759f6471-908e-4e08-9e0b-a56c4d0ece96
d2c(Gh) # Transform to continuous time

# ╔═╡ 2b05b91e-abd1-44d1-8016-32d78f02b1c7
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ControlSystemIdentification = "3abffc1c-5106-53b7-b354-a47bfc086282"
ControlSystems = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
ControlSystemIdentification = "~2.2.3"
ControlSystems = "~0.12.4"
MonteCarloMeasurements = "~1.0.7"
Optim = "~1.6.1"
Plots = "~1.25.11"
PlutoUI = "~0.7.35"
Revise = "~3.3.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "56c347caf09ad8acb3e261fe75f8e09652b7b05b"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.10"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "ce68f8c2162062733f9b4c9e3700d5efc4a8ec47"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.16.11"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "5e98d6a6aa92e5758c4d58501b7bf23732699fa3"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "849799453de85b55e78550fc7b0c8f442eb497ab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.8"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "03dc838350fbd448fca0b99285ed4d60fc229b72"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.5"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "759a12cefe1cd1bb49e477bc3702287521797483"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.7"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ComponentArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "Requires"]
git-tree-sha1 = "470e3dd0295dc7838c60e2f59272b8e998dde916"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.11.9"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.ControlSystemIdentification]]
deps = ["BandedMatrices", "ComponentArrays", "ControlSystems", "DSP", "DelimitedFiles", "FFTW", "FillArrays", "LinearAlgebra", "LowLevelParticleFilters", "MatrixEquations", "MonteCarloMeasurements", "Optim", "Parameters", "QuadGK", "Random", "RecipesBase", "Roots", "Statistics", "StatsBase", "TotalLeastSquares"]
git-tree-sha1 = "90c6eb9f15c9d2aaa5aa72feb4a521cc06d5da42"
uuid = "3abffc1c-5106-53b7-b354-a47bfc086282"
version = "2.2.3"

[[deps.ControlSystems]]
deps = ["Colors", "DSP", "DelayDiffEq", "DiffEqCallbacks", "IterTools", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MatrixEquations", "MatrixPencils", "OrdinaryDiffEq", "Polynomials", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1d19487e223ad375dccd4b5af08d67124bb22354"
uuid = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
version = "0.12.4"

[[deps.DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "3e03979d16275ed5d9078d50327332c546e24e68"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.5"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "ceb3463f2913eec2f0af5f0d8e1386fb546fdd32"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.34.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "433291c9e63dcfc1a0e42c6aeb6bb5d3e5ab1789"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.81.4"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "e57ecaf9f7875714c164ccca3c802711589127cf"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.20.1"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3e1289d9a6a54791c1ee60da0850f4fd71188da6"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.11.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "0f8ef5dcb040dbb9edd98b1763ac10882ee1ff03"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.12"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ec299fdc8f49ae450807b0cb1d161c6b76fd2b60"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.1"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "3965a3216446a6b020f0d48f1ba94ef9ec01720d"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.6"

[[deps.Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "0a815f0060ab182f6c484b281107bfcd5bbb58dc"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.7"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "cae5e3dfd89b209e01bcd65b3a25e74462c67ee0"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.3.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "6333cc5b848295895f3b23eb763d020fc8e05867"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.7.12"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "97e2adfcbe7ac07112ca79f03e34fc88cac6b9e7"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.2"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "6dd77ee76188b0365f7d882d674b95796076fa2c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.5"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "dbb14c604fc47aa4f2e19d0ebb7b6416f3cfa5f5"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.5.1"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "f27bb8e4eabdb93ed3703c55025b111e045ffe81"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "534aa24fae56f5f0956134d8789ab30d6fe2f615"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.102"

[[deps.LowLevelParticleFilters]]
deps = ["Distributions", "ForwardDiff", "Lazy", "LinearAlgebra", "LoopVectorization", "PDMats", "Parameters", "Printf", "Random", "RecipesBase", "Requires", "StaticArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "f36f172ca61a201d65ed83ade9190898d3ba18a4"
uuid = "d9d29d28-c116-5dba-9239-57a5fe23875b"
version = "2.0.1"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "6b0440822974cab904c8b14d79743565140567f6"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.1"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MatrixEquations]]
deps = ["LinearAlgebra", "LinearMaps"]
git-tree-sha1 = "49b86537cefbfbf4b86afa3cc17226f1781d4389"
uuid = "99c1a7ee-ab34-5fd5-8076-27c950a045f4"
version = "2.1.0"

[[deps.MatrixPencils]]
deps = ["LinearAlgebra", "Polynomials", "Random"]
git-tree-sha1 = "caf22dd55c613095638a113b39f38e20016e2c36"
uuid = "48965c70-4690-11ea-1f13-43a2532b2fa8"
version = "1.6.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[deps.MonteCarloMeasurements]]
deps = ["Distributed", "Distributions", "LinearAlgebra", "MacroTools", "Random", "RecipesBase", "Requires", "SLEEFPirates", "StaticArrays", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "bce0a32d6c64165388eec2573b68000ed06f39c1"
uuid = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
version = "1.0.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "b61c51cd5b9d8b197dfcbbf2077a0a4e1505278d"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.14"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "045d10789f5daff18deb454d5923c6996017c2f3"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.6.1"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "df82fa0f9f90f669cc3cf9e3f0400e431e0704ac"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.6.6"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "95a4038d1011dfdbde7cecd2ad0ac411e53ab1bc"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.10.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5c907bdee5966a9adb8a106807b7c387e51e4d6c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.11"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "2232d3865bc9a098e664f69cbe340b960d48217f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.6"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "dc11fa882240c43a875b48e21e6423704927d12f"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.4"

[[deps.Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "a1f7f4e41404bed760213ca01d7f384319f717a5"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.25"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "e4cb8d4a2edf9b3804c1fb2c2de57d634ff3f36e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "736699f42935a2b19b37a6c790e2355ca52a12ee"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.2"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "7ad4c2ef15b7aecd767b3921c0d255d39b3603ea"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.9"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "606ddc4d3d098447a09c9337864c73d017476424"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.3.2"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "0abe7fc220977da88ad86d339335a4517944fea2"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.14"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "61a96d8b89083a53fb2b745f3b59a05359651bbe"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.30"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "8ff1bf96965b3878ca5d235752ff1daf519e7a26"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "87efd1676d87706f4079e8e717a7a5f02b6ea1ad"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.2"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "e0a02838565c4600ecd1d8874db8cfe263aaa6c7"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.12"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "0f1017f68dc25f1a0cb99f4988f78fe4f2e7955f"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.1"

[[deps.TotalLeastSquares]]
deps = ["FillArrays", "LinearAlgebra", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "61ab1b838aa6df7773bf593ba1281b60e35903e6"
uuid = "028f657a-7ace-5159-a694-8cfd97933b0c"
version = "1.7.1"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "5cbc1a4551fcf8afe8f80bb4f1f13e3271ee2656"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "e9a35d501b24c127af57ca5228bcfb806eda7507"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.24"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─1569ee96-9567-11ec-2f28-419213ea3fab
# ╠═8fe0da54-e41f-45a2-9d79-3df73932ee40
# ╟─37d05643-4906-4df7-9905-bd6a88f82a21
# ╟─7dbef7ff-437d-4f1f-ad3e-c03f1394107e
# ╟─a26ef44f-b7d7-4811-b235-2054b190214b
# ╟─13965d46-319a-4921-89b8-a2f27304ca69
# ╟─b53f9a9e-366e-422a-bf1d-124b77c501f4
# ╟─8fa558af-0b87-4cd1-93d3-b1f5710a594b
# ╟─7a871898-e00b-4f10-ab55-58058758c3cf
# ╟─fbac90ff-299b-4512-9fe5-2294baeab12f
# ╟─9fa1ed83-6463-46ef-9b5d-db828b99cfae
# ╟─53e19605-57b7-4e3c-8285-5c3808685408
# ╠═9e372c3b-90e1-4c59-81d2-44a39d8165e6
# ╟─347f3a17-dd66-48e7-b97b-b6585d643a40
# ╟─601b389c-4bad-45d9-b079-1ddfdbc9b9a2
# ╟─484dc785-59b3-4e5d-918f-6b71aa2d057c
# ╟─cca92baf-988a-4f17-a472-9e0870d4239d
# ╟─9d90fae3-5cfc-4ff1-a6aa-112e25c98ef0
# ╟─42b8af96-db21-4f2c-bb79-5a8570a710fa
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
# ╟─e8e80d7d-3476-4825-9799-9a083bad8aeb
# ╟─4dbad986-5975-40f2-aad0-8795e53d0a92
# ╠═c7d9f4f3-e01a-425f-9037-4807cea0ab08
# ╟─6476dc63-52e0-4304-a849-190e82f55a8f
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
