### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 562d445c-f12e-11ec-245c-e5a7d7c682d1
begin
	using Revise, Pkg
	Pkg.activate()
	using ControlSystems, Plots, PlutoUI, Random, LinearAlgebra
	default(label="")
end

# ╔═╡ 80d45c7f-f914-4a68-a4a4-c90e0e8ac5c8
using Symbolics

# ╔═╡ 8f43ecab-2a47-4e3c-b55b-3d84101a1c37
begin
	using MonteCarloMeasurements
	k  = 1..1.2   # Uniform uncertainty
	ω  = 1 ± 0.1  # Gaussian uncertainty
	ζ  = 0.1..0.2
	Pᵤ = tf(k*ω^2, [1 , 2ζ*ω, ω^2])
	plot( step(Pᵤ, 50),
		ri=false, N=1000, α=0.1, title="Step response of \$ P_u(s) = k\\dfrac{ω^2}{s^2 + 2ζω s + ω^2}\$")
end

# ╔═╡ b46fc3fe-050c-4f16-8e4a-ff5f2240c6cc
begin
	using ControlSystemIdentification

	# Create a system
	w = 2pi .* exp10.(LinRange(-3, log10(0.5), 500)) # Frequency vector for plots
	G0 = tf(1, [10, 1]) # The true system, 10ẋ = -x + u
	Gtemp = c2d(G0, 1)  # discretize with a sample time of 1s

	# Generate data
	u = sign.(sin.((0:0.01:20) .^ 2))' 	# sample a control input for identification
	y, t, x = lsim(ss(Gtemp), u) 		# Simulate the true system to get test data
	yn = y .+ 0.2 .* randn.() 			# add measurement noise
	data = iddata(yn, u, t[2] - t[1]) 	# create a data object
	md"Generate data for system identification"
end

# ╔═╡ ee89ca86-f692-40b8-85db-f5fb20203f17
html"<button onclick='present()'>present</button>"

# ╔═╡ e22bc962-67e7-4d14-b98d-9424d47365e9
md"""
# JuliaControl
## JuliaControl -- Outline

- Who am I?
- What does JuliaControl offer?
    - [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
    - [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl)
    - [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl)
    - [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)
- Interaction with the Julia ecosystem
    - [ModelingToolkit.jl](https://mtk.sciml.ai/dev/)
    - JuliaSim

![logo](https://avatars.githubusercontent.com/u/10605979?s=400&u=7b2efdd404c4db3b3f067f04c305d40c025a8961&v=4)
"""

# ╔═╡ 399c861f-c478-4097-a66e-4ac9aa4ecca0
md"""
## Who am I?

My name is **Fredrik Bagge Carlson**, I live in Lund in southern **Sweden**.

My background is within **robotics, control and system identification**.

I enjoy writing **software toolboxes** 😄 [baggepinnen @ github](https://github.com/baggepinnen)

I work with simulation and control at **Julia Computing** (I represent myself here, all opinions are my own).
"""

# ╔═╡ 337f4c13-4458-4fdc-9039-b370a22b9606
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

# ╔═╡ 90e060d8-cb27-4f82-bdba-62f8ac8d8df5
md"""
# Why do we need control?
### Uncertainty.

Our knowledge of the world is inaccurate, models are inaccurate, sensors are noisy, the world is changing.

**Uncertainty is the fundamental motivation for feedback**.

```
 u  ┌───┐ y 
───►│ P ├──►   Open loop
    └───┘
```

### Feedback
Feedback implies making and acting on **measurements**.
```
r   ┌─────┐     ┌─────┐
───►│     │  u  │     │ y
    │  C  ├────►│  P  ├─┬─►  Closed loop
  ┌►│     │     │     │ │
  │ └─────┘     └─────┘ │
  │                     │
  └─────────────────────┘
```	
Feedback allows you to make *accurate systems* out of *inaccurate components*.
"""

# ╔═╡ d07debe4-bb74-4842-8cb1-24ea985241ea
md"""
# What is JuliaControl?
[ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) is a toolbox for analyzing and designing control systems
- Similar in spirit and scope to Control systems toolbox
- Linear systems (some nonlinear simulation)
- Transfer functions, statespace systems, frequency response, LQR, LQG, Kalman filter, Bode and Nyquist plots, root locus, PID controllers
- [Paper describing the toolbox](https://portal.research.lu.se/en/publications/controlsystemsjl-a-control-toolbox-in-julia)
> FBC., Fält, Heimerson & Troeng, (2021). ControlSystems.jl: A Control Toolbox in Julia. 60th IEEE Conference on Decision and Control (CDC)https://doi.org/10.1109/CDC45484.2021.9683403

[RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl) contains facilities for robust analysis and design
[ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl) contains LTI identification
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

# ╔═╡ 09d17cf3-543c-4088-adcf-2cc89a637bd5
md"""
## Do we offer something unique?
System types are generic w.r.t. element types
- Systems with symbolic coefficients, uncertain parameters, `BigFloat` etc.
- Delay systems
- GPU support
- Automatic differentiation
"""

# ╔═╡ ddc3e3f7-5d98-40ee-a8b1-6dfb4876e3fa
md"""
# Look and feel
"""

# ╔═╡ 7262ac30-2153-4ad1-94e8-81264d1f42c4
md"""
## Transfer functions:
- Linear ODE in time domain ``\quad \xrightarrow{\text{Laplace transform}}\quad`` Transfer function in frequency domain  
- Differential equation ``\quad \xrightarrow{\text{Laplace transform}}\quad`` Algebraic equation
```math
P(s) = \dfrac{B(s)}{A(s)}
```
"""

# ╔═╡ 1c0a6cc3-4f3c-41aa-9390-ffeda5c3ad0f
tf(1, [1, 2, 1])

# ╔═╡ b6ac5888-8ba9-4d73-9ec2-6385d9567878
md"""
## Discrete-time systems:
```math
P(z) = \dfrac{B(z)}{A(z)}
```
"""

# ╔═╡ e8e38785-4ccb-4ef5-8412-255961b3df79
c2d(tf(1, [1,2,1]), 0.1)

# ╔═╡ d091e596-4d71-4d7b-b84b-8e38bd8dec15
md"## Zero-pole-gain representation:"

# ╔═╡ 9d08530e-7bc1-4883-9715-31c1f56ee50e
zpk(tf(1, [1,2,1]))

# ╔═╡ 28d5f01d-d49d-49c6-b6a7-26cc9922770e
md"""
## State-space systems:
```math
\begin{aligned}
\dot x &= Ax + Bu\\
y &= Cx + Du
\end{aligned}
```
"""

# ╔═╡ aa19f719-29f5-49b4-a1bd-4c1e84fec309
sys = let
	A = [1 0.1; 0 0.7]
	B = [0; 1]
	C = [1 0]
	ss(A,B,C,0, 0.1)
end

# ╔═╡ 443ab987-f5b4-403c-ab72-a865dafac72f
md"## Standard plots"

# ╔═╡ 0d6c654b-0d42-4508-b8d9-6e88b1288cd1
plot(
	bodeplot(sys, title=["Bode plot" ""]),
	
	plot(impulse(sys, 10), title="Impulse response"),
	
	nyquistplot(sys, ylim=(-2, 1), xlims=(-2, 1), Ms_circles=1, xlab="Re", ylab="Im"),
	
	layout=(1,3),
	size=(1000, 400),
	leftmargin = Plots.PlotMeasures.mm * 5,
	bottommargin = Plots.PlotMeasures.mm * 5,
)

# ╔═╡ 328ce7d8-f182-4eeb-bb97-f6ebbdd969c9
md"## Support for time delays:"

# ╔═╡ 0ea92c94-2cd5-4850-ba8b-cdb48ee6de7f
Pdelay = tf(1, [1,2,1]) * delay(1.5)

# ╔═╡ 36384e19-5a8b-4bd8-adb0-93af32e6a5de
plot(
	      plot(step(Pdelay, 12)), 	               nyquistplot(Pdelay)
)

# ╔═╡ f52521a2-2c5d-4680-994f-885430520762
md"## Support for nonlinearities:"

# ╔═╡ 4008c41e-bce0-4aaa-9693-0a50a427f914
md"""
saturation: = $(@bind sat Slider(0.5:0.5:8, default=5, show_value=true))
"""

# ╔═╡ 574a9970-5718-4d3e-a7c3-4e25c38b4e50
let
	using ControlSystems: saturation
	P  = tf(1, [1,2,1])
	C  = saturation(sat) * tf(10) 					# Saturated P controller
	T  = feedback(P*C)
	CS = feedback(C, P)
	stepsim = step([T; CS], 5) 						# Step-response simulation
	plot(stepsim, ylim=[(0, 1.3) (-4,8)], ylabel=["y" "u"], framestyle=:zerolines, title=["Closed-loop step response" ""])
end

# ╔═╡ 0a3f9306-493a-4ea1-95fd-862c02aa1c34
md"## Basic design:"

# ╔═╡ 28ff405a-fe91-42e2-abc7-70fec9372f3c
md"""
``T_s`` = $(@bind Ts Slider(round.(exp10.(-2:0.1:0), sigdigits=2), default=0.1, show_value=true))
``R`` = $(@bind r Slider(round.(exp10.(-3:0.2:3), sigdigits=1), default=10, show_value=true))
"""

# ╔═╡ 8245d095-d730-4534-a477-3a370cdc000d
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
	res    = lsim(sys,u,t,x0=[1, 0])
	plot(res, lab=["Output" "Input + disturbance" "Position" "Velocity"], xlabel="Time [s]", plotx = true, plotu=true, ploty=false, size=(400, 500))
end

# ╔═╡ 85d116d5-f393-46c3-99d6-c3041f63fb36
md"## Symbolics"

# ╔═╡ 9af80959-b4f2-4896-8ca0-c17f262bf00c
let
	@variables δ k2 γ
	s  = tf("s")
	Gc = -k2 * (s + δ) / (s + γ)
	H  = tf(1, [1, 1])
	sys_cl = feedback(H, Gc)
end

# ╔═╡ f700d0db-aef9-490b-a316-8d9146713e6b
md""" ## GPU
```julia
using CUDA
A   = cu([-1.2 0; 0.1 -0.4]);  B = cu([1.2; 0])
C   = cu([1 0]);               D = cu([0])
sys = HeteroStateSpace(A, B, C, D)
```
"""

# ╔═╡ 6903f5c8-ef41-446a-855b-5e5adbe9e7d6
md"## Automatic differentiation"

# ╔═╡ 7b77f641-4a01-4630-9961-1a68285a8348
function cost_function(params)
	gain, pole = params
	P = tf(1, [1, 1])
	C = tf(gain, [1, pole])
	closed_loop = feedback(P*C)
	y = step(closed_loop, 0:0.1:5).y
	sum(abs2, y .- 1)
end

# ╔═╡ 2a6e7eac-49f7-4363-9c83-2e22827ac344
let
	using ForwardDiff
	params = [1.0, 1.0]
	∇cost = ForwardDiff.gradient(cost_function, params)
end

# ╔═╡ a03461fc-9910-4a7c-a7e7-1280dc49610d
md"## Uncertain systems"

# ╔═╡ 815217df-b5a8-45ae-9b59-703e572c7440
md"""
# [RobustAndOptimalControl.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl)
An extension to ControlSystems.jl containing
- Uncertainty modeling and $\mu$-analysis
- Named systems
- Extra LQG functionality
-  $\mathcal{H}_\infty$ and $\mathcal{H}_2$ optimal design.
"""

# ╔═╡ b3ef6e4d-bbad-41ba-92c2-d9dc7c417282
md"""
# System identification
The art and science of estimating models of dynamical systems.

[ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl)
- Black-box linear models.
- Time and frequency domain.
- Subspace-based methods (N4SID etc.)
- Prediction-error method (PEM)
"""

# ╔═╡ 1388489b-e749-4dc3-80de-47c36db6236b
md"## Example"

# ╔═╡ 0ffa64ff-247e-48fb-add1-aede23b5c871
data

# ╔═╡ 0496f355-687c-4171-adf7-a55b736839a1
plot(data)

# ╔═╡ 151cf22f-604e-4326-8418-89cd372ccb5c
P̂ = subspaceid(data, 1)

# ╔═╡ 5ef8e90c-7ba3-4804-83e4-16020adad29d
predplot(P̂, data[1:500])

# ╔═╡ bea93409-5d31-4050-920e-bdc99d2cc2ee
md"""
# [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)
- Tools for working with systems with symbolic coefficients.
- C-code generation of controllers and filters.
"""

# ╔═╡ 349c30e2-cfa4-49de-a691-40cca1f9ef3a
md"""
# ControlSystems ❤ ModelingToolkit
"""

# ╔═╡ dd961e0b-6c19-476f-a9a6-123e8ae77587
md"""
# JuliaSim
"""

# ╔═╡ bdac9817-e65e-4a31-b85f-bd5c5e559757
md"""
# Summary
- [JuliaControl](https://github.com/JuliaControl/) tries to offer a fully-featured environment for analysis, design and estimation of linear control systems.
- Don't hesitate to reach out!
    - [Issue](https://github.com/JuliaControl/ControlSystems.jl/issues)
    - [Discourse](https://discourse.julialang.org/tag/control)
    - [`#control`](https://julialang.slack.com/archives/CKH1UTZT9) on slack
    - PRs are welcome, including tutorials and examples.
- If you use the software, be sure to say hi and mention what you like and what can be improved.
- This notebook is available in [ControlExamples.jl/presentations](https://github.com/JuliaControl/ControlExamples.jl/tree/master/presentations)

##### Thank you 😃
"""

# ╔═╡ 0d212db3-2ac1-40e6-ac77-617e10cda336
html"""<style>
main {
    max-width: 60%;
}
"""

# ╔═╡ 13355e2a-8c46-4892-9300-51439ff0bb17
PlutoUI.TableOfContents()

# ╔═╡ Cell order:
# ╟─ee89ca86-f692-40b8-85db-f5fb20203f17
# ╟─e22bc962-67e7-4d14-b98d-9424d47365e9
# ╟─399c861f-c478-4097-a66e-4ac9aa4ecca0
# ╟─337f4c13-4458-4fdc-9039-b370a22b9606
# ╟─90e060d8-cb27-4f82-bdba-62f8ac8d8df5
# ╟─d07debe4-bb74-4842-8cb1-24ea985241ea
# ╟─09d17cf3-543c-4088-adcf-2cc89a637bd5
# ╟─ddc3e3f7-5d98-40ee-a8b1-6dfb4876e3fa
# ╟─7262ac30-2153-4ad1-94e8-81264d1f42c4
# ╠═1c0a6cc3-4f3c-41aa-9390-ffeda5c3ad0f
# ╟─b6ac5888-8ba9-4d73-9ec2-6385d9567878
# ╠═e8e38785-4ccb-4ef5-8412-255961b3df79
# ╟─d091e596-4d71-4d7b-b84b-8e38bd8dec15
# ╠═9d08530e-7bc1-4883-9715-31c1f56ee50e
# ╟─28d5f01d-d49d-49c6-b6a7-26cc9922770e
# ╠═aa19f719-29f5-49b4-a1bd-4c1e84fec309
# ╟─443ab987-f5b4-403c-ab72-a865dafac72f
# ╠═0d6c654b-0d42-4508-b8d9-6e88b1288cd1
# ╟─328ce7d8-f182-4eeb-bb97-f6ebbdd969c9
# ╠═0ea92c94-2cd5-4850-ba8b-cdb48ee6de7f
# ╠═36384e19-5a8b-4bd8-adb0-93af32e6a5de
# ╟─f52521a2-2c5d-4680-994f-885430520762
# ╠═574a9970-5718-4d3e-a7c3-4e25c38b4e50
# ╟─4008c41e-bce0-4aaa-9693-0a50a427f914
# ╟─0a3f9306-493a-4ea1-95fd-862c02aa1c34
# ╠═8245d095-d730-4534-a477-3a370cdc000d
# ╟─28ff405a-fe91-42e2-abc7-70fec9372f3c
# ╟─85d116d5-f393-46c3-99d6-c3041f63fb36
# ╠═80d45c7f-f914-4a68-a4a4-c90e0e8ac5c8
# ╠═9af80959-b4f2-4896-8ca0-c17f262bf00c
# ╟─f700d0db-aef9-490b-a316-8d9146713e6b
# ╟─6903f5c8-ef41-446a-855b-5e5adbe9e7d6
# ╠═7b77f641-4a01-4630-9961-1a68285a8348
# ╠═2a6e7eac-49f7-4363-9c83-2e22827ac344
# ╟─a03461fc-9910-4a7c-a7e7-1280dc49610d
# ╠═8f43ecab-2a47-4e3c-b55b-3d84101a1c37
# ╟─815217df-b5a8-45ae-9b59-703e572c7440
# ╟─b3ef6e4d-bbad-41ba-92c2-d9dc7c417282
# ╟─1388489b-e749-4dc3-80de-47c36db6236b
# ╠═0ffa64ff-247e-48fb-add1-aede23b5c871
# ╠═0496f355-687c-4171-adf7-a55b736839a1
# ╠═151cf22f-604e-4326-8418-89cd372ccb5c
# ╠═5ef8e90c-7ba3-4804-83e4-16020adad29d
# ╟─bea93409-5d31-4050-920e-bdc99d2cc2ee
# ╟─349c30e2-cfa4-49de-a691-40cca1f9ef3a
# ╟─dd961e0b-6c19-476f-a9a6-123e8ae77587
# ╟─bdac9817-e65e-4a31-b85f-bd5c5e559757
# ╠═0d212db3-2ac1-40e6-ac77-617e10cda336
# ╠═562d445c-f12e-11ec-245c-e5a7d7c682d1
# ╠═13355e2a-8c46-4892-9300-51439ff0bb17
# ╟─b46fc3fe-050c-4f16-8e4a-ff5f2240c6cc
