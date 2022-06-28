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

# ╔═╡ 562d445c-f12e-11ec-245c-e5a7d7c682d1
begin
		using Revise, Pkg
		Pkg.activate()
		using ControlSystems, Plots, PlutoUI, Random, LinearAlgebra
end

# ╔═╡ 80d45c7f-f914-4a68-a4a4-c90e0e8ac5c8
using Symbolics

# ╔═╡ b46fc3fe-050c-4f16-8e4a-ff5f2240c6cc
begin
	using ControlSystemIdentification
	w = 2pi .* exp10.(LinRange(-3, log10(0.5), 500)) # Frequency vector for plots
	G0 = tf(1, [10, 1]) # The true system, 10ẋ = -x + u
	Gtemp = c2d(G0, 1)  # discretize with a sample time of 1s
	G0
end

# ╔═╡ ee89ca86-f692-40b8-85db-f5fb20203f17
html"<button onclick='present()'>present</button>"

# ╔═╡ e22bc962-67e7-4d14-b98d-9424d47365e9
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

# ╔═╡ 399c861f-c478-4097-a66e-4ac9aa4ecca0
md"""
## Who am I?

My name is **Fredrik Bagge Carlson**, I live in Lund in southern **Sweden**.

I got my MSc and PhD from Dept. Automatic Control at **Lund University**.

My background is within **robotics, control and system identification**.

I enjoy writing **software toolboxes** 😄 [baggepinnen@github](https://github.com/baggepinnen)

I work with simulation and control at **Julia Computing** (I represent myself here, all opinions are my own).
"""

# ╔═╡ d07debe4-bb74-4842-8cb1-24ea985241ea
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

# ╔═╡ ddc3e3f7-5d98-40ee-a8b1-6dfb4876e3fa
md"""
## Look and feel
"""

# ╔═╡ 7262ac30-2153-4ad1-94e8-81264d1f42c4
md"""
#### Transfer functions:
```math
P(s) = \dfrac{B(s)}{A(s)}
```
"""

# ╔═╡ 1c0a6cc3-4f3c-41aa-9390-ffeda5c3ad0f
tf(1, [1,2,1])

# ╔═╡ b6ac5888-8ba9-4d73-9ec2-6385d9567878
md"""
#### Discrete-time systems:
```math
P(z) = \dfrac{B(z)}{A(z)}
```
"""

# ╔═╡ e8e38785-4ccb-4ef5-8412-255961b3df79
c2d(tf(1, [1,2,1]), 0.1)

# ╔═╡ d091e596-4d71-4d7b-b84b-8e38bd8dec15
md"#### Zero-pole-gain representation:"

# ╔═╡ 9d08530e-7bc1-4883-9715-31c1f56ee50e
zpk(tf(1, [1,2,1]))

# ╔═╡ 28d5f01d-d49d-49c6-b6a7-26cc9922770e
md"""
#### State-space systems:
```math
\begin{aligned}
\dot x &= Ax + Bu\\
y &= Cx + Du
\end{aligned}
```
"""

# ╔═╡ aa19f719-29f5-49b4-a1bd-4c1e84fec309
let
	A = [1 0.1; 0 1]
	B = [0; 1]
	C = [1 0]
	ss(A,B,C,0, 0.1)
end

# ╔═╡ 328ce7d8-f182-4eeb-bb97-f6ebbdd969c9
md"#### Support for time delays:"

# ╔═╡ 0ea92c94-2cd5-4850-ba8b-cdb48ee6de7f
Pdelay = tf(1, [1,2,1]) * delay(1.5)

# ╔═╡ 36384e19-5a8b-4bd8-adb0-93af32e6a5de
plot(step(Pdelay, 12))

# ╔═╡ f52521a2-2c5d-4680-994f-885430520762
md"#### Support for nonlinearities:"

# ╔═╡ 4008c41e-bce0-4aaa-9693-0a50a427f914
md"""
saturation: = $(@bind sat Slider(0.5:0.5:5, default=1, show_value=true))
"""

# ╔═╡ 574a9970-5718-4d3e-a7c3-4e25c38b4e50
let
	using ControlSystems: saturation
	P  = tf(1, [1,2,1])
	C  = saturation(sat) * tf(10) 					# Saturated P controller
	T  = feedback(P*C)
	CS = feedback(C, P)
	stepsim = step([T; CS], 5) 						# Step-response simulation
	plot(stepsim, ylim=[(0, 1.3) (-5,5)], ylabel=["y" "u"])
end

# ╔═╡ 0a3f9306-493a-4ea1-95fd-862c02aa1c34
md"### Basic design:"

# ╔═╡ 28ff405a-fe91-42e2-abc7-70fec9372f3c
md"""
``T_s`` = $(@bind Ts Slider(exp10.(-2:0.1:0), default=0.1, show_value=true))
``R`` = $(@bind r Slider(exp10.(-3:0.2:3), default=10, show_value=true))
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
	plot(res, lab=["Input + disturbance" "Position" "Velocity"], xlabel="Time [s]", plotx = true, plotu=true, ploty=false)
end

# ╔═╡ 09d17cf3-543c-4088-adcf-2cc89a637bd5
md"""
## Do we offer something unique?
System types are generic w.r.t. element types
- Systems with symbolic coefficients, uncertain parameters, `BigFloat` etc.
- Delay systems
- GPU support
- Automatic differentiation
"""

# ╔═╡ 85d116d5-f393-46c3-99d6-c3041f63fb36
md"#### Symbolics"

# ╔═╡ 9af80959-b4f2-4896-8ca0-c17f262bf00c
let
	@variables δ k2 γ
	s  = tf("s")
	Gc = -k2 * (s + δ) / (s + γ)
	H  = tf(1, [1, 1])
	sys_cl = feedback(H, Gc)
end

# ╔═╡ f700d0db-aef9-490b-a316-8d9146713e6b
md""" #### GPU
```julia
using CUDA
A   = cu([-1.2 0; 0.1 -0.4]);  B = cu([1.2; 0])
C   = cu([1 0]);               D = cu([0])
sys = HeteroStateSpace(A, B, C, D)
```
"""

# ╔═╡ 6903f5c8-ef41-446a-855b-5e5adbe9e7d6
md"#### Automatic differentiation"

# ╔═╡ 7b77f641-4a01-4630-9961-1a68285a8348
function cost_function(params)
	gain, pole = params
	P = tf(1, [1,1])
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

# ╔═╡ 0496f355-687c-4171-adf7-a55b736839a1
begin
	u = sign.(sin.((0:0.01:20) .^ 2))' 	# sample a control input for identification
	y, t, x = lsim(ss(Gtemp), u) 		# Simulate the true system to get test data
	yn = y .+ 0.2 .* randn.() 			# add measurement noise
	data = iddata(yn, u, t[2] - t[1]) 	# create a data object
	plot(data)
end

# ╔═╡ cbe9bb03-eed7-4a74-a17f-2dec9b1ae749
md"""
We see that the data we're going to use for identification is a chirp input. Chirps are excellent for identification as they have a well defined and easily controllable interval of frequencies for identification. We start by inspecting the coherence plot to ensure that the data is suitable for identification of a linear system
"""

# ╔═╡ 4b075f42-43d2-4979-98d1-ccf3868894dc
coherenceplot(data, hz=true)

# ╔═╡ 74c78155-33ec-499a-9adc-f260931ff1f8
md"""
The coherence is high for all frequencies spanned by the chirp, after which it drops significantly. This implies that we can only ever trust the identified model to be accurate up to the highest frequency that was present in the chirp input.

Next we set the parameters for the estimation, the numerator and denominator have one parameter each, so we set $n_a = n_b = 1$ and estimate two models.
"""

# ╔═╡ b889e691-32d9-4f33-b9c9-52d18e0113a2

Gh = subspaceid(data, 1, r=6, zeroD=true).sys

# ╔═╡ 18d66398-b853-4837-ae3a-1eed48951b63
md"""
We can plot the results in several different ways:
"""

# ╔═╡ b37914b7-5f7b-446d-a336-d5370db69a19
d2c(Gh) # Transform to continuous time

# ╔═╡ ecc4efc7-1d89-4b5a-8cdc-57055e2ac5bb
begin
	bp = bodeplot(Gtemp, w, lab = "G (true)", hz = true, l = 3)
	bodeplot!(Gh, w, lab = "subspace", hz = true)
	
	sp = plot(step(Gtemp, 150), lab="G (true)")
	plot!(step(Gh, 150), lab = "subspace", ticks = :default)
	hline!([1], primary = false, l = (:black, :dash))
	
	lp = plot(lsim(ss(Gtemp), u), lab="G (true)")
	plot!(lsim(ss(Gh), u), lab = "subspace")
	plot!(data.t, yn[:], lab = "Estimation data", alpha=0.5)
	
	plot(bp, sp, lp, layout = @layout([[a b]; c]), size=(1000,800))
end

# ╔═╡ bea93409-5d31-4050-920e-bdc99d2cc2ee
md"""
# [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl)
- Tools for working with systems with symbolic coefficients.
- C-code generation of controllers and filters.
"""

# ╔═╡ 349c30e2-cfa4-49de-a691-40cca1f9ef3a
md"""
# ModelingToolkit
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
# ╠═562d445c-f12e-11ec-245c-e5a7d7c682d1
# ╟─e22bc962-67e7-4d14-b98d-9424d47365e9
# ╟─399c861f-c478-4097-a66e-4ac9aa4ecca0
# ╟─d07debe4-bb74-4842-8cb1-24ea985241ea
# ╟─ddc3e3f7-5d98-40ee-a8b1-6dfb4876e3fa
# ╟─7262ac30-2153-4ad1-94e8-81264d1f42c4
# ╠═1c0a6cc3-4f3c-41aa-9390-ffeda5c3ad0f
# ╟─b6ac5888-8ba9-4d73-9ec2-6385d9567878
# ╠═e8e38785-4ccb-4ef5-8412-255961b3df79
# ╟─d091e596-4d71-4d7b-b84b-8e38bd8dec15
# ╠═9d08530e-7bc1-4883-9715-31c1f56ee50e
# ╟─28d5f01d-d49d-49c6-b6a7-26cc9922770e
# ╠═aa19f719-29f5-49b4-a1bd-4c1e84fec309
# ╠═328ce7d8-f182-4eeb-bb97-f6ebbdd969c9
# ╠═0ea92c94-2cd5-4850-ba8b-cdb48ee6de7f
# ╠═36384e19-5a8b-4bd8-adb0-93af32e6a5de
# ╟─f52521a2-2c5d-4680-994f-885430520762
# ╠═574a9970-5718-4d3e-a7c3-4e25c38b4e50
# ╟─4008c41e-bce0-4aaa-9693-0a50a427f914
# ╟─0a3f9306-493a-4ea1-95fd-862c02aa1c34
# ╠═8245d095-d730-4534-a477-3a370cdc000d
# ╟─28ff405a-fe91-42e2-abc7-70fec9372f3c
# ╟─09d17cf3-543c-4088-adcf-2cc89a637bd5
# ╟─85d116d5-f393-46c3-99d6-c3041f63fb36
# ╠═80d45c7f-f914-4a68-a4a4-c90e0e8ac5c8
# ╠═9af80959-b4f2-4896-8ca0-c17f262bf00c
# ╟─f700d0db-aef9-490b-a316-8d9146713e6b
# ╟─6903f5c8-ef41-446a-855b-5e5adbe9e7d6
# ╠═7b77f641-4a01-4630-9961-1a68285a8348
# ╠═2a6e7eac-49f7-4363-9c83-2e22827ac344
# ╟─815217df-b5a8-45ae-9b59-703e572c7440
# ╟─b3ef6e4d-bbad-41ba-92c2-d9dc7c417282
# ╠═b46fc3fe-050c-4f16-8e4a-ff5f2240c6cc
# ╠═0496f355-687c-4171-adf7-a55b736839a1
# ╟─cbe9bb03-eed7-4a74-a17f-2dec9b1ae749
# ╠═4b075f42-43d2-4979-98d1-ccf3868894dc
# ╟─74c78155-33ec-499a-9adc-f260931ff1f8
# ╠═b889e691-32d9-4f33-b9c9-52d18e0113a2
# ╟─18d66398-b853-4837-ae3a-1eed48951b63
# ╠═b37914b7-5f7b-446d-a336-d5370db69a19
# ╠═ecc4efc7-1d89-4b5a-8cdc-57055e2ac5bb
# ╟─bea93409-5d31-4050-920e-bdc99d2cc2ee
# ╠═349c30e2-cfa4-49de-a691-40cca1f9ef3a
# ╟─bdac9817-e65e-4a31-b85f-bd5c5e559757
# ╠═0d212db3-2ac1-40e6-ac77-617e10cda336
# ╠═13355e2a-8c46-4892-9300-51439ff0bb17
