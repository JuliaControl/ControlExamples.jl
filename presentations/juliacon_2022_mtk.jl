### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 39fde46a-f13d-11ec-13f9-1f07466e8243
begin
	using Pkg
	Pkg.activate()
	using Revise, ModelingToolkit, ControlSystems, PlutoUI, Plots, OrdinaryDiffEq, RobustAndOptimalControl, LinearAlgebra
end

# ╔═╡ 0877bb10-7b50-4804-b1f9-0232383f279b
using SymbolicControlSystems # Writing C-code requires SymPy

# ╔═╡ cfbf9a99-2baf-4a99-9995-815636e015af
using Kroki # For the block diagram

# ╔═╡ 0c90bc5f-9a0b-4e10-ae34-f7460805c8f0
html"<button onclick='present()'>present</button>"

# ╔═╡ 0abd1ff0-d13c-4850-a16d-40a89dc4fc06
md"""
#### TODO:
- Change call to jacobian calculation to `linearize` when available. 
- Incorporate operating point in StateSpace (PR pending)
"""

# ╔═╡ 9ea9b12b-a3e3-4f17-ab01-ac22a8034257
md"""
# Complete control design workflow--ControlSystems meets ModelingToolkit
How to go from a nonlinear system model to a nonlinear closed-loop system through linear control design
"""

# ╔═╡ 638e8fad-ac00-44b4-8daf-b594ad1801c8
Kroki.Diagram(:Mermaid, "
graph LR
    MTK[ModelingToolkit<br>model] -->|Linearize| CS[ControlSystems<br>model]
    CS --> D{Control<br>design}
    D --> V1{Linear<br>Verification}
    V1 -->|ODESys| MTK2[ModelingToolkit<br>model]
    MTK2 --> S{Nonlinear<br>simulation}
    V1 -->|Refine| D
    S --> V2{Nonlinear<br>Verification}
    V2 -->|Refine| D
    V2 --> C-code
    classDef cMTK fill:#f66,stroke:#333;
    classDef cCS fill:#66f,stroke:#333;
    class MTK,MTK2,S,V2 cMTK;
    class CS,D,V1 cCS;
")

# ╔═╡ 5da7140f-9042-4250-b3e1-dbc8673c2502
md"""
# Control design for a quadruple-tank system

In this example, we will implement an LQG controller with integral action for a nonlinear MIMO process.

The process we will consider is a quadruple tank, where two upper tanks feed into two lower tanks, depicted in the schematics below. The quad-tank process is a well-studied example in many multivariable control courses, this particular instance of the process is borrowed from the Lund University [introductory course on automatic control](https://control.lth.se/education/engineering-program/frtf05-automatic-control-basic-course-for-fipi/).

 $(PlutoUI.LocalResource("/home/fredrikb/.julia/dev/JuliaSimControls/docs/src/examples/quadtank.png"))

The process has a *cross coupling* between the tanks, governed by a parameters $\gamma_i$: The flows from the pumps are
divided according to the two parameters $γ_1 , γ_2 ∈ [0, 1]$. The flow to tank 1
is $γ_1 k_1u_1$ and the flow to tank 4 is $(1 - γ_1 )k_1u_1$. Tanks 2 and 3 behave symmetrically.

The dynamics are given by
```math
\begin{aligned}
\dot{h}_1 &= \dfrac{-a_1}{A_1   \sqrt{2g h_1}} + \dfrac{a_3}{A_1 \sqrt{2g h_3}} +     \dfrac{γ_1 k_1}{A_1   u_1} \\
\dot{h}_2 &= \dfrac{-a_2}{A_2   \sqrt{2g h_2}} + \dfrac{a_4}{A_2 \sqrt{2g h_4}} +     \dfrac{γ_2 k_2}{A_2   u_2} \\
\dot{h}_3 &= \dfrac{-a_3}{A_3 \sqrt{2g h_3}}                         + \dfrac{(1-γ_2) k_2}{A_3   u_2} \\
\dot{h}_4 &= \dfrac{-a_4}{A_4 \sqrt{2g h_4}}                          + \dfrac{(1-γ_1) k_1}{A_4   u_1}
\end{aligned}
```
where $h_i$ are the tank levels and $a_i, A_i$ are the cross-sectional areas of outlets and tanks respectively. For this system, if $0 \leq \gamma_1 + \gamma_2 < 1$ the system is *non-minimum phase*, i.e., it has a zero in the right half plane. 

In the examples below, we assume that we can only measure the levels of the two lower tanks, and need to use a state observer to estimate the levels of the upper tanks.

The interested reader can find more details on the quadruple-tank process from the manual provided [here (see "lab 2")](https://canvas.education.lu.se/courses/16044/pages/course-materials?module_item_id=486293), from where the example is taken.

We start by defining the dynamics.
"""

# ╔═╡ 786cc34f-7cc1-417c-abf6-eb111c2ec9de
md"Create an `ODESystem` from ModelingToolkit:"

# ╔═╡ 57f236b7-8476-4805-a881-b861674bfdc0
md"""
## Trim the system
Find a steady-state operating point by simulating for a long time:
"""

# ╔═╡ d58306fc-97c1-4e77-9b8f-3ede99e1c8c3
md"## Linearize the system:"

# ╔═╡ 01a93b57-c826-430c-b99a-aaea609b8386
md"Contruct a ControlSystems.jl `StateSpace` object"

# ╔═╡ ca22eddf-d003-4ca6-9e7c-ea8913a3e776
md"""
## LQG design with integral action
We add integral action by augmenting the system model with a model of a low-frequency disturbance
"""

# ╔═╡ 1b2be79e-7379-4e2b-8d2c-86fc36e1b7f6
md"""
We design the LQG controller by choosing cost (``Q``) and covariance (``R``) matrices and then calling `LQGProblem`:
"""

# ╔═╡ 3acdbd8f-4f8c-46bd-b5fa-ce878c269b8f
md"## Verify robustness
Verify robustness by inspecting the sensitivity functions. `gangoffourplot` is a convenience function that plots 4 common sensitivity functions, ``S`` denotes the standard sensitivity function and ``T`` denotes the complimentary sensitivity function."

# ╔═╡ b661521d-53fb-4c2b-9e8c-9c7e32329db4
md"## Linear simulation
Perform a time-domain simulation, the bottom two plots are the control signals"

# ╔═╡ a9170eff-7aa8-408b-9b86-af59ec4d2866
md"## Nonlinear simulation"

# ╔═╡ c9dabf55-22f9-49e7-80a2-9b68b1345de9
md"Extract the controller from the `LQGProblem`:"

# ╔═╡ d5d83e2f-6a92-4e76-9e50-926b027dcaa2
md"Create an `ODESystem` that can be used to simulate the controller with ModelingToolkit:"

# ╔═╡ 21967b68-5571-4812-b544-4c3d4413d14a
md"Set up the closed-loop simulation with a step load disturbance like above:"

# ╔═╡ 682ffca1-cb58-4da3-b8ea-5166e009336c
md"Perform the nonlinear simulation and compare with the linear simulation performed above:"

# ╔═╡ 60065f60-bfd0-4353-9a75-d4662d99ca0f
md"Check that the simulation did not cause negative tank levels, the model does not handle this case!"

# ╔═╡ 8766bf2e-75a4-484f-8a69-80b27363b235
md"""
### Write C-code for the controller
We make use of [SymbolicControlSystems.jl](https://github.com/JuliaControl/SymbolicControlSystems.jl) for C-code generation of the controller represented by a linear system. The code generation requires SymPy.jl
"""

# ╔═╡ d7fb4856-9dd7-4bd4-9e9e-faf5d4a5c3c6
md"We first discretize the continuous-time controller using the Tustin method:"

# ╔═╡ 47de7694-4245-47b4-a4fe-c2e98b286897
html"""<style>
main {
    max-width: 65%;
}
"""

# ╔═╡ db810395-88a6-4767-8f8c-2229a49a6372
PlutoUI.TableOfContents()

# ╔═╡ 9ec318db-9388-4c20-bbbe-34e9bd4d2e00
import ModelingToolkitStandardLibrary.Blocks

# ╔═╡ f64e62a1-7d14-4d50-ace5-8ad8c889b398
function round_coefficients(G)
	rounder(X) = round.(X, sigdigits=2)
	ss(rounder.(ssdata(G))..., G.timeevol)
end

# ╔═╡ 06451abe-6a31-4150-9afa-c6e1f11b9d21
function quadtank(h, u)
    k1, k2, g = 1.6, 1.6, 9.81
    A1 = A3 = A2 = A4 = 4.9
    a1, a3, a2, a4 = 0.03, 0.03, 0.03, 0.03
    γ1, γ2 = 0.2, 0.2

    ssqrt(x) = √(max(x, zero(x)) + 1e-5) # For numerical robustness at x = 0

    xd = [
        -a1 / A1 * ssqrt(2g * h[1]) + a3 / A1 * ssqrt(2g * h[3]) + γ1 * k1 / A1 * u[1]
        -a2 / A2 * ssqrt(2g * h[2]) + a4 / A2 * ssqrt(2g * h[4]) + γ2 * k2 / A2 * u[2]
        -a3 / A3 * ssqrt(2g * h[3]) + (1 - γ2) * k2 / A3 * u[2]
        -a4 / A4 * ssqrt(2g * h[4]) + (1 - γ1) * k1 / A4 * u[1]
    ]
end

# ╔═╡ 4a8ce5c0-a157-4c3a-b66a-351d34533084
begin
	u0 = 0.2 # Input at operating point
	@parameters t
	D = Differential(t)
	@variables h[1:4](t)=1
	@variables u[1:2](t)=u0 [input=true]
	@variables y[1:2](t)=0  [output=true]
	u = collect(u)
	h = collect(h)
	y = collect(y)
end;

# ╔═╡ b3df1faa-ba5e-447d-999a-f101a85da390
begin
	dh = quadtank(h, u);
	eqs = [
		collect(D.(h) .~ dh);
		y[1] ~ h[1]
		y[2] ~ h[2]
	]
	@named quadsys = ODESystem(eqs, t)
end;

# ╔═╡ 08416248-e523-4717-83b9-d8ebae22f21b
quadsys

# ╔═╡ 6dd1e043-49dd-4ec1-a748-9d17fa8108eb
begin
	ssys = structural_simplify(quadsys)
	prob = ODEProblem(ssys, Pair[], (0, 1500))
	sol = solve(prob, Rodas5());
end;

# ╔═╡ 840b037e-742b-4ff7-aecb-220688b53c58
h0 = sol[end][1:4] # Steady state

# ╔═╡ 3e4986f9-c4f8-4f1f-ae32-756cba3dd485
plot(sol, title="Steady state solution")

# ╔═╡ 71d096b4-1e02-49c7-8cc0-b2af871172da
states(quadsys)

# ╔═╡ ea39b8cf-1fed-41ec-962a-ea1571e3be91
AB = calculate_jacobian(quadsys) # TODO: replace this block with a call to linearize

# ╔═╡ d156bd7a-1a49-4aba-90aa-b6336789eb16
begin
	# TODO: replace this block with a call to linearize
	ABnum = substitute.(AB, Ref(Dict(collect(h .=> h0)..., collect(u .=> u0)...)))
	ABnum = ABnum .|> ModelingToolkit.value .|> Float64
	A = ABnum[1:4, 1:4] 
	B = ABnum[1:4, 5:6]
	(; A, B)
end

# ╔═╡ 5f784148-e381-4e1c-851c-361e6bcd4246
begin
	C = [1 0 0 0; 0 1 0 0] # We can measure the levels in tanks 1 and 2
	sys = ss(A,B,C,0)
end;

# ╔═╡ 7b7b031e-3c09-499d-be80-e9bbf227fbf6
augsys = add_low_frequency_disturbance(sys; ϵ=1e-6);

# ╔═╡ c3fcd41f-ffd8-4480-9215-3d52220954cd
lqg = let
	C = augsys.C
	Q1 = C'*I(2)*C                      # We penalize the output only
	Q2 = I(2)
	R1 = diagm([ones(4); 0.002ones(2)]) # The last two state correspond to the integral states.
	R2 = I(2)
	LQGProblem(augsys, Q1, Q2, R1, R2)
end;

# ╔═╡ c0c911ae-43f2-4ae4-8ffe-bca989127746
gangoffourplot(lqg, ylabel="", xlabel="", lab="", legend=:bottomright, titlefontsize=12, plot_titlefontsize=12, size=(800,800))

# ╔═╡ f3e93258-32e4-4ab1-8097-be612c7f784d
begin
	disturbance = (x, t) -> [-0.5(t ≥ 10), 0]
	Gcl  = [G_PS(lqg); -comp_sensitivity(lqg)] # -comp_sensitivity(prob) is the same as the transfer function from load disturbance to control signal
	res  = lsim(Gcl, disturbance, 0:5:1000, dtmin=1e-1, force_dtmin=true)
	plot(res, framestyle=:zerolines, ylab=["h1" "h2" "u1" "u2"], layout=4)
end

# ╔═╡ c4f7d1ed-4a72-40d6-9eab-d67f65774c48
controller,_ = balance_statespace(-observer_controller(lqg))

# ╔═╡ 8d7ce3b5-7a19-40a4-8630-ca6a6c4a06d1
# ╠═╡ show_logs = false
@named lqg_controller = Blocks.StateSpace(ssdata(controller)...);

# ╔═╡ c1c623e2-fa18-4c9d-8b0a-50ec82fff840
discrete_controller = c2d(controller, 1.0, :tustin);

# ╔═╡ 6d5381a4-5786-40c6-8dfc-f9fbaed139aa
# ╠═╡ show_logs = false
begin
	code = SymbolicControlSystems.ccode(discrete_controller
		|> modal_form |> first |> round_coefficients  # Rounding and modal form for nicer printing, don't do this in practice
	)
	Markdown.parse("```\n" * code * "\n```") # Make the code display nicely
end

# ╔═╡ 132f7bf8-20bf-4067-a091-c4f1b36d4a0c
# ╠═╡ show_logs = false
begin
	@named dist = Blocks.Step(start_time=10, height=-0.5, smooth=true)
	
	connections = [
		# TODO: add operating point to Blocks.StateSpace
		lqg_controller.output.u[1] + u0 + dist.output.u ~ quadsys.u[1]
		lqg_controller.output.u[2] + u0 ~ quadsys.u[2]
		collect(lqg_controller.input.u  .~ quadsys.y .- h0[1:2])
	]

	@named closed_loop = ODESystem(connections, t, systems=[lqg_controller, quadsys, dist])
end;

# ╔═╡ 771c9cb9-9dce-4800-9983-01348ed40548
begin
	prob2 = ODEProblem(structural_simplify(closed_loop), Pair[collect(quadsys.h .=> h0);], (0, 1000))
	sol2 = solve(prob2, Rodas5(), dtmax=0.1)
	plot(sol2, vars=[quadsys.h[1:2]..., lqg_controller.output.u...], layout=(2,2), ylab=["h1" "h2" "u1" "u2"], lab="nonlinear")
	hline!(h0[1:2]', l=(:black, :dash), primary=false)
	plot!(res.t, (res.y .+ [h0[1:2]; 0; 0])', sp=(1:4)', lab="linear", framestyle=:zerolines, legend=:bottomright)
end

# ╔═╡ 02b39d84-da82-4949-9baa-463c46a9c243
any(<(0), reduce(hcat, sol2[quadsys.h])) && @error "Negative tank heights during simulation!"

# ╔═╡ Cell order:
# ╟─0c90bc5f-9a0b-4e10-ae34-f7460805c8f0
# ╟─0abd1ff0-d13c-4850-a16d-40a89dc4fc06
# ╟─9ea9b12b-a3e3-4f17-ab01-ac22a8034257
# ╟─638e8fad-ac00-44b4-8daf-b594ad1801c8
# ╟─5da7140f-9042-4250-b3e1-dbc8673c2502
# ╟─786cc34f-7cc1-417c-abf6-eb111c2ec9de
# ╠═08416248-e523-4717-83b9-d8ebae22f21b
# ╟─57f236b7-8476-4805-a881-b861674bfdc0
# ╠═6dd1e043-49dd-4ec1-a748-9d17fa8108eb
# ╠═840b037e-742b-4ff7-aecb-220688b53c58
# ╠═3e4986f9-c4f8-4f1f-ae32-756cba3dd485
# ╟─d58306fc-97c1-4e77-9b8f-3ede99e1c8c3
# ╠═71d096b4-1e02-49c7-8cc0-b2af871172da
# ╠═ea39b8cf-1fed-41ec-962a-ea1571e3be91
# ╠═d156bd7a-1a49-4aba-90aa-b6336789eb16
# ╟─01a93b57-c826-430c-b99a-aaea609b8386
# ╠═5f784148-e381-4e1c-851c-361e6bcd4246
# ╟─ca22eddf-d003-4ca6-9e7c-ea8913a3e776
# ╠═7b7b031e-3c09-499d-be80-e9bbf227fbf6
# ╟─1b2be79e-7379-4e2b-8d2c-86fc36e1b7f6
# ╠═c3fcd41f-ffd8-4480-9215-3d52220954cd
# ╟─3acdbd8f-4f8c-46bd-b5fa-ce878c269b8f
# ╠═c0c911ae-43f2-4ae4-8ffe-bca989127746
# ╟─b661521d-53fb-4c2b-9e8c-9c7e32329db4
# ╠═f3e93258-32e4-4ab1-8097-be612c7f784d
# ╟─a9170eff-7aa8-408b-9b86-af59ec4d2866
# ╟─c9dabf55-22f9-49e7-80a2-9b68b1345de9
# ╠═c4f7d1ed-4a72-40d6-9eab-d67f65774c48
# ╟─d5d83e2f-6a92-4e76-9e50-926b027dcaa2
# ╠═8d7ce3b5-7a19-40a4-8630-ca6a6c4a06d1
# ╟─21967b68-5571-4812-b544-4c3d4413d14a
# ╠═132f7bf8-20bf-4067-a091-c4f1b36d4a0c
# ╟─682ffca1-cb58-4da3-b8ea-5166e009336c
# ╠═771c9cb9-9dce-4800-9983-01348ed40548
# ╟─60065f60-bfd0-4353-9a75-d4662d99ca0f
# ╠═02b39d84-da82-4949-9baa-463c46a9c243
# ╟─8766bf2e-75a4-484f-8a69-80b27363b235
# ╠═0877bb10-7b50-4804-b1f9-0232383f279b
# ╟─d7fb4856-9dd7-4bd4-9e9e-faf5d4a5c3c6
# ╠═c1c623e2-fa18-4c9d-8b0a-50ec82fff840
# ╠═6d5381a4-5786-40c6-8dfc-f9fbaed139aa
# ╠═47de7694-4245-47b4-a4fe-c2e98b286897
# ╠═db810395-88a6-4767-8f8c-2229a49a6372
# ╠═39fde46a-f13d-11ec-13f9-1f07466e8243
# ╠═9ec318db-9388-4c20-bbbe-34e9bd4d2e00
# ╠═cfbf9a99-2baf-4a99-9995-815636e015af
# ╠═f64e62a1-7d14-4d50-ace5-8ad8c889b398
# ╠═06451abe-6a31-4150-9afa-c6e1f11b9d21
# ╠═4a8ce5c0-a157-4c3a-b66a-351d34533084
# ╠═b3df1faa-ba5e-447d-999a-f101a85da390
