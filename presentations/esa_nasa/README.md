The notebooks that were presented were
- [`nasa_esa.jl`](https://github.com/JuliaControl/ControlExamples.jl/blob/master/presentations/esa_nasa/nasa_esa.jl)
- [`juliacon_2022_mtk.jl`](https://github.com/JuliaControl/ControlExamples.jl/blob/master/presentations/juliacon_2022_mtk.jl) (one level up)
- [ReachabilityAnalysis.jl Quadrotor example](https://juliareach.github.io/ReachabilityAnalysis.jl/dev/generated_examples/Quadrotor/)

The notebooks are available as rendered html files [here](https://drive.google.com/drive/folders/1U5T6-KQs_bhtmzSpTh_cMMSIdy0cXei6?usp=sharing)

To run a notebook, install Julia and the [Pluto.jl package](https://github.com/fonsp/Pluto.jl). Pluto is responsible for running the notebooks. Once you start one of the notebooks inside of pluto, the rest of the required packages will be automatically installed. The notebook [`juliacon_2022_mtk.jl`](https://github.com/JuliaControl/ControlExamples.jl/blob/master/presentations/juliacon_2022_mtk.jl) is an exception to this automatic package installation by Pluto, since it relies on a development version of [Kroki.jl](https://github.com/bauglir/Kroki.jl) that draws block diagrams. This notebook requires you to manually install all packages that are loaded inside the notebook. Once Kroki.jl makes a new release with the required updates, this procedure will be simplified.