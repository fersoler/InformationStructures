# Informational Structures

Code to generate and measure Information Structures following the algorithm in Pablo Almaraz, Piotr Kalita, José A. Langa, Fernando Soler-Toscano, "Structural stability of invasion graphs for generalized Lotka--Volterra systems", [arXiv:2209.09802](https://arxiv.org/abs/2209.09802).

Some relevant files:

-   *R/ISbuild.R* is the file to build *Information Structures* (IS).
-   *R/ISgraph.R* is the file with functions for graphical representation of IS.
-   *R/ISmeasures.R* contains function to compute several measures on IS that have been used to [analyse brain dynamics](https://doi.org/10.1371/journal.pcbi.1010412).
-   *R/IGmeasures.R* contains functions to measure [*Invasion Graphs* (IG)](https://doi.org/10.1007/s00285-022-01815-2).

There are several examples of IS and IG creation and measurement. They can be found in the 'examples' directory:

-   *ISandIGexamples.R* contains several examples of creation and visualization of both IS and IG. IG measures are applied to examples. IS and IG are compared with the same connectivity parameters. Under certain stability conditions both structures coincide.
-   *exampleFiguresR* contains several visualizations of Information Fields, phase space, biodiversity cones, etc.

Several ShinyApps are provided to interact with IS. Specifically:

-   *InfoStructureDemo2D* is an app to interact with a system of 2 species. Information Structure, Information Field, biodiversity cones and phase space can be created with desired parameters.
-   *CooperativeSphere*, *CompetitiveSphere* and *CoopCompSphere* are apps to look at the inner structure of biodiversity cones of systems with 3 species.

A performance comparison between IS and IG creation is provided in *tests/invasionGraphs*. It is introduced a modification to Schreiber code when building communities.
