neurotic
========

Project to infer 3D structure from Open Connectome data

Requirements:

 - Minka's [Lightspeed](http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/) MATLAB toolbox
 - [vlfeat](http://www.vlfeat.org/)
 - [PGF and TikZ](http://sourceforge.net/projects/pgf/)


Files:

 - latex/ -- see submission.tex. 
 - src/ -- code for this project. cd to this directory and run `startup` to prepare the environment. In particular,
   - experiment/ -- main functions; entry point for running code
   - neal3_iter.m, neal8_iter.m -- Implementations of MCMC sampling based on Radford Neal's papers (details in comments)
   - MDP.m, MDPLight.m, NormalWishart.m, GammaGamma.m, etc. -- Class implementations of conjugate distributions supporting online inference.
   - GraphDist.m, GraphDistDistribution.m -- Support for Geodesic distances (details in paper)
 - src/scrape - scripts (mostly curl) to download data from openconnecto.me and compatible servers. Some scripts placeholders for secret keys. Please email the repository collaborators for details.
 - ddpcode/ - Dependence Dirichlet Process tracking module due to Willie Neiswanger.
 - data/ - Code expects data to live here. Scraping code puts it here.

