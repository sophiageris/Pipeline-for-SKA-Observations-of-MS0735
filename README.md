# Pipeline-for-SKA-Observations-of-MS0735
A pipeline to simulate Square Kilometer Array Sunyaev Zel'dovich observations of the AGN cavities of galaxy cluster MS0735. 

A brief idea of the process: 1) make a model map of the SZ signal being recieved from the cluster, 2) run the SKA simulation using 'Profile' software (Grainge et al., 2002), 3) based on the output uv fits files, create an input file of the model to use in McAdam where we will run Bayesian inference — specifically nested sampling via MultiNest (Feroz and Hobson 2008, Feroz et al. 2009, Feroz et al. 2019.) for parameter estimate and model comparison.

## Creating a simulated SZ signal map of cluster MS0735
The purpose of this investigation is to observe the galaxy cluster cavities which are produced by the AGN jets, and are described by the suppression factor 'f' defined in Abdulla et al. (2019). For an in-depth description of the calculations behind the creation of a simulated y-map and subsequent signal map of a galaxy cluster, see section 4.2.1 of the Master's thesis which can be found here [Sophia Geris Master's Thesis](https://vuw-my.sharepoint.com/:b:/g/personal/gerisso_staff_vuw_ac_nz1/Ed5ZLI0h3r1DmTRMXnmXU3wBA6ukzjePt2zDtwWFgKZh9g?e=fAg5I5). This process is carried out using the Python script Create_SZ_signal_map_for_MS0735.py. the the 


## References
Z. Abdulla, J. E. Carlstrom, A. B. Mantz, D. P. Marrone, C. H. Greer, J. W. Lamb, E. M. Leitch, S. Muchovej, C. O’Donnell, T. J. Plagge, and D. Woody. Constraints on the Thermal Contents of the X-Ray Cavities of Cluster MS 0735.6+7421 with Sunyaev- Zel’dovich Effect Observations. The Astrophysical Journal, 871(2):195, 2019.

F. Feroz and M. P. Hobson. Multimodal nested sampling: an efficient and robust alter- native to markov chain monte carlo methods for astronomical data analyses. Monthly Notices of the Royal Astronomical Society, 384(2):449–463, 2008.

F. Feroz, M. P. Hobson, and M. Bridges. MultiNest: an efficient and robust bayesian inference tool for cosmology and particle physics. Monthly Notices of the Royal Astro- nomical Society, 398(4):1601–1614, 2009.

F. Feroz, M. P. Hobson, E. Cameron, and A. N. Pettitt. Importance nested sampling and the MultiNest algorithm. The Open Journal of Astrophysics, 2(1), 2019.

K. Grainge, M. E. Jones, G. Pooley, R. Saunders, A. Edge, W. F. Grainger, and R. Kneissl. Measuring the Hubble constant from Ryle Telescope and X-ray observations, with application to Abell 1413. Monthly Notices of the Royal Astronomical Society, 333(2): 318–326, June 2002.
