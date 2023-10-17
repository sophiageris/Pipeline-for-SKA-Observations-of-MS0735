# Pipeline-for-SKA-Observations-of-MS0735
A pipeline to simulate Square Kilometer Array Sunyaev Zel'dovich observations of the AGN cavities of galaxy cluster MS0735. 

A brief idea of the process is to make a model map of the SZ signal being recieved from the cluster; run the SKA simulation using 'Profile' software; then based on the output uv fits files create an input file of the model to use in McAdam where we will run Bayesian inference — specifically nested sampling via MultiNest (Feroz and Hobson 2008, Feroz et al. 2009, Feroz et al. 2019.) for parameter estimate and model comparison.

## Creating a simulated SZ signal map of cluster MS0735


## References
FEROZ, F., AND HOBSON, M. P. Multimodal nested sampling: an efficient and robust alternative to markov chain monte carlo methods for astronomical data analyses. Monthly Notices of the Royal Astronomical Society 384, 2 (jan 2008), 449– 463.

FEROZ, F., HOBSON, M. P., AND BRIDGES, M. MultiNest: an efficient and robust bayesian inference tool for cosmology and particle physics. Monthly Notices of the Royal Astronomical Society 398, 4 (oct 2009), 1601–1614.

FEROZ, F., HOBSON, M. P., CAMERON, E., AND PETTITT, A. N. Importance nested sampling and the MultiNest algorithm. The Open Journal of Astrophysics 2, 1 (nov 2019).
