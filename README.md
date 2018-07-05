# Stability

For a high-dimensional dataset 𝐗 with 𝑁 samples and 𝑃 dimensions (traits), 
dimensionality reduction techniques aim:
1. to provide a meaningful low-dimensional representation 𝐙 of 𝐾 dimensions while only losing minor amounts of information,
1. to use only a small number of free parameters
1. to preserve the quantities of interest in the data. 

There are a variety of approaches for dimensionality reduction with different underlying mathematical concepts and parameters
and choosing the most appropriate method for a given dataset is not trivial. Fundamentally, the problem is finding an 
objective criterion of what a good dimensionality reduction method is.

Here, I define a stability criterion that measures the number of components that can be reliably recovered 
in cross-validation and thus helps to determine the stable dimensions of the low-dimensional space.

