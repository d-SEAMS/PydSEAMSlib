# Pyseams: python bindings for D-SEAMS

D-SEAMS is a tool for the analysis of molecular dynamics trajectories which is specifically able to qualitatively classify ice structures in both strong-confinement and bulk systems. Pyseams makes D-SEAMS user friendly.

# Installation for windows/ubuntu
Assuming access has been granted to the Github repository:

```{code-block} bash
git clone https://github.com/d-SEAMS/pyseams.git
cd pyseams
```
We provide a `conda` environment:

```{code-block} bash
micromamba create -f environment.yml
micromamba activate pyseamsdev
```

# Installation with dockerfile

```{code-block} bash
git clone https://github.com/RuhiRG/D-seams-docker-recipes.git
docker run -it dseamsdockerrecipes:latest /usr/bin/bash
git clone https://github.com/d-SEAMS/pyseams.git
cd pyseams
```



## Additional Topics

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/library_root






## Indices and tables

- [](genindex)
- [](modindex)
- [](search)