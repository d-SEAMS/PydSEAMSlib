# PydSEAMSlib: python bindings for D-SEAMS

D-SEAMS is a tool for the analysis of molecular dynamics trajectories which is specifically able to qualitatively classify ice structures in both strong-confinement and bulk systems. PydSEAMSlib makes D-SEAMS user friendly.

# Installation for windows/ubuntu
Assuming access has been granted to the Github repository:

```{code-block} bash
git clone https://github.com/d-SEAMS/PydSEAMSlib.git
cd PydSEAMSlib
```
We provide a `conda-lock` environment:

```{code-block} bash
micromamba create -f conda-lock.yml -n pyseamsdev
micromamba activate pyseamsdev
```

# Installation with dockerfile

```{code-block} bash
git clone https://github.com/RuhiRG/D-seams-docker-recipes.git
docker run -it dseamsdockerrecipes:latest /usr/bin/bash
git clone https://github.com/d-SEAMS/PydSEAMSlib.git
cd PydSEAMSlib
```



## documentation for PydSEAMSlib-bindings 

```{toctree}
:maxdepth: 2
:caption: Contents

```

```{eval-rst}
.. automodule:: pydseamslib.cyoda
```

<!-- ```{autodoc2-summary}
pydseamslib.cyoda
```

```{eval-rst}

``` -->

## Indices and tables

- [](genindex)
- [](modindex)
- [](search)
