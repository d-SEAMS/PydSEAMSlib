# PydSEAMSlib: python bindings for d-SEAMS

[d-SEAMS](https://dseams.info/) (Deferred Structural Elucidation Analysis for
Molecular Simulations) is a tool for the analysis of molecular dynamics
trajectories which is specifically able to qualitatively classify ice structures
in both strong-confinement and bulk systems {cite:t}`idx-goswamiDSEAMSDeferredStructural2020`. ``PydSEAMSlib`` makes d-SEAMS user
friendly and accessible to the larger Scientific Python ecosystem.

```{note}
Historically, in 2023, the project was initiated under the name ``PySeams``, which could not be registered on PyPI due to naming similarities, so the project in 2024 switched to its current name, ``PydSEAMSlib``.
```

# Installation 

```{code-block} sh
git clone https://github.com/d-SEAMS/PydSEAMSlib.git
cd PydSEAMSlib
```

We provide a ``conda-lock`` environment:

```{code-block} sh
micromamba create -f conda-lock.yml -n pyseamsdev
micromamba activate pyseamsdev
```

Which can then be used directly for a ``pip`` install:

```{code-block} sh
# Pure bindings
pip install .
# For ASE integration
pip install .[adapters]
# With tests
pip install .[testing]
# Everything
pip install .[testing,adapters,docs]
```

Development local builds can be prepared via ``meson``:

```{code-block} sh
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib
meson compile install bbdir
```


## ``docker`` installation for emulated hardware

For ``arm64v8`` machines, it may be easier to use the following ``docker`` image.


```{code-block} sh
git clone https://github.com/RuhiRG/D-seams-docker-recipes.git
docker run -it dseamsdockerrecipes:latest /usr/bin/bash
git clone https://github.com/d-SEAMS/PydSEAMSlib.git
cd PydSEAMSlib
```

Subsequently, the earlier steps can be followed.

# License

MIT. For more details, including contribution guidelines visit the GitHub repo.


## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: IDX_
keyprefix: idx-
---
```


## Documentation

```{toctree}
:maxdepth: 2
:caption: API Documentation

api
history
```

## Indices and tables

- [](genindex)
- [](modindex)
- [](search)
