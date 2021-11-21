# GeodesicTracer.jl

<a href="http://www.star.bris.ac.uk/fbaker/docs/GeodesicTracer.jl/">
<img alt="Docs" src="https://img.shields.io/badge/docs-dev-blue.svg"/>
</a>

## Setup

Install with

```julia
julia>] add "https://github.com/astro-group-bristol/GeodesicTracer.jl"
```

The package currently supports serial and parallel threaded CPU computation, and is designed to take advantage of as many resources as the Julia environment is permitted.

To change the number of threads, simply start a new Julia environment with the desired number, e.g.
```bash
julia --threads=8
```

## Documentation

Documentation hosted on [the Bristol Astrophysics Server](http://www.star.bris.ac.uk/fbaker/docs/GeodesicTracer.jl/).

To build the documentation locally
```bash
julia --project=. docs/make.jl
```
which will create and populate the directory `docs/build`, containing a static webpage with auto documentation.

Note, this can take a while, since `Documenter.jl` needs to load full Julia environment to properly inspect the code.