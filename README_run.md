## start julia and activate the project environment
```
]activate .  # use the current directoty as the project environment
]st
instantiate
```
In this environment the files Project.toml and Manifest.toml live that specify what packages and what versions of those packages are available in your environment. However, the default environment (@v1.9) is still active.

```
using Silico
```

## Add the package by running the command add PackageName.
```
]add Silico
```

## Running the demos
```
include("examples/activate.jl")
include("examples/demo/3d_polytope_drop.jl")
```