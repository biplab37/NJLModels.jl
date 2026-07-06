# Parameters
```@meta
CurrentModule = NJLModels
```

The parameters are specified via `Parameters` abstract type. Some specific parameter sets are already predifined in the module. Some examples are given below.

## Vanilla NJL model
```@example params
using NJLModels, NJLModels.Zablocki
Zablocki.Parameters1()
```

The parameters are mutable and one can use the constructor to create other parameter sets as well,
```@example params
Zablocki.Parameters1(eta_d=1.0)
```

The parameters can directly be passed into almost all the functions as an argument.

For example, the massgap in the NJL model with vector interaction can be calculated as,
```@example params
T = 0.01
mu = 0.0
param = Zablocki.Parameters1()
m, ω0 = massgap(T, mu, param).zero
```

## Polyakov Loop parameters
```@example params
dump(Zablocki.PolyakovParameters())
```

## Parameters for NJL model with RDF and Polyakov loop
The parameter set used for my PhD thesis,
```@example params
dump(Zablocki.ParametersRDF4())
```
