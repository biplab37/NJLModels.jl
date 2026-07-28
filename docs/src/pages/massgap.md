# Massgap Calculations

The solution of the gap equations are provided via the function `massgap`.

## Mass gap NJL

First setup the parameters,
```@example massgap
using NJLModels, NJLModels.Zablocki
param = Zablocki.Parameters1()
```

The massgap can be calculated as,
```@example massgap
T = 0.01
mu = 0.0
m, ω₀ = massgap(T, mu, param).zero
```
which gives the constituent mass $m$, and the shift of the chemical potential $\omega_0$ coming from the vector interaction.

The `massgap` has multiple dispatch to also take an array as input

```@example massgap
trange = 0.05:0.001:0.3
mlist, wlist = massgap(trange, mu, param);
```

```@example massgap
using Plots
plot(trange, mlist, lab="", xlabel="T [GeV]", ylabel="m [GeV]")
```

## Massgap RDF
For the model used in my thesis,
```@example massgap
param2 = Zablocki.ParametersRDF4()
plp = Zablocki.PolyakovParameters(T0=0.208)
dump(plp)
```

```@example massgap
condensate_list, phi_list, phibar_list = Zablocki.massgap_rdf(trange, mu, param2, plp)
mlist_rdf = [Zablocki._get_mass(condensate, param2) for condensate in condensate_list];
dT = trange[2] - trange[1]
plot(
    plot(trange, [mlist_rdf/mlist_rdf[1] phi_list], lab=["m/m(T=0)" "Phi"], xlabel="T [GeV]", legend=:left),
    plot(trange[2:end], [-diff(mlist_rdf)/dT diff(phi_list)/dT], lab=["-dm/dT" "dPhi/dT"], xlabel="T [GeV]")
)
```
