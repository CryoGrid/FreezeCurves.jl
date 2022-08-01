# Solvers

The `Solvers` sub-module provides an interface for performing non-linear root-finding, e.g. solving implicit freeze curve formulations or recovering temperature given enthalpy, on the soil freeze characteristic curves implemented in this package. Note that this is *not* to be confused with parameter inference which is handled by a separate sub-module.

```@meta
DocTestSetup = quote
    using FreezeCurves
    using FreezeCurves.Solvers
end
```

```@autodocs
Modules = [FreezeCurves.Solvers]
Private = false
Order = [:type, :function, :macro]
```
