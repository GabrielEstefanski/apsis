# apsis-py-core

Capsule transport and duck-type extractors for apsis Python bindings.

Used by the `apsis` Python distribution (`crates/apsis-python`) and by
any external `apsis-plugin-X` cdylib that publishes operators
consumable by `apsis.System.add_hamiltonian_perturbation`.

## API

| function                  | purpose                                                                 |
|---------------------------|-------------------------------------------------------------------------|
| `box_into_capsule`        | Wrap `Box<dyn HamiltonianOperator>` into a `PyCapsule` (single-consume).|
| `take_box_from_capsule`   | Extract the boxed operator at the `add_hamiltonian_perturbation` boundary.|
| `extract_unit_system`     | Read an `apsis::UnitSystem` from a Python object's scale attributes.    |

## Plugin author template

External plugins publish operators by building a `Box<dyn HamiltonianOperator>`,
wrapping it in a capsule, and presenting it as an `apsis.Perturbation`
to the user. The wrap dance imports the `apsis` Python package by name,
which `apsis-py-core` does not do — copy this 8-line template into your
own cdylib:

```rust
use apsis::physics::integrator::HamiltonianOperator;
use apsis_py_core::box_into_capsule;
use pyo3::prelude::*;

pub fn wrap_in_apsis_perturbation(
    py: Python<'_>,
    inner: Box<dyn HamiltonianOperator>,
    label: &str,
) -> PyResult<PyObject> {
    let capsule = box_into_capsule(py, inner)?;
    let apsis = py.import("apsis")?;
    let perturbation_cls = apsis.getattr("Perturbation")?;
    Ok(perturbation_cls.call1((capsule, label))?.into())
}
```

Plugin's PyO3 factory then calls it:

```rust
#[pymethods]
impl PyMyForce {
    #[staticmethod]
    fn for_units(py: Python<'_>, units: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        let u = apsis_py_core::extract_unit_system(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(my_crate::MyForce::for_units(u)),
            "MyForce(for_units)",
        )
    }
}
```

User attaches the result to a `System` exactly like an internal
operator:

```python
import apsis
import apsis_plugin_myforce

sys = apsis.System(bodies=[...], units=apsis.units.SOLAR_CANONICAL, ...)
sys.add_hamiltonian_perturbation(
    apsis_plugin_myforce.MyForce.for_units(units=apsis.units.SOLAR_CANONICAL),
)
```

## Capsule version

`apsis_perturbation_box_v3`. Bumping the suffix is the breaking-change
marker; consumers must recompile against the new tag.

| version | change                                                                  |
|---------|-------------------------------------------------------------------------|
| `_v2`   | `accumulate` migrated from `&mut [(f64, f64)]` to `&mut [Vec3]` (3D).   |
| `_v3`   | Payload narrowed to `HamiltonianOperator`; non-conservative reserved.   |
