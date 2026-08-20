# Nonlinear FEM: Large-Strain Analysis of Solid and Porous Hyperelastic (Gent) Specimens

A Total Lagrangian, Newton–Raphson nonlinear finite element (FE) code, implemented from scratch in
MATLAB, for 2D plane-strain finite elasticity of a compressible **Gent hyperelastic** material. The
code is used to study the large-strain tensile response of solid and porous (five-hole) dog-bone
elastomeric specimens, using 4-node quadrilateral and 6-node triangular elements with both full and
selective-reduced Gauss integration, and is validated against the classical **Cook's membrane**
benchmark.

📄 **Full report:** [`report/report.pdf`](report/report.pdf) (built from [`report/report.tex`](report/report.tex))

---

## Overview

Elastomers are nearly incompressible and undergo large deformations, which makes their finite element
analysis prone to two classical difficulties: geometric/material nonlinearity (requiring an
incremental-iterative Newton–Raphson solution) and **volumetric locking** in low-order displacement
elements. This project:

- Derives the second Piola–Kirchhoff stress **S** and consistent material tangent **D** for the Gent
  strain-energy function, and implements them in a Total Lagrangian FE code.
- Implements two isoparametric element types: a 4-node bilinear quadrilateral (**4NR**) and a 6-node
  quadratic triangle (**6NT**), plus a 3-node linear triangle used only for validation.
- Implements both **full** Gauss integration and a **selective reduced integration (SRI)** scheme
  (under-integrating the volumetric/bulk part of the stress) to study and mitigate volumetric locking.
- Validates the code against the **Cook's membrane** benchmark (Total and Updated Lagrangian
  formulations), comparing against the mixed FE results of Brink & Stein (1996).
- Applies the validated code to a **solid** and a **porous** (five circular voids) dog-bone elastomer
  specimen under uniaxial tension, comparing stress fields, load–displacement response, and locking
  behaviour across all 8 combinations of geometry × element × integration scheme.

## Key Results

- The porous specimen extends more than the solid specimen under the same applied traction, and
  develops pronounced local stress concentrations at the hole boundaries.
- Volumetric locking is clearly observed for the **fully-integrated 4-node quadrilateral** element
  (most severely in the solid specimen), and is removed by selective reduced integration.
- The **6-node triangular** element is intrinsically much less prone to volumetric locking, with little
  difference between full and reduced integration.

## Repository Structure

```
.
├── Solid_4NR_FullIntegration/         # Solid specimen, 4-node quad, full integration
├── Solid_4NR_SelectiveReduced/        # Solid specimen, 4-node quad, selective reduced integration
├── Solid_6NT_FullIntegration/         # Solid specimen, 6-node triangle, full integration
├── Solid_6NT_SelectiveReduced/        # Solid specimen, 6-node triangle, selective reduced integration
├── Porous_4NR_FullIntegration/        # Porous specimen, 4-node quad, full integration
├── Porous_4NR_SelectiveReduced/       # Porous specimen, 4-node quad, selective reduced integration
├── Porous_6NT_FullIntegration/        # Porous specimen, 6-node triangle, full integration
├── Porous_6NT_SelectiveReduced/       # Porous specimen, 6-node triangle, selective reduced integration
├── Validation_CooksMembrane/          # Cook's membrane benchmark (3NT/4NR/6NT, TL & UL formulations)
│   └── Reference_article_for_validating_results.pdf   # Brink & Stein (1996)
├── Tapp_vs_Displacement/              # Combined load–displacement comparison plot (all 8 cases)
└── README.md
```

Each element/case folder contains:
| File | Description |
|---|---|
| `run_*.m` | Main MATLAB driver script (mesh input, material params, Newton–Raphson loop, post-processing) |
| `nodes.txt` | Nodal coordinates |
| `elements.txt` | Element connectivity |
| `fixed_nodes.txt` | Clamped (fixed-displacement) node list |
| `force_nodes.txt` | Nodes on the loaded (traction) edge |
| `sigma11.jpg`, `sigma22.jpg` | Stress contour plots on the deformed configuration |
| `mesh4n` / `mesh6n` | Mesh visualization file |

## Theory Summary

**Material model — compressible Gent hyperelasticity:**

$$
\Omega(\mathbf{F}) = -\frac{\mu\alpha}{2}\ln\left[1-\frac{I_1-3}{\alpha}\right] - \mu\ln J + \frac{G'}{2}(J-1)^2
$$

**Second Piola–Kirchhoff stress** (derived and implemented):

$$
\mathbf{S} = \frac{\mu\alpha}{\alpha - I_1 + 3}\,\mathbf{F} - \mu\,\mathbf{F}^{-T} + G'(J-1)J\,\mathbf{F}^{-T}
$$

with the consistent material tangent $\mathbb{D} = \partial\mathbf{S}/\partial\mathbf{F}$ used for the
Newton–Raphson linearization. Full derivations, the weak form, FE discretization, and the selective
reduced integration scheme are in the [report](report/report.pdf).

## Material Parameters Used

| Parameter | Value |
|---|---|
| Shear modulus, μ | 0.6667 MPa |
| Volumetric penalty modulus, G′ | 666.6667 MPa |
| Gent locking parameter, α | 53.3028 |
| Thickness | 1 (plane strain) |

## Running the Code

Each case is self-contained. To run one, e.g. the solid specimen with the 4-node element and full
integration:

```matlab
cd Solid_4NR_FullIntegration
run_solid_4NR_full   % or: run('run_solid_4NR_full.m')
```

This reads the mesh files in that folder, runs the load-stepped Newton–Raphson solution, and produces
the stress contour and load–displacement plots. Requires MATLAB (no additional toolboxes).
```

## Validation

The code was validated against the Cook's membrane benchmark (tapered, clamped, shear-loaded panel)
using a compressible neo-Hookean limit of the Gent model, comparing Total and Updated Lagrangian
formulations and 3-, 4-, and 6-node elements against the mixed finite-element results of:

> U. Brink and E. Stein, "On some mixed finite element methods for incompressible and nearly
> incompressible finite elasticity," *Computational Mechanics*, vol. 19, pp. 105–119, 1996.

## License

This project is released under the MIT License (see `LICENSE`).
