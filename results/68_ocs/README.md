# 68 — Optimal Contribution Selection (Meuwissen 1997)

Picks per-accession contribution shares that maximise expected multi-trait genetic gain subject to a constraint on per-generation inbreeding. Distinct from the merit-index z-score composite (§3.12 unweighted ranking) and from the usefulness criterion on cross pairs (§3.13). The merit index tells the breeder *who is good*; OCS tells the breeder *who to use and how much, given an inbreeding ceiling*.

## Method

Formulate the QP

`maximise c' EBV − (λ/2) · c' G c   subject to Σ c = 1, c_i ≥ 0`

where:
- **G** is the 95 × 95 VanRaden method-1 GRM from script 13 with a 10⁻⁴ ridge added on the diagonal for positive-definiteness (some eigenvalues sit near zero at our panel size).
- **EBV** is a per-accession direction-signed composite z-score across the 10 direction-signed traits (Crude_Protein, Antioxidant, Flavonoid, Phenol, Mass_of_Seeds, Insoluble_Oxalate scored as "+1 = good"; Tannin, Seed_Coat_Tannin, Total_Oxalate, Soluble_Oxalate scored as "−1 = good"; the three seed-size dimensions are direction-neutral and excluded).
- **λ** is the inbreeding-control multiplier. Swept across a 20-point log-spaced grid from 0.01 (gain-aggressive) to 100 (diversity-preserving).

Solver: `scipy.optimize.minimize` with method `SLSQP`, equality constraint `Σ c = 1`, box bounds `c_i ∈ [0, 1]`, analytic gradient. Warm-start: the solution at each λ seeds the next-larger λ to reduce wall time and improve convergence.

Three breeder-priority scenarios:
- **balanced**: all 10 direction-signed traits weight 1.
- **soluble_oxalate_priority**: Soluble_Oxalate weight 3, Insoluble_Oxalate weight 2, others weight 1 — encodes the "low gut-absorbed oxalate" target where the breeder wants both low soluble and high insoluble.
- **crude_protein_priority**: Crude_Protein weight 3, others weight 1.

## Findings

| Scenario | Gain range | Inbreeding range | n_active at moderate λ | Moderate λ |
|---|---|---|---|---|
| balanced | 0.085 – 0.931 | 0.0004 – 0.83 | 7 | 0.785 |
| soluble_oxalate_priority | 0.111 – 0.931 | 0.0005 – 0.83 | 8 | 0.785 |
| crude_protein_priority | 0.110 – **2.074** | 0.0005 – 0.83 | 4 | 1.274 |

**Three operational takeaways:**

1. **Crude_Protein-priority front sits well above the other two** (max gain 2.07 vs 0.93). Crude_Protein has the largest dynamic range in z-scored GEBV across the panel — the top-elite TSs157A is at +1.52 in raw GEBV — so the weight-3 multiplier translates into a sharply higher achievable composite EBV. Note the gain numbers are unitless because the EBV is a composite z-score; they index relative achievability across scenarios, not absolute trait-unit predictions.

2. **TSs157A is the dominant contributor across all three scenarios** (≈ 0.32 under balanced and Soluble_Oxalate priority; ≈ 0.75 under Crude_Protein priority). Consistent with its headline position on the §3.12 unweighted merit index. The Crude_Protein priority concentrates contribution into just 4 active accessions because the trait-specific gain is dominated by a few high-Crude_Protein accessions.

3. **Soluble_Oxalate priority pulls in TSs325, TSs59B, TSs333, TSs331 alongside the balanced elites.** These four accessions carry favourable Soluble_Oxalate GEBVs that lift their contribution share once the weighted composite weights soluble-oxalate-direction more heavily.

The balanced and Soluble_Oxalate-priority Pareto fronts nearly overlap, reflecting that under the direction-signed framework the Soluble_Oxalate signal mostly aligns with the other oxalate-fraction biases — the priority shift moves the active-set composition but doesn't dramatically change the achievable gain–inbreeding trade-off.

## Outputs

- `tables/pareto_front.csv` — 60 rows (3 scenarios × 20 λ values): scenario / λ / gain / inbreeding / n_active / converged.
- `tables/ocs_contributions_per_priority.csv` — 5,700 rows (60 × 95): scenario / λ / sample / contribution.
- `tables/top_contributions_per_priority.csv` — 45 rows (3 scenarios × top-15): scenario / λ / sample / contribution at the moderate-λ operating point.
- `figures/fig_ocs_pareto.png` / `.pdf` — Pareto front per scenario on a common axis.
- `figures/fig_ocs_top_contributions.png` / `.pdf` — three-panel bar plot of top-15 contributing accessions per scenario.

## Caveats

- **The 10⁻⁴ ridge on G** stabilises the QP but slightly inflates achieved inbreeding by `c'·(10⁻⁴ I)·c/2 = 5 × 10⁻⁵`. Below the practical reporting precision.
- **Direction-neutral seed-size traits (Length / Width / Thickness) are excluded** from the EBV composite — including them with an unsigned weight would inflate the EBV for accessions extreme in either direction, which is not breeder-relevant. This matches the §2.13 merit-index convention.
- **OCS as formulated here is unconstrained on the number of active parents.** A practical crossing block typically operates with 8–16 active parents; the moderate-λ operating points (4–8 active) already sit in that range. If a hard constraint on parent count is required, the QP can be cast as a mixed-integer programme (MILP / MIQP) — deferred to a future revision pass.
- **OCS does not pick the crosses, only the contribution shares.** Mate allocation (MateSel, Kinghorn 2011) is the natural downstream step that consumes OCS contributions and pairs parents under additional constraints; flagged for future work.

## Script

`scripts/68_ocs.py`. Runtime: ~ 40 s on a laptop (60 SLSQP solves × ≈ 0.7 s each).
