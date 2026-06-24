# β-drift projection test — does the testbed bias control landfall cross-track?

Status: 2026-06-23. This is the **item-0 gating analysis** for the structural-init pivot.
It is the cheap, no-new-code test the six-storm note flagged as "the single most valuable next
analysis" ([SIX_STORM_VALIDATION_NOTE.md], lines ~100–106) and that
[STRUCTURAL_INIT_RECIPE.md] walks past into a code-heavy experiment. It must be answered
**before** any of the red-team hold items, because it decides whether structural-init is worth
building at all.

## Question

If the characterized testbed β-gyre aim error were the **dominant** source of landfall track
error, then a single fixed bias vector, projected through each storm's landfall geometry, should
reproduce the six observed along/cross splits. Does it?

## Method

Take the model's documented mature β-drift aim-error vector. For a storm moving with heading θ
(clockwise from due north), decompose the error velocity into:

- **cross-track** (positive = right of track): `C_vel = E_east·cosθ − E_north·sinθ`
- **along-track** (positive = ahead of obs): `A_vel = E_east·sinθ + E_north·cosθ`

Accumulated displacement over transit time T: `displacement_km = velocity_m/s × T_h × 3.6`.

## Inputs

**Bias vector — firm.** The cheatsheet's own characterization: the model's β-drift is a spurious
~NE drift, **≈ +1.0 m/s north + 1.0 m/s east** relative to canonical NW (equivalently, a westward
component of ~0.4 m/s vs a canonical ~1.4 m/s — a ~1.0 m/s eastward deficit). So
**E = (+1.0 east, +1.0 north) m/s**. The eastward component (the cross-track driver for poleward
movers) is robustly positive ~0.6–1.0 across reasonable canonical-magnitude assumptions, so the
qualitative result does not hinge on the exact value.

**Observed along/cross — firm.** Control-baseline decompositions from [STRUCTURAL_INIT_RECIPE.md].

**Headings and transit times — confirmed by Justin (2026-06-23).** Headings from the track-character
column in [SIX_STORM_VALIDATION_NOTE.md]; transit times firm for Hugo/Katrina/Ivan, confirmed for
Fran/Michael/Laura.

| Storm   | heading θ        | transit T (h) | observed cross | observed along |
| ------- | ---------------- | ------------- | -------------- | -------------- |
| Hugo    | ~325° (NNW–NW)   | 28            | +110           | +23            |
| Katrina | ~350° (N)        | 35            | +125           | +77            |
| Ivan    | ~340° (NNW)      | 42            | +126           | +249           |
| Fran    | ~335° (NNW,fast) | 28            | +8             | −45            |
| Michael | ~010° (N→NE)     | 36            | −99            | +124           |
| Laura   | ~350° (→N)       | 30            | −32            | +37            |

## Result

Predicted from the single bias vector E = (+1, +1) m/s vs. observed (km):

| Storm   | **pred cross** | obs cross | cross sign | **pred along** | obs along | along sign |
| ------- | -------------- | --------- | ---------- | -------------- | --------- | ---------- |
| Hugo    | **+140**       | +110      | ✓          | **+25**        | +23       | ✓          |
| Katrina | **+146**       | +125      | ✓          | **+102**       | +77       | ✓          |
| Ivan    | **+194**       | +126      | ✓          | **+90**        | +249      | ✓          |
| Fran    | **+134**       | +8        | ✗ (~0)     | **+49**        | −45       | ✗          |
| Michael | **+105**       | −99       | ✗ (flip)   | **+150**       | +124      | ✓          |
| Laura   | **+125**       | −32       | ✗ (flip)   | **+88**        | +37       | ✓          |

## Findings

**1. The single bias vector predicts large EAST cross-track for all six storms (+105 to +194 km).**
All six head roughly poleward, and a NE bias projects rightward (east) for every poleward heading.
It matches the in-sample three within ~25%. It is catastrophically wrong for all three
out-of-sample: Fran (~0 observed), Michael and Laura (observed west).

**2. The clean falsifier — Katrina vs. Laura.** Both make landfall in/near Louisiana moving ~north
(θ ≈ 350°). A β-drift bias projected through the *same heading* must give the *same cross-track
sign*. Observed: Katrina **+125 east**, Laura **−32 west**. **No single geometry-projected vector
can produce opposite cross-track signs for the same heading.** This conclusion needs nothing from
the heading/transit estimates beyond "both moving roughly north," which both are. The factor that
flips the sign between them is their *different steering*, not their β-drift.
**⇒ Steering controls landfall cross-track; the β-drift bias is a subdominant term.**

**3. The bias is real — the test localizes where it lives.** It accounts for the in-sample three
cleanly, and it predicts the *along-track* of the poleward movers well: Hugo +25 vs +23,
Michael +150 vs +124. The Michael match confirms the six-storm note's reading that **Michael's −99
"west" is a metric artifact** — the poleward bias dumps into along-track overshoot, and the
same-latitude scalar reprojects it as spurious west. The bias is detectable, but as *along-track on
northward storms*, not as a universal eastward cross-track. (Independent argument for retiring the
same-latitude scalar in favor of along/cross.)

## The structural-init consequence (no new runs needed)

Structural-init broadens every storm's vortex (the six recipe fits give taper 238–314 km, all
broader than the 200 km default). The `gate-beta-renv` sweep already in hand shows what a broader
vortex does to the model's β-drift: R_env 400→800 km gives 2.36→2.95 m/s with heading 344°→353° —
i.e. westward component 0.65→0.36 (**less** west) and northward 2.27→2.93 (**more** poleward). So
broadening pushes every storm **more poleward / more eastward**:

- Hugo/Katrina/Ivan are already +110/+125/+126 too far east → structural-init makes them **worse**.
- It cannot manufacture Laura's or Michael's westward errors — those are steering, not size.
- The structural-init "success scenario" (broader vortex → more NW drift → recovers accuracy)
  assumes the model behaves like canonical β-drift. **The model's own sweeps say it does not** —
  broader is more poleward, not more westward.

This also largely discharges red-team hold item #1: the component-decomposed size response already
exists (gate-beta-renv + gate-beta-taper at Vmax 64; f-plane/structure probes give W≈const,
N∝Vmax). Decomposed, it says **broader → more east** — the wrong direction to fix the residual.

## Verdict

The bridge from the testbed β-gyre bias to landfall cross-track is **falsified as a dominant
mechanism**, and structural-init is **contraindicated by data already on the shelf**. Building the
EBTRK parser and the two-arm treatment would, in the best case, make the in-sample storms worse and
still miss the out-of-sample ones.

**Recommended pivot:** the two-contribution paper.
1. The over-rotating β-gyre as a characterized **model property** (testbed; survives f-plane null,
   intensity-independence, 2× grid convergence, the formulation probe).
2. **Out-of-sample track accuracy** with honest along/cross reporting — the bias described as
   *detectable in along-track on poleward movers, subdominant to steering for cross-track*, which is
   exactly what this test shows. No landfall-attribution bridge.

## Caveats

- The bias vector is a quiescent-environment testbed property; the projection assumes it accumulates
  linearly over transit with steering otherwise correct. That is the **strong** form of the
  hypothesis — which is the point: the strong form is what structural-init needs, and it fails.
- Findings 1 and 2 (predicts-east-for-all-six; Katrina/Laura same-heading-opposite-sign) are robust
  to the heading/transit inputs. The table *magnitudes* are not — they scale with |E| and T.
- Michael's −99 is largely a same-latitude-scalar artifact (see finding 3); Laura's −32 and Fran's
  ~0 are the cleaner falsifiers (a real westward cross-track, and a fully nulled bias, respectively).
