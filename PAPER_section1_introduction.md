# Draft section — §1 Introduction

*Drop-in for §1 per [PAPER_outline.md]. Manuscript register, BAMS narrative voice; author notes
at the end, plus a bonus draft of the required BAMS Significance Statement (the outline lists it
but it wasn't yet drafted — it falls out of the same material). Drafted 2026-07-01, synthesizing
the outline's abstract beats, §3 (PAPER_section3_clean_init_cascade.md), §2.1
(PAPER_section2_1_dynamical_core.md), the finished Results subsection, and
SIX_STORM_VALIDATION_NOTE.md. Main text ≈ 950 words (outline target 800–1000).*

---

## 1. Introduction

This paper reports two facts about one tropical-cyclone model — a characterized error and a
blind success — and the test that stopped us from claiming the obvious bridge between them.

The facts concern track. A tropical cyclone's motion is, to leading order, a superposition of
two influences: *environmental steering*, the storm carried along by the surrounding deep-layer
flow, and *self-propagation*, the storm pushing itself. The classical self-propagation mechanism
is β-drift: a vortex on the rotating Earth organizes a pair of counter-rotating gyres out of the
planetary-vorticity gradient, and the flow between those gyres propels the storm — canonically
toward the northwest, at 1–3 m s⁻¹, for a Northern Hemisphere cyclone [cite: Chan 2005; Holland
1983; Chan and Williams 1987; Fiorino and Elsberry 1989]. Neither influence is directly
observable in a track; what is observable, at the moment of landfall, is their sum. Throughout
this paper we decompose track error into components along and across the observed track: the
along-track component is a timing error — the storm arrives early or late — while the
cross-track component displaces landfall itself, and with it every downstream judgment that
depends on where a storm comes ashore.

That decomposition names an attribution problem it cannot by itself solve. A landfall error is a
convolution of self-propagation error, steering error, and whatever track features the model
cannot represent at all; landfall data alone cannot separate them, and a compensating pair of
errors in any two produces a small landfall miss for free. The standard development practice —
adjust the model until the hindcast matches — is, seen from this angle, a machine for
manufacturing exactly such cancellations. Each tuning step is accepted because it reduces the
error against the cases in hand; nothing in that procedure distinguishes a model that is right
from a model whose errors have been arranged to cancel, and the arrangement fails silently on
the first storm that breaks its geometry. Section 3 argues that this hazard is not hypothetical:
we document a cascade in our own development in which seven successive errors — in the physics,
the configuration, the input data, and finally our measurement instruments — each hid behind a
result that looked good.

Our response is a discipline we call *clean initialization*. Every storm is initialized directly
from the HURDAT2 best track, with the initialization values verified against the primary source
file. Every storm runs under one storm-agnostic configuration. The configuration's single
structural free parameter — the radius at which the initial vortex's outer wind profile begins
to taper — is calibrated in a quiescent-environment testbed against the canonical β-drift range
from the published literature, never against any landfall. And residuals are reported as
findings rather than absorbed into the configuration. The cost of this discipline is that the
model misses visibly, including in ways a tuned model would not. The payoff is that its misses
mean something — and so do its hits.

The two facts are as follows. First, clean initialization exposed a structural error worth
characterizing in its own right: the model's β-gyres *over-rotate*. In the isolated testbed the
simulated β-drift has canonical magnitude but is aimed nearly due poleward — its westward
component is a fraction of the canonical value — and the deficit survives every lever that could
plausibly control it: outer-vortex structure, diffusion, a doubling of horizontal resolution,
and an f-plane null test that rules out numerical and tracker artifacts (Section 4.1). It is a
bounded, reproducible, mechanistically located bias in self-propagation aim, and we report it as
a model property rather than removing it. Second, the same configuration exhibits blind track
skill: six historical Atlantic landfalls — Hugo (1989), Katrina (2005), and Ivan (2004), on
which the residual was first characterized, and Fran (1996), Michael (2018), and Laura (2020),
which the configuration had never run — are reproduced with cross-track errors ranging from a
few tens of kilometers to roughly 125 km, with no track-fitted parameter anywhere in the system
(Section 4.2).

The obvious paper connects these facts: the discovery storms all landed east of their observed
tracks by a strikingly consistent ~120 km, and an over-rotating β-gyre pushes poleward-moving
storms east. We drafted that paper. It was wrong, and the manner of its failure is the third
thing this paper reports. The three test storms broke the cluster — two of them missed *west* —
and a projection test showed that no single self-propagation bias vector, pushed through each
storm's landfall geometry, can reproduce the six observed errors: two storms approaching on the
same heading missed on opposite sides. Steering controls landfall cross-track in these cases;
the β-gyre bias is real, detectable — it predicts one test storm's along-track error almost
exactly — but subdominant (Section 4.2). We therefore present the two results as two results,
with the honest relation between them, and we present the methodology that forced that honesty —
the compensating-errors cascade, the registered predictions, the claims retired by our own
probes before review — as a contribution in its own right (Section 3), of particular relevance
as model development workflows become increasingly automated.

The remainder of the paper is organized as follows. Section 2 describes the model — a
nonhydrostatic anelastic core run, deliberately, in a reduced barotropic configuration — and the
experimental design: initialization, the single calibrated parameter, steering, and track
verification. Section 3 presents the clean-initialization methodology and the cascade of
compensating errors it exposed. Section 4 presents the two results: the over-rotating β-gyre as
a characterized model property (4.1) and the six-storm landfall record with its along/cross
attribution (4.2). Section 5 discusses what the over-rotation is and is not, the limits of
landfall attribution, and what blind skill under a no-tuning discipline does and does not
demonstrate. Section 6 concludes.

---

### BONUS — Significance Statement draft (BAMS requirement, ~115 words)

> Weather models are usually adjusted until they reproduce past events — a practice that can
> hide errors that happen to cancel. We built a hurricane track model that is never tuned to
> observed tracks: it is initialized from the official hurricane database, and its one
> structural setting is calibrated against textbook vortex physics. Run on six historical
> hurricanes — three of them never used during development — it reproduced landfall positions
> without any storm-specific adjustment. Just as important, the discipline exposed a specific,
> measurable error in how the model's storms steer themselves, and a chain of hidden,
> compensating errors in our own physics, data, and measurement tools that conventional tuning
> would have buried. We offer the approach as a template for trustworthy model evaluation.

---

### Author notes (delete before submission)

- **The opening line** ("...and the test that stopped us from claiming the obvious bridge
  between them") is a deliberate hook — BAMS tolerates and rewards this register, but it sets a
  self-critical tone for the whole paper. If you want a conventional opening, swap paragraphs 1→2
  (start with the superposition) and let the two-facts framing arrive in paragraph 5. I
  recommend keeping it: it is the paper's thesis in one sentence, and it matches title
  candidate #1 ("Catching Our Own Errors...").
- **"We drafted that paper. It was wrong."** — strongest sentence in the section; verify you're
  comfortable owning it in print. It is accurate (the eastward-cluster + bridge claim was the
  manuscript's centerpiece before the six-storm validation and projection test).
- **Numbers used:** cross-track range "few tens of km to ~125 km" spans the six same-latitude
  values (+103/+114/+68/+23/−75/−39) and the landfall-fix cluster (~+120). "Two of them missed
  west" = Michael (−75, timing-dominated per the metric note) and Laura (−39). "Predicts one
  test storm's along-track error almost exactly" = Michael, +123 predicted vs +124 observed
  (BETA_DRIFT_PROJECTION_TEST.md). "Same heading, opposite sides" = Katrina +125 / Laura −32 at
  θ ≈ 350°. All consistent with §3 and §4.2 usage.
- **"increasingly automated" clause** (end of the two-results paragraph) is the light-touch
  version of the AI-assistance point — same dial as §3's chapter-6 author note. Delete the
  clause if you want §1 agnostic; expand it if you decide §3 names the workflow explicitly. The
  two sections should move together.
- **[cite] block** in paragraph 2: same four references §4.1 already carries (Chan 2005 review +
  Holland 1983, Chan & Williams 1987, Fiorino & Elsberry 1989) — vol/pages still to verify per
  the outline's open-items list. No new citations introduced.
- **Deliberately NOT in §1:** any statement of the ~8–10° aim numbers, the 0.4 m s⁻¹ westward
  component, or per-storm tables — §4 owns the numbers; §1 speaks in mechanism and magnitude
  class only ("canonical magnitude, nearly poleward aim," "a fraction of the canonical value").
  This keeps the intro broad-audience and avoids double-anchoring figures that §4 must own.
- **Roadmap paragraph** mentions the barotropic reduction ("deliberately") to pre-empt the
  methods-section surprise for expert readers — coordinated with §2.1's reviewer-critical
  paragraph.
- **Length:** ≈ 950 words, top of the 800–1000 target. Trim candidates: the roadmap paragraph
  can lose the §5–6 sentence detail; paragraph 2's final sentence (the §3 forward-reference)
  duplicates a beat the two-results paragraph also hits.
