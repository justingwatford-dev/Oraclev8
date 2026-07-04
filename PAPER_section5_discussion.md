# Draft section — §5 Discussion

*Drop-in for §5 per [PAPER_outline.md] (target 1000–1200 words; this draft ≈ 1150). Drafted
2026-07-01 from the §4 results subsection, SIX_STORM_VALIDATION_NOTE.md (open questions),
CHEATSHEET (what's-next + known issues, incl. the DLM double-counting question), and the
red-team audit (hold item 2 — the annulus caveat §2.3 explicitly defers here, now discharged).
Register per the approved voice target: calm, calibrated, drama in the facts only. Author notes
at the end.*

---

## 5. Discussion

### 5.1 What the over-rotation is — and is not

The testbed diagnosis of Section 4.1 can be stated in one sentence: the model's β-gyres reach
the right amplitude and the wrong orientation, and the orientation never locks. Drift speed
saturates near the canonical magnitude while the drift heading precesses steadily poleward
through — and past — due north, at a rate insensitive to vortex intensity. In the canonical
picture, the gyre pair equilibrates where the vortex circulation's advection of the gyres
balances their generation by the β-effect, and the storm settles onto a steady northwestward
bearing [cite: Fiorino and Elsberry 1989; Chan and Williams 1987]. In our model the amplitude
side of that balance is achieved and the phase side is not: whatever arrests the cyclonic
winding of the gyre pair in nature — and in the models that reproduce it — is too weak here.
We know three things it is *not*. It is not an advection or measurement artifact (the f-plane
null and the validated instrument); it is not diffusive smearing (reducing ν₄ degrades the
solution without rotating the aim); and it is not under-resolution in any simple sense (the
heading is flat across a doubling of resolution). Three candidates remained, untested at the
time of this analysis, all cheap to test: the imposed outer-wind cutoff, which sits inside the
radial annulus where the gyres live, and whose sweeps move the aim only to a floor; the reduced
barotropic configuration, since a vortex with vertical structure disperses β-Rossby energy
differently than a barotropic one; and the specific shape — rather than the radius — of the
taper. A subsequent profile-family test has since implicated the first: replacing the
compact-support taper with a smooth outer envelope recovers canonical, phase-locked β-drift in
the testbed, and the mechanism and its six-storm validation are reported in a companion paper
[ref: in preparation].

### 5.2 What a landfall can attribute — and a test worth exporting

The six-storm record enforces a discipline on attribution that we suspect generalizes. A
landfall error convolves self-propagation error, steering error, and unrepresented track
features; landfall data alone cannot deconvolve them, and Section 4.2 shows how persuasive the
resulting illusions can be — three discovery storms produced a 16-km-tight eastward cluster
that was, nonetheless, not a property of the model. The instrument that caught this cost one
afternoon: take the candidate bias vector, project it through each storm's landfall geometry,
and demand it reproduce the held-out storms' along/cross splits. We commend the projection test
as a general gate. Any proposed model change that is motivated by a landfall residual — a
structural initialization scheme, in our case — can be evaluated on paper against held-out
geometry before a line of it is built. Ours failed the gate; the sweeps then showed the fix
would have moved every discovery storm the wrong way. The expensive experiment was cancelled by
the cheap one, which is the correct order.

### 5.3 How much of the skill belongs to the steering?

An honest accounting of the blind track skill must start with what the configuration is given:
ERA5 is a reanalysis, so the environmental flow the model relaxes toward is close to the truth.
The skill demonstrated here is therefore hindcast skill under near-true steering — evidence
about the *model's dynamics and the method*, not a claim of operational forecast value, where
steering must itself be forecast. Within that frame, the attribution of Section 4.2 says the
steering carries the cross-track outcome and the model's contribution is to integrate it
faithfully while adding a bounded, characterized self-propagation. Two caveats on the steering
pathway remain open, deferred here from Section 2.3. First, the 3–7° annulus average is not a
storm-removed environmental field; at those radii it retains some of the storm's own outer
circulation, and the sensitivity of track to annulus width — or to storm-removed or
time-lagged steering — has not been tested. Second, and more subtly, a deep-layer mean drawn
from reanalysis in a storm-centered annulus may already contain part of the storm's
*propagation*, not just its environment; if so, the model's own β-drift partially
double-counts it, and the symptom would be systematic poleward over-translation on storms
whose DLM tracks their motion — a signature at least consistent with our largest timing miss
(Ivan, −8 h). Both are bracketing experiments, not rebuilds, and both belong ahead of any
larger claim for the skill.

### 5.4 One model, or a class?

We do not know whether the over-rotating β-gyre is a property of this model or of a family of
models that share its ingredients — soundproof cores, bounded vortices, hyperdiffusive
closures, or simply β-gyres asked to equilibrate at 15-km resolution. The diagnostic, however,
is portable and cheap: a quiescent β-plane, a balanced vortex, a mature-window drift vector,
and an f-plane null. Any track model can be asked the same question in an afternoon of compute,
and the answer is a vector, not an interpretation. We would be glad to learn whether other
cores over-rotate; either answer is informative, and the comparison would locate the mechanism
faster than we can from inside one model.

### 5.5 What blind skill demonstrates, and what would strengthen it

Six storms — three of them never run during development — are enough to *retire* a false claim,
which is the use we made of them; they are thin support for a strong positive claim, and we
state the skill accordingly. Strengthening it is now mechanical rather than conceptual: the
unified pipeline reduces a new storm to a registry entry and a reanalysis download, and the
storms that would stress the method most are known from the present six — more fast
straight-movers and slow recurvers, and a second basin with different climatological steering.
Beyond track, the core's thermodynamic pathway — validated but switched off throughout this
paper — is the natural frontier: with buoyancy active, intensity becomes a prognostic quantity
rather than a capped one, and the clean-initialization discipline would face a harder test,
since intensity error and track error feed each other through exactly the kind of compensation
this paper is about. We intend to arrive at that test with the same rule we arrived here with:
the model is never granted the benefit of the doubt.

---

### Author notes (delete before submission)

- **§2.3's deferred caveat is discharged in 5.3** (annulus ≠ storm-removed field + the
  double-counting question) — red-team hold item 2 is thereby represented in the paper as an
  open, named limitation rather than omitted. The Ivan −8 h linkage is phrased as "at least
  consistent with," not attribution — keep it that soft; the cheatsheet lists the DLM
  double-counting as an untested hypothesis.
- **5.1 candidate list** is deliberately three items. The taper-inside-
  the-gyre-annulus observation: gyres live roughly 150–450 km from center (gyre instrumentation
  read 75–450 km m=1), taper onset is 200 km — verify the overlap phrasing against the gyre
  figure before print. The barotropic/baroclinic dispersion candidate is stated without
  citation; if you want one, the vertical-structure β-drift literature (e.g., Wang & Holland
  mid-90s JAS) is the place, but verify before citing — I have not confirmed those papers'
  availability or exact claims.
- **Forward reference added 2026-07-03** (closing sentence of 5.1): the over-rotation campaign
  (OVERROTATION_CANDIDATES.md, all stages committed) implicated the compact-support cutoff and
  validated the envelope fix across six storms — the companion paper is outlined in
  PAPER2_outline.md. Replace "[ref: in preparation]" with the real citation when it exists;
  if paper 2's status changes, this sentence and the abstract/§6 need no edits (they never
  claimed the candidates were open — only §5.1 did).
- **"16-km-tight"** = the landfall-fix cluster span (110/125/126). §4.2 owns the table; this is
  the only number 5.2 repeats — drop it if it feels double-anchored.
- **5.3's ERA5-reanalysis honesty paragraph** is new framing not present elsewhere in the
  paper — check you're happy with "evidence about the model's dynamics and the method, not a
  claim of operational forecast value." It is the kind of sentence reviewers reward and
  press releases hate.
- **Deferred pressure-Holland arm** (outline's future-work list): NOT included — I couldn't
  ground what it refers to from the current repo docs. If it's the pressure-based Vmax/structure
  initialization variant, add one clause to 5.5's "storms that would stress the method" sentence;
  needs your one-line description.
- **5.5's closing sentence** intentionally reprises the paper's rule ("never granted the benefit
  of the doubt") — final echo before §6 uses its own variant. If §6's "arranged to cancel"
  family of echoes ever feels heavy in the assembled read, this is the first echo to cut.
- **Length** ≈ 1150 words (target 1000–1200). Trim order if needed: 5.4 (can compress to three
  sentences), then the 5.2 final two sentences.
