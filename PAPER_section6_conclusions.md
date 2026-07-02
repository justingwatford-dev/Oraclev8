# Draft section — §6 Conclusions

*Drop-in for §6 per [PAPER_outline.md] (~400-word target; this draft ≈ 430). Drafted 2026-07-01.
The outline's beat: "Two separate facts + no bridge; the methodological lesson." Deliberately
echoes §1's opening frame and closes the loop §3 opened. Author notes at the end.*

---

## 6. Conclusions

We set out to characterize a tropical-cyclone model's track error without giving the model any
opportunity to flatter us, and we report what survived. Two facts stand. First, the model's
β-gyres over-rotate: in isolation its storms propel themselves with canonical β-drift magnitude
but nearly poleward aim, a westward deficit that survives structure, diffusion, and resolution
levers, an f-plane null, and a validated measurement instrument. It is a bounded, reproducible,
mechanistically located model property, and we have left it in the model, characterized rather
than concealed. Second, the same configuration — one set of coefficients, one structural
parameter calibrated to idealized-vortex physics, nothing anywhere adjusted against a landfall —
reproduces six historical landfalls, including three storms it had never run, with cross-track
errors from a few tens of kilometers to roughly 125 km.

The relation between these facts is the caution at the center of this paper. The obvious
manuscript connects them — a poleward-biased self-propagation pushes poleward-moving storms
east, and our first three storms all missed east by a strikingly consistent margin. That
manuscript was drafted, and it was wrong: the three test storms broke the cluster, and a
single projection calculation showed that no fixed self-propagation bias can reproduce six
landfall geometries in which two storms on the same heading miss on opposite sides. Steering
controls landfall cross-track in these cases; the characterized bias is detectable — it
predicts the along-track error of the slowest poleward mover almost exactly — but subdominant.
We publish the two facts separately because that is what the evidence supports, and we publish
the falsification because it is the most instructive result we have.

That instruction is the third contribution. Under clean initialization, every flattering number
we ever produced eventually dissolved into a pair of canceling errors — seven times, in strata
descending from model physics through configuration and input data to the measurement
instruments themselves — and every dissolution was forced by a cheap, prespecified test:
registered prediction bands, primary-source verification, null controls, a Galilean instrument
check. None of this required new theory or new computing; it required only that the model never
be granted the benefit of the doubt. A model whose errors are characterized, bounded, and
honestly placed is more useful — to forecasters weighing its guidance, to developers hunting its
next error, and to anyone deciding how far to trust it — than a model whose errors have been
arranged to cancel. Six hurricanes suggest that standard is attainable without sacrificing
skill. We offer the discipline, as much as the model, as this paper's contribution.

---

### Author notes (delete before submission)

- **Echoes, by design:** "characterized rather than concealed" ↔ §4's closing; "That manuscript
  was drafted, and it was wrong" ↔ §1's "We drafted that paper. It was wrong." (same admission,
  past-tense callback — if §1's sentence gets softened, soften this one identically);
  "arranged to cancel" ↔ §1 ¶3. The repetitions are rhetorical bookends, not accidents.
- **"the slowest poleward mover"** = Michael (+123 predicted vs +124 observed along-track,
  BETA_DRIFT_PROJECTION_TEST.md). Kept nameless here to stay at conclusion altitude; §4.2 owns
  the storm-level detail.
- **"a few tens of kilometers to roughly 125 km"** — same span quoted in §1; the two must stay
  in sync if either is revised.
- **No future-work list here** — the outline puts more storms / second basin / deferred
  pressure-Holland arm in §5. If §5 ends up not carrying a future-work paragraph, one sentence
  could be appended here ("The obvious next tests — more storms, a second basin — are now
  mechanical."), but resist more than that; the 400-word cap is what gives this section its
  force.
- **Length:** ≈ 430 words vs the outline's ~400 — close enough; the trim candidate is the
  penultimate sentence ("None of this required...") if BAMS asks for cuts.
