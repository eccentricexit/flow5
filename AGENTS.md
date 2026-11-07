# Developing this fork

This fork of flow5 serves a headless use: a potential-flow solver behind a command line, driven by an agent that designs RC
aircraft. The GUI stays buildable but is not the product; what matters is that the solver returns derived numbers or
refuses with a reason. Work lands by pull request on this repository.

## Principles

- Correct, never matched. A model is the physics. AVL, a wind tunnel, a lattice or another method is an indicator that
  something may be wrong, never a target, and not assumed right itself. A constant chosen so that a case matches a
  reference hides the cause; where no principled form is within reach, the number is left out and the output says so.
- Converge, then agree. A result converges with refinement. Two methods that model the same physics agree where neither
  leaves a term out; a difference is the term one has and the other lacks, in size and in sign, or a defect to find.
- A fix is for a wrong number returned without a word, two of the solver's own quantities that disagree where they
  should not, or a result that does not converge. A change of what the solver models says so.
- A PR says what it changes and what it measured before and after, on its own: the physics, the measurement, the
  decision, with no reference to any other repository. A pointer, where one is needed, is to this repository's own PR.
- Say less. A file says what is true now; the commit message says why it changed. A comment only where the code cannot
  say it; no histories in files.
- Breaking changes are welcome when they serve the headless use: upstream compatibility, the GUI's conveniences and the
  `.fl5` export yield to a correct headless result. Say what breaks.
- Keep each file's line endings: the sources mix CRLF and LF.
