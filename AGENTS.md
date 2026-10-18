# Developing this fork

This fork of flow5 serves a headless use: a potential-flow solver behind a command line that an agent drives to design RC
aircraft. The GUI stays buildable but is not the product; what matters is that the solver libraries return derived numbers
or refuse with a reason. Work lands by pull request on this repository; each PR says what it changes, on what basis, and
what it measured before and after.

## Rules

- **Derived, never fitted.** A model is the physics with its basis cited where it is used, and a reference checks it, never
  tunes it. A constant chosen so that a case matches a reference hides the cause. Where no principled form is within reach,
  the method leaves the number out and says so in its output.
- **Methods differ.** Hold a method to a derivation or a measurement, not to another method: two methods agreeing proves
  nothing, and a term kept so that they agree is a fit. Each method's limits are stated at the function that carries them.
- **Breaking changes are welcome** when they serve the headless use: upstream compatibility, the GUI's conveniences and the
  `.fl5` export yield to a correct headless result. Say in the PR what breaks.
- **Every physical model is documented at its function**: what it computes, its basis (the equation, the paper, the AVL
  routine), its limit with the measurement that sets it. Comments say what the code cannot; never what the code does.
- **Self-contained.** Everything here (code, comments, commit messages, PR and issue text) is readable on its own: the
  physics, the measurement and the decision, with no reference to any other repository or its issues. A pointer, where one
  is needed, is to this repository's own issue or PR.
- **Keep each file's line endings** (the sources mix CRLF and LF; `git diff --ignore-cr-at-eol` shows a real diff).
- **Measure before and after.** A PR that moves a number says which, on which fixtures, and why the new value is the
  derived one.
