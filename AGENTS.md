# Developing this fork

This is the flow5 that [eccentricexit/uav](https://github.com/eccentricexit/uav) ships: a potential-flow solver behind a command line an
agent uses to design RC aircraft. The GUI stays buildable but is not the product; what matters is that `planerun` and `foilrun` in uav
return derived numbers or refuse with a reason. Work arrives as issues on uav, each with a reproducer; a change here lands by PR on this
repository and a pin bump there, where the tests live (`uav/tests/`, `pixi run test`, `pixi run validate`).

## Rules

- **Derived, never fitted.** A model is the physics with its basis cited where it is used, and a reference checks it, never tunes it. A
  constant chosen so that a case matches a reference hides the cause. Where no principled form is within reach, the method leaves the
  number out and says so in its output, and the issue records why.
- **Methods differ.** Hold a method to a derivation or a measurement, not to another method: two methods agreeing proves nothing, and
  a term kept so that they agree is a fit. Each method's limits are stated at the function that carries them.
- **Breaking changes are welcome** when they serve uav's vision: upstream compatibility, the GUI's conveniences and the `.fl5` export
  yield to a correct headless result. Say in the PR what breaks.
- **Every physical model is documented at its function**: what it computes, its basis (the equation, the paper, the AVL routine), its
  limit with the measurement that sets it, and the uav issue. Comments say what the code cannot; never what the code does.
- **Keep each file's line endings** (the sources mix CRLF and LF; `git diff --ignore-cr-at-eol` shows a real diff).
- **Measure before and after.** A PR that moves a number says which, on which fixtures, and why the new value is the derived one.
