# Figure workflow (week 2)

Skills and agents for making and reviewing data figures, built from the rules, good and bad examples and guidelines supplied by the user.

- `skills/figure-rubric/` - shared numbered rules, 46 in all, covering chart style fit, unbiased presentation and typography (font sizes, titles, legends) (each with source slide, severity, and origin tag NOTES / SLIDE / ADDED; rules tagged ADDED await the user's confirmation) plus `kernel.py` helpers: `lint_script`, `audit_figure`, color, scale and gap helpers.
- `skills/figure-generation/` - workflow: raw data -> week-1 analysis checkpoints -> figure spec -> build -> self-lint -> automatic review (max 3 rounds).
- `agents/FIGURE_MAKER.md` - builds figures; sends every figure to the reviewer.
- `agents/FIGURE_REVIEWER.md` - read-only critic; returns a PASS/WARN/FAIL card per rule.

These are exported copies of the skills and agent profiles in the Claude Science account (the account versions are the live ones). Rules tagged ADDED are additions not stated in class and need student confirmation.

## How to run

Start the session as `FIGURE_MAKER` (not as a sub-agent of another agent) and give it the data and the question. It runs the analysis checkpoints, writes the figure spec, builds the figure, and sends each saved figure to `FIGURE_REVIEWER`. A sub-agent cannot dispatch another agent, so if the maker runs as a sub-agent it returns its figures labeled UNREVIEWED DRAFTS and the parent session must send them to the reviewer.

Tested on a simulated raw counts table (empty samples, a batch effect, three hidden subpopulations and a treatment effect): the checkpoints removed the empty samples, found the three groups and the treated genes, and the reviewer returned full cards with no FAIL and many specific WARNs.
