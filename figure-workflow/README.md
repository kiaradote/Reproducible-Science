# Figure workflow (week 2)

Skills and agents for making and reviewing data figures, built from the rules, good and bad examples and guidelines supplied by the user.

- `skills/figure-rubric/` - shared numbered rules, 42 in all, covering chart style fit and unbiased presentation (each with source slide, severity, and origin tag NOTES / SLIDE / ADDED; rules tagged ADDED await the user's confirmation) plus `kernel.py` helpers: `lint_script`, `audit_figure`, color, scale and gap helpers.
- `skills/figure-generation/` - workflow: raw data -> week-1 analysis checkpoints -> figure spec -> build -> self-lint -> automatic review (max 3 rounds).
- `agents/FIGURE_MAKER.md` - builds figures; sends every figure to the reviewer.
- `agents/FIGURE_REVIEWER.md` - read-only critic; returns a PASS/WARN/FAIL card per rule.

These are exported copies of the skills and agent profiles in the Claude Science account (the account versions are the live ones). Rules tagged ADDED are additions not stated in class and need student confirmation.
