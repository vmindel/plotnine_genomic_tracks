# Claude Code skill

plotnine-genomic-tracks ships an
[Agent Skill](https://docs.claude.com/en/docs/claude-code/skills) for
[Claude Code](https://claude.com/claude-code). It teaches the agent the
track functions, how to assemble multi-track figures, and the mistakes
that give a wrong figure without an error.

## Install

Claude Code does not look inside installed Python packages for skills, so
copy it once:

```bash
pgt-install-skills
```

This installs to `~/.claude/skills/plotnine-genomic-tracks/`, available in
every project from the next session on.

| command | effect |
|---|---|
| `pgt-install-skills` | install to `~/.claude/skills/plotnine-genomic-tracks` |
| `pgt-install-skills --dest .claude/skills/plotnine-genomic-tracks` | install into the current project only |
| `pgt-install-skills --force` | overwrite, e.g. after upgrading |
| `pgt-install-skills --print-path` | print where the bundled copy is |

Then ask in plain language: "plot Msn2 ChEC and the CREs around Gata2 in
mm10 with plotnine", "facet this signal by condition", "stack these three
tracks with the genes at the bottom", "why is my gene track empty?".

It complements the pygenomic_tracks skill (`pygt-install-skills`), which
covers the matplotlib API; the descriptions tell Claude which one a task
needs.

## What's inside

A short router, `SKILL.md`, sends the agent to one reference file; the
agent reads it only when the task needs it.

| file | covers |
|---|---|
| `references/tracks.md` | `scale_x_genome`, `geom_signal`, `geom_intervals`, `geom_genes`: signatures, defaults, colour, overlays, facets, resolution |
| `references/layout.md` | `stack_tracks`, heights, saving, `theme_track`, shared y-ranges, highlights, several loci |
| `references/recipes.md` | pymisha to tracks, sample overlays, condition facets, coloured CREs and motifs, rasterising dense signal |
| `references/pitfalls.md` | empty tracks, mismatched regions, `ylim` on feature tracks, replaced coords, facets without `by=`, saving |

The test suite checks that every reference file is routed, every
`pgt.`/`pygt.` name it mentions exists, every argument in its signature
blocks is a real parameter, and that the install command works.
