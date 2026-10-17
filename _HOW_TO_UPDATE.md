# Updating gcoqueret.com

Edit, then run `quarto render` in this folder (or click Render in RStudio), then commit and push.

| To change…                               | Edit this file                                  |
|------------------------------------------|-------------------------------------------------|
| A job, position or degree                | `data/cv.yml`                                   |
| A published article                      | `data/articles.yml`                             |
| A working paper                          | `data/working-papers.yml`                       |
| Practitioner press, chapters, talks…     | `data/practitioner.yml`, `data/chapters.yml`, `data/permanent-papers.yml`, `data/presentations.yml` |
| A course on the teaching page            | `data/courses.yml`                              |
| Intro texts, books, profiles, contact    | `index.qmd`, `research.qmd`, `teaching.qmd`     |
| A course or research material page       | `ml_factor.qmd`, `macro.qmd`, `rstats.qmd`, `shiny.qmd`, `PRs.qmd`, `horizons.qmd`, `tot.qmd` |
| "Updated October 2026" and footer links  | `_partials/footer.html`                         |
| Colors, fonts, spacing                   | `assets/site.css` (colors are at the top)       |

Rules of thumb for the `.yml` files:

- Copy an existing block and change the values. Order in the file = order on the page.
- Put text in "quotes" when it contains a colon.
- Markdown works inside values: `*italic*`, `[label](https://...)`.
- `links:` takes one `Label: https://...` per line.
- In `articles.yml`, `topic:` decides the filter button and `selected: true` shows the paper in the default view.

On the material pages, a numbered session looks like:

```markdown
1. **Title**: short description [[Notebook](files/...html) [Slides](files/...pdf)]{.links}
```

To add a new page: create `newpage.qmd` (copy `tot.qmd` as a start) and add it to the `render:` list in `_quarto.yml`.
