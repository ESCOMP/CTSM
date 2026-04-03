(tips-for-working-with-markdown)=

# Tips for working with Markdown

Markdown is great for very simple documentation files—it's much easier to write and read Markdown source than reStructuredText source.

Note that our documentation build system uses the MyST Markdown parser. MyST is a special "flavor" of Markdown that has a lot of reStructuredText-like features added to it; see the guide to MyST Markdown syntax [here](https://myst-parser.readthedocs.io/en/v5.0.0/syntax/typography.html), although note that we don't include a lot of the referenced extensions. 

(md-cross-references)=

## Markdown: Cross-references

You can [link to Markdown section headers](#md-cross-references) by adding a label before the header:
```markdown
(md-cross-references)=

## Markdown: Cross-references
```

and then referring to it like `[link to section headers](#md-cross-references)`. You can also have the display text just be the title of the relevant section—e.g., [](#md-cross-references)—by doing `[](#md-cross-references)`.

That linking syntax also works for links to [section labels in reStructuredText files](#rst-cross-references): `[reStructuredText section labels](#rst-cross-references)`.

(md-math)=

## Markdown: Math

Inline math can be achieved with the typical Markdown syntax of just surrounding your expression with dollar signs. E.g., $y = mx + b$:
```
E.g., $y = mx + b$:
```

You can also use Markdown's math block syntax for big equations on their own lines, with an optional label after the second `$$` that will give it a number and make it cross-referenceable:
```
$$
y = mx + b
$$ (my-equation-label)
```
$$
y = mx + b
$$ (my-equation-label)

Then you can get the equation number {eq}`my-equation-label` like so:
```
Then you can get the equation number {eq}`my-equation-label` like so:
```

## Markdown: Comments
If you want to add some text that's only visible in the documentation source file, you can put a `%` at the beginning of a line. E.g.:

```markdown
% This will not appear on the webpage.
```

% This will not appear on the webpage.

## Markdown: Tables

Markdown tables are supported. See [GitHub's "Organizing information with tables"](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/organizing-information-with-tables) for more info.
