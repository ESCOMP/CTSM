.. _caveats-for-working-with-markdown:

# Caveats for working with Markdown

.. _cross-references:

## Cross-references
You can [link to section headings](#cross-references) with the `[link to section headings](#cross-references)` format, but the Sphinx compiler will complain. To fix that, you will need to also add a reStructuredText label like `.. _cross-references:` above the heading, with blank lines before and after it. This rST label will appear if you render the .md file as pure Markdown (e.g., in VSCode's Markdown preview pane), but it will be invisible on the final generated webpage.

If you forget to surround the label with blank lines, you will get errors like "Explicit markup ends without a blank line; unexpected unindent [docutils]" that often point to lines far away from the actual problem.

## Inline math
(Note that parts of this section will be rendered incorrectly by Markdown parsers!)

Inline math can't be achieved with the typical Markdown syntax of just surrounding your expression with dollar signs. Instead, you need to surround THAT with backticks. So to render `$y = mx + b$`, we can't do
```
So to render $y = mx + b$, ...
```
because we'd just see $y = mx + b$ on the generated webpage. Instead, we do
```
So to render `$y = mx + b$`, ...
```
We could also use rST's syntax like so:
```
So to render :math:`y = mx + b`, ...
```