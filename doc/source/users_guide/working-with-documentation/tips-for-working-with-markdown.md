.. _tips-for-working-with-markdown:

# Tips for working with Markdown

Markdown is great for very simple documentation files—it's much easier to write and read Markdown source than reStructuredText source. However, there are some compromises that you should be aware of, and you may find yourself needing to mix in some reStructuredText.

.. _md-cross-references:

## Markdown: Cross-references
You can [link to section headings](#md-cross-references) with the `[link to section headings](#md-cross-references)` format, but the Sphinx compiler will complain. Instead, use the :ref:`reStructuredText cross-reference and label<rst-cross-references>` syntax.

## Markdown: Math
(Note that parts of this section will be rendered incorrectly by Markdown parsers!)

Inline math can't be achieved with the typical Markdown syntax of just surrounding your expression with dollar signs. Instead, you need to surround THAT with backticks. So to render `$y = mx + b$`, we can't do
```
So to render $y = mx + b$, ...
```
because we'd just see $y = mx + b$ on the generated webpage. Instead, we do
```
So to render `$y = mx + b$`, ...
```
We could also use :ref:`rST's syntax<rst-math>` like so:
```
So to render :math:`y = mx + b`, ...
```

You can also use Markdown's math block syntax for big equations on their own lines:
```
$$
y = mx + b
$$
```
$$
y = mx + b
$$

However, you won't get the equation numbering or labeling that you would with the :ref:`reStructuredText math format<rst-math>`.

## Markdown: Comments
If you want to add some text that's only visible in the documentation source file, there's not really a way to do that in Markdown. However, you can use the :ref:`reStructuredText comment syntax<rst-comments>` in a Markdown document.

## Markdown: Tables

Markdown tables are supported. See [GitHub's "Organizing information with tables"](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/organizing-information-with-tables) for more info.