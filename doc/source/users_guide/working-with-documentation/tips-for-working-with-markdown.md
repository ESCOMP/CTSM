(tips-for-working-with-markdown)=

# Tips for working with Markdown

Markdown is great for very simple documentation files—it's much easier to write and read Markdown source than reStructuredText source.

Note that our documentation build system uses the MyST Markdown parser. MyST is a special "flavor" of Markdown that has a lot of reStructuredText-like features added to it; see the guide to MyST Markdown syntax [here](https://myst-parser.readthedocs.io/en/v5.0.0/syntax/typography.html), although note that we don't include a lot of the referenced extensions.

If you use VS Code, you may want to install the [MyST-Markdown extension](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) to get syntax highlighting, auto-complete, and improved rendering in the VS Code previewer.

```{contents}
:depth: 1
:local:
```

(md-cross-references)=

## Markdown: Cross-references

You can [link to Markdown section headers](#md-cross-references) by adding a label before the header:
```markdown
(md-cross-references)=

## Markdown: Cross-references
```

and then referring to it like `[link to section headers](#md-cross-references)`. You can also have the display text just be the title of the relevant section—e.g., [](#md-cross-references)—by doing `[](#md-cross-references)`.

That linking syntax also works for links to [section labels in reStructuredText files](#rst-cross-references): `[section labels in reStructuredText files](#rst-cross-references)`.

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

## Markdown: Admonitions
[Admonitions](https://myst-parser.readthedocs.io/en/v5.0.0/syntax/admonitions.html) are rendered as special "call-out" boxes. The general syntax is:
~~~
```{admonition} This is the title of a generic admonition
  It needs a title specified. Synonyms you can put in the `{}` instead of `admonition` include `note` and `seealso`; if you use one of those, you can't specify a title. In fact, `admonition` is the only type of admonition that you can specify a custom title for.
```
~~~
```{admonition} This is the title of a generic admonition
  It needs a title specified. Synonyms you can put in the `{}` instead of `admonition` include `note` and `seealso`; if you use one of those, you can't specify a title. In fact, `admonition` is the only type of admonition that you can specify a custom title for.
```

There are also a number of built-in admonition types that get their own special rendering:

~~~
```{attention}
The reader should pay special attention to this. Synonyms you can put in the `{}` instead of `attention` include `caution` and `warning`.
```
~~~
```{attention}
The reader should pay special attention to this. Synonyms you can put in the `{}` instead of `attention` include `caution` and `warning`.
```

~~~
```{danger}
This tells the reader about something dangerous. You can also put `error` in the `{}` instead of `danger`.
```
~~~
```{danger}
This tells the reader about something dangerous. You can also put `error` in the `{}` instead of `danger`.
```

~~~
```{hint}
Here's a hint. Synonyms you can put in the `{}` instead of `hint` include `important` and `tip`.
```
~~~
```{hint}
Here's a hint. Synonyms you can put in the `{}` instead of `hint` include `important` and `tip`.
```
