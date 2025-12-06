.. _tips-for-working-with-rst:

# Tips for working with reStructuredText

If you've never used reStructuredText before, you should be aware that its syntax is pretty different from anything you've ever used before. We recommend the following resources as references for the syntax:
- [Sphinx's reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- The [Quick reStructuredText](https://docutils.sourceforge.io/docs/user/rst/quickref.html) cheat sheet

Some especially useful bits:
- [Section headers](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#sections)
- [Hyperlinks](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#hyperlinks)
- [Callout blocks (e.g., warning, tip)](https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#admonitions-messages-and-warnings)

On this page, we've compiled some supplemental information that might be helpful, including a list of common errors and their causes.

.. contents::
   :depth: 1
   :backlinks: top
   :local:

.. _rst-math:

## reStructuredText: Math
You can write inline math like ``:math:`y = mx + b``` → :math:`y = mx + b`. You can also write bigger equations on their own line that will automatically be numbered:

```reStructuredText
.. math::
  :label: equation for a line

  y = mx + b
```
.. math::
  :label: equation for a line

  y = mx + b

Note (a) the leading spaces for each line after `.. math::` and (b) the empty line after the label. If you don't include the `:label:` line, the equation will not be numbered.

reStructuredText math largely follows LaTeX syntax.

.. _rst-cross-references:

## reStructuredText: Cross-references
reStructuredText lets you define labels that can be cross-referenced as links elsewhere in the documentation. A label looks like ``.. _this-is-my-label:`` or ``.. _This is my label with CAPS and spaces:``, on its own line surrounded by blank lines. The leading ``.. _`` and trailing ``:`` are what tell rST "this line is a label." E.g.:

```
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.

.. _this-is-my-label:

My cool information
^^^^^^^^^^^^^^^^^^^
Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
```

You could then refer to that section with ``:XYZ:`this-is-my-label``` (leaving off the leading `.. _`), where `XYZ` can be `ref`, `numref`, or `eq` (see examples below). This will create a link that, when clicked, takes the reader to the "My cool information" section.

Here are some examples. Note that the displayed link text will update automatically as needed (e.g., a section number or figure caption gets changed).
- Section headings, text: ``:ref:`rst-cross-references``` → :ref:`rst-cross-references`
- Section headings, number: ``:numref:`rst-cross-references``` → :numref:`rst-cross-references`. Note that, unlike `numref` for other things mentioned here, "Section" is not automatically prepended to the section number in the link text.
- Table, text: ``:ref:`Table Crop plant functional types``` → :ref:`Table Crop plant functional types`
- Table, number: ``:numref:`Table Crop plant functional types``` → :numref:`Table Crop plant functional types`
- Figure, text (uses entire caption): ``:ref:`Figure CLM subgrid hierarchy``` → :ref:`Figure CLM subgrid hierarchy`
- Figure, number: ``:numref:`Figure CLM subgrid hierarchy``` → :numref:`Figure CLM subgrid hierarchy`
- Equation, number: ``:eq:`equation for a line``` → :eq:`equation for a line`. The parentheses in the link text seem unavoidable, and there seems to be no way to refer to have the link show the label text or anything else aside from the number.

You can have any link (except for equations) show custom text by putting the referenced label at the end in `<angled brackets>`. E.g., ``:ref:`Diagram of CLM subgrid hierarchy<Figure CLM subgrid hierarchy>``` → :ref:`Diagram of CLM subgrid hierarchy<Figure CLM subgrid hierarchy>`. 

Note that this is necessary for labels that aren't immediately followed by a section heading, a table with a caption, or a figure with a caption. For instance, to refer to labels in our bibliography, you could do ``:ref:`(Bonan, 1996)<Bonan1996>``` → :ref:`(Bonan, 1996)<Bonan1996>`.

.. _rst-comments:

## reStructuredText: Comments
If you want to add some text that's only visible in the documentation source file, you can use the reStructuredText comment syntax:

```
..
  This will not appear on the webpage or even anywhere in the generated HTML.

```

Make sure to include at least one empty line after the comment text.


## reStructuredText: Tables
Tables defined with the [:table: directive](https://docutils.sourceforge.io/docs/ref/rst/directives.html#table) can be annoying because they're very sensitive to the cells inside them being precisely the right widths, as defined by the first `====` strings. If you don't get the widths right, you'll see "Text in column margin" errors. Instead, define your tables using the [list-table](https://docutils.sourceforge.io/docs/ref/rst/directives.html#list-table) directive.

If you already have a table in some other format, like comma-separated values (CSV), you may want to check out the R package [knitr](https://cran.r-project.org/web/packages/knitr/index.html). Its [kable](https://bookdown.org/yihui/rmarkdown-cookbook/kable.html) command allows automatic conversion of R dataframes to tables in reStructuredText and other formats.


## reStructuredText: Common error messages and how to handle them

.. _error-unexpected-unindent:

### "ERROR: Unexpected indentation"

Like Python, reStructuredText is very particular about how lines are indented. Indentation is used, for example, to denote [code ("literal") blocks](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#literal-blocks) and [quote blocks](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#lists-and-quote-like-blocks). An error like
```
/path/to/file.rst:102: ERROR: Unexpected indentation. [docutils]
```
indicates that line 102 is indented but not in a way that reStructuredText expects.

### "WARNING: Block quote ends without a blank line; unexpected unindent"

This is essentially the inverse of :ref:`error-unexpected-unindent`: The above line was indented but this one isn't. reStructuredText tried to interpret the indented line as a [block quote](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#lists-and-quote-like-blocks), but block quotes require a blank line after them.

.. _inline-literal-start-without-end:

### "WARNING: Inline literal start-string without end-string"

An "inline literal" is when you want to mix code into a normal line of text (as opposed to in its own code block) ``like this``. This is accomplished with double-backticks:
```reStructuredText
An "inline literal" is when you want to mix code into a normal line of
text (as opposed to in its own code block) ``like this``.
```
(A backtick is what you get if you press the key to the left of 1 on a standard US English keyboard.)

If you have a double-backtick on a line, reStructuredText will think, "They want to start an inline literal here," then look for another double-backtick to end the literal. The "WARNING: Inline literal start-string without end-string" means it can't find one on that line.

This might happen, for example, if you try to put a [Markdown code block](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/creating-and-highlighting-code-blocks) in a .rst file. In that case, use the [reStructuredText code block syntax](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#literal-blocks) instead (optionally with [syntax highlighting](https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-highlight)).

### "WARNING: Inline interpreted text or phrase reference start-string without end-string"

Like :ref:`inline-literal-start-without-end`, this is probably related to having one double-backtick without another on the same line. As with that other error, it could be the result of a Markdown code block in a .rst file.

### "ERROR: Error in "code" directive: maximum 1 argument(s) allowed, 19 supplied"

This error might show something other than "code," like "highlight" or "sourcecode". It also will probably show a second number that's not 19. The problem is that you tried to write a [reStructuredText code block with syntax highlighting](https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-highlight) but didn't include a blank line after the first one:

```reStructuredText
.. code:: shell
 # How to list all the available grids
 cd cime/scripts
 ./query_config --grids
```

Fix this by adding a blank line:
```reStructuredText
.. code:: shell
 
 # How to list all the available grids
 cd cime/scripts
 ./query_config --grids
```

### 'ERROR: Error in "math" directive: invalid option block'

You might have forgotten the empty line after an equation label.

### "WARNING: Explicit markup ends without a blank line; unexpected unindent"

You might have forgotten the leading spaces for every line after `.. math::`. As a reminder, you need at least one leading space on each line.

You can also get this error if you forget to surround a :ref:`cross-reference label<rst-cross-references>` with blank lines. In this case, the error message might point to lines far away from the actual problem.

### "WARNING: Failed to create a cross reference: A title or caption not found"
This probably means you tried to `:ref:` a label that's not immediately followed by (a) a table/figure with a caption or (b) a section.

### "WARNING: undefined label"

If you're sure the label you referenced actually exists, this probably means you tried to ``:numref:`` a label that's not immediately followed by a table, figure, or section (see above). Alternatively, you might have tried to ``:ref:`` an :ref:`equation<rst-math>`; in that case, use ``:eq:`` instead.

### "WARNING: malformed hyperlink target"

You may have forgotten the trailing `:` on a label line.