.. _common-docs-errors:

# Common doc build errors and how to handle them

## Common docs errors: doc-builder

.. _common-rst-errors:

## Common docs errors: reStructuredText

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