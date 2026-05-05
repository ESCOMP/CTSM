(converting-between-md-and-rst)=

# Converting between Markdown and reStructuredText

You may find yourself wanting to convert a Markdown file to a reStructuredText file or vice versa. [`pandoc`](https://pandoc.org/) is a command-line utility that can do this for you, or at least get you something that's most of the way there. Installation instructions can be found [here](https://pandoc.org/installing.html), but it's pretty complicated. An easier alternative is to use the online demo page [here](https://pandoc.org/try/):

1. Paste the contents of the file into the text box on the left.
1. Select the correct file types in the "from" and "to" menus (just use plain "Markdown", not one of the other flavors, unless you've copied something in one of those other flavors—e.g., if you've copied Markdown from a GitHub wiki page).
1. **Make sure "Preserve breaks" is selected, not "Auto wrap"**—see [style guide](#ctsm-docs-style-guide).
1. Click "Convert". You will be given the results of the conversion on the right side of the page.
1. Copy and paste that into a new file with the new extension.

Note that `pandoc` doesn't support [MyST Markdown](https://mystmd.org/guide/quickstart-myst-markdown), the special flavor of Markdown we support in our docs. Thus, if you're converting a Markdown file from our docs to reStructuredText, there may be some MyST Markdown left over. E.g.:

```
(common-docs-errors)=

# Common doc build errors and how to handle them

(common-doc-builder-errors)=

## Common docs errors: doc-builder

### "RuntimeError: No compatible container software found: docker, podman"

You tried to build the documentation using our container (`./build_docs ... -d`) but didn't have any software running that could handle the container. Try again after starting up container software according to the instructions for your platform:
- {ref}`bld-prev-docs-casper` **(recommended)**
- {ref}`bld-prev-docs-mac`
- {ref}`bld-prev-docs-windows`
```

becomes:

```reStructuredText
(common-docs-errors)=

Common doc build errors and how to handle them
==============================================

(common-doc-builder-errors)=

Common docs errors: doc-builder
-------------------------------

“RuntimeError: No compatible container software found: docker, podman”
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You tried to build the documentation using our container (``./build_docs ... -d``) but didn’t have any software running that could handle the container. Try again after starting up container software according to the instructions for your platform:
- {ref}\ ``bld-prev-docs-casper`` **(recommended)**
- {ref}\ ``bld-prev-docs-mac``
- {ref}\ ``bld-prev-docs-windows``
```

The labels `(common-docs-errors)=` and `(common-doc-builder-errors)=` will need to be manually converted to their [reStructuredText equivalents](#rst-cross-references) `.. _common-docs-errors:` and `.. _common-doc-builder-errors:`, respectively.

This can also happen when converting in the other direction, from reStructuredText to Markdown. In that case, however, there are tools that can convert from reStructuredText to MyST Markdown specifically. See the command-line tool [`rst-to-myst`](https://rst-to-myst.readthedocs.io/en/latest/) and the [MySTyc website](https://astrojuanlu.github.io/mystyc/), although at the time of writing the latter wasn't working.
