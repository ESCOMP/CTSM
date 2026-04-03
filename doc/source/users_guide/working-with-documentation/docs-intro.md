.. _docs-intro:

# Working with the CTSM documentation

If you're reading this page, you're probably interested in making an improvement to this CTSM documentation website—thank you! Please read this page thoroughly before starting to make it as smooth a process as possible.

## How the documentation works
The documentation starts life as files in the `doc/source/tech_note/` and `doc/source/users_guide/` directories. These files are written in a mixture of what are called "markup languages." You may already be familiar—for better or for worse—with the LaTeX markup language. Fortunately, our documentation is simpler than that. It was originally written entirely in [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html), and it still mostly is, as you can tell by the predominance of .rst files. However, it's also possible to write Markdown documents (.md), which is nice because it's a much simpler and more widespread format (although see :ref:`tips-for-working-with-markdown`). If you've formatted text on GitHub, for instance, you've used Markdown.

During the "build" process, these files are converted to HTML webpages using a tool called Sphinx. This creates a directory with HTML files that can then be "previewed" to make sure they look right. When you submit a pull request (PR; see [Contribution guidelines](#docs-contribution-guidelines) below), tests will automatically run to make sure your changed documentation builds with no errors. Then, when your PR is merged to the main CTSM branch, it will be rebuilt and published automatically.

.. _docs-contribution-guidelines:

## Contribution guidelines
We use the [CTSM GitHub repo](https://github.com/ESCOMP/CTSM) to track issues with the documentation and bring in changes. Please have a look at our ["How to contribute"](https://github.com/ESCOMP/CTSM/blob/master/CONTRIBUTING.md) readme for some general guidelines.

If you have found a problem with the documentation but aren't able to fix it immediately, or if you're not sure whether something is truly a problem or how to fix it, please [file an issue](https://github.com/ESCOMP/CTSM/issues/new?template=03_documentation.md).

If you've made changes that you'd like us to bring in, you can file a [pull request](https://github.com/ESCOMP/ctsm/pulls) (PR; "New pull request" button at that link). Please try to avoid filing many small documentation PRs in a short time, as it's easier for the maintenance team if you combine similar edits into one larger PR. For example, several typo fixes in a file or across files would be a good single PR. You can also mark a PR as "draft" if you think you may be adding more to it, or if you think it's otherwise not ready for review.

Whenever you submit a documentation PR or commit new changes to one, automated testing will check the updated documentation for any errors. If you get failures, please try to diagnose and fix them; see :ref:`common-docs-errors` for tips. If you resolve the errors, add a comment on your PR saying so, and one of the CTSM software engineers will have a look.

## Editing the documentation
First, you will need a clone of CTSM to get all the documentation files and infrastructure. Once you have that, editing the documentation is as simple as opening the source file for the page you want to edit, then changing text. Make sure to use either reStructuredText or Markdown syntax, depending on the file's extension (.rst or .md, respectively). You can mix some reStructuredText into Markdown files, but generally not the other way around. For more information, see:
- :ref:`tips-for-working-with-markdown`
- :ref:`tips-for-working-with-rst`

We strongly recommend building and previewing the docs yourself as you're editing them, or at least before requesting review on your PR. You can find instructions for doing so here:
- :ref:`bld-prev-docs-casper` **(recommended)**
- :ref:`bld-prev-docs-mac`
- :ref:`bld-prev-docs-windows`
