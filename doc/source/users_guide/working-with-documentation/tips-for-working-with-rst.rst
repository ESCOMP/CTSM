.. _tips-for-working-with-rst:

Tips for working with reStructuredText
========================================

If you've never used reStructuredText before, you should be aware that its syntax is pretty different from anything you've ever used before. We recommend the following resources as references for the syntax:

- `Sphinx's reStructuredText Primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
- The `Quick reStructuredText <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_ cheat sheet

Some especially useful bits:

- `Section headers <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#sections>`_
- `Hyperlinks <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#hyperlinks>`_

On this page, we've compiled some supplemental information that might be helpful, including a list of common errors and their causes.

.. contents::
   :depth: 1
   :backlinks: top
   :local:

.. _rst-math:

reStructuredText: Math
----------------------

You can write inline math like ``:math:`y = mx + b``` → :math:`y = mx + b`. You can also write bigger equations on their own line that will automatically be numbered:

.. code-block:: reStructuredText

  .. math::
    :label: equation for a line

    y = mx + b

.. math::
  :label: equation for a line

  y = mx + b

Note (a) the leading spaces for each line after ``.. math::`` and (b) the empty line after the label. If you don't include the ``:label:`` line, the equation will not be numbered.

reStructuredText math largely follows LaTeX syntax.

.. _rst-cross-references:

reStructuredText: Cross-references
----------------------------------

reStructuredText lets you define labels that can be cross-referenced as links elsewhere in the documentation. A label looks like ``.. _this-is-my-label:`` or ``.. _This is my label with CAPS and spaces:``, on its own line surrounded by blank lines. The leading ``.. _`` and trailing ``:`` are what tell rST "this line is a label." E.g.:

.. code-block:: reStructuredText

  Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore   magna aliqua.
  
  .. _this-is-my-label:
  
  My cool information
  ^^^^^^^^^^^^^^^^^^^
  Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.

You could then refer to that section with ``:XYZ:`this-is-my-label``` (leaving off the leading ``.. _``), where ``XYZ`` can be ``ref``, ``numref``, or ``eq`` (see examples below). This will create a link that, when clicked, takes the reader to the "My cool information" section.

Here are some examples. Note that the displayed link text will update automatically as needed (e.g., a section number or figure caption gets changed).

- Section headings, text: ``:ref:`rst-cross-references``` → :ref:`rst-cross-references`
- Section headings, number: ``:numref:`rst-cross-references``` → :numref:`rst-cross-references`. Note that, unlike ``numref`` for other things mentioned here, "Section" is not automatically prepended to the section number in the link text.
- Table, text: ``:ref:`Table Crop plant functional types``` → :ref:`Table Crop plant functional types`
- Table, number: ``:numref:`Table Crop plant functional types``` → :numref:`Table Crop plant functional types`
- Figure, text (uses entire caption): ``:ref:`Figure CLM subgrid hierarchy``` → :ref:`Figure CLM subgrid hierarchy`
- Figure, number: ``:numref:`Figure CLM subgrid hierarchy``` → :numref:`Figure CLM subgrid hierarchy`
- Equation, number: ``:eq:`equation for a line``` → :eq:`equation for a line`. The parentheses in the link text seem unavoidable, and there seems to be no way to refer to have the link show the label text or anything else aside from the number.

You can have any link (except for equations) show custom text by putting the referenced label at the end in ``<angled brackets>``. E.g., ``:ref:`Diagram of CLM subgrid hierarchy<Figure CLM subgrid hierarchy>``` → :ref:`Diagram of CLM subgrid hierarchy<Figure CLM subgrid hierarchy>`. 

Note that this is necessary for labels that aren't immediately followed by a section heading, a table with a caption, or a figure with a caption. For instance, to refer to labels in our bibliography, you could do ``:ref:`(Bonan, 1996)<Bonan1996>``` → :ref:`(Bonan, 1996)<Bonan1996>`.

.. _rst-comments:

reStructuredText: Comments
--------------------------

If you want to add some text that's only visible in the documentation source file, you can use the reStructuredText comment syntax:

.. code-block:: reStructuredText

  ..
    This will not appear on the webpage or even anywhere in the generated HTML.


Make sure to include at least one empty line after the comment text.


reStructuredText: Tables
------------------------

Tables defined with the `:table: directive <https://docutils.sourceforge.io/docs/ref/rst/directives.html#table>`_ can be annoying because they're very sensitive to the cells inside them being precisely the right widths, as defined by the first ``====`` strings. If you don't get the widths right, you'll see "Text in column margin" errors. Instead, define your tables using the `:list-table: <https://docutils.sourceforge.io/docs/ref/rst/directives.html#list-table>`_ directive.

If you already have a table in some other format, like comma-separated values (CSV), you may want to check out the R package `knitr <https://cran.r-project.org/web/packages/knitr/index.html>`_. Its `kable <https://bookdown.org/yihui/rmarkdown-cookbook/kable.html>`_ command allows automatic conversion of R dataframes to tables in reStructuredText and other formats.

reStructuredText: Admonitions (e.g., warning, tip)
--------------------------------------------------

`Admonitions <https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#admonitions-messages-and-warnings>`_ are rendered as special "call-out" boxes. The general syntax is:

.. code-block:: reStructuredText

  .. admonition:: This is the title of a generic admonition

    It needs a title specified. Synonyms you can put in the ``{}`` instead of ``admonition`` include ``note`` and   ``seealso``; if you use one of those, you don't need to specify a title.
  

.. admonition:: This is the title of a generic admonition

  It needs a title specified. Synonyms you can put in the ``{}`` instead of ``admonition`` include ``note`` and ``seealso``; if you use one of those, you don't need to specify a title.


There are also a number of built-in admonition types that get their own special rendering:

.. code-block:: reStructuredText

  .. attention::

    The reader should pay special attention to this. Synonyms you can put in the ``{}`` instead of ``attention`` include ``caution`` and ``warning``.

.. attention::

  The reader should pay special attention to this. Synonyms you can put in the ``{}`` instead of ``attention`` include ``caution`` and ``warning``.

.. code-block:: reStructuredText

  .. danger::

    This tells the reader about something dangerous. You can also put ``error`` in the ``{}`` instead of ``danger``.

.. danger::

  This tells the reader about something dangerous. You can also put ``error`` in the ``{}`` instead of ``danger``.

.. code-block:: reStructuredText

  .. hint::

    Here's a hint. Synonyms you can put in the ``{}`` instead of ``hint`` include ``important`` and ``tip``.

.. hint::

  Here's a hint. Synonyms you can put in the ``{}`` instead of ``hint`` include ``important`` and ``tip``.



reStructuredText: Text subscripts
--------------------------------------------------------------

Note that these instructions apply to text (i.e., not in inline or block math).

**Superscripts**

- Incorrect: ``Fe:sup:`3+``` becomes Fe:sup:`3+`
- Incorrect: ``Fe :sup:`3+``` becomes Fe :sup:`3+` (note extraneous space)
- Correct: ``Fe\ :sup:`3+``` becomes Fe\ :sup:`3+`

**Subscripts**

- Incorrect: ``CO:sub:`2``` becomes CO:sub:`2`
- Incorrect: ``CO :sub:`2``` becomes CO :sub:`2` (note extraneous space)
- Correct: ``CO\ :sub:`2``` becomes CO\ :sub:`2`



reStructuredText: Common error messages and how to handle them
--------------------------------------------------------------

See :ref:`common-rst-errors`.