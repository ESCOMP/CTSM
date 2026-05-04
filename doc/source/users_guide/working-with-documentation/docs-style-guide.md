(docs-style-guide)=

# Documentation style guide

- Please don't add manual line breaks when writing text, as this harms searchability. (Note that it's fine to do this in multi-line math blocks.) For more information, see [here](https://github.com/ESCOMP/CTSM/issues/2135#issuecomment-1764999337).
- You can write the degree symbol ° with Opt-Shift-8 on Mac or Alt+0176 on Windows, or you can copy it from here. Note that this is different from the [masculine ordinal indicator](https://en.wikipedia.org/wiki/Ordinal_indicator) (typed with Opt-0 on Mac). This is much cleaner and more searchable than using RestructuredText syntax to write a superscript-o!
- Whenever possible, please give equations meaningful labels. E.g., for [eq. 2.26.2](https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/Crop_Irrigation/CLM50_Tech_Note_Crop_Irrigation.html#equation-25-2), `:label: gdds_for_cfts` would be better than `:label: 25.2`. Numeric labels become obsolete—as you can see in that example!—whenever new equations and/or sections are added. (Note that the equation numbering in the rendered HTML is automatic.)
- Tables defined with the [:table: directive](https://docutils.sourceforge.io/docs/ref/rst/directives.html#table) can be annoying because they're very sensitive to the cells inside them being precisely the right widths, as defined by the first `====` strings. If you don't get the widths right, you'll see "Text in column margin" errors. Instead, define your tables using the `:list-table:` directive.
