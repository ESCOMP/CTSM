<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_definition">
  <head>
    <title>Configuration Definition</title>
  </head>
  <body>
    <h2>Configuration Definition</h2>

    <h3>Physics Configurations</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Valid Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='physics']"/>
    </table>

    <h3>Biogeochemistry Configurations</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Valid Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='bgc']"/>
    </table>

    <h3>Directories</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Valid Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='directories']"/>
    </table>

    <h3>Machine Options</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Valid Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='mach_options']"/>
    </table>

    <h3>Glacier model Options</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Valid Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='glc']"/>
    </table>

    <h3>Standalone CLM Testing Options (NOT used by normal CCSM scripts)</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Valid Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='standalone_test']"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td><xsl:value-of select="@value"/></td>
    <td><xsl:value-of select="@valid_values"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
</xsl:template>


</xsl:stylesheet>
