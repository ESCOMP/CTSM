<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="namelist_definition">
  <head>
    <title>CLM Namelist Definition</title>
  </head>
  <body>
    <h2>CLM namelist variables</h2>
    <p>Included in the table are the following pieces of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Variable type (<code>char</code>, <code>integer</code>,
    <code>real</code>, or <code>logical</code>).  The type
    <code>char</code> has the length appended
    following an asterisk, e.g., <code>char*256</code>.  Variables that are
    arrays have their dimension specifier appended inside parentheses.  For
    example <code>char*1(6)</code> denotes a array of six
    <code>char*1</code> values.
    </li>
    <li>Valid values (if restricted).</li>
    <li>Variable description (includes information on defaults).</li>
    </ul>

    <h3>CLM Physics Options</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='clm_physics']"/>
    </table>

    <h3>Solar, orbital and radiation options</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='radiation']"/>
    </table>

    <h3>Data atmosphere model settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='datm']"/>
    </table>

    <h3>Driver</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='driver']"/>
    </table>

    <h3>Time Manager</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='time_mgr']"/>
    </table>

    <h3>Datasets</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='datasets']"/>
    </table>

    <h3>History output settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='history']"/>
    </table>

    <h3>Restart (Continuation and Branch) Runs</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='restart']"/>
    </table>

    <h3>Performance Tuning and Profiling</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='performance']"/>
    </table>

    <h3>Single Column Atmosphere Model (part of CAM)</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='scam']"/>
    </table>

    <h2>Command Line Options to Build-namelist</h2>
    <p>Variables that are entered as options to build-namelist (but NOT used by 
        namelists in code). Included in the table are the following pieces 
        of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Type.</li>
    <li>Valid values.</li>
    <li>Variable description.</li>
    </ul>


    <h3>Default Settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='default_settings']"/>
    </table>

    <h3>Data atmosphere model internal settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Type</th>
      <th>Valid values</th>
      <th>Description</th>
      <xsl:apply-templates select="entry[@category='datm_internal']"/>
    </table>


  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td><xsl:value-of select="@type"/></td>
    <td><xsl:value-of select="@valid_values"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
</xsl:template>

<xsl:template match="entry/varname">
  <font color="#ff0000"><xsl:apply-templates/></font>
</xsl:template>

<xsl:template match="entry/default">
  <p><xsl:apply-templates/></p>
</xsl:template>

<xsl:template match="entry/listelm">
  <li><xsl:apply-templates/></li>
</xsl:template>

<xsl:template match="entry/unlist">
  <ul>
  <xsl:apply-templates/>
  </ul>
</xsl:template>


</xsl:stylesheet>
