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
<p>
</p>
<hr></hr>
<p>
</p>
    <h2>CLM namelist variables</h2>
    <p>Note, these all would go into the clm.buildnml.csh file:</p>
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
    <li>Variable description (includes information on defaults).</li>
    <li>Valid values (if restricted).</li>
    </ul>

    <h3>CLM Physics Options</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_physics']"/>
    </table>

    <h3>CLM Datasets</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datasets']"/>
    </table>

    <h3>CLM History output settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='history']"/>
    </table>

    <h3>CLM Restart settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_restart']"/>
    </table>

    <h3>CLM Performance Tuning</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_performance']"/>
    </table>

    <h3>CASA Namelist items in clm_inparm namelist (EXPERIMENTAL NOT SUPPORTED!)</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='casa']"/>
    </table>

<!--
    <h3>PIO Namelist items in clm_inparm namelist (EXPERIMENTAL NOT WORKING!)</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='pio']"/>
    </table>
-->
<p>
</p>
<hr></hr>
<p>
</p>
    <h2>Command Line Options to CLM Build-namelist</h2>
    <p>Variables that are entered as options to build-namelist (but NOT used by 
        namelists in code). Included in the table are the following pieces 
        of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Type.</li>
    <li>Valid values.</li>
    <li>Variable description.</li>
    </ul>

    <h3>CLM Default Settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='default_settings']"/>
    </table>
<p>
</p>
<hr></hr>
<p>
</p>
    <h2>Driver and DATM namelist variables</h2>
    <p>Note, these all would go into the cpl.buildnml.csh or datm.buildnml.csh files:</p>
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
    <li>Variable description (includes information on defaults).</li>
    <li>Valid values (if restricted).</li>
    </ul>

    <h3>Driver Solar, orbital and radiation options</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='radiation']"/>
    </table>

    <h3>Data atmosphere model settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm']"/>
    </table>

    <h3>Driver</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='driver']"/>
    </table>

    <h3>Driver Time Manager</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='time_mgr']"/>
    </table>

    <h3>Driver and DATM Restart (Continuation and Branch) Runs</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='restart']"/>
    </table>

    <h3>Driver and DATM History Options</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='drv_history']"/>
    </table>

    <h3>Driver Performance Tuning and Profiling</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='performance']"/>
    </table>

    <h3>Driver Single Column Atmosphere Model (part of CAM)</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='scam']"/>
    </table>

    <h3>Data atmosphere model internal settings</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm_internal']"/>
    </table>

    <h3>Driver Dry-Deposition Options</h3>
    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='dry-deposition']"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td rowspan="2"><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td rowspan="2"><xsl:value-of select="@type"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
  <tr>
    <td><xsl:value-of select="@valid_values"/></td>
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
