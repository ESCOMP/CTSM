<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
  <head>
    <title>CLM Namelist Defaults</title>
  </head>
  <body>
    <h2>Default values for namelist variables</h2>
    <p>Included in the table are the following pieces of information:</p>
    <ul>
       <li>Name of variable</li>
       <li>Horizontal grid resolution</li>
       <li>Land ocean mask type</li>
       <li>Simulation year</li>
       <li>Simulation year range (for transient datasets)</li>
       <li>Biogeochemistry (BGC) type (none, CN, etc.)</li>
       <li>Initial condition date (ymd - year month day)</li>
       <li>Initial condition time of day (tod) (sec)</li>
       <li>River Transport Model (RTM) on or off</li>
       <li>Maximum number of Plant Function Types (maxpft)</li>
       <li>Number of glacier multiple elevation classes (glc_nec)</li>
       <li>Glacier model grid size (glc_grid)</li>
    </ul>

    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Horz. Grid</th>
      <th>Mask</th>
      <th>Sim year</th>
      <th>Sim year range</th>
      <th>Miscellaneous</th>
      <xsl:for-each select="namelist_defaults/*">
      <xsl:sort select="name()"/>
      <tr>
        <td rowspan="2"><font color="#ff0000">
        <xsl:value-of select="name()"/>
        </font></td>
        <td><xsl:value-of select="@hgrid"/></td>
        <td><xsl:value-of select="@mask"/></td>
        <td><xsl:value-of select="@sim_year"/></td>
        <td><xsl:value-of select="@sim_year_range"/></td>
        <td>
        <xsl:if test="string-length(@bgc)>0">
        bgc=<xsl:value-of select="@bgc"/>
        </xsl:if>
        <xsl:if test="string-length(@ic_ymd)>0">
        ymd=<xsl:value-of select="@ic_ymd"/>
        </xsl:if>
        <xsl:if test="string-length(@ic_tod)>0">
        tod=<xsl:value-of select="@ic_tod"/>
        </xsl:if>
        <xsl:if test="string-length(@rtm)>0">
        rtm=<xsl:value-of select="@rtm"/>
        </xsl:if>
        <xsl:if test="string-length(@maxpft)>0">
        maxpft=<xsl:value-of select="@maxpft"/>
        </xsl:if>
        <xsl:if test="string-length(@glc_nec)>0">
        glc_nec=<xsl:value-of select="@glc_nec"/>
        </xsl:if>
        <xsl:if test="string-length(@glc_grid)>0">
        glc_grid=<xsl:value-of select="@glc_grid"/>
        </xsl:if>
        <xsl:if test="string-length(@datm_presaero)>0">
        datm_presaero=<xsl:value-of select="@datm_presaero"/>
        </xsl:if>
        </td>
      </tr>
      <tr>
        <td colspan="9"><b>Content: </b><xsl:value-of select="."      /></td>
      </tr>
      </xsl:for-each>
    </table>

  </body>

  </html>
</xsl:template>

<xsl:template match="name">
  <tr>
    <td><font color="#ff0000">----</font></td>
    <td><xsl:value-of select="@hgrid"/></td>
    <td><xsl:value-of select="@sim_year"/></td>
    <td><xsl:value-of select="@sim_year_range"/></td>
    <td><xsl:value-of select="@maxpft"/></td>
    <td><xsl:value-of select="@mask"/></td>
    <td><xsl:value-of select="@bgc"/></td>
    <td><xsl:value-of select="@glc_nec"/></td>
    <td><xsl:value-of select="@rtm"/></td>
    <td><xsl:value-of select="@ic_ymd"/></td>
    <td><xsl:value-of select="@ic_tod"/></td>
    <td><xsl:value-of select="@content"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
</xsl:template>

<xsl:template match="name/varname">
  <font color="#ff0000"><xsl:apply-templates/></font>
</xsl:template>

<xsl:template match="name/default">
  <p><xsl:apply-templates/></p>
</xsl:template>

<xsl:template match="name/listelm">
  <li><xsl:apply-templates/></li>
</xsl:template>

<xsl:template match="name/unlist">
  <ul>
  <xsl:apply-templates/>
  </ul>
</xsl:template>


</xsl:stylesheet>
