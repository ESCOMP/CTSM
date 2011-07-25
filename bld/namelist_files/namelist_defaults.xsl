<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
  <head>
    <title>CLM Namelist Defaults</title>
  </head>
  <body>
    <h2>Default Values for Namelist Variables</h2>
    <p>Included in the table are the following pieces of information:</p>
    <h3>Table headers include:</h3>
    <ul>
       <li>Name of variable</li>
       <li>Horizontal grid resolution</li>
       <li>Land ocean mask type</li>
       <li>Simulation year</li>
       <li>Simulation year range (for transient datasets)</li>
    </ul>
    <h3>Miscellaneous items include:</h3>
    <ol>
       <li>Biogeochemistry (BGC) type (none, CN, etc.)</li>
       <li>Initial condition date (ymd - year month day)</li>
       <li>Initial condition time of day (tod) (sec)</li>
       <li>River Transport Model (RTM) on or off</li>
       <li>Maximum number of Plant Function Types (maxpft)</li>
       <li>Number of glacier multiple elevation classes (glc_nec)</li>
       <li>Site specific point name (sitespf_pt)</li>
       <li>Glacier model grid size (glc_grid)</li>
       <li>Crop model (crop)</li>
       <li>Irrigation model (irrig)</li>
       <li>Data model forcing source (forcing)</li>
       <li>Accelerated Decomposition spin-up mode (ad_spin-up)</li>
    </ol>

    <table border="1" cellpadding="10">
    <caption>Namelist Defaults</caption>
    <tr>
      <th rowspan="2">Name</th>
      <th>Horz. Grid</th>
      <th>Mask</th>
      <th>Sim year</th>
      <th>Sim year range</th>
      <th>Miscellaneous</th>
    </tr>
    <tr>
      <th colspan="5">Default Value for this Configuration</th>
    </tr>
      <xsl:for-each select="namelist_defaults/*">
      <xsl:sort select="name()"/>
      <tr>
        <td rowspan="2"><font color="#ff0000">
        <xsl:value-of select="name()"/>
        </font></td>
        <td>
        <xsl:choose>
        <xsl:when test="string-length(@hgrid)>0">
             <xsl:value-of select="@hgrid"/>
        </xsl:when>
        <xsl:otherwise>
             All res
        </xsl:otherwise>
        </xsl:choose>
        </td>
        <td>
        <xsl:choose>
        <xsl:when test="string-length(@mask)>0">
             <xsl:value-of select="@mask"/>
        </xsl:when>
        <xsl:otherwise>
             All masks
        </xsl:otherwise>
        </xsl:choose>
        </td>
        <td>
        <xsl:choose>
        <xsl:when test="string-length(@sim_year)>0">
             <xsl:value-of select="@sim_year"/>
        </xsl:when>
        <xsl:otherwise>
             All yrs
        </xsl:otherwise>
        </xsl:choose>
        </td>
        <td>
        <xsl:choose>
        <xsl:when test="string-length(@sim_year_range)>0">
             <xsl:value-of select="@sim_year_range"/>
        </xsl:when>
        <xsl:otherwise>
             All sim-yr-rng
        </xsl:otherwise>
        </xsl:choose>
        </td>
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
        <xsl:if test="string-length(@sitespf_pt)>0">
        sitespf_pt=<xsl:value-of select="@sitespf_pt"/>
        </xsl:if>
        <xsl:if test="string-length(@glc_grid)>0">
        glc_grid=<xsl:value-of select="@glc_grid"/>
        </xsl:if>
        <xsl:if test="string-length(@datm_presaero)>0">
        datm_presaero=<xsl:value-of select="@datm_presaero"/>
        </xsl:if>
        <xsl:if test="string-length(@crop)>0">
        crop=<xsl:value-of select="@crop"/>
        </xsl:if>
        <xsl:if test="string-length(@irrig)>0">
        irrig=<xsl:value-of select="@irrig"/>
        </xsl:if>
        <xsl:if test="string-length(@ad_spinup)>0">
        ad_spinup=<xsl:value-of select="@ad_spinup"/>
        </xsl:if>
        <xsl:if test="string-length(@source)>0">
        forcing=<xsl:value-of select="@source"/>
        </xsl:if>
        </td>
      </tr>
      <tr>
        <td colspan="5"><b>Value: </b><xsl:value-of select="."      /></td>
      </tr>
      </xsl:for-each>
    </table>

  </body>

  </html>
</xsl:template>

</xsl:stylesheet>
