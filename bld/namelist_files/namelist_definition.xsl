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
<hr/>
<p>
</p>
    <h1>Definition of CESM namelist variables</h1>
    <p>We list all of the relevant namelist variables for CLM I cases. This includes
    CLM Namelist items as well as CLM build-namelist settings, and datm and driver
    namelists.</p>
<hr/>
    <h2>Definition of CLM namelist variables</h2>
    <p>Note, these all would go into the clm.buildnml.csh file (or into user_nl_clm
    before configure):</p>
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

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Physics Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_physics']"/>
    </table>

    <p>
      <b>NOTE:</b>
      Below "atm" and "lnd" refer to the resolutions for the fine-mesh.
      Datasets with "atm" in the name or "atm grid" in the description are
      at the course resolution the model is coupled to the atmosphere with.
      Datasets with "lnd" in the name or "lnd grid" are at the fine mesh
      resolution. By default datasets are all on the same resolution, and the
      fine-mesh is NOT activated. The fine-mesh is considered experimental and
      should NOT normally be used.
    </p>
    <table border="1" cellpadding="10">
    <caption>CLM Namelist Datasets</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description
      </th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datasets']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist History output settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='history']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Restart settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_restart']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Performance Tuning</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_performance']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CASA Namelist items in clm_inparm namelist (EXPERIMENTAL NOT SUPPORTED!)</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='casa']"/>
    </table>

<p>
</p>
<hr/>
<p>
</p>
    <h2>Command Line Options to CLM Build-namelist</h2>
    <p>Variables that are entered as options to build-namelist (but NOT used by 
        namelists in code). Most of these are options that could be added to
        CLM_BLDNML_OPTS. Included in the table are the following pieces 
        of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Type.</li>
    <li>Valid values.</li>
    <li>Variable description.</li>
    </ul>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Default Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='default_settings']"/>
    </table>
<p>
</p>
<hr/>
<p>
</p>
    <h2>Namelist items for CLM Tools</h2>
    <p>These are namelist items that appear in the CLM Tools under models/lnd/clm/tools.
    </p>
    <table border="1" cellpadding="10">
    <caption>CLM mksurfdata</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='mksurfdata']"/>
    </table>
    <table border="1" cellpadding="10">
    <caption>CLM mkgriddata</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='mkgriddata']"/>
    </table>
    <table border="1" cellpadding="10">
    <caption>Miscellaneous CLM tools</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='tools']"/>
    </table>

<hr/>
<p>
</p>
    <h2>Driver namelist variables</h2>
    <p>Note, these all would go into the cpl.buildnml.csh file.
    Many of these items are set with env settings in env_conf.xml. In
    most cases of those that do, the env name is similar enough to make
    the mapping obvious. In some cases where it is not, we give the env
    variable setting that is used to set it.
    </p>
    <p>
    Also note that many of the driver namelist items do NOT apply to "I" 
    compsets, but to fully coupled cases with active atmosphere, ocean
    and/or sea-ice.
    </p>
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

    <table border="1" cellpadding="10">
    <caption>Driver Namelist Solar, Orbital and Radiation Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='radiation']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Driver Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='driver']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Driver Time Manager Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='time_mgr']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Driver Restart (Continuation and Branch) Run Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='drv_restart']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Driver History Namelist Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='drv_history']"/>
    </table>

    <p>
    <b>NOTE:</b> Most of these physics settings have to do with fully coupled cases
    with active ocean and atmosphere. Most do NOT have to do with CLM.
    </p>
    <table border="1" cellpadding="10">
    <caption>Driver Physics Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='drv_physics']"/>
    </table>
    <p>
    <b>Note:</b> Most of the following performance items set the tasks
    and threads and are in the env_mach_pes.xml file. The pio settings
    are all in cpl.buildnml.csh.
    </p>
    <table border="1" cellpadding="10">
    <caption>Driver Performance Tuning Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='performance']"/>
    </table>
    <table border="1" cellpadding="10">
    <caption>Driver Performance Profiling Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='profiling']"/>
    </table>


    <table border="1" cellpadding="10">
    <caption>Driver Single Column Model Namelist Options  (PTS_MODE)</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='scam']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Driver Dry-Deposition Namelist Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='dry-deposition']"/>
    </table>

<hr/>
<p>
</p>
    <h2>DATM namelist variables</h2>
    <p>Note, these all would go into the datm.buildnml.csh file:</p>
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

    <table border="1" cellpadding="10">
    <caption>Data Atmosphere Model Physics Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm_physics']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Data Atmosphere Model Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm']"/>
    </table>

    <p>
    <b>NOTE:</b> The following settings are NOT namelist settings, but
     settings used by datm-build-namelist to figure out actual namelist
     settings.
    </p>
    <table border="1" cellpadding="10">
    <caption>Data Atmosphere Model Internal Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm_internal']"/>
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
    <td colspan="1"><xsl:if test="string-length(@valid_values)>0"><b>Valid Values: </b>
         <xsl:value-of select="@valid_values"/></xsl:if></td>
  </tr>
</xsl:template>

</xsl:stylesheet>
