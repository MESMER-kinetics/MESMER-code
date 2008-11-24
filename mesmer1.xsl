<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:me="http://www.chem.leeds.ac.uk/mesmer">

  <xsl:include href="mesmerDiag.xsl"/>
  <xsl:include href="switchcontent.xsl"/>
  
<xsl:key name="molrefs" match="cml:molecule" use="@id"/>
  <xsl:variable name="title">
    <xsl:choose>
      <xsl:when test="//cml:metadataList">
        <xsl:value-of select="cml:metadata[@name='dc:title']/@content"/>
      </xsl:when>
      <xsl:otherwise>
        Mesmer datafile
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>
<xsl:template match="me:mesmer">
  <html>
    <head>
      <title>
        <xsl:value-of select="$title"/>
      </title>
      <script type="text/javascript">
        <xsl:value-of select="$importedjavascript"/>
        <!--If we had used src="switchcontent.js" it would have been relative to the
        position of the xml file. The href in xsl:include is relative to the xsl file.-->
      </script>

      <style>
        <![CDATA[
        body{margin:20px;padding:0;}
        table.mol{border-spacing:10px;}
        .name{font-weight:bold;}
        .tableheader
        {
          font-family: Arial, Helvetica, sans-serif;
          font-weight:bold;
          text-align:center;
        }
        h3{color:teal;font-family: Arial, Helvetica, sans-serif;}
        .normal{color:black; font-size:smaller;}
        .handcursor{cursor:hand; cursor:pointer;}
        #header{color:black;font-family: Arial, Helvetica, sans-serif;font-weight:bold;}
        #title{font-size:larger;font-weight:bold;}
        #metadata{color:teal;font-size:smaller;}
        ]]>
      </style>
    </head>
    <body>
      <div id="header">
        <xsl:apply-templates select="//cml:metadataList"/>
      </div>
      <h3 id="mols-title" class="handcursor">Molecules</h3>
      <div id="mols" class="switchgroup3">
        <table class="mol">
          <tr class="tableheader">
            <td>Name</td>
            <td>Energy/kJ mol<sup>-1</sup>
          </td>
            <td>Rotational constants/cm<sup>-1</sup>
          </td>
            <td>Vibrational frequencies/cm<sup>-1</sup>
          </td>
          </tr>
          <xsl:apply-templates select="cml:moleculeList"/>
        </table>
      </div>
      <h3 id="reactions-title" class="handcursor">Reactions</h3>
      <table id="reactions" class="switchgroup4">
      <xsl:apply-templates select="cml:reactionList"/>
    </table>

    <!--Show the "results"-->
    <h3 id="densityOfStates-title" class="handcursor">Density of States</h3>
    <div id="densityOfStates" class="switchgroup1">
      <!--<xsl:apply-templates select="//*[@calculated]"/>-->
      <xsl:apply-templates select="//me:densityOfStatesList"/>
    </div>

    <h3 id="microRates-title" class="handcursor">Microcanonical Rate Coefficients</h3>
    <div id="microRates" class="switchgroup2">
      <xsl:apply-templates select="//me:microRateList"/>
    </div>

    <!--Script for expanding an contracting sections-->
    <script type="text/javascript">
      <![CDATA[
        for(var i=1; i <=4; i++)
        {
          var mc=new switchcontent("switchgroup" + i)
          mc.setStatus('- ','+ ')
          mc.setPersist(true)
          mc.init()
        }
      ]]>
    </script>

    <xsl:call-template name="drawDiag"/>
  </body>
    </html>
</xsl:template>

  
  <xsl:template match="cml:molecule">
    <tr>
      <td class="name">
        <xsl:value-of select="@id"/>
      </td>
      <td>
        <xsl:choose>
          <xsl:when test=".//cml:property[@dictRef='me:ZPE']/cml:scalar">
            <xsl:value-of select=".//cml:property[@dictRef='me:ZPE']/cml:scalar"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:if test=".//cml:property[@dictRef='me:Hf298']/cml:scalar">
              <xsl:value-of select=".//cml:property[@dictRef='me:Hf298']/cml:scalar"/>
              <xsl:value-of select="'(Hf298)'"/>
            </xsl:if>
          </xsl:otherwise>
        </xsl:choose>
      </td>
      <td>
        <xsl:value-of select=".//cml:property[@dictRef='me:rotConsts']/cml:array"/>
      </td>
      <td>
        <xsl:value-of select=".//cml:property[@dictRef='me:vibFreqs']/cml:array"/>
      </td>
    </tr>
  </xsl:template>
  
  <xsl:template match="cml:reaction">
    <tr>
      <td>
        <xsl:value-of select="@id"/>
      </td>
      <td>
        <xsl:for-each select=".//cml:reactant/cml:molecule/@ref">
          <xsl:if test="position()!=1"> + </xsl:if>
          <xsl:value-of select="."/>
        </xsl:for-each>
      </td>
      <td>
        <xsl:if test="@reversible='true'">&lt;</xsl:if>=></td>
      <td>
        <xsl:for-each select=".//cml:product/cml:molecule/@ref">
          <xsl:if test="position()!=1"> + </xsl:if>
          <xsl:value-of select="."/>
        </xsl:for-each>
      </td>
      <td>
        (
        <xsl:if test="not(.//me:transitionState/cml:molecule/@ref)">No </xsl:if>
        Transition State <xsl:value-of select=".//me:transitionState/cml:molecule/@ref"/>)
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="me:densityOfStatesList">
    <h4>
      Density of states for
      <xsl:value-of select="../@id"/>
      <span class="normal">
        (Calculated 
        <xsl:value-of select="@calculated"/>
        )
      </span>
    </h4>
    <table>
      <tr class="tableheader">
        <td>T</td>
        <td>qtot</td>
        <td>sumc</td>
        <td>sumg</td>
      </tr>
      <xsl:for-each select="me:densityOfStates">
        <tr>
          <td>
            <xsl:value-of select="me:T"/>
          </td>
          <td>
            <xsl:value-of select="me:qtot"/>
          </td>
          <td>
            <xsl:value-of select="me:sumc"/>
          </td>
          <td>
            <xsl:value-of select="me:sumg"/>
          </td>
        </tr>
      </xsl:for-each>
    </table>
  </xsl:template>

    <xsl:template match="me:microRateList">
    <h4>
      Microcanonical rate coefficients for
      <xsl:value-of select="../@id"/>
      <span class="normal">
        (Calculated 
        <xsl:value-of select="@calculated"/>
        )
      </span>
    </h4>
    <table>
      <tr class="tableheader">
        <td>T/K</td>
        <td>Rate/s-1</td>
      </tr>
      <xsl:for-each select="me:microRate">
        <tr>
          <td>
            <xsl:value-of select="me:T"/>
          </td>
          <td>
            <xsl:value-of select="me:val"/>
          </td>
        </tr>
      </xsl:for-each>
    </table>
  </xsl:template>

  <xsl:template match="cml:metadataList">
      <p id="title"><xsl:value-of select="cml:metadata[@name='dc:title']/@content"/></p>
    <div id="metadata">
        <xsl:value-of select="cml:metadata[@name='dc:creator']/@content"/>:
        <xsl:value-of select="cml:metadata[@name='dc:date']/@content"/>,
        <xsl:value-of select="cml:metadata[@name='dc:contributor']/@content"/>
    </div>
  </xsl:template>

  </xsl:stylesheet> 
