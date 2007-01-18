<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:me="http://www.chem.leeds.ac.uk/mesmer">

<xsl:key name="molrefs" match="cml:molecule" use="@id"/>
<xsl:template match="me:mesmer">
  <html>
    <head>
      <title>Mesmer datafile</title>
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
        h3{color:red;font-family: Arial, Helvetica, sans-serif;}
        .normal{color:black; font-size:smaller;}
        ]]>
      </style>
    </head>
    <body>
      <h3>Molecules</h3>
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
    <h3>Reactions</h3>
    <table>
      <xsl:apply-templates select="cml:reactionList"/>
    </table>
    
    <!--<xsl:apply-templates select="cml:reactionList" mode="diagram"/>-->
      
    <!--Show the "results"-->
    <xsl:apply-templates select="//*[@calculated]"/>
      
    </body>
    </html>
</xsl:template>

  
  <xsl:template match="cml:molecule">
    <tr>
      <td class="name">
        <xsl:value-of select="@id"/>
      </td>
      <td>
        <xsl:value-of select=".//cml:property[@dictRef='me:ZPE']/cml:scalar"/>
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
    </tr>
  </xsl:template>

  <xsl:template match="cml:reactionList" mode="diagram">
    <svg:svg xmlns:svg="http://www.w3.org/2000/svg" version="1.1"
      width="700px" height="400px">
      <svg:g  stroke="black" stroke-width="1">
        <xsl:for-each select="cml:reaction">
          <xsl:variable name="rmols" select="key('molrefs','cml:reactant/cml:molecule/@ref')"/>
          <xsl:variable name="E" select="$rmols/cml:property[@dictRef='me:Energy']/cml:scalar"/>
          <svg:line>
            <xsl:attribute name="x1">
              <xsl:number value="(position()-1)*100" />
            </xsl:attribute>
            <xsl:attribute name="y1">
              <xsl:number value="200"/>
            </xsl:attribute>
            <xsl:attribute name="x2">
              <xsl:number value="(position()-1)*100+ 20"/>
            </xsl:attribute>
            <xsl:attribute name="y2">
              <xsl:number value="200"/>
            </xsl:attribute>
          </svg:line>
          <svg:text  y="250" font-family="Verdana" font-size="15" fill="none">
            <xsl:attribute name="x">
              <xsl:number value="(position()-1)*100" />
            </xsl:attribute>
            <xsl:for-each select="cml:reactant/cml:molecule/@ref">
              <xsl:if test="position()!=1"> + </xsl:if>
              <xsl:value-of select="."/>
            </xsl:for-each>

          </svg:text>
        </xsl:for-each>
      </svg:g>
    </svg:svg>
  </xsl:template>

  <xsl:template match="me:densityOfStatesList">
    <h3>
      Density of states for
      <xsl:value-of select="../@id"/>
      <span class="normal">
        (Calculated 
        <xsl:value-of select="@calculated"/>
        )
      </span>
    </h3>
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
    <h3>
      Microcanonical rate coefficients for
      <xsl:value-of select="../@id"/>
      <span class="normal">
        (Calculated 
        <xsl:value-of select="@calculated"/>
        )
      </span>
    </h3>
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

</xsl:stylesheet> 
