<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
  xmlns:svg="http://www.w3.org/2000/svg"
  xmlns:exsl="http://exslt.org/common" >
 
  <xsl:variable name="colors">
    <c>red</c>
    <c>green</c>
    <c>blue</c>
    <c>orange</c>
    <c>black</c>
    <c>teal</c>
  </xsl:variable>
  <xsl:variable name="popspecies" select="//me:analysis[1]/me:populationList[1]/me:population[1]/me:pop/@ref"/>

<!--  <xsl:template match="/">
        <xsl:apply-templates select="//me:analysis" mode="diagram"/>
      </xsl:template>
-->  
  <xsl:template match="//me:analysis" mode="diagram">
    <p class="paramheader">
      <xsl:for-each select="me:parameters/@*">
        <xsl:value-of select="concat(name(),'=',.,'  ')"/>
      </xsl:for-each>
    </p>
    <xsl:apply-templates select="me:populationList" mode="diagram"/>
  </xsl:template>
  <xsl:template name="populationDiagram" match="me:populationList" mode="diagram">
    <xsl:variable name="p" select="me:population"/> <!--Save because context node will change in for-eachs-->
    <xsl:if test="position()=1">
      <xsl:call-template name="legend"/>
    </xsl:if>
  
    <!--Frame of graph-->
    <svg:svg version="1.1" width="200px" height="220px">
    <svg:text x="30" y="14" font-family="Verdana" font-size="13">
      <xsl:value-of select="concat(@T,'K ',@conc)"/>
    </svg:text>
    <svg:rect x="1" y="21" width="198" height="158" fill="white" stroke="black" stroke-width="1"/>
    <xsl:variable name="startval" select="me:population[1]/@logTime"/>
    <xsl:variable name="endval" select="me:population[last()]/@logTime"/>
    <xsl:call-template name="xaxis">
      <xsl:with-param name="val" select="$startval + 1.0"/><!--label from -->
      <xsl:with-param name="maxval" select="$endval"/>
      <xsl:with-param name="pxperstep" select="(200 * 2) div ($endval - $startval)"/>
      <xsl:with-param name="xpx" select="200 div ($endval - $startval)"/>
      <xsl:with-param name="ypx" select="180"/>
    </xsl:call-template>
    <svg:text text-anchor="middle"  font-family="Verdana" font-size="10">
      <xsl:attribute name="x">
        <xsl:value-of select="100"/>
      </xsl:attribute>
      <xsl:attribute name="y">
        <xsl:value-of select="208"/>
      </xsl:attribute>
      log10(time/secs)
    </svg:text>

    <!--Contents of graph-->    
    <svg:svg y="20" height="160" preserveAspectRatio="none">
        <xsl:attribute name="viewBox">
          <xsl:value-of select="concat($startval, ' -1.0 ', $endval - $startval, ' 1.0')"/>
        </xsl:attribute>
        <xsl:for-each select="me:population[1]/me:pop/@ref"> <!--for each species-->
          <svg:path fill="none" stroke-width="0.02" >
            <xsl:attribute name="stroke">
              <xsl:variable name="pos" select="position()"/>
              <xsl:value-of select="exsl:node-set($colors)/c[$pos]"/>
            </xsl:attribute>
            <xsl:attribute name="d">
              <xsl:variable name="cursp" select="."/>
              <xsl:value-of select="concat('M ', $p/@logTime ,' -1.0 L ')"/>
              <xsl:for-each select="$p/me:pop[@ref=$cursp]">
                <xsl:value-of select="concat(../@logTime,' -', ., ' ')"/><!--see note at end-->
              </xsl:for-each>
            </xsl:attribute>
          </svg:path>
        </xsl:for-each>
      </svg:svg>
    </svg:svg>
  </xsl:template>
  
  <xsl:template name="legend">
    <svg:svg width="160" height="200">
      <xsl:for-each select="$popspecies">
        <xsl:variable name="pos" select="(position()-1)"/>
        <svg:line stroke-width="5" x1="10" x2="40">
          <xsl:attribute name="y1">
            <xsl:value-of select="15 + $pos*20"/>
          </xsl:attribute>
          <xsl:attribute name="y2">
            <xsl:value-of select="15 + $pos*20"/>
          </xsl:attribute>
          <xsl:attribute name="stroke">
            <xsl:value-of select="exsl:node-set($colors)/c[$pos+1]"/>
          </xsl:attribute>
        </svg:line>
        <svg:text x="50">
          <xsl:attribute name="y">
            <xsl:value-of select="20 + $pos*20"/>
          </xsl:attribute>
          <xsl:value-of select="."/>
        </svg:text>
      </xsl:for-each>

      <svg:g font-family="Verdana" font-size="10">
        <svg:text x="140" y="8">1.0</svg:text>
        <svg:text x="140" y="160">0.0</svg:text>
        <svg:text x="-136" y="150" transform="rotate(-90)">fractional population</svg:text>
      </svg:g>   
   </svg:svg>
  </xsl:template>

  <!--Called recursively to add the tick marks and labels on the x-axis-->
  <xsl:template name="xaxis">
    <xsl:param name="val"/>
    <xsl:param name="maxval"/>
    <xsl:param name="xpx"/>
    <xsl:param name="ypx"/>   
    <xsl:param name="pxperstep"/>
    <xsl:if test="$val &lt; $maxval">
      <!--tick-->
      <svg:path stroke="black">
        <xsl:attribute name="d">
          <xsl:value-of select="concat(' M ', $xpx, ' ', $ypx)"/>
          <xsl:value-of select="concat(' l 0 ', 5)"/>
        </xsl:attribute>
      </svg:path>
      <!--tick label-->
      <svg:text text-anchor="middle"  font-family="Verdana" font-size="10">
        <xsl:attribute name="x">
          <xsl:value-of select="$xpx"/>
        </xsl:attribute>
        <xsl:attribute name="y">
          <xsl:value-of select="$ypx + 15"/>
        </xsl:attribute>
        <xsl:value-of select="$val"/>
      </svg:text>
      <xsl:call-template name="xaxis">
        <xsl:with-param name="val" select="$val + 2"/>
        <xsl:with-param name="maxval" select="$maxval"/>
        <xsl:with-param name="pxperstep" select="$pxperstep"/>
        <xsl:with-param name="xpx" select="$xpx + $pxperstep"/>
        <xsl:with-param name="ypx" select="$ypx"/>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>
  
  <!-- Probably because XPath 1.0 doesn't use scientific numbers, some numbers give NaN when  
       arithmetic is done on them in XSLT before passing on to SVG. Consequently, the minus
       sign needed to invert the y axis is provided separately and the number is presumably
       handled as a string.-->
</xsl:stylesheet>
