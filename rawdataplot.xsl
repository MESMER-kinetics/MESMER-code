<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
  xmlns:svg="http://www.w3.org/2000/svg"
  xmlns:xlink="http://www.w3.org/1999/xlink"
  xmlns:exsl="http://exslt.org/common"
  xmlns:math="http://exslt.org/math"
  xmlns:set="http://exslt.org/sets"
  xmlns:str="http://exslt.org/strings"
  extension-element-prefixes="exsl math set str"
>
<!--This version requires that the data values are not in scientific format-->
  <xsl:template name="rawDataPlot" match="me:rawData" >
    <xsl:variable name="splittimes" select="str:split(me:times)"/>
    <xsl:variable name="splitsignals" select="str:split(me:signals)"/>

    <xsl:variable name="splitcalcsignals" select="str:split(me:calcSignals)"/>
    <xsl:variable name="maxval" select ="math:max($splitsignals)"/>
    <xsl:variable name="minval" select ="math:min($splitsignals)"/>
    <xsl:variable name="maxtime" select ="math:max($splittimes)"/>
    <xsl:variable name="mintime" select ="math:min($splittimes)"/>

    <svg:svg version="1.1" width="340px" height="220px" >
      <!--title of graph-->
      <svg:text x="0" y="14"  font-size="11">
        
        <xsl:value-of select="concat(@name,' T=',../@T,'K P=',../@P,../@units,' ',../@excessReactantConc)"/>
      </svg:text>
      
      <!--frame of graph-->
      <svg:rect x="1" y="21" width="338" height="180" fill="none" stroke="black"/>
      <svg:text x="5" y="212" font-size="10" text-anchor="middle">
        <xsl:value-of select="$mintime"></xsl:value-of>
      </svg:text>
      <svg:text x="330" y="212" font-size="10" text-anchor="middle">
        <xsl:value-of select="$maxtime"></xsl:value-of>
      </svg:text>
      <svg:text x="160" y="212" font-size="10" text-anchor="middle">
        <xsl:value-of select ="../@timeUnits"/>
      </svg:text>
      
      <!--plot calculated points-->
        <xsl:variable name="pathtext">
          <xsl:value-of select="'M '"/>
          <xsl:for-each select="$splittimes">
            <xsl:variable name="pos" select="position()"/>
            <xsl:value-of select="concat(., ' ',-$splitsignals[$pos], ' L ')"/>
          </xsl:for-each>
        </xsl:variable>
      
      <!--Viewbox puts -ymax at top and -ymin at bottom. y values are negated before plotting.-->
      <!--viewBox(minx, -miny, width, height) for a normal graph-->
      <svg:svg x="1" y="21" width="338" height="180"  preserveAspectRatio="none" >

        <xsl:attribute name="viewBox">
          <!--xsl:value-of select="concat(0, -$maxval, $maxtime, $maxval - $minval)"/>-->
          <xsl:value-of select="concat($mintime, ' -',$maxval, ' ',$maxtime - $mintime,' ', $maxval - $minval)"/>
        </xsl:attribute>

        <svg:path stroke="red" fill="none" stroke-width="1" vector-effect="non-scaling-stroke" >
          <xsl:attribute name="d">
            <xsl:value-of select="$pathtext"/>
          </xsl:attribute>
        </svg:path>

        <xsl:if test="me:calcSignals">
          <svg:path stroke="blue" fill="none" stroke-width="2" vector-effect="non-scaling-stroke">
            <xsl:attribute name="d">
              <xsl:value-of select="'M '"/>
              <xsl:for-each select="$splittimes">
                <xsl:variable name="pos" select="position()"/>
                <xsl:value-of select="concat(., ' -',$splitcalcsignals[$pos], ' L ')"/>
              </xsl:for-each>
            </xsl:attribute>
          </svg:path>
        </xsl:if>
      </svg:svg>
      
   </svg:svg>
  
  </xsl:template>

</xsl:stylesheet>
