<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
  xmlns:svg="http://www.w3.org/2000/svg">

  <xsl:variable name ="ybase" select="260"/>
  <xsl:variable name ="xleft" select="20"/>
  <xsl:variable name="xspacing" select="100"/>  <!--x spacing of species-->
  <xsl:variable name="spwidth" select="40"/>    <!--width of species bars-->
  <xsl:variable name="tswidth" select="30"/>    <!--width of TS bars-->
  <xsl:variable name="dysptitle" select="15"/>  <!--downward displacement of species titles-->
  <xsl:variable name="dyspenergy" select="-8"/> <!--downward displacement of species energies-->
  <xsl:variable name="dytstitle" select="15"/>  <!--upward displacement of transition state titles-->
  <xsl:variable name="dytsenergy" select="5"/>  <!--upward displacement of transition state energies-->
  <xsl:variable name="dxtstitle" select="10"/>  <!--leftward dx from mid point of reactant/product of ts title-->

  <xsl:variable name="yscale">
    <xsl:for-each select="//cml:molecule[@me:type='modelled' or @me:type='transitionState' 
      or @me:type='displayed']/cml:propertyList/cml:property[@dictRef='me:ZPE']">
      <xsl:sort select="cml:scalar" data-type="number"/>
      <xsl:if test="position()=last()">
        <xsl:value-of select="($ybase - 40) div cml:scalar "/>
      </xsl:if>
    </xsl:for-each>
  </xsl:variable>

  <xsl:key name="molrefs" match="cml:molecule" use="@id"/>
  
  <xsl:template name="drawDiag" match="me:mesmer">
    <svg:svg  version="1.1">
      <xsl:attribute name="viewBox">
        <xsl:value-of select="concat('0 0 ', 
        100+80*(count(//cml:molecule[@me:type='modelled' or @me:type='displayed']) -1), ' 300')"/>
      </xsl:attribute>
      <!--<xsl:call-template name="drawYaxis"/>-->
      
      <!--Draw the energy levels of the modelled molecules-->
      <xsl:apply-templates select="cml:moleculeList" mode="diagram"/>
      
      <!--Draw the energy levels of the transition states and the lines to reactant and product-->
      <xsl:apply-templates select="cml:reactionList" mode="diagram"/>
    </svg:svg>
  </xsl:template>

  <xsl:template match="cml:moleculeList" mode="diagram">
    <xsl:variable name="mols" select="cml:molecule[@me:type='modelled' or @me:type='displayed']"/>

    <svg:g style="stroke:teal">
      <xsl:for-each select="$mols">
        <xsl:variable name="yval"
           select="$ybase - $yscale * cml:propertyList/cml:property[@dictRef='me:ZPE']/cml:scalar"/>
        <svg:text font-family="Verdana" font-size="9"  >
          <xsl:attribute name="x" >
            <xsl:number value="(position()-1)* $xspacing+ $xleft" />
          </xsl:attribute>
          <xsl:attribute name="y" >
            <xsl:number value="$yval+$dysptitle" />
          </xsl:attribute>
          <xsl:value-of select="@id"/>
        </svg:text>
        <svg:text font-family="Verdana" font-size="8"  >
          <xsl:attribute name="x" >
            <xsl:number value="(position()-1)* $xspacing+ $xleft + 0.1 * $xspacing"/>
          </xsl:attribute>
          <xsl:attribute name="y" >
            <xsl:number value="$yval+$dyspenergy" />
          </xsl:attribute>
          <xsl:value-of select="cml:propertyList/cml:property[@dictRef='me:ZPE']/cml:scalar"/>
        </svg:text>
        <svg:line style="stroke-width:3">
          <xsl:attribute name="x1">
            <xsl:number value="(position()-1)* $xspacing+ $xleft" />
          </xsl:attribute>
          <xsl:attribute name="y1">
            <xsl:value-of select="$yval"/>
          </xsl:attribute>
          <xsl:attribute name="x2">
            <xsl:number value="(position()-1)* $xspacing+ $spwidth + $xleft"/>
          </xsl:attribute>
          <xsl:attribute name="y2">
            <xsl:value-of select="$yval"/>
          </xsl:attribute>
        </svg:line>
      </xsl:for-each>
    </svg:g>

  </xsl:template>


  <xsl:template match="cml:reactionList" mode="diagram">
    <xsl:variable name="mols" select="//cml:molecule[@me:type='modelled' or @me:type='displayed']"/>
    <xsl:for-each select="cml:reaction">

      <xsl:variable name="reactantmol" select="key('molrefs', cml:reactant/cml:molecule/@ref)" />
      <xsl:variable name="productmol" select="key('molrefs', cml:product/cml:molecule/@ref)" />
      <xsl:variable name="TSmol" select="key('molrefs', me:transitionState/cml:molecule/@ref)" />
      <xsl:variable name="EnergyTS" select="$TSmol/cml:propertyList/cml:property[@dictRef='me:ZPE']/cml:scalar"/>
      <xsl:variable name="yTS" select="$ybase - $yscale * $TSmol/cml:propertyList/cml:property[@dictRef='me:ZPE']/cml:scalar"/>

      <!--There must be a better way of obtaining the position in mols!-->
      <xsl:variable name="reactantpos"><!-- middle of bar-->
        <xsl:for-each select="$mols">
          <xsl:if test="@id=$reactantmol/@id">
            <xsl:value-of select="(position()-1)* $xspacing+ $xleft + $spwidth*0.5"/>
          </xsl:if>
        </xsl:for-each>
      </xsl:variable>

      <xsl:variable name="productpos"><!-- middle of bar-->
        <xsl:for-each select="$mols">
          <xsl:if test="@id=$productmol/@id">
            <xsl:value-of select="(position()-1)* $xspacing+ $xleft + $spwidth*0.5"/>
          </xsl:if>
        </xsl:for-each>
      </xsl:variable>

      <!-- offsets -ve when reactants are to the left of the products-->
      <xsl:variable name="spoffset">
        <xsl:choose>
          <xsl:when test="$productpos - $reactantpos &lt; 0">
            <xsl:value-of select="-$spwidth*0.5"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$spwidth*0.5"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="tsoffset">
        <xsl:choose>
          <xsl:when test="$productpos - $reactantpos &lt; 0">
            <xsl:value-of select="-$tswidth*0.5"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$tswidth*0.5"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>

      <!--Draw line from reactant, through TS to product-->
      <svg:path stroke ="black" fill="none">

        <!--Make line dashed if reactant and product are not adjacent-->
        <xsl:if test = "$productpos - $reactantpos &gt;  $xspacing
                     or $productpos - $reactantpos &lt; - $xspacing">
          <xsl:attribute name="stroke-dasharray">
            <xsl:value-of select="'4,2'"/>
          </xsl:attribute>
        </xsl:if>

        <xsl:attribute name="d">
          <xsl:value-of select="'M'" />
          <xsl:value-of select="concat($reactantpos  + $spoffset,' ')"/>
          <xsl:value-of select="$ybase - $yscale * 
            $reactantmol[@me:type='modelled' or @me:type='displayed']/cml:propertyList/cml:property[@dictRef='me:ZPE']/cml:scalar"/>

          <xsl:if test="$yTS">
          <!--Only if there is a valid TS-->
            <!--Draw lines and TS level-->
            <xsl:value-of select="concat('L ', ($reactantpos+$productpos)*0.5 - $tsoffset, ' ')"/>
            <xsl:value-of select="$yTS"/>
            <xsl:value-of select="concat(' ',  ($reactantpos+$productpos)*0.5 + $tsoffset, ' ')"/>
            <xsl:value-of select="$yTS"/>
          </xsl:if>
          
          <xsl:value-of select="concat(' L',$productpos - $spoffset,' ')"/>
          <xsl:value-of select="$ybase - $yscale * $productmol[@me:type='modelled' or @me:type='displayed']/cml:propertyList/cml:property[@dictRef='me:ZPE']/cml:scalar"/>
        </xsl:attribute>
      </svg:path>
      
     <xsl:if test="$yTS  and string-length($TSmol/@id) &lt; 40">
     <!--Label TS if less than 40 characters-->
      <svg:text font-family="Verdana" font-size="9"  >
        <xsl:attribute name="x" >
          <xsl:number value="($reactantpos+$productpos -$tswidth * 0.7)*0.5" />
        </xsl:attribute>
        <xsl:attribute name="y" >
          <xsl:number value="$yTS - $dytstitle" />
        </xsl:attribute>
        <xsl:value-of select="$TSmol/@id"/>
      </svg:text>
       <svg:text font-family="Verdana" font-size="8"  >
         <xsl:attribute name="x" >
           <xsl:number value="($reactantpos+$productpos -$tswidth * 0.7)*0.5" />
         </xsl:attribute>
         <xsl:attribute name="y" >
           <xsl:number value="$yTS - $dytsenergy" />
         </xsl:attribute>
         <xsl:value-of select="$EnergyTS"/>
       </svg:text>
     </xsl:if>

    </xsl:for-each>
  </xsl:template>

  <xsl:template name="drawYaxis">
    <svg:line stroke="black" x1="40" x2="40" y1="20">
      <xsl:attribute name="y2">
        <xsl:value-of select="$ybase"/>
      </xsl:attribute>
    </svg:line>
    <xsl:call-template name="drawTicks">
      <xsl:with-param name="yabove" select="$ybase"/>
    </xsl:call-template>
  </xsl:template>

  <xsl:template name="drawTicks">
    <xsl:param name="yabove"/>
  </xsl:template>  
  </xsl:stylesheet> 
