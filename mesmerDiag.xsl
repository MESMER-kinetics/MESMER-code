<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
  xmlns:svg="http://www.w3.org/2000/svg"
  xmlns:set="http://exslt.org/sets"
  xmlns:exsl="http://exslt.org/common"
  xmlns:math="http://exslt.org/math"
  extension-element-prefixes="set math"
  exclude-result-prefixes="exsl">

  <xsl:variable name ="diagheight" select="300"/>
  <xsl:variable name ="xleft" select="20"/>
  <xsl:variable name="xspacing" select="100"/>  <!--x spacing of species-->
  <xsl:variable name="spwidth" select="40"/>    <!--width of species bars-->
  <xsl:variable name="tswidth" select="30"/>    <!--width of TS bars-->
  <xsl:variable name="dysptitle" select="15"/>  <!--downward displacement of species titles-->
  <xsl:variable name="dyspenergy" select="-8"/> <!--downward displacement of species energies-->
  <xsl:variable name="dytstitle" select="15"/>  <!--upward displacement of transition state titles-->
  <xsl:variable name="dytsenergy" select="5"/>  <!--upward displacement of transition state energies-->
  <xsl:variable name="dxtstitle" select="10"/>  <!--leftward dx from mid point of reactant/product of ts title-->
  
  <xsl:variable name ="debug">
       <xsl:value-of select="count(//me:debugXSLT)"/>
  </xsl:variable>

  <xsl:variable name="TSEnergies"
    select="key('molrefs',//me:transitionState/cml:molecule/@ref)//cml:property[@dictRef='me:ZPE']/cml:scalar" />
  
  <!--The notional energy of a product set without an assigned energy. 
  This is just a magic number, with a value unlikely to actually occur-->
  <xsl:variable name="dummyE" select="99999999"/>
  
  <xsl:key name="molrefs" match="cml:molecule" use="@id"/>

  <!--This contains the ref attributes of the first reactant and the first product molecule
      for all reactions. Duplicates have been removed.-->
  <xsl:variable name="distinctRandP" 
    select="set:distinct(//cml:reaction//cml:reactant[1]/cml:molecule/@ref
     | //cml:reaction//cml:product[1]/cml:molecule/@ref)"/>

  <!--Make a temporary tree which has <Energy> elements with the energies of each node in distinctRandP-->
  <xsl:variable name="RandPEnergies">
    <xsl:for-each select="$distinctRandP">
      <!--The type of the current set- 'reactant' or 'product'-->
      <xsl:variable name="RorP" select="local-name(../..)"/>
      <!--Node set with one or two reactant or product elements-->
      <xsl:variable name="bothSpecies" 
        select="../../../*[local-name()=$RorP]"/>

      <xsl:variable name="bothEnergies"
        select="key('molrefs',$bothSpecies/cml:molecule/@ref)//cml:property[@dictRef='me:ZPE']/cml:scalar" />
      <!--The // in the above expression means that surrounding <cml:property> by a <cml:propertyList> is optional-->
      <xsl:element name="Energy">
        <xsl:choose>
          <!--If the only or both species have energies specified, use their sum.
          If not, use assume 0 for an unspecified reactant, and a special dummy value for the sum of products.-->
          <xsl:when test="
          boolean($RorP = 'reactant')
          or
          (boolean(key('molrefs',$bothSpecies[1]/cml:molecule/@ref)//cml:property[@dictRef='me:ZPE']/cml:scalar)
          and 
          boolean(key('molrefs',$bothSpecies[count($bothSpecies)]/cml:molecule/@ref)//cml:property[@dictRef='me:ZPE']/cml:scalar))
          ">
            <xsl:value-of select="sum($bothEnergies)"/>                    
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$dummyE"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:element>
     </xsl:for-each>
  </xsl:variable>

  <xsl:variable name="okEnergies" select="exsl:node-set($RandPEnergies)/Energy[. != $dummyE] | $TSEnergies"/>
  <xsl:variable name="Emin" select="math:min($okEnergies)"/>
  <xsl:variable name="Emax" select="math:max($okEnergies)"/>

  <!-- energyOffset (subtracted from each non-zero energy)
  If the value of me:diagramEnergyOffset is 0, the lowest energy in the system will be labelled 0.0.
  If a species whose energy is now labelled E needs to be set to 0.0, change the value of
  me:diagramEnergyOffset to E.-->
  <xsl:variable name="energyOffset">
    <xsl:choose>
      <xsl:when test="//me:diagramEnergyOffset">
        <xsl:choose>
          <!--An enhancement to set aspecified reaction's reactants
          as zero energy is not yet working-->
          <xsl:when test="//me:diagramEnergyOffset/@ref">
            <xsl:variable name="posn">
              <xsl:for-each select="//cml:reaction">
                <xsl:if test="@id=//me:diagramEnergyOffset/@ref">
                  <xsl:value-of select="position()"/>
                </xsl:if>
              </xsl:for-each>
            </xsl:variable>
            <xsl:value-of select=" //me:diagramEnergyOffset  
                          + exsl:node-set($RandPEnergies)/Energy[number($posn)]"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$Emin + //me:diagramEnergyOffset"/>
          </xsl:otherwise>
        </xsl:choose>       
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="0"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>

  <!--Calculate an appropriate  vertical scaling factor-->
  <xsl:variable name="yscale">
    <xsl:value-of select="($diagheight - 90) div ($Emax - $Emin)"/>
  </xsl:variable>

  <xsl:variable name ="ybase" select="$diagheight+($Emin * $yscale) - 60"/>

  <xsl:template name="drawDiag" match="me:mesmer">
    <xsl:variable name="wellsRatio">
      <xsl:choose>
        <xsl:when test="(count($distinctRandP) div 10) &gt; 1.0">
          <xsl:value-of select="count($distinctRandP) div 10"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="1"/>
        </xsl:otherwise>                   
      </xsl:choose>
    </xsl:variable>
    <svg:svg  version="1.1" font-family="Verdana" width="100%" preserveAspectRatio="none">
      <xsl:attribute name="height">
        <xsl:value-of select="$diagheight * $wellsRatio + 60" />
      </xsl:attribute>
      
      <!--Have viewBox only with large number of wells--> 
      <xsl:if test="$wellsRatio &gt; 1.0">
        <xsl:attribute name="viewBox">
          <xsl:value-of select="
            concat(0, ' ', 0, ' ', 
            count($distinctRandP)*$xspacing + 40, ' ', $diagheight + 60)" />
        </xsl:attribute>
      </xsl:if>
      
      <xsl:if test="$debug!=0">
        <xsl:call-template name="test3"/>
      </xsl:if>
      
      <!--Draw the energy levels of the modelled molecules-->
      <xsl:call-template name="drawWells"/>
      
      <!--Draw the energy levels of the transition states and the lines to reactant and product-->
      <xsl:apply-templates select="cml:reactionList" mode="diagram"/>

      <xsl:if test="cml:reactionList">
        <svg:text x="20" font-size="9">
          <xsl:attribute name="y">
            <xsl:value-of select="$diagheight"/>
          </xsl:attribute>
          <xsl:value-of select="'Energies in '"/>
          <xsl:value-of select="//cml:property[@dictRef='me:ZPE']/cml:scalar/@units"/>
        </svg:text>
      </xsl:if>
      
    </svg:svg>
  </xsl:template>

  <!--===================================================================-->
  <xsl:template name="drawWells">
    <svg:g style="stroke:teal;"> <!--All wells-->
      <xsl:for-each select="$distinctRandP">
        <xsl:variable name="pos" select="position()"/>
        <!--The type of the current set- 'reactant' or 'product'-->
        <xsl:variable name="RorP" select="local-name(../..)"/>
        <!--Node set with one or two reactant or product elements-->
        <xsl:variable name="bothSpecies" 
          select="../../../*[local-name()=$RorP]"/>
        <xsl:variable name="bothEnergies"
          select="key('molrefs',$bothSpecies/cml:molecule/@ref)//cml:property[@dictRef='me:ZPE']/cml:scalar" />

        <xsl:variable name="yval">
          <xsl:choose>
            <xsl:when test="exsl:node-set($RandPEnergies)/Energy[number($pos)] != $dummyE">
              <xsl:value-of select="$ybase - $yscale * exsl:node-set($RandPEnergies)/Energy[number($pos)]"/>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="$diagheight - 100"/>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:variable>
        
        <svg:g> <!--This well-->
          <xsl:if test="key('molrefs', .)/@active='false'">
            <xsl:attribute name="class">inactive</xsl:attribute>
          </xsl:if>
          <!--Write names of species-->
          <svg:text font-size="9"  >
            <xsl:attribute name="x" >
              <xsl:number value="(position()-1)* $xspacing+ $xleft" />
            </xsl:attribute>
            <xsl:attribute name="y" >
              <xsl:number value="$yval+$dysptitle" />
            </xsl:attribute>

            <xsl:for-each select="$bothSpecies">
              <xsl:value-of select="cml:molecule/@ref"/>
              <xsl:if test="position()!=last()">
                <xsl:text> + </xsl:text>
              </xsl:if>
            </xsl:for-each>
          </svg:text>

          <!--Write energies-->
          <xsl:if test="exsl:node-set($RandPEnergies)/Energy[number($pos)] != $dummyE">
            <xsl:variable name="EWell" select="sum($bothEnergies)" />
            <svg:text font-size="8"  >
              <xsl:attribute name="x" >
                <xsl:number value="(position()-1)* $xspacing+ $xleft + 0.1 * $xspacing"/>
              </xsl:attribute>
              <xsl:attribute name="y" >
                <xsl:number value="$yval+$dyspenergy" />
              </xsl:attribute>
              <xsl:value-of select="format-number($EWell - $energyOffset,'#.0')" />
            </svg:text>

            <!--Draw horizontal bar at the appropriate energy-->
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
          </xsl:if>
        </svg:g>
       </xsl:for-each>
    </svg:g>
  </xsl:template>

  <!--===================================================================--> 
  <xsl:template match="cml:reactionList" mode="diagram">
    <xsl:for-each select="cml:reaction">

      <xsl:variable name="reactantmol" select="key('molrefs', .//cml:reactant/cml:molecule/@ref)" />
      <xsl:variable name="productmol" select="key('molrefs', .//cml:product/cml:molecule/@ref)" />
      <xsl:variable name="TSmol" select="key('molrefs', me:transitionState/cml:molecule/@ref)" />
      <xsl:variable name="EnergyTStemp" select="$TSmol//cml:property[@dictRef='me:ZPE']/cml:scalar"/>
      
      <xsl:variable name="reactantIndex">
        <xsl:for-each select="$distinctRandP">
          <xsl:if test=".=$reactantmol/@id">
            <xsl:value-of select="position()"/>
          </xsl:if>
        </xsl:for-each>
      </xsl:variable>
      
      <xsl:variable name="reactantpos"><!-- middle of bar-->
            <xsl:value-of select="($reactantIndex - 1)* $xspacing+ $xleft + $spwidth*0.5"/>
      </xsl:variable>

      <xsl:variable name="productIndex">
        <xsl:for-each select="$distinctRandP">
          <xsl:if test=".=$productmol/@id">
            <xsl:value-of select="position()"/>
          </xsl:if>
        </xsl:for-each>
      </xsl:variable>
      
      <xsl:variable name="productpos"><!-- middle of bar-->
            <xsl:value-of select="($productIndex - 1)* $xspacing+ $xleft + $spwidth*0.5"/>
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

      <xsl:variable name="yTStemp" select="$ybase - $yscale * $TSmol//cml:property[@dictRef='me:ZPE']/cml:scalar"/>
      <xsl:variable name="isILT" select="me:MCRCMethod='MesmerILT' or me:MCRCMethod='SimpleILT'"/>
      <xsl:variable name="EnergyTS">
        <xsl:choose>
          <xsl:when test="$TSmol">
            <xsl:value-of select="$EnergyTStemp" />
          </xsl:when>
          <!--For ILT source reactions use the activationEnergy + energy of reactants-->
          <!--Look first for alternative format(from OB)-->
          <xsl:when test="$isILT and cml:rateParameters/cml:E[not(@reverse)]">
            <xsl:value-of select="(cml:rateParameters/cml:E + 
                 exsl:node-set($RandPEnergies)/*[number($reactantIndex)])"/>
          </xsl:when>
          <xsl:when test="$isILT and cml:rateParameters/cml:E[@reverse]">
            <xsl:value-of select="(cml:rateParameters/cml:E + 
                 exsl:node-set($RandPEnergies)/*[number($productIndex)])"/>
          </xsl:when>
          <xsl:when test="$isILT and me:activationEnergy[not(@reverse)]">
            <xsl:value-of select="(me:activationEnergy + 
                 exsl:node-set($RandPEnergies)/*[number($reactantIndex)])"/>
          </xsl:when>
          <xsl:when test="$isILT and me:activationEnergy[@reverse]">
            <xsl:value-of select="(me:activationEnergy + 
                 exsl:node-set($RandPEnergies)/*[number($productIndex)])"/>            
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$EnergyTStemp" />
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="yTS" select="$ybase - $yscale * $EnergyTS"/>

      <!--Draw line from reactant, through TS to product-->
      <svg:path stroke ="black" fill="none">
        <xsl:if test="@active='false'">
          <xsl:attribute name="class">inactive</xsl:attribute>
        </xsl:if>
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
          <xsl:value-of select="$ybase - $yscale * exsl:node-set($RandPEnergies)/Energy[number($reactantIndex)]"/>
            

          <xsl:if test="$yTS">
           <!--Only if there is a valid TS-->
            <!--Draw lines and TS level-->
            <xsl:value-of select="concat('L ', ($reactantpos+$productpos)*0.5 - $tsoffset, ' ')"/>
            <xsl:value-of select="$yTS"/>
            <xsl:value-of select="concat(' ',  ($reactantpos+$productpos)*0.5 + $tsoffset, ' ')"/>
            <xsl:value-of select="$yTS"/>
          </xsl:if>
          
          <xsl:value-of select="concat(' L',$productpos - $spoffset,' ')"/>
          <xsl:choose>
            <xsl:when test="exsl:node-set($RandPEnergies)/Energy[number($productIndex)] != $dummyE">
              <xsl:value-of select="$ybase - $yscale * exsl:node-set($RandPEnergies)/Energy[number($productIndex)]"/>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="$diagheight - 100"/>
            </xsl:otherwise>
          </xsl:choose>

        </xsl:attribute>
      </svg:path>
      
     <xsl:if test="$yTS  and string-length($TSmol/@id) &lt; 40">
     <!--Label TS if less than 40 characters-->
       <svg:g>
         <xsl:if test="@active='false'">
           <xsl:attribute name="class">inactive</xsl:attribute>
         </xsl:if>
         <svg:text font-size="9"  >
          <xsl:attribute name="x" >
            <xsl:number value="($reactantpos+$productpos -$tswidth * 0.7)*0.5" />
          </xsl:attribute>
          <xsl:attribute name="y" >
            <xsl:number value="$yTS - $dytstitle" />
          </xsl:attribute>
          <xsl:value-of select="$TSmol/@id"/>
        </svg:text>
         <svg:text font-size="8"  >
           <xsl:attribute name="x" >
             <xsl:number value="($reactantpos+$productpos -$tswidth * 0.7)*0.5" />
           </xsl:attribute>
           <xsl:attribute name="y" >
             <xsl:number value="$yTS - $dytsenergy" />
           </xsl:attribute>
           <xsl:value-of select="format-number($EnergyTS - $energyOffset,'#.0')"/>
         </svg:text>
       </svg:g>
     </xsl:if>

    </xsl:for-each>
  </xsl:template>

  <xsl:template name="test">
    <!--Displays the unique reaction/product sets with the names of both species if present.-->
    <xsl:for-each select="$distinctRandP">
      <svg:text x="160">
        <xsl:attribute name="y">
          <xsl:value-of select="40+position()*20"/>
        </xsl:attribute>
        <xsl:value-of select="."/>
        <xsl:choose>
          <xsl:when test=". = ../../../cml:reactant[1]/cml:molecule/@ref">
            <!--Above: better way to test whether the parent node is a cml:reactant or a cml:product?-->
                <xsl:if test="count(../../../cml:reactant)>1">
                  <xsl:text> + </xsl:text>
                </xsl:if>
            <xsl:value-of select="../../../cml:reactant[2]/cml:molecule/@ref"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:if test="count(../../../cml:product)>1">
              <xsl:text> + </xsl:text>
            </xsl:if>
            <xsl:value-of select="../../../cml:product[2]/cml:molecule/@ref"/>
          </xsl:otherwise>
        </xsl:choose>
        
      </svg:text>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="test2">
    <xsl:for-each select="//cml:reaction">
      <svg:text x="160">
        <xsl:attribute name="y">
          <xsl:value-of select="40+position()*20"/>
        </xsl:attribute>
        <xsl:value-of select="cml:reactant[1]/cml:molecule/@ref"/>
        <xsl:if test="count(cml:reactant)">
          <xsl:text> + </xsl:text>
        </xsl:if>
        <xsl:value-of select="cml:reactant[2]/cml:molecule/@ref"/>
      </svg:text>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="test3">
    <svg:text x="0">
      <xsl:attribute name="y">
        <xsl:value-of select="$diagheight"/>
      </xsl:attribute>
      $energyOffset=
      <xsl:value-of select="$energyOffset"/>
      $Emin=
      <xsl:value-of select="$Emin"/>
      Count RandP sets =
      <xsl:value-of select="count($distinctRandP)"/>
      Count TS =
      <xsl:value-of select="count(key('molrefs',//me:transitionState/cml:molecule/@ref))"/>
      Count TSEnergies =
      <xsl:value-of select="count($TSEnergies)"/>
    </svg:text>
    <svg:text x="0">
      <xsl:attribute name="y">
        <xsl:value-of select="$diagheight + 20"/>
      </xsl:attribute>
      yscale =
      <xsl:value-of select="$yscale"/>
      Emax =
      <xsl:value-of select="$Emax"/>
      Emin =
      <xsl:value-of select="$Emin"/>
      Count energy levels =
      <xsl:value-of select="count(exsl:node-set($RandPEnergies)/Energy | $TSEnergies)"/>
    </svg:text>
    <svg:text x="0">
      <xsl:attribute name="y">
        <xsl:value-of select="$diagheight + 40"/>
      </xsl:attribute>
      <xsl:for-each select="exsl:node-set($RandPEnergies)/Energy">
        <xsl:value-of select="current()"/>,
      </xsl:for-each>
    </svg:text>
    
  </xsl:template>

 
</xsl:stylesheet>
