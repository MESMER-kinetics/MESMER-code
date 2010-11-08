<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
  xmlns:dc="http://purl.org/dc/elements/1.1/">

  <xsl:include href="mesmerDiag.xsl"/>
  <xsl:include href="switchcontent.xsl"/>
  <xsl:include href="popDiag.xsl"/>

  <xsl:key name="molrefs" match="cml:molecule" use="@id"/>
  
<xsl:variable name="title">
  <xsl:choose>
    <xsl:when test="//cml:metadataList">
      <xsl:value-of select="//cml:metadataList/dc:title"/>
    </xsl:when>
    <xsl:when test="//cml:title">
      <xsl:value-of select="//cml:title"/>
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
      <script type="text/javascript">
        <![CDATA[
          function toggle()
          {
            var oRule = document.getElementsByTagName('style')[1].sheet.cssRules[0];
            if(oRule.style.visibility=='hidden')
            {
            oRule.style.visibility = 'visible';
            }
            else
            {
            oRule.style.visibility = 'hidden';
            }
          }
          ]]>
      </script>
      <style>
        <![CDATA[
        body{margin:20px;padding:0;}
        table.mol{border-spacing:10px;}
        .name{font-weight:bold;}
        .tableheader
        {
          font-family: Arial, Helvetica, sans-serif;
          font-size:smaller;
          
          text-align:center;
        }
        .tablehead1
        {
          font-family: Arial, Helvetica, sans-serif;
          color:black;
          font-weight:bold;
          
        }
        .tablehead2{text-decoration:underline;padding-top:10px;}
        .tablehead3{font-weight:bold;padding-top:10px;text-align:center}
        .paramheader{
          font-family: Arial, Helvetica, sans-serif;
          color:black;
          font-weight:bold;
          text-decoration: underline;
        }
        table{background-color:#e0f8f8; margin-bottom:12px;}
        td{padding:0px 4px;}
        h3{color:teal;font-family: Arial, Helvetica, sans-serif;}
        hh5{color:black;font-family: Arial, Helvetica, sans-serif;font-weight:bold;font-size:smaller;}
        .normal{color:black; font-size:smaller;}
        .handcursor{cursor:hand; cursor:pointer;}
        .inactive{color:silver;stroke:silver;}
        .error{font-weight:bold;font-size:large;background-color:red;padding:20px;}
        #header{color:black;font-family: Arial, Helvetica, sans-serif;font-weight:bold;}
        #title{font-size:larger;font-weight:bold;}
        #description{color:black;font-size:60%;}
        #metadata{color:teal;font-size:60%;}
        #hide{font-size:small;text-decoration:underline;color:blue;cursor:pointer;}
        #Punchout{font-size:12px;color:gray;}
        ]]>
      </style>
      <xsl:variable name="inactval">
        <xsl:choose>
          <xsl:when test="//me:hideInactive">
            <xsl:value-of select="'hidden'"/>
          </xsl:when>         
          <xsl:otherwise>
            <xsl:value-of select="'visible'"/>
          </xsl:otherwise>               
        </xsl:choose>      
      </xsl:variable>
      <style>
        <xsl:value-of select="concat('.inactive{visibility:', $inactval, ';}')"/>
      </style>
    </head>
    <body>
      <div id="header">
        <p id="title">
          <xsl:value-of select="$title"/>
        </p>
        <p id="description">
          <xsl:value-of select="/me:mesmer/cml:description"/>
        </p>
        <xsl:apply-templates select="//cml:metadataList"/>
      </div>
      <xsl:variable name ="Eunits">
        <xsl:call-template name="ZPEunits"/>
      </xsl:variable>
      <xsl:for-each select=
         "//cml:property[@dictRef='me:ZPE']/cml:scalar/@units | //me:activationEnergy/@units">
        <xsl:if test=".!=$Eunits">
          <div class="error">
            <xsl:value-of select=
             "concat(.,' is different energy unit. All energy units need to be the same.')"/>
          </div>
        </xsl:if>
      </xsl:for-each>
      <h3 id="mols-title" class="handcursor">Molecules</h3>
      <div id="mols" class="switchgroup3">
        <table class="mol">
          <tr class="tableheader">
            <td>Name</td>
            <td>Energy<br /><xsl:value-of select="$Eunits"/></td>
            <td>Rotational constants<br />cm<sup>-1</sup></td>
            <td>Vibrational frequencies<br />cm<sup>-1</sup></td>
          </tr>
          <xsl:apply-templates select="cml:moleculeList"/>
        </table>
      </div>
      <h3 id="reactions-title" class="handcursor">Reactions</h3>
      <table id="reactions" class="switchgroup4">
      <xsl:apply-templates select="cml:reactionList"/>
    </table>

    <!--Show the "results"-->
    <xsl:if test="//me:densityOfStatesList">
      <h3 id="densityOfStates-title" class="handcursor">Partition Functions</h3>
      <div id="densityOfStates" class="switchgroup1">
        <!--<xsl:apply-templates select="//*[@calculated]"/>-->
      <xsl:variable name="txt1" select="substring-after(//me:densityOfStatesList/me:description,'.')"/>
      <xsl:variable name="txt2" select="substring-after($txt1,'.')"/>
      <xsl:value-of select="substring-before($txt1,'.')"/>
      <br/>
      <xsl:value-of select="substring-before($txt2,'.')"/>
      <br/>
      <xsl:value-of select="substring-after($txt2,'.')"/>    
      <p></p>
        <xsl:apply-templates select="//me:densityOfStatesList"/>
      </div>
    </xsl:if>

    <xsl:if test="//me:microRateList">
      <h3 id="microRates-title" class="handcursor">Microcanonical Rate Coefficients</h3>
      <div id="microRates" class="switchgroup2">
        <xsl:apply-templates select="//me:microRateList"/>
      </div>
    </xsl:if>

    <xsl:if test="//me:rateList">
      <h3 id="BWrates-title" class="handcursor">Bartis-Widom Phenomenological Rate Coefficients</h3>
      <div id="BWrates" class="switchgroup5">
        <hh5 id="Punchout-title" class="handcursor">Copy and paste for spreadsheets, etc.</hh5>
        <div id="Punchout" class="switchgroup8">
          <xsl:call-template name="punchheader"/>
          <xsl:call-template name="punchoutput"/>
        </div>
        <xsl:apply-templates select="//me:rateList"/>
      </div>
    </xsl:if>

    <xsl:if test="//me:populationList">
      <h3 id="Populations-title" class="handcursor">Species / time profiles</h3>
      <div id="Populations" class="switchgroup6">
        <xsl:apply-templates select="//me:populationList"/>
      </div>
    </xsl:if>

    <xsl:if test="//@fitted">
      <h3 id="FittedParams-title" class="handcursor">Fitted Parameters</h3>
      <div id="FittedParams" class="switchgroup7">
        <table>
          <xsl:apply-templates select="//@fitted"/>
        </table>
      </div>
    </xsl:if>
      
      <!--Show toggle for inactive display if the file contains any-->
    <xsl:if test="//*[@active='false']">
      <p id="hide" onclick="toggle()">Hide/show inactive</p>
    </xsl:if>
      
    <!--Script for expanding an contracting sections-->
    <script type="text/javascript">
      <![CDATA[
        for(var i=1; i <=8; i++)
        {
          var mc=new switchcontent("switchgroup" + i)
          mc.setStatus('- ','+ ')
          mc.setPersist(true)
          mc.init()
        }
      ]]>
    </script>

    <xsl:call-template name="drawDiag"/>
    <xsl:apply-templates select="//me:analysis" mode="diagram"/>

    </body>
 </html>
</xsl:template>

  
  <xsl:template match="cml:molecule">
    <tr>
      <xsl:if test="@active='false'">
        <xsl:attribute name="class">inactive</xsl:attribute>
      </xsl:if>
      <td class="name">
        <xsl:value-of select="@id"/>
      </td>
      <td align="center">
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
      <xsl:if test="@active='false'">
        <xsl:attribute name="class">inactive</xsl:attribute>
      </xsl:if>
      <td>
        <xsl:value-of select="@id"/>
      </td>
      <td class="name">
        <xsl:for-each select=".//cml:reactant/cml:molecule/@ref">
          <xsl:if test="position()!=1"> + </xsl:if>
          <xsl:value-of select="."/>
        </xsl:for-each>
      </td>
      <td>
        <xsl:if test="@reversible='true'">&lt;</xsl:if>=>
      </td>
      <td class="name">
        <xsl:for-each select=".//cml:product/cml:molecule/@ref">
          <xsl:if test="position()!=1"> + </xsl:if>
          <xsl:value-of select="."/>
        </xsl:for-each>
      </td>
      <td>
        <xsl:if test=".//me:transitionState/cml:molecule/@ref">
          (Transition State <xsl:value-of select=".//me:transitionState/cml:molecule/@ref"/>)
        </xsl:if>
      </td>
      <td>
        <xsl:value-of select="me:MCRCMethod"/>
      </td>
      <td>
        <xsl:if test="me:preExponential">
          <xsl:if test="me:activationEnergy/@reverse">
            <xsl:value-of select="'(reverse) '"/>
          </xsl:if>
          <xsl:value-of select="concat('A = ', me:preExponential, 
              ' E = ', me:activationEnergy, me:activationEnergy/@units)"/>
        </xsl:if>
        <xsl:if test="me:preExponential/@stepsize | me:activationEnergy/@stepsize">
          <xsl:value-of select="' with range of parameters'"/>
        </xsl:if>
                  
       <xsl:if test="cml:rateParameters">
          <xsl:value-of select="concat( 'A = ', cml:rateParameters/cml:A,
              ' E = ', cml:rateParameters/cml:E, cml:rateParameters/cml:E/@units)"/>
        </xsl:if>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="me:densityOfStatesList">
    <table>
      <tr>
        <td class="tablehead1" colspan="5" align="center">
          Partition functions for <xsl:value-of select="../@id"/>
        </td>
      </tr>
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
      Canonical rate coefficients for
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

  <xsl:template match="me:rateList">
    <xsl:if test="not(//@fitted)"><!--don't show params when fitting-->
      <p class="paramheader">
        <xsl:for-each select="../me:parameters/@*">
          <xsl:value-of select="concat(name(),'=',.,', ')"/>
        </xsl:for-each>
      </p>
    </xsl:if>
    <table>
     <tr>
       <td class="tablehead1" colspan="5" align="center">
         At <xsl:value-of select="concat(@T,' K, ', @conc, ' molecules cm')"/><sup>-3</sup>
       </td>
     </tr>
      <xsl:if test="me:firstOrderRate">
       <tr>
         <td class="tablehead2" colspan="5" align="center">Conversion Rate Coefficients</td>
       </tr>
       <xsl:for-each select="me:firstOrderRate">
          <tr>
            <td><xsl:value-of select="@fromRef"/></td>
            <td>&#8195;&#8594;&#8195;</td>
            <td><xsl:value-of select="@toRef"/></td>
            <td> <xsl:value-of select="."/></td>
            <td>s<sup>-1</sup></td>
          </tr>
        </xsl:for-each>
      </xsl:if>
      <tr>
        <td class="tablehead2" colspan="5" align="center">Loss Rate Coefficients</td>
      </tr>
      <xsl:for-each select="me:firstOrderLoss">
        <tr>
          <td><xsl:value-of select="@ref"/> </td>
          <td></td>
          <td></td>
          <td><xsl:value-of select="."/></td>
          <td>s<sup>-1</sup></td>
        </tr>
       </xsl:for-each>
     </table>
    </xsl:template>

  <xsl:template match="//me:populationList">
    <xsl:variable name="speciesNames" select="me:population[1]/me:pop/@ref"/>
    <p class="paramheader">
      <xsl:for-each select="../me:parameters/@*">
        <xsl:value-of select="concat(name(),'=',.,', ')"/>
      </xsl:for-each>
    </p>
    <table>
      <tr><td class="tablehead1" colspan="5" align="center">
        <xsl:value-of select="concat('Populations (mole fractions) at ',@T,'K ',
                      @conc, ' molecules cm')"/><sup>-3</sup>
      </td></tr>
      <tr><td class="tablehead3">Time, sec</td>
        <xsl:for-each select="$speciesNames">
          <td class="tablehead3">
            <xsl:value-of select="."/>
          </td>
        </xsl:for-each>
      </tr>
      <xsl:for-each select="me:population">
        <tr>
          <td><xsl:value-of select="@time"/></td>
          <xsl:for-each select="me:pop">
            <td>
              <xsl:value-of select="."/>
            </td>
          </xsl:for-each>
        </tr>
      </xsl:for-each>
    </table>
  </xsl:template>

  <xsl:template match="//@fitted">
    <tr>
      <td> <xsl:value-of select="ancestor::*[@id]/@id"/> </td>
      <td>
        <xsl:choose>
          <xsl:when test="local-name(../..)='property'">
            <xsl:value-of select="../../@dictRef"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="name(..)"/>
          </xsl:otherwise>
      </xsl:choose>
      </td>
      <td> <xsl:value-of select=".."/> </td>
      <td> <xsl:value-of select="../@units"/> </td>
      <td> <xsl:value-of select="concat('(fitted ', ., ')')"/></td>
    </tr>
  </xsl:template>


    <xsl:template match="cml:metadataList">
    <div id="metadata">
        <xsl:value-of select="dc:creator"/>:
        <xsl:value-of select="dc:date"/>,
        <xsl:value-of select="dc:contributor"/>
    </div>
  </xsl:template>

  <xsl:template name="ZPEunits">
    <xsl:value-of select=
      "//cml:property[@dictRef='me:ZPE']/cml:scalar/@units | //me:activationEnergy/@units"/>
  </xsl:template>

  <xsl:template name="punchheader">
    <xsl:value-of select="concat(/me:mesmer/cml:title,' calculated ',//me:analysis[1]/@calculated)"/>
    <br/>
    <xsl:for-each select="//me:analysis[1]/me:parameters/@*">
      <xsl:value-of select="concat(name(),',')"/>
    </xsl:for-each>
    <xsl:value-of select="concat('T',',','conc',',')"/>
    <xsl:for-each select="//me:analysis[1]/me:rateList[1]/me:firstOrderRate">
      <xsl:value-of select="concat(@fromRef,'->',@toRef,',')"/>
    </xsl:for-each>
  </xsl:template>
  
  <xsl:template name="punchoutput">
    <xsl:for-each select ="//me:rateList">
      <br/>
      <xsl:for-each select="../me:parameters/@*">
        <xsl:value-of select="concat(.,',')"/>
      </xsl:for-each>
      <xsl:value-of select="concat(@T,',',@conc,',')"/>
      <xsl:for-each select="me:firstOrderRate">
        <xsl:value-of select="concat(.,',')"/>
      </xsl:for-each>
    </xsl:for-each>
  </xsl:template>

</xsl:stylesheet> 
