<?xml version="1.0" encoding="utf-8"?>

<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:str="http://exslt.org/strings"
  xmlns:set="http://exslt.org/sets">

  <xsl:include href="chemStrucEL.xsl"/>
  <xsl:include href="mesmerDiag.xsl"/>
  <xsl:include href="switchcontent.xsl"/>
  <xsl:include href="popDiag.xsl"/>
  <xsl:include href="rawdataplot.xsl"/>
  
  <xsl:key name="molrefs" match="cml:molecule" use="@id"/>
  
<xsl:variable name="title">
  <xsl:choose>
    <xsl:when test="/me:mesmer/cml:metadataList/dc:title">
      <xsl:value-of select="/me:mesmer/cml:metadataList/dc:title"/>
    </xsl:when>
    <xsl:when test="//me:title|//cml:title">
      <xsl:value-of select="//me:title|//cml:title"/>
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
        body{margin:20px;padding:4;}
        caption{font-weight:bold;}
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
		.warn{color:red;}
        .tablehead2{text-decoration:underline;padding-top:10px;}
        .tablehead3{font-weight:bold;padding-top:10px;text-align:center}
        .tablehead4{
          padding-top:10px;text-align:center;
          
          border-bottom-width:2;}
         .tablehead5{font-weight:bold;}
        .paramheader{
          font-family: Arial, Helvetica, sans-serif;
          color:black;
          font-weight:bold;
          text-decoration: underline;
        }
        .poplabel{color:black; font-weight:bold;} 

        thead{border-style:solid;border-width:2;}
        table{background-color:#e0f8f8; margin-bottom:12px;}
        table.thermo{ text-align:center;width:450px;}
        .highlight{background-color:#f8f8f8;}

        td{padding:0px 2px;}
        td.SAid{font-size:smaller;font-weight:bold;}
        .SA2{text-align:center;}
        th{font-size:smaller;border-bottom:1px dashed black;}
        h3{color:teal;font-family: Arial, Helvetica, sans-serif;font-weight:bold;}
        hh5{color:black;font-family: Arial, Helvetica, sans-serif;font-weight:bold;font-size:smaller;}
        .normal{color:black; font-size:smaller;}
        .normal2{color:black; background-color:#e0f8f8}
        .handcursor{cursor:hand; cursor:pointer;}
        .inactive{color:silver;stroke:silver;}
        .error{font-weight:bold;font-size:large;background-color:red;padding:20px;}
        #header{color:black;font-family: Arial, Helvetica, sans-serif;font-weight:bold;}
        #title{font-size:larger;font-weight:bold;}
        #description{color:black;font-size:60%;}
        #metadata{color:teal;font-size:60%;}
        #hide{font-size:small;text-decoration:underline;color:blue;cursor:pointer;}
        #Punchout{font-size:12px;color:gray;}
        .thermo tr td:nth-of-type(odd){ border-right:1px dashed;padding-right:10px;}
        .thermo tr th:nth-of-type(odd){ border-right:1px dashed;padding-right:10px;}
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
          <xsl:value-of select="/me:mesmer/cml:description|me:description"/>
        </p>
        <xsl:if test="cml:metadataList">
          <xsl:apply-templates select="cml:metadataList"/>
        </xsl:if>
      </div>
      
      <xsl:variable name ="Econvention" select="//cml:moleculeList/@convention"/>
      <xsl:variable name ="Etext">
        <xsl:choose>
          <xsl:when test="$Econvention='thermodynamic'">&#916;Hf at 0K</xsl:when>
          <xsl:otherwise>Energy</xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      
      <xsl:variable name ="Eunits">
        <xsl:call-template name="ZPEunits"/>
      </xsl:variable>
      <h3 id="mols-title" class="handcursor">Molecules</h3>
      <div id="mols" class="switchgroup3">
        <xsl:if test="$Econvention">
          <xsl:value-of select="concat('Energy convention is ', $Econvention)"/>
        </xsl:if>
        <table class="mol">
          <tr class="tableheader">
            <td>Name</td>
            <td><xsl:value-of select="$Etext"/><br /><xsl:value-of select="$Eunits"/></td>
            <td>Rotational constants<br />cm<sup>-1</sup></td>
            <td>Vibrational frequencies<br />cm<sup>-1</sup></td>
          </tr>
          <xsl:apply-templates select="cml:moleculeList"/>
        </table>
      </div>

      <!--Get molecule structure names in a variable so section heading
      can be hidden if there are none. Drawing structure twice is very inefficient!-->
      <xsl:variable name="structText">
        <xsl:apply-templates select="//cml:molecule" mode="chemStructure"/>
      </xsl:variable>
      <xsl:if test="string-length($structText)>0">
        <h3 id="structs-title" class="handcursor">Chemical Structures</h3>
        <div id="structs" class="switchgroup10">
         <xsl:apply-templates select="//cml:molecule" mode="chemStructure"/>
        </div>
      </xsl:if>

      <xsl:if test="//cml:reactionList">
        <h3 id="reactions-title" class="handcursor">Reactions</h3>
        <table id="reactions" class="switchgroup4">
          <xsl:apply-templates select="cml:reactionList"/>
        </table>
      </xsl:if>

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

    <xsl:if test="//me:microRateList | //me:canonicalRateList">
      <h3 id="microRates-title" class="handcursor">Canonical Rate Coefficients</h3>
      <div id="microRates" class="switchgroup2">
        <xsl:apply-templates select="//me:microRateList | //me:canonicalRateList"/>
      </div>
    </xsl:if>

    <xsl:if test="//me:rateList">
      <h3 id="BWrates-title" class="handcursor">Bartis-Widom Phenomenological Rate Coefficients</h3>
      <div id="BWrates" class="switchgroup5">
        <div class="warn">
          <xsl:value-of select="//me:rateList/me:warning"/>
        </div>
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

    <xsl:if test="//me:grainPopulationList">
      <h3 id="GrainPopulations-title" class="handcursor">Grain populations at selected times</h3>
      <div id="GrainPopulations" class="switchgroup11">
        <xsl:apply-templates select="//me:grainPopulationList"/>
      </div>
    </xsl:if>
      
      <xsl:if test="//me:eigenvalueList">
        <h3 id="Eigenvalues-title" class="handcursor">Eigenvalues</h3>
        <div id="Eigenvalues" class="switchgroup9">
          <xsl:apply-templates select="//me:eigenvalueList"/>
        </div>
      </xsl:if>

      <xsl:if test="//@fitted">
      <h3 id="FittedParams-title" class="handcursor">Fitted Parameters</h3>
      <div id="FittedParams" class="switchgroup7">
        <table>
          <xsl:apply-templates select="//@fitted"/>
        </table>
        <xsl:if test="//me:errorPropagationTable">
          <table>
            <tr><td class="tablehead1" colspan="4">Rate error estimates (2 standard deviations)</td></tr>
            <xsl:apply-templates select="//me:errorPropagationTable"/>
          </table>
        </xsl:if>
      </div>
    </xsl:if>

      <xsl:if test="//me:representation">
        <h3 id="ChebyshevRepresentation-title" class="handcursor">Chebyshev representation of rates</h3>
        <div id="ChebyshevRepresentation" class="switchgroup14">
          <xsl:apply-templates select="//me:representation"/>
        </div>
      </xsl:if>
      
      <xsl:if test="//me:thermoTable">
        <h3 id="ThermodynamicData-title" class="handcursor">Thermodynamics for molecules</h3>
        <div id="ThermodynamicData" class="switchgroup12">
          <xsl:apply-templates select="//me:thermoTable"/>
            <pre> <!--Fixed format NASA polynomials-->
              <xsl:for-each select="//cml:scalar[@dictRef='NasaPolynomial']">
                <xsl:value-of select="."/>
              </xsl:for-each>
            </pre>
          <!--Cantera form of NASA polynomials
          species(name = "O2", atoms = " O:2 ",
          thermo = (
            NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
            NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )
      Indentation is significant in the following block.-->
      <xsl:for-each select="//cml:property[@dictRef='NasaPolynomial']">
        <pre>
          <div>
            <xsl:value-of select="concat('name = &quot;',ancestor::cml:molecule/@id,'&quot;, ',
                          'atoms = &quot; ')"/>
            <xsl:variable name="atoms" select="ancestor::cml:molecule//cml:atom"/>
            <xsl:variable name="elements" select="set:distinct($atoms/@elementType)"/>
            <xsl:for-each select="$elements">
              <xsl:value-of select="concat(., ':', count($atoms[@elementType=current()]),' ')"/>
            </xsl:for-each>
            <xsl:text>&quot;, </xsl:text>
          </div>
          <xsl:variable name="coeffs" select="str:tokenize(cml:array[@dictRef='NasaCoeffs'])"/>
          <xsl:variable name="lowercoeffs" select="$coeffs[position()&lt;=7]"/>
          <xsl:variable name="uppercoeffs" select="$coeffs[position()>7 and position()!=15]"/>
          <xsl:choose>
            <xsl:when test="cml:scalar[@dictRef='NasaMidT']!=0"> 
              <xsl:value-of select="concat('thermo = (',
              'NASA([',format-number(cml:scalar[@dictRef='NasaLowT'], '#.##'),',',
                        format-number(cml:scalar[@dictRef='NasaMidT'], '#.##'),'],[')"/>
              <xsl:for-each select="$uppercoeffs">
                <xsl:value-of select="."/>
                <xsl:if test="position() != last()">,</xsl:if>
              </xsl:for-each>]),
          <xsl:value-of select="concat(
          'NASA([',format-number(cml:scalar[@dictRef='NasaMidT'], '#.##'),',',
                    format-number(cml:scalar[@dictRef='NasaHighT'], '#.##'),'],[')"/>
          <xsl:for-each select="$lowercoeffs">
            <xsl:value-of select="."/>
            <xsl:if test="position() != last()">,</xsl:if>
          </xsl:for-each>]))
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="concat('thermo = (',
              'NASA([',format-number(cml:scalar[@dictRef='NasaLowT'], '#.##'),',',
                        format-number(cml:scalar[@dictRef='NasaHighT'], '#.##'),'],[')"/>
              <xsl:for-each select="$lowercoeffs">
                <xsl:value-of select="."/>
                <xsl:if test="position() != last()">,</xsl:if>
              </xsl:for-each>]))
            </xsl:otherwise>
          </xsl:choose>
        </pre>
      </xsl:for-each>
          
    </div> 
      </xsl:if>
      
    <xsl:if test="//me:experimentalRate | //me:experimentalYield | //experimentalEigenvalue">
      <h3 id="ExperimentalData-title" class="handcursor">Comparison with Experimental Data</h3>
      <div id="ExperimentalData" class="switchgroup13">
        <xsl:apply-templates select="//me:PTs"/>
      </div>
    </xsl:if>

    <xsl:if test="//me:sensitivityAnalysisTables">
      <h3 id="sensitivityAnalysis-title" class="handcursor">Sensitivity Analysis</h3>
      <div id="sensitivityAnalysis" class="switchgroup15">
        <xsl:apply-templates select="//me:sensitivityAnalysisTables/me:sensitivityAnalysisTable"/>
      </div>
    </xsl:if>

      <!--Show toggle for inactive display if the file contains any-->
    <xsl:if test="//*[@active='false']">
      <p id="hide" onclick="toggle()">Hide/show inactive</p>
    </xsl:if>
      
    <!--Script for expanding an contracting sections-->
    <script type="text/javascript">
      <![CDATA[
        for(var i=1; i <=15; i++)
        {
          var mc=new switchcontent("switchgroup" + i)
          mc.setStatus('- ','+ ')
          mc.setPersist(true)
          mc.init()
        }
      ]]>
    </script>

    <xsl:call-template name="drawDiag"/>
    <xsl:apply-templates select="//me:rawData"/>
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

    <xsl:template match="me:microRateList | me:canonicalRateList">
    <h4>
      Canonical (high pressure) rate coefficients for
      <xsl:value-of select="../@id"/>
      <span class="normal">
        (calculated from microcanonical rates 
        <xsl:value-of select="@calculated"/>
        )
      </span>
    </h4>
    <table>
      <tr class="tableheader">
        <td>T/K</td>
        <td>
          <xsl:value-of select="concat('kf/', me:kinf/me:val/@units)"/>
        </td>
        <xsl:if test="me:kinf/me:rev">
          <td>
            <xsl:value-of select="concat('krev/', me:kinf/me:rev/@units)"/>
          </td>
          <td>
            <xsl:value-of select="concat('Keq/', me:kinf/me:Keq/@units)"/>
          </td>
        </xsl:if>
      </tr>
      <xsl:for-each select="me:kinf">
        <tr>
          <td>
            <xsl:value-of select="me:T"/>
          </td>
          <td>
            <xsl:value-of select="me:val"/>
          </td>
          <td>
            <xsl:value-of select="me:rev"/>
          </td>
          <td>
            <xsl:value-of select="me:Keq"/>
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
         <xsl:value-of select="concat(' in ', @bathGas)"/>
       </td>
     </tr>
      <xsl:if test="me:firstOrderRate">
       <tr>
         <td class="tablehead2" colspan="5" align="center">First Order Conversion Rate Coefficients</td>
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
      <xsl:if test="me:secondOrderRate">
        <tr>
          <td class="tablehead2" colspan="5" align="center">Second Order Conversion Rate Coefficients</td>
        </tr>
        <xsl:for-each select="me:secondOrderRate">
          <tr>
            <td>
              <xsl:value-of select="@fromRef"/>
            </td>
            <td>&#8195;&#8594;&#8195;</td>
            <td>
              <xsl:value-of select="@toRef"/>
            </td>
            <td>
              <xsl:value-of select="."/>
            </td>
            <td>
              cm<sup>3</sup>molecule<sup>-1</sup>s<sup>-1</sup>
            </td>
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

<!--===========================================================-->
  <xsl:template match="//me:grainPopulation">
    <xsl:if test="position()=2"><!--why 2 not 1?-->
      <p>See graphs below. This section is mainly for copying and pasting to a spreadsheet.</p>
    </xsl:if>
    <p class="paramheader">
      <xsl:for-each select="../../me:parameters/@*">
        <xsl:value-of select="concat(name(),'=',.,', ')"/>
      </xsl:for-each>
    </p>
    <table>
      <tr>
        <td class="tablehead1" colspan="3" align="center">
          <xsl:value-of select="concat(@ref,' at ',@time,' ',../@T,'K ', ../@conc,' molecules cm')"/>
          <sup>-3</sup>
        </td>
      </tr>
      <tr><td class="tablehead3">Grain Energy,cm-1</td>
          <td class="tablehead3">Normalised Pop</td>
          <td class="tablehead3">Log Pop</td>
      </tr>
      <xsl:apply-templates select="me:grain"/>
    </table>
  </xsl:template>
  
  <!--===========================================================-->
  <xsl:template match="me:grain">
    <!--Writes a row in the grain population table--> 
    <tr>
      <td>
        <xsl:value-of select="@energy"/>
      </td>
      <td>
        <xsl:value-of select="@normpop"/>
      </td>
      <td>
        <xsl:value-of select="@logpop"/>
      </td>
    </tr>
  </xsl:template>

  <!--===========================================================-->
  <xsl:template match="//me:eigenvalueList">
    <p class="normal2">
      <span class="tableheader">
        <xsl:choose>
          <xsl:when test="@selection&lt;=0">
            <xsl:text>All of </xsl:text>
          </xsl:when>
          <xsl:otherwise>
        <xsl:value-of select="@selection"/>
            <xsl:text> out of </xsl:text>
          </xsl:otherwise>
        </xsl:choose><xsl:value-of select="@number"/>        
         eigenvalues (sec<sup>-1</sup>) are shown
      </span>
      <br/>
      <xsl:for-each select="me:eigenvalue">
        <xsl:value-of select="concat(.,', ')"/>
      </xsl:for-each>
    </p>
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
    <td> <xsl:value-of select="concat('ChiSquared=',../@chiSquared)"/> </td>
      <td> <xsl:value-of select="concat('(fitted ', ., ')')"/></td>
    </tr>
  </xsl:template>

  <xsl:template match ="//me:errorPropagationTable">
    <xsl:variable name="reactants" 
                  select="substring-before(me:propagatedErrors/@reaction, '=>')"/>
    <!--Count number of '+' in the reactants to determine the reaction order-->
    <xsl:variable name="coeffunits">
      <xsl:choose>
        <xsl:when 
          test="string-length($reactants) - string-length(translate($reactants,'+',''))=0">
          s⁻¹
        </xsl:when>
        <xsl:otherwise>cm³ molecule⁻¹ s⁻¹</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>      
      <tr class="tablehead2">
        <td><xsl:value-of select="concat('At ', @temperature,'K , ', @concentration, ' molecule cm⁻³')"/></td>
        <td colspan="2" align="center">
          <xsl:value-of select="concat('Rate  coeff, ', $coeffunits)"/>
        </td>
      </tr>
      <tr>
        <td> <xsl:value-of select="me:propagatedErrors/@reaction"/></td>
        <td>
          <xsl:value-of select="me:propagatedErrors/me:rateCoefficient"></xsl:value-of>
        </td>
        <td>
          <xsl:value-of select="concat('&#177;', number(me:propagatedErrors/me:standardDeviation)*2)"></xsl:value-of>
        </td>
      </tr>

  </xsl:template>
  
  <xsl:template match="//me:PTs">
  <div class="tablehead1">
    <!--"R1 => R2" or "R1 loss"(when ref1=ref2) or "R yield"-->
    <xsl:choose>
      <xsl:when test="me:PTpair/*/@ref1[1]=me:PTpair/*/@ref2[1]">
        <xsl:value-of select="concat(me:PTpair/*/@ref1[1], ' loss')" />
      </xsl:when>
      <xsl:when test="me:PTpair/me:experimentalYield/@ref[1]">
        <xsl:value-of select="concat(me:PTpair/me:experimentalYield/@ref[1], ' yield')" />
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="concat(me:PTpair/*/@ref1[1], ' => ', me:PTpair/*/@ref2[1])"/>
      </xsl:otherwise>
    </xsl:choose>	  
  </div>
    <xsl:value-of select="concat('  calculated ', me:PTpair/*/@calculated[1])"/>
    <table>
    <tr class="tablehead1">
      <td>Temperature(K)</td>
      <td>
        <xsl:value-of select="concat('Concentration(', me:PTpair/@units|@me:units, ')')"/>
      </td>
      <td>bathGas</td>
      <xsl:for-each select="me:PTpair[1]/me:experimentalRate|me:PTpair[1]/me:experimentalYield|me:PTpair[1]/me:experimentalEigenvalue">
        <!--just the first node is selected-->
        <td> <xsl:value-of select="local-name()"/> </td>
        <td> <xsl:value-of select="concat('calculated', substring-after(local-name(),'experimental'))"/> </td>
      </xsl:for-each>

    </tr>
    <xsl:for-each select="me:PTpair/me:experimentalRate|me:PTpair/me:experimentalYield|me:PTpair/me:experimentalEigenvalue">
      <tr style="text-align:center">
        <td> <xsl:value-of select="../@T|@me:T"/></td>
        <td> <xsl:value-of select="../@P|@me:P"/></td>
        <td> <xsl:value-of select="../@bathGas"/></td>
        <td> <xsl:value-of select="."/></td>
        <td> <xsl:value-of select="@calcVal"/></td>
      </tr>
    </xsl:for-each>
      </table>
  </xsl:template>

  <xsl:template match="//me:representation">
    <pre><xsl:value-of select="."/></pre>
  </xsl:template>

  <xsl:template match="//me:thermoTable">
    <h4><xsl:value-of select="../@id"/></h4>
    <xsl:choose>
      <xsl:when test="me:thermoValue/@cellH">
      <p>Data is derived in two ways: from analytical expressions and from cell averages</p>
      <table class="thermo">
        <tr class="tablehead1">
          <td>T</td> <td colspan="2">H(T)-H(0)</td>
          <td colspan="2">S(T)</td>
          <td colspan="2">G(T)</td>
          <td colspan="2">Cp(T)</td>
          <xsl:if test="@unitsHf">
            <td>&#916;Hf</td>
          </xsl:if>
        </tr>
        <tr>
          <td></td>
          <td>analytic</td><td>cell</td>
          <td>analytic</td><td>cell</td>
          <td>analytic</td><td>cell</td>
          <td>analytic</td><td>cell</td>
        </tr>
          <tr>
          <th><xsl:value-of select="@unitsT"/></th>
          <th><xsl:value-of select="@unitsH"/></th>
          <th><xsl:value-of select="@unitsH"/></th>
          <th><xsl:value-of select="@unitsS"/></th>
          <th><xsl:value-of select="@unitsS"/></th>
          <th><xsl:value-of select="@unitsG"/></th>
          <th><xsl:value-of select="@unitsG"/></th>
          <th><xsl:value-of select="@unitsCp"/></th>
          <th><xsl:value-of select="@unitsCp"/></th>
          <th><xsl:value-of select="@unitsHf"/></th>
          </tr>
        <xsl:for-each select="me:thermoValue">
          <tr>
            <xsl:if test="@T=298.15">
              <xsl:attribute name="class">highlight</xsl:attribute>
            </xsl:if> 
            <td> <xsl:value-of select="@T"/></td>
            <td> <xsl:value-of select="@H"/></td>
            <td> <xsl:value-of select="@cellH"/></td>
            <td> <xsl:value-of select="@S"/></td>
            <td> <xsl:value-of select="@cellS"/></td>
            <td> <xsl:value-of select="@G"/></td>
            <td> <xsl:value-of select="@cellG"/></td>
            <td> <xsl:value-of select="@Cp"/></td>
            <td> <xsl:value-of select="@cellCp"/></td>
            <td> <xsl:value-of select="@Hf"/></td>
          </tr>
        </xsl:for-each>
      </table>
    </xsl:when>
      <xsl:otherwise>
        <!--Values only from analytical version-->
        <table>
          <tr class="tablehead3">
            <td>T</td>
            <td> &#160;H(T)-H(0)  &#160; </td>
            <td> &#160;S(T)  &#160; </td>
            <td> &#160;G(T)  &#160; </td>
            <td> &#160;Cp(T)  &#160; </td>
            <xsl:if test="@unitsHf">
              <td>&#916;Hf</td>
            </xsl:if>
          </tr>
          <tr>
            <th><xsl:value-of select="@unitsT"/></th>
            <th><xsl:value-of select="@unitsH"/></th>
            <th><xsl:value-of select="@unitsS"/></th>
            <th><xsl:value-of select="@unitsG"/></th>
            <th><xsl:value-of select="@unitsCp"/></th>
            <th><xsl:value-of select="@unitsHf"/></th>
          </tr>
          <xsl:for-each select="me:thermoValue">
            <tr>
              <xsl:if test="@T=298.15">
                <xsl:attribute name="class">highlight</xsl:attribute>
              </xsl:if>
              <td><xsl:value-of select="@T"/></td>
              <td><xsl:value-of select="@H"/></td>
              <td><xsl:value-of select="@S"/></td>
              <td><xsl:value-of select="@G"/></td>
              <td><xsl:value-of select="@Cp"/></td>
              <td><xsl:value-of select="@Hf"/></td>
            </tr>
          </xsl:for-each>
        </table>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template match="//me:sensitivityAnalysisTables/me:sensitivityAnalysisTable">
    <p class="tablehead1">
      At <xsl:value-of select="concat(@temperature,' K, ', @concentration, ' molecules cm')"/><sup>-3</sup>
    </p>
    <xsl:for-each select="me:sensitivityIndices">
      <xsl:variable name="params" select="me:firstOrderIndex/@key"/>
      <div class="tablehead1">
        <xsl:value-of select="concat('For ',@reaction)"/>
      </div>
      <div>
        <xsl:value-of select="concat('Standard deviation: ', me:standardDeviation)"/>
        <br/>
        <xsl:value-of select="concat('R-Squared statistic: ', me:R_Squared)"/>
      </div>
      <table>
        <caption class="tablehead5">First order indices</caption>
        <tr>
          <th>ID</th><th>Input variable</th><th>First order index</th>
        </tr>
        <xsl:for-each select="me:firstOrderIndex">
          <tr>
            <td class="SAid">
              <xsl:value-of select="concat('(',position(),')')"/>
            </td>
            <td>
              <xsl:value-of select="@key"/>
            </td>
            <td>
              <xsl:value-of select="."/>
            </td>
          </tr>
        </xsl:for-each>
      </table>
      <xsl:variable name="secondorders" select="me:secondOrderIndex"/>
      <table>
        <caption class="tablehead5">Second order indices</caption>
        <tr>
          <th> IDs </th>
          <xsl:for-each select="$params">
            <th>
              <xsl:value-of select="concat('(',position(),')')"/>
            </th>
          </xsl:for-each>
        </tr>

        <xsl:for-each select="$params">
          <xsl:variable name="curparam" select="."/>
          <tr>
            <td class="SAid">
              <xsl:value-of select="concat('(',position(),')')"/>
            </td>
            <xsl:for-each select="$params">
              <xsl:variable name="val2" select="$secondorders[(@key1=$curparam)and(@key2=current())]"/>
              <td class="SA2">
                <xsl:choose>
                  <xsl:when test ="$val2">
                    <xsl:value-of select="$val2"/>
                  </xsl:when>
                  <xsl:otherwise>----</xsl:otherwise>
                </xsl:choose>
              </td>
            </xsl:for-each>
          </tr>
        </xsl:for-each>
      </table>
      <br/>
    </xsl:for-each>
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
    <xsl:for-each select="//me:analysis[1]/me:rateList[1]/me:firstOrderLoss">
      <xsl:value-of select="concat('Loss of ',@ref,',')"/>
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
      <xsl:for-each select="me:secondOrderRate">
        <xsl:value-of select="concat(.,',')"/>
      </xsl:for-each>
      <xsl:for-each select="me:firstOrderLoss">
        <xsl:value-of select="concat(.,',')"/>
      </xsl:for-each>
    </xsl:for-each>
  </xsl:template>

</xsl:stylesheet> 
