<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:cml="http://www.xml-cml.org/schema"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:me="http://www.chem.leeds.ac.uk/mesmer"
                xmlns:exsl="http://exslt.org/common"
                xmlns:set="http://exslt.org/sets"
                xmlns:svg="http://www.w3.org/2000/svg"
                xmlns:math="http://exslt.org/math"
                extension-element-prefixes="exsl set">
  
  <!--The tranformation from CML to SVG by XSLT needs two stages to allow
  z sorting to get proper hiding of objects by others.
  In the first stage the CML input is read to make a variable aatoms which contains
  centered atom coordinates and a variable bbonds with the bond information.
  The coordinates are modified from the original by a rotation matrix.
  These variables contain several <object>s :
    <object x="1.1" y="2.2" eltype="C"> 3.3 </object> for atoms
    <object x1="1.1" y1="2.2" x2="4.4" y2="5.5"> 6.6 </object> for bonds
  For bonds the value of z  is derived the average of its atoms.
  All the <object>s are sorted by their z coordinates and output in a second
  stage template, which writes a svg:circle for atoms and a svg:line for bonds
  to the output.
  -->

  <!--Global variables-->
  <xsl:variable name="atomRadius">0.45</xsl:variable>
  <xsl:variable name="bondWidth">0.01</xsl:variable>

  <!--Element colors, base and highlight-->
  <xsl:variable name="atomcolors">
    <el e="H" c="#808080" h="#C0C0C0"/>
    <el e="C" c="#000000" h="#808080"/>
    <el e="O" c="#0000E0" h="#6060FF"/>
    <el e="N" c="#D00000" h="#FF6060"/>
    <el e="S" c="#D0B000" h="#EEDEAA"/>
  </xsl:variable>
  
  <xsl:variable name="sincos">
    <sin>0.000</sin><cos>1.000</cos><sin>0.259</sin><cos>0.966</cos>
    <sin>0.500</sin><cos>0.866</cos><sin>0.707</sin><cos>0.707</cos>
    <sin>0.866</sin><cos>0.500</cos><sin>0.966</sin><cos>0.259</cos>
    <sin>1.000</sin><cos>0.000</cos><sin>0.966</sin><cos>-0.259</cos>
    <sin>0.866</sin><cos>-0.500</cos><sin>0.707</sin><cos>-0.707</cos>
    <sin>0.500</sin><cos>-0.866</cos><sin>0.259</sin><cos>-0.966</cos>
    <sin>0.000</sin><cos>-1.000</cos><sin>-0.259</sin><cos>-0.966</cos>
    <sin>-0.500</sin><cos>-0.866</cos><sin>-0.707</sin><cos>-0.707</cos>
    <sin>-0.866</sin><cos>-0.500</cos><sin>-0.966</sin><cos>-0.259</cos>
    <sin>-1.000</sin><cos>0.000</cos><sin>-0.966</sin><cos>0.259</cos>
    <sin>-0.866</sin><cos>0.500</cos><sin>-0.707</sin><cos>0.707</cos>
    <sin>-0.500</sin><cos>0.866</cos><sin>-0.259</sin><cos>0.966</cos>
  </xsl:variable>
  
  <!--Rotation about y-axis-->
  <xsl:variable name="nsc" select="2"/>
  <xsl:variable name="cos" select="exsl:node-set($sincos)/cos[$nsc]"/>
  <xsl:variable name="sin" select="exsl:node-set($sincos)/sin[$nsc]"/>

  <xsl:template name="chemStructure" mode="chemStructure" match="cml:molecule">
    <!--Do not output molecules without atoms or missing any 3D coordinates-->
    <xsl:if test="(count(cml:atomArray/cml:atom) and
            (count(cml:atomArray/cml:atom) * 3) = count(cml:atomArray/cml:atom/@x3
            | cml:atomArray/cml:atom/@y3 | cml:atomArray/cml:atom/@z3))">
      <!--Calculate the dimension scaling factor, $range-->
      <xsl:variable name="xmax" select="math:max(//cml:atom/@x3)"/>
      <xsl:variable name="xmin" select="math:min(//cml:atom/@x3)"/>
      <xsl:variable name="ymax" select="math:max(//cml:atom/@y3)"/>
      <xsl:variable name="ymin" select="math:min(//cml:atom/@y3)"/>
      <xsl:variable name="zmax" select="math:max(//cml:atom/@z3)"/>
      <xsl:variable name="zmin" select="math:min(//cml:atom/@z3)"/>
      <xsl:variable name="xspan" select="$xmax - $xmin"/>
      <xsl:variable name="yspan" select="$ymax - $ymin"/>
      <xsl:variable name="zspan" select="$zmax - $zmin"/>
      <xsl:variable name="range">
        <xsl:choose>
          <xsl:when test="($xspan>=$yspan) and ($xspan>$zspan)">
            <xsl:value-of select="2*$atomRadius + $xspan"/>
          </xsl:when>
          <xsl:when test="($yspan>$zspan) and ($yspan>=$xspan)">
            <xsl:value-of select="2*$atomRadius + $yspan"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="2*$atomRadius + $zspan"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <!--Populate a variable with atom and bonds having rotated coordinates-->
      <xsl:variable name="aatoms">
        <xsl:apply-templates select="cml:atomArray/cml:atom">
          <xsl:with-param name="xc" select="0.5*($xmax + $xmin)"/>
          <xsl:with-param name="yc" select="0.5*($ymax + $ymin)"/>
          <xsl:with-param name="zc" select="0.5*($zmax + $ymin)"/>
        </xsl:apply-templates>

      </xsl:variable>
      <xsl:variable name="bbonds">
        <xsl:apply-templates select="cml:bondArray/cml:bond">
          <xsl:with-param name="ats" select="exsl:node-set($aatoms)"/>
        </xsl:apply-templates>
      </xsl:variable>

      <!--Do drawing in second stage-->
      <svg:svg viewBox="{concat(-0.5*$range,' ', -0.5*$range,' ',$range,' ',$range)}" width="200"
               stroke-width="{$bondWidth*$range}" >
        <!--Define coloring gradient for each element present in molecule-->
        <svg:defs>
          <xsl:for-each select="set:distinct(//cml:atom/@elementType)">
            <svg:radialGradient id="{.}" fx="30%" fy="30%">
              <svg:stop offset="0%" stop-color="{exsl:node-set($atomcolors)/el[@e=current()]/@h}"/>
              <svg:stop offset="100%" stop-color="{exsl:node-set($atomcolors)/el[@e=current()]/@c}"/>
            </svg:radialGradient>
          </xsl:for-each>
        </svg:defs>

        <!--Draw each atom and bond in z order-->
        <xsl:apply-templates select="exsl:node-set($bbonds)/object | exsl:node-set($aatoms)/object">
          <xsl:sort select="@z" data-type="number"/>
        </xsl:apply-templates>
        <svg:text x="-15%" y="47%" stroke-width=".005" font-size="0.5" >
          <xsl:value-of select="@id"/>
        </svg:text>
      </svg:svg>
    </xsl:if>
  </xsl:template>

  <!--First stage atom template-->
  <xsl:template name="atomstempl1" match="cml:atomArray/cml:atom">
    <xsl:param name="xc"/>
    <xsl:param name="yc"/>
    <xsl:param name="zc"/>
    
    <xsl:variable name="eltype" select="@elementType"/>
    <object x="{(@x3 - $xc)* $cos + (@z3 - $zc)* $sin}"
            y="{@y3 -$yc}"
            z="{(@z3 -$zc) * ($cos - $sin)}"
            eltype="{@elementType}" id="{@id}"/>
  </xsl:template>

  <!--First stage bond template-->
  <xsl:template name="bondstempl1" match="cml:bondArray/cml:bond">
    <xsl:param name="ats"/>
    <xsl:variable name="rref1" select="substring-before(@atomRefs2,' ')"/>
    <xsl:variable name="rref2" select="substring-after(@atomRefs2,' ')"/>
    <xsl:variable name="firstatom" select="$ats/object[@id=$rref1]"/>
    <xsl:variable name="secondatom" select="$ats/object[@id=$rref2]"/>
    <xsl:variable name="bondcolor">
      <xsl:choose>
        <xsl:when test="@color">
          <xsl:value-of select="@color"/>
        </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="'darkslategray'"/>
      </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <object x1="{$firstatom/@x}"  y1="{$firstatom/@y}"
            x2="{$secondatom/@x}" y2="{$secondatom/@y}"
            z ="{($firstatom/@z + $secondatom/@z) div 2}"
            order="{@order}" color="{$bondcolor}" id="{@id}"/>
  </xsl:template>

  <!--Second stage template to draw atoms-->
  <xsl:template match="*[@eltype]">
    <xsl:variable name="et" select="./@eltype"/>
    <svg:circle cx="{@x}" cy="{@y}" r="{$atomRadius}" z="{@z}" id="{@id}"
                fill="{concat('url(#', $et, ')')}" stroke="black"
                stroke-width="{$atomRadius * 0.03}"/>
  </xsl:template>

  <!--Second stage template to draw bonds-->
  <xsl:template match="*[@x1]">
    <svg:line x1="{@x1}" x2="{@x2}" y1="{@y1}" y2="{@y2}" z="{@z}" id="{@id}"
              stroke="{@color}"/>
  </xsl:template>
</xsl:stylesheet>
