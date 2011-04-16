<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0"  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:me="http://www.chem.leeds.ac.uk/mesmer">
  
  <xsl:output method="text"/>
  <xsl:template match="/">
    <xsl:value-of select="concat(me:mesmer/cml:title,' calculated ',//me:analysis[1]/@calculated,'&#xa;')"/>
    <xsl:call-template name="header"/>
    <xsl:apply-templates select="//me:analysis/me:rateList"/>
  </xsl:template>

  <xsl:template name="header">
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
    <xsl:text>&#xa;</xsl:text>
  </xsl:template>
  
  <xsl:template match="me:analysis/me:rateList">
    <xsl:for-each select="../me:parameters/@*">
      <xsl:value-of select="concat(.,',')"/>    
    </xsl:for-each>
    <xsl:value-of select="concat(@T,',',@conc,',')"/>
    <xsl:for-each select="me:firstOrderRate">
      <xsl:value-of select="concat(.,',')"/>
    </xsl:for-each>
    <xsl:for-each select="me:firstOrderLoss">
      <xsl:value-of select="concat(.,',')"/>
    </xsl:for-each>
    <xsl:text>&#xa;</xsl:text>
  </xsl:template>
</xsl:stylesheet>
