@echo off
rem Makes a "punch file" - a csv file for input to spreadsheet, etc.
rem punchout outfile.xml
rem will produce outfile.csv from a Mesmer output file outfile.xml
if exist punch.xsl (
msxsl %1 punch.xsl -o "%~dpn1.csv"
) else (
msxsl %1 ..\..\punch.xsl -o "%~dpn1.csv")
