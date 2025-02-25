#!/bin/sh

# qashell.sh
# mesmer
#
# Created by Chi-Hsiu Liang on 13/05/2010.

# module switch intel gnu/7.2.0

starttime=`date`
directive=
otfn=mesmer.test

# To execute the 
if [ -f "../bin/mesmer" ]; then
  executable="../../bin/mesmer"
else
  echo "----------------------------------------------------------------------"
  echo "The executable mesmer cannot be found, the QA test stopped."
  echo "Please make sure the executable mesmer is located in the \"bin\" folder."
  if [ "$TERM_PROGRAM" == "Apple_Terminal" ]; then 
    echo "If you use XCode to compile mesmer, make sure you set the \"Copy Files\""
    echo "build phase of target mesmer to have Destination \"Absolute Path\", copying the full"
    echo "path of mesmer\bin to the \"Full Path\" box. This will make sure you get a fresh executable"
    echo "each time you compile mesmer."
    echo "----------------------------------------------------------------------"
    exit
  fi
fi

tfn=mesmer.test
lfn=mesmer.log
bline=baselines/Linux64/
outf=out.xml

case $1 in
  -u)
    directive=-q
    otfn=test.test
    echo "----------------------------------------------------------------------"
    echo "User QA check mode: output will copy to test.test in the baselines folder."
    echo "Use your own \"diff\" program to check changes between test.test and "
    echo "mesmer.test in baseline folders. Please ensure the original files in "
    echo "the baseline folders were not previously modified by user."
    echo "----------------------------------------------------------------------"
    ;;
  -o)
    echo "----------------------------------------------------------------------"
    echo "Developer QA check mode: output will overwrite the baselines. Use"
    echo "\"SVN check for modifications\" to check the changes compared with the baselines."
    echo "----------------------------------------------------------------------"
    ;;
  "")
    echo ""
    echo "----------------------------------------------------------------------"
    echo "This file execute mesmer executable on several quality assessment files."
    echo "There are two options to execute this file:"
    echo ""
    echo "1. user test mode; simply to confirm that if the files (mesmer.test) "
    echo "generated by current machine is identical to the mesmer baselines."
    echo ""
    echo "Syntax:"
    echo ""
    echo "    \MESMER_PATH\examples>./Examples.sh -u"
    echo ""
    echo "2. developer SVN test mode; it overwrites mesmer.test files in the baseline"
    echo "   folders, thus the developer can use SVN check for modifications from the "
    echo "   examples folder to see what files are changed and the detail of the "
    echo "   modifications. To use it in the developer mode simply give a -o directive"
    echo "   after QA command."
    echo ""
    echo "Syntax:"
    echo ""
    echo "    \MESMER_PATH\examples>./Examples.sh -o"
    echo "----------------------------------------------------------------------"
    echo ""
    ;;
esac

pwd
$executable -V
mpicommand="mpirun -n 12 $executable" 
cmdline=$mpicommand

cd AcetylO2
$cmdline Acetyl_O2_associationEx.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ]; then
  cp ./$lfn ./$bline$lfn 
fi
cd ..

cd AcetylPrior
$cmdline AcetylPrior.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn
fi
cd ..

cd Butyl_H_to_Butane
$cmdline Butyl_H_to_Butane.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn 
fi
cd ..

cd Ethyl_H_to_Ethane
testName="Ethyl_H_to_Ethane"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi

testName="Ethyl_H_to_Ethane_inversion"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd i-propyl
$cmdline ipropyl_LM.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn 
fi
cd ..

cd Methyl_H_to_Methane
$cmdline Methyl_H_to_Methane.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn 
fi

$cmdline -N Methyl_H_to_Methane_FTST.xml $directive
cp ./Methyl_H_to_Methane_FTST.test ./$bline/Methyl_H_to_Methane_FTST.test
if [ "$1" == "-o" ] ; then
  cp ./Methyl_H_to_Methane_FTST.log ./$bline/Methyl_H_to_Methane_FTST.log 
fi
cd ..

cd reservoirSink
$cmdline reservoirSinkAcetylO2.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn 
fi
cd ..

cd spin_forbidden_kinetics
$cmdline -N HCCH_methylene.xml  $directive
cp ./HCCH_methylene.test ./$bline/HCCH_methylene.test
if [ "$1" == "-o" ] ; then
  cp ./HCCH_methylene.log ./$bline/HCCH_methylene.log 
fi

$cmdline -N LZ_test.xml $directive
cp ./LZ_test.test ./$bline/LZ_test.test
if [ "$1" == "-o" ] ; then
  cp ./LZ_test.log ./$bline/LZ_test.log
fi

$cmdline -N WKB_test.xml $directive
cp ./WKB_test.test ./$bline/WKB_test.test
if [ "$1" == "-o" ] ; then
  cp ./WKB_test.log ./$bline/WKB_test.log
fi

$cmdline -N ZNtest.xml $directive
cp ./ZNtest.test ./$bline/ZNtest.test
if [ "$1" == "-o" ] ; then
  cp ./ZNtest.log ./$bline/ZNtest.log
fi
cd ..

cd OH_NO_to\ HONO
$cmdline OH_NO_HONO.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn 
fi
cd ..

cd OH-acetylene
$cmdline OH_HCCH-irreversibleBim-publish.xml -o $outf $directive
cp ./$tfn ./$bline$otfn
if [ "$1" == "-o" ] ; then
  cp ./$lfn ./$bline$lfn 
fi
cd ..

cd "AnalyticalRepresentation"
testName=Chebyshev
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName=ChebyshevCK
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName=Plog
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName=PlogCK
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd SensitivityAnalysis
testName=ipropyl_SA
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName=pentyl_isomerization_SA
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "cis_to_trans_But-2-ene"
testName=cis_to_trans_But-2-ene
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName=cis_to_trans_But-2-eneEx
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "C4H9O2_NO2_to_C4H9O2NO2"
testName=C4H9O2_NO2_association
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "CH3O2_NO2_to_CH3O2NO2"
testName=CH3O2_NO2_associationEx2
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "ErrorPropagation"
testName=ipropyl_EP
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "2Methyl_to_Ethane"
testName=CH3_CH3_to_C2H6
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "Glyoxyl"
testName=Glyoxal
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi

testName=Glyoxal_Exchange
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "DefinedTunellingCoefficients"
testName=OH+methanol
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "diamond"
testName=diamond-gaussianModel
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "Tunnelling"
testName="H+H2,T+T2"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName="Acetyl_O2_associationEx1"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
testName="Acetyl_O2_associationEx2"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "CoupledRotors"
testName="ipropyl_UCR"
extension="_out"
$cmdline -N $testName.xml
if [ "$1" == "-o" ] ; then
  cp ./mesmer_out.xml ./$bline$testName$extension.xml
fi
testName="ipropyl_UQR"
$cmdline -N $testName.xml
if [ "$1" == "-o" ] ; then
  cp ./mesmer_out.xml ./$bline$testName$extension.xml
fi
testName="ipropyl_CCR"
$cmdline -N $testName.xml
if [ "$1" == "-o" ] ; then
  cp ./mesmer_out.xml ./$bline$testName$extension.xml
fi
testName="Octyl_H_to_Octane_UCR"
$cmdline -N $testName.xml
if [ "$1" == "-o" ] ; then
  cp ./mesmer_out.xml ./$bline$testName$extension.xml
fi
testName="Octyl_H_to_Octane_UQR"
$cmdline -N $testName.xml
if [ "$1" == "-o" ] ; then
  cp ./mesmer_out.xml ./$bline$testName$extension.xml
fi
testName="Octyl_H_to_Octane_CCR"
$cmdline -N $testName.xml
if [ "$1" == "-o" ] ; then
  cp ./mesmer_out.xml ./$bline$testName$extension.xml
fi
cd ..

cd "OH_Ethene"
testName="Ethylene_abstraction"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "Methyl_O_to_Methoxy"
testName="Methyl_O_to_Methoxy"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "POandO2"
testName="PO+O2"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "H2Ominimal"
testName="H2Ominimal"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi

testName="H2OBend"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "Butan-1-ol"
testName="Butoxy_H_to_Butan-1-ol"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi

testName="Butan-1-peroxy_H_to_Butan-1-hydroperoxy"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "InternalRotor"
testName="InternalRotor"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

cd "H2O_SiO2_to_H2SiO3"
testName="H2O_SiO2_to_H2SiO3_Aij"
$cmdline -N $testName.xml
cp ./$testName.test ./$bline$testName.test
if [ "$1" == "-o" ] ; then
  cp ./$testName.log ./$bline$testName.log
fi
cd ..

echo Start Time=$starttime - End Time=`date`

