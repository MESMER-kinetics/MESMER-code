echo off
setlocal
set starttime=%time%

cd pentyl 
"../../Windows VC8/Mesmer/Mesmer.exe" "pentyl_isomerization_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
copy "./mesmer.log" "./baselines/Win32/mesmer.log"
cd ..

cd HSO2
"../../Windows VC8/Mesmer/Mesmer.exe" "HSO2_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
copy "./mesmer.log" "./baselines/Win32/mesmer.log"
cd ..

cd "cyclopropene isomerization"
"../../Windows VC8/Mesmer/Mesmer.exe" "Cyclopropene_isomerization_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
copy "./mesmer.log" "./baselines/Win32/mesmer.log"
cd ..

cd "OH acetylene association" 
"../../Windows VC8/Mesmer/Mesmer.exe" "OH_acetylene_association_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
copy "./mesmer.log" "./baselines/Win32/mesmer.log"
cd ..

cd "Acetyl O2 association" 
"../../Windows VC8/Mesmer/Mesmer.exe" "Acetyl_O2_association.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
copy "./mesmer.log" "./baselines/Win32/mesmer.log"
cd ..

echo Start Time=%starttime% - End Time=%time%

:: This line below makes DOS to create a system Beep (err yup)
echo 
echo 
echo 

