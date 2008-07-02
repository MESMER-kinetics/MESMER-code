echo off

cd HSO2
"../../Windows VC8/Mesmer/Mesmer.exe" "HSO2_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
cd ..

cd "cyclopropene isomerization"
"../../Windows VC8/Mesmer/Mesmer.exe" "Cyclopropene_isomerization_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
cd ..

cd pentyl 
"../../Windows VC8/Mesmer/Mesmer.exe" "pentyl_isomerization_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
cd ..

cd "OH acetylene association" 
"../../Windows VC8/Mesmer/Mesmer.exe" "OH_acetylene_association_test.xml"
copy "./mesmer.test" "./baselines/Win32/mesmer.test"
cd ..
