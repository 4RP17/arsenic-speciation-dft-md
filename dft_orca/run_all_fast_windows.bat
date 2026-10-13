@echo off
echo Running ORCA fast Opt/Freq jobs with full ORCA path...

C:\ORCA_6.1.1\orca.exe 01_H3AsO3_AsIII_fast_optfreq.inp > 01_H3AsO3_AsIII_fast_optfreq.out
C:\ORCA_6.1.1\orca.exe 02_H2AsO4_AsV_fast_optfreq.inp > 02_H2AsO4_AsV_fast_optfreq.out
C:\ORCA_6.1.1\orca.exe 03_HAsO4_AsV_fast_optfreq.inp > 03_HAsO4_AsV_fast_optfreq.out
C:\ORCA_6.1.1\orca.exe 04_H2AsO3_AsIII_anion_optional_fast_optfreq.inp > 04_H2AsO3_AsIII_anion_optional_fast_optfreq.out

echo Finished.
findstr /C:"ORCA TERMINATED NORMALLY" *.out
pause