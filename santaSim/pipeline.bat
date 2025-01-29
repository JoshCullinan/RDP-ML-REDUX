@ECHO OFF

@REM Preferably run this script as pipeline.bat [number of iterations] e.g. pipeline.bat 100 
@REM Otherwise set EPOCH in this script and do not enter a number after pipeline.bat
set /A EPOCH=1

@REM Looks at the input after pipeline.bat to see if iterations was given.
if "%~1"=="" goto iterations 
:iterations
set EPOCH=%1
GOTO END
:END

@REM Creating a "unique" identifier based on the current time for the alignment files. If you know how to do it better please do - this code is crap. I don't know BAT.
SET T=%date%-%time%
for /f "tokens=1-3 delims=/" %%I in ("%T%") do @set T=%%I-%%J-%%K
for /f "tokens=1-3 delims=:" %%I in ("%T%") do @set T=%%I-%%J-%%K
for /f "tokens=1-2 delims=." %%I in ("%T%") do @set T=%%I-%%J
call conda activate RDP

FOR /L %%A IN (1,1,%EPOCH%) DO (

    echo Starting EPOCH: %%A out of %EPOCH%

    @REM Run simulation.
    java -jar santa.jar low_recomb_rate.xml 1>>out_%%A.txt 2>>&1

    @REM Housekeeping to allow for easy file management.
    Ren "alignment_0.fa" "alignment_%%A.fa"
    Ren "recombination_events.txt" "%%A_recombination_events.txt"
    Ren "sequence_events_map.txt" "%%_Asequence_events_map.txt" 
    
    @REM Run RDP scan on the simulation
    RDP\RDP5CL.exe -f ../alignment_%%A.fa -ds >>out_%%A.txt 2>>&1

    @REM Parse the output files if they exist else make note
    if exist alignment_%%A.faRecIDTests.csv (python output_parser.py --rdpcsv "alignment_%%A.faRecombIdentifyStats.csv" --seq "sequence_events_map_%%A.txt" --rec "recombination_events_%%A.txt" --IDTests "alignment_%%A.faRecIDTests.csv") else (echo RDP did not detect any recombination on run: %%A)
    
    @REM
    md "output\alignments" >nul 2>&1
    md "output\alignments\%T%" >nul 2>&1

    del "alignment_%%A.faRecombIdentifyStats.csv", "alignment_%%A.faRecIDTests.csv", "alignment_%%A.fa.csv", "alignment_%%A.fa.rdp5" >nul 2>&1 

    @REM Move all of the output files generated into the output folder.
    move "*.txt" "..\RDP-ML\output\alignments\%T%" >nul 2>&1
    move "*.fa" "..\RDP-ML\output\alignments\%T%" >nul 2>&1
)

@REM Tidy up some extra RDP files
cd RDP
del RDP5Redolist* >nul 2>&1
cd ..

echo Fin