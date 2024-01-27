@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="ORION" (taskkill /f /pid 10184)
if /i "%LOCALHOST%"=="ORION" (taskkill /f /pid 18756)

del /F cleanup-ansys-ORION-18756.bat
