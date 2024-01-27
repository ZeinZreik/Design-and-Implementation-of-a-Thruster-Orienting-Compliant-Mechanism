@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="ORION" (taskkill /f /pid 15848)
if /i "%LOCALHOST%"=="ORION" (taskkill /f /pid 18040)

del /F cleanup-ansys-ORION-18040.bat
