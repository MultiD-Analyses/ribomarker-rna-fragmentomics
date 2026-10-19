@echo off
rem ribomarker: the RiboMarker pipeline in Docker.
rem Starts ribomarker_windows.ps1 without changing the PowerShell script setting (execution policy).
rem Usage: ribomarker SHEET OUTPUT [options], ribomarker --demo OUTPUT, ribomarker --help
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0ribomarker_windows.ps1" %*
exit /b %ERRORLEVEL%
