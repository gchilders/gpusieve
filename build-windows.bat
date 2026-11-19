@echo off
rem Windows build of gpusieve as a BOINC application, from a "x64 Native Tools Command Prompt for VS" (MSVC) with the
rem CUDA Toolkit. BOINC is a source tree whose win_build solution has been built (x64 Release): set BOINC below.
rem MSVC has no __int128, so gs_port.h uses its portable 128-bit code there (the same code as -DGS_NO_INT128 on Linux).

if "%BOINC%"=="" set BOINC=C:\src\boinc
if "%CUDA_PATH%"=="" set CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8
set BOINC_LIB=%BOINC%\win_build\Build\x64\Release

cl /nologo /O2 /MT /EHsc /std:c++17 /DHAVE_BOINC /D_CRT_SECURE_NO_WARNINGS ^
  /I"%BOINC%" /I"%BOINC%\api" /I"%BOINC%\lib" /c boinc_glue.cpp /Foboinc_glue.obj || exit /b 1

nvcc ^
  -gencode arch=compute_120,code=sm_120 ^
  -gencode arch=compute_100,code=sm_100 ^
  -gencode arch=compute_90,code=sm_90 ^
  -gencode arch=compute_89,code=sm_89 ^
  -gencode arch=compute_86,code=sm_86 ^
  -gencode arch=compute_80,code=sm_80 ^
  -gencode arch=compute_75,code=sm_75 ^
  -gencode arch=compute_70,code=sm_70 ^
  -gencode arch=compute_61,code=sm_61 ^
  -gencode arch=compute_60,code=sm_60 ^
  -gencode arch=compute_52,code=sm_52 ^
  -gencode arch=compute_50,code=sm_50 ^
  -gencode arch=compute_80,code=compute_80 ^
  -O3 -std=c++17 --cudart static -Wno-deprecated-gpu-targets ^
  -Xcompiler "/O2 /MT /EHsc" -D_CRT_SECURE_NO_WARNINGS ^
  -o gpusieve.exe gpusieve.cu boinc_glue.obj ^
  "%BOINC_LIB%\libboincapi.lib" "%BOINC_LIB%\libboinc.lib" user32.lib || exit /b 1

echo built gpusieve.exe
