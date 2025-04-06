@echo off
cd out
del *.ppm
cd ..
renderRose.exe
set count=0
for %%A in ("out/*.ppm") do set /a count+=1
for /F %%a in ('echo prompt $E ^| cmd') do @set "ESC=%%a"
if %count% == 1 (
ffmpeg -i "out/renderRose_0.ppm" "renderRose.png" -y
echo %ESC%[32mexport image to renderRose.png%ESC%[0m
) else (
ffmpeg -framerate 30 -start_number 0 -i "out/renderRose_%%d.ppm" -c:v libx264 -pix_fmt yuv420p "renderRose.mp4" -y
echo %ESC%[36mexport video to renderRose.mp4%ESC%[0m
)
pause