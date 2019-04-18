function [] = CompileMexFiles()
    %Compiles mex files that accelerate B-spline evaluation
mex -v CFLAGS="$CFLAGS -std=c99" evalBin.c
mex -v CFLAGS="$CFLAGS -std=c99" evalBSpline.c
mex -v CFLAGS="$CFLAGS -std=c99" evalBinTimesY.c
end

