if not(libisloaded('libtest_matlab'))
    loadlibrary('libtest_matlab', 'test_matlab.h')
end


xp = libpointer('doublePtr', 1.0)
yp = libpointer('doublePtr', 0.0)

calllib('libtest_matlab', 'test_matlab_MOD_test', yp, xp)

xp.value
yp.value

clear xp yp
