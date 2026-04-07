# Windows
In Windows, although clang-cl is installed (in visual studio) and
.pythranrc/pythran-win32.cfg is set, pythran still use cl.exe and then cause
compile failed.

Directly run `pythran _gdca.cpp` would use clang-cl.exe and pass the compile,
but link step is failed due to missing link path info.

# MacOS
Not tested
