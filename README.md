## Nori Version 2

Nori is a simple ray tracer written in C++. It runs on Windows, Linux, and
Mac OS.

# Installation
## Linux / Mac OS X

```
git --clone recursive https://github.com/alexus37/NoriV2.git
cd NoriV2
mkdir build
cd build
cmake-gui ..
```
(You can use `cmake` instead of `cmake-gui`)

After the Makefiles are generated, simply run make to compile all dependencies and Nori itself.

```
make -j 4
```

This can take quite a while; the above command compiles with four processors at the same time. Note that you will probably see many warning messages while the dependencies are compiled you can ignore them.

## Windows

Begin by installing Visual Studio 2013 (older versions won't do) and a reasonably recent (≥ 3.x) version of [CMake][cmake]. Download the version 1.2 [libpng][libpngURL] package and the [zlib][zlibURL] package and extent the PATH variable with the installation path of the libaries.   

`
setx path "%path%;C:\Program Files (x86)\GnuWin32\bin\"
`   

Start CMake and navigate to the location where you cloned the Nori repository.
After setting up the project, click the Configure and Generate button. This will create a file called ***nori.sln*** —double-click it to open Visual Studio.

The Build->Build Solution menu item will automatically compile all dependency libraries and Nori itself; the resulting executable is written to the Release or Debug subfolder of your chosen build directory. Note that you will probably see many warning messages while the dependencies are compiled—you can ignore them.


## Arch

Download http://ala.seblu.net/packages/l/libpng12 and install with pacman.

[cmake]: http://www.cmake.org/download/
[libpngURL]: http://gnuwin32.sourceforge.net/packages/libpng.htm
[zlibURL]: http://gnuwin32.sourceforge.net/packages/zlib.htm

