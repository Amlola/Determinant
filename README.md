# Matrix
Matrix class for working with standard linear algebra methods

## Dependencies
1. g++
2. CMake 3.10 version (or higher)
3. GTest (for testing)

## Compiling 

### Matrix

To compile matrix:

``` cmd
$ mkdir build
$ cmake -S ./ -B build/Release
$ cmake --build build/Release
```

### Test

To compile unit tests you need the ```gtest``` library:

``` cmd
$ cd test
$ mkdir build
$ cd build
$ сmake ..
$ make
```

## Run the program:

### Matrix
``` cmd
$ ./build/Debug/discr
```

### Test:
``` cmd
$ ./UnitTest
```
