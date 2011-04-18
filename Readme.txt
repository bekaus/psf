Peak Shape Function Framework
-----------------------------

Compiles on Windows Server 2008. Tested on Ubuntu 8.04, OpenSUSE 11.1 and Fedora
11.

Dependencies:
A current development snapshot of vigra.

How to build:
./cmake .
./ccmake .. (set include path to vigra, if not /usr/include) (depends only on vigra headers!)
./make
./make test
