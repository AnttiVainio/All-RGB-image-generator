Images generated by the program and Windows and Linux binaries compiled by Antti Vainio are available on his site at:
http://www.anttivainio.net/allrgb

A video of this program generating images can be seen here:
https://www.youtube.com/watch?v=ZelU28SUB_k

This program creates random images in which each pixel has a unique color. The program also tries to use as wide variety of colors from the RGB color space as possible. The resulting image and some intermediate steps are saved into BMP-files.

The program only has a command line interface and doesn't have any interactivity. Making changes to the settings to generate different kind of images requires recompiling the program. The program can be simply started without any arguments and the generated images are automatically placed in the images-folder.

The goal of this program was to generate these images very quickly while still producing reasonably nice quality. Creating a massive 4K image takes only around half a minute. It could be possible to create better images that have less noise but it would take much more time.

**NOTES:**

Any changes to the settings such as the image size, shape, color settings, etc. have to be made in imagestuff.cpp after which the program has to be recompiled. The settings are at the beginning of the file with comments explaining everything.

The code also contains OpenMP pragmas that will make the program multithreaded when compiled with OpenMP. If using the Makefile, you can enable OpenMP by compiling with "make openmp".
