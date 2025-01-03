# sharpen
Simple parallel image sharpening exercise for training purposes

This is a simple program that uses smoothing and edge detection to
sharpen up a fuzzy image. As written, it produces a Portable Grey Map
(PGM) file "sharpened.pgm" from an input "fuzzy.pgm". The PGM files
are a relatively standard text format for greyscale images which can
be viewed using, for example, the "display" program from the
ImageMagick suite. If your desktop machine doesn't have "display" then
you can try converting the file on the HPC platform, e.g. to produce a
JPEG file: "convert sharpened.pgm sharpened.jpg".

It is deliberately written to be very simple and easy to understand so
it can be used as a teaching example. It uses a very simple
parallelisation strategy, designed to illustrate Amdahl's law with a
purely serial part (file IO) and a parallel part (embarrassingly
parallel calculation). The program could easily be made orders of
magnitude faster with some simple modifications. It is deliberately
written in a very naive manner so that it runs for a measurable amount
of time. Typical timings for the compiled serial program will be a
few seconds, with around 0.1 to 0.2 seconds for IO.

Here are the results before and after the sharpening process - thanks
to Visa Korkiakoski for permission to use this image.

<p align="center">
    <img src="doc/fuzzy.png" height="300">
    <img src="doc/sharpened.png" height="300">
</p>


