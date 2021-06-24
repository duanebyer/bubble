# Bubble

`bubble` is a library for producing random events according to an arbitrary
multi-dimensional distribution. It works best in fewer than 10 dimensions. The
approach used is very similar to the [FOAM][1] library: the distribution is
approximated by an n-d tree, with events being uniformly produced within each
cell of the tree. Rejection sampling is then applied to improve the quality of
the events further.

The `bubble` library is still in development.

[1]: http://jadach.web.cern.ch/jadach/Foam/Index.html
