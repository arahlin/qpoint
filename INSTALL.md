Some users may want to build the `C` library for linking to their own projects,
or have more control over the build process.

To just build and install the `C` library without OpenMP support:

```
make
make install
```

To enable OpenMP support in the `C` library:

```
make ENABLE_OMP=y
```

To build a shared library instead of a static one:

```
make ENABLE_SHARED=y
```

To build a "lite" version of the library without support for OpenMP or any of the healpix mapmaking backend:

```
make ENABLE_LITE=y
```

To install to a user directory (by default this is `$HOME/.local`):

```
make install-user
```

To install to a specific prefix (the directory that will contain the `lib/`
and `include/` subfolders):

```
make PREFIX=/your/install/prefix install
```
