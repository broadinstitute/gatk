Docker for GATK4-protected
--------------------------

In this directory lives a simple, *unsupported* docker script.

This doc does not give usage info on Docker and assumes you are already versed in its use.

Please see the excellent docker documentation if you need info about docker:  https://docs.docker.com/

Notes:
- Image is built on ``ubuntu:14.04``
- Once image is built, the hellbender-protected.jar (symlink) is found in ``/root`` (i.e. ``$HOME``).
- HDF5 jni library is in ``/usr/lib/jni``, since this is the default for Ubuntu 14.04.  This has been added to the ``JAVA_LIBRARY_PATH`` environment variable.

#### Create docker image

From this directory, run: ``sudo docker build -t my_tag_name/my_tag_version .`` to create your image

#### Run gradle test

Currently, this is a bit manual.

```
# bash is the default command.
sudo docker run -i -t <image> 

# On the docker prompt
cd hellbender-protected
gradle test
```
