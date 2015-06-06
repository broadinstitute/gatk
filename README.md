[![Build Status](https://magnum.travis-ci.com/broadinstitute/hellbender-protected.svg?token=yyvS66phxJv6t5Zzsp8M&branch=master)](https://magnum.travis-ci.com/broadinstitute/hellbender-protected)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/hellbender-protected/badge.svg?branch=master&t=fjUaFR)](https://coveralls.io/r/broadinstitute/hellbender-protected?branch=master)

Hellbender Protected
====================

GATK4 development of the license-protected part of the toolkit

Requirements
------------
* Java 8

* Gradle 2.2.1

* HDF5-Java JNI Libraries Release 2.6


Read Hellbender's README
------------------------

Please refer to Hellbender's public repo readme file for general guidelines and how to set-up your development enviroment:

https://github.com/broadinstitute/hellbender/blob/master/README.md

Get HDF5-Java JNI Libraries Set-up
----------------------------------

### Library acquisition

The Java portion of the HDF5 IO library is automatically pulled as a Maven repository dependency (org.hdfgroup:hdf-java:2.6.1).

The native JNI C portion however depends on your development environment.

First you need to download a binary from http://www.hdfgroup.org or build it from source.

#### Ubuntu (Linux) 

You should be able to get it from the standard packages. 

```
   sudo apt-get install hdfview
```

At least this works with release 12.04 LTS (Travis-CI images release).

#### MacOSX developers:

Unfortunatelly for MacOSX, their build seems to be broken so one has to rely on the pre-compiled binary:

```
  curl https://www.hdfgroup.org/ftp/HDF5/releases/HDF-JAVA/HDF-JAVA-2.6/bin/macintel64/hdf-java/lib/macosx/libjhdf5.jnilib > wherever-i-downloaded-the-jnilib/libjhdf5.jnilib
``` 

### Get ```gradle test``` to work.

The VM will search for the JNI library in the path indicated by the ```java.library.path``` system property which
by default is set to the system "standard" library locations (e.g. /usr/lib will be included in most Unix flavors).

If you cannot or you don't want to install libjhdf5.jnilib in one of those standard locations and your VM's default ```java.library.path``` does not include your choice, 
you will need to tell explicitly to the build process where to look for it. 

Here you have a couple of options:

1. set the ```JAVA_LIBRARY_PATH``` enviroment variable to the value you want to set ```java.library.path``` to during testing,

2. or define the gradle project property ```testJavaLibraryPath=wherever-i-downloaded-the-jnilib``` in ```~/.gradle/gradle.properties```.

Please refrain from using gradle.properties in the root project directory for this as you don't 
want to share your set-up with other developers thru the source repo; .gitignore should prevent this from happening for now
but is best to avoid it all together as in the future we might want to use gradle.properties for common set-up.

### Get ```java -jar hellbender.jar``` to work.

If you didn't need to indicate the location of ```libjhdf5.jnilib``` explicitly to get the testing working then you are all set.

Otherwise you will need to tell the VM each time as well like so:

```
java -Djava.library.path=wherever-i-downloaded-the-jnilib -jar hellbender-prorected.jar ...
```
