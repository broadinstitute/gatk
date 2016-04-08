[![Build Status](https://travis-ci.org/broadinstitute/gatk-protected.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk-protected)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gatk-protected/badge.svg?branch=master&service=github)](https://coveralls.io/github/broadinstitute/gatk-protected?branch=master)

GATK4-Protected (codename Hellbender-protected)
===============================================

GATK4 development of the license-protected part of the toolkit

This README is aimed at developers.  For user information, please see the [GATK 4 forum](http://http://gatkforums.broadinstitute.org/gatk/categories/gatk-4-alpha)

Requirements
------------
* R 3.1.3 see additional requirements below: [R package requirements](#r-required-packages)

* Java 8

* HDF5 1.8.13 

* HDF5-Java JNI Libraries Release 2.9 (2.11 for Macs)

* Gradle 2.12 is needed for building the GATK. We recommend using the `./gradlew` script which will 
download and use an appropriate gradle version automatically


Read GATK 4 README
------------------------

Please refer to the GATK 4 public repo readme file for general guidelines and how to set-up your development environment:

https://github.com/broadinstitute/hellbender/blob/master/README.md


#### R Required Packages
R packages can be installed using the install_R_packages.R script inside the scripts directory.


Get HDF5-Java JNI Libraries Set-up
----------------------------------

There are two external libraries needed for HDF5 support in GATK4-protected:

1. hdf -- native code only.
2. hdf-java -- includes both Java (JAR files) and native JNI code (.so/.dynlib files). 

*The Maven repository (org.hdfgroup:hdf-java:2.6.1) for the HDF5 IO library is out of date and should not be used.*

For more information about HDF:  https://www.hdfgroup.org/

Developer note:

The gradle build will handle the Java (hdf-java) dependency, without any intervention from the user, if instructions for your platform are followed (see below).  If IntelliJ is configured correctly, it will
  automatically create the dependency to the JARs in your project.  You will still need to note the location of the JNI native files
  from the description below for your platform.

#### Ubuntu (Linux) 12.10 and above (requires sudo)

*We do not guarantee that this will work with versions of Ubuntu after 15.04.*

You simply need to install the hdfview package, which includes hdf and hdf-java:

```
   sudo apt-get install hdfview
```

This will install all of the required HDF5 libraries (hdf and hdf-java).  

By default:
- The jnilib native files will be installed at: ``/usr/lib/jni/``  This location can be used in the instructions below.
- The JAR files will be installed in ``/usr/share/java/``.
  

#### MacOSX (requires admin):

You simply need to install hdfview:

1. Download the binary (https://www.hdfgroup.org/ftp/HDF5/hdf-java/current/bin/).  Select the darwin dmg file.
2. Launch the installer and follow any instructions.

This will install all of the required HDF5 libraries (hdf and hdf-java).

By default:
- The jnilib native files will be installed at: ``/Applications/HDFView.app/Contents//Resources/lib/``  This location can be used in the instructions below.
- The JAR files will be installed in ``/Applications/HDFView.app/Contents/Java/``.


#### Travis CI (Ubuntu (Linux) 12.04 LTS)

``.travis.yml`` implements the installation from binaries as per instructions above.


#### Broad VMs

Building on a Broad VM is similar to Ubuntu, except that you use a dotkit, instead of installing HDFview.

```
use .hdfview-2.9
```

- The jnilib native files will be in:  `` /broad/software/free/Linux/redhat_6_x86_64/pkgs/hdfview_2.9/HDFView/lib/linux/``
- The JAR files will be in: `` /broad/software/free/Linux/redhat_6_x86_64/pkgs/hdfview_2.9/HDFView/lib/``

The gradle build is already configured to search the JAR directory.


#### Other platforms (from binaries)

These limited and *mostly untested* instructions must be used when applications and libraries cannot be installed into default locations.

1. Download the hdfview binary from: https://www.hdfgroup.org/ftp/HDF5/releases/HDF-JAVA/hdf-java-2.9/hdfview/
2. Locate the jar files that were installed.  Use this location in the gradle build command in step 5.
3. Locate the JNI native files (.so/.dylib) that were installed.  Use this location in the instructions below. 
4. Rebuild the hellbender-protected.jar and provide the directory that contains jhdf.jar.

See example:

```
# Linux 64-bit example
HDF_DIR=/opt/hdf/
wget "https://www.hdfgroup.org/ftp/HDF5/releases/HDF-JAVA/hdf-java-2.9/hdfview/hdfview_install_linux64.bin" -O ${HDF5_DIR}/hdfview_install_linux64.bin
chmod +x ${HDF5_DIR}/hdfview_install_linux64.bin
${HDF5_DIR}/hdfview_install_linux64.bin

# Jar files will be in ${HDF5_DIR}/HDFView/lib
# JNI lib files will be in ${HDF5_DIR}/HDFView/lib/linux

# When calling gradle commands, you must add: -Pcustom.jar.dir=${HDF5_DIR}/HDFView/lib
#  or put the jar location in your ~/.gradle/gradle.properties
#  custom.jar.dir=${HDF5_DIR}/HDFView/lib

./gradlew -Pcustom.jar.dir=${HDF5_DIR}/HDFView/lib build.gradle shadowJar
```

### Get ```./gradlew test``` to work.

The VM will search for the JNI library in the path indicated by the ```java.library.path``` system property which
by default is set to the system "standard" library locations (e.g. /usr/lib will be included in most Unix flavors).

If you cannot or you don't want to install libjhdf5.jnilib in one of those standard locations and your VM's default ```java.library.path``` does not include your choice, 
you will need to tell explicitly to the build process where to look for it. 

Here you have a couple of options:

1. set the ```JAVA_LIBRARY_PATH``` environment variable to the value you want to set ```java.library.path``` to during testing,

2. or define the gradle project property ```testJavaLibraryPath=wherever-i-downloaded-the-jnilib``` in ```~/.gradle/gradle.properties```.

Please refrain from using gradle.properties in the root project directory for this as you don't 
want to share your set-up with other developers through the source repo; .gitignore should prevent this from happening for now
but is best to avoid it all together as in the future we might want to use gradle.properties for common set-up.

### Create the jar file

`` ./gradlew shadowJar ``

A jar file will appear in ``build/libs``.


### Get ```java -jar hellbender-protected.jar``` to work.

If you didn't need to indicate the location of ```libjhdf5.jnilib``` explicitly to get the testing working then you are all set.

Otherwise you will need to tell the VM each time as well like so:

```
java -Djava.library.path=wherever-i-downloaded-the-jnilib -jar hellbender-protected.jar ...
```

The CNV case and PoN workflows (description and examples)
---------------------------------------------------------

This can be found [here](http://gatkforums.broadinstitute.org/gatk/discussion/6791/description-and-examples-of-the-steps-in-the-cnv-case-and-cnv-pon-creation-workflows)


Running the CNV case and PoN creation Workflows with premade Queue scripts
--------------------------------------------------------------------------

For [Broad Internal instructions](http://gatkforums.broadinstitute.org/gatk/discussion/6786/howto-run-gatk-cnv-using-premade-queue-scripts-broad-internal)

