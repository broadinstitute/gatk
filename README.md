[![Build Status](https://travis-ci.org/broadinstitute/hellbender-protected.svg?branch=master)](https://travis-ci.org/broadinstitute/hellbender-protected)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/hellbender-protected/badge.svg?branch=master&t=fjUaFR)](https://coveralls.io/r/broadinstitute/hellbender-protected?branch=master)

GATK4-Protected (codename Hellbender-protected)
===============================================

GATK4 development of the license-protected part of the toolkit

Requirements
------------
* R 3.1.3 see additional requirements below: [R package requirements](#r-required-packages)

* Java 8

* Gradle 2.7

* HDF5 1.8.13 

* HDF5-Java JNI Libraries Release 2.9 (2.11 for Macs)


Read GATK 4 README
------------------------

Please refer to the GATK 4 public repo readme file for general guidelines and how to set-up your development environment:

https://github.com/broadinstitute/hellbender/blob/master/README.md


#### R Required Packages
R packages can be installed using the install_R_packages.R script inside the scripts directory. Reproduced below:

```
source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
#Make sure to use http not https as this will give an "unsupported URL scheme" error
getoptUrl="http://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz"
if (!("getopt" %in% rownames(installed.packages))) {
  install.packages(getoptUrl, repos=NULL, type="source")
}
optparseUrl="http://cran.r-project.org/src/contrib/optparse_1.3.2.tar.gz"
if (!("optparse" %in% rownames(installed.packages))) {
  install.packages(optparseUrl, repos=NULL, type="source")
}
dependencies = c("naturalsort")
if (!all(dependencies %in% rownames(installed.packages()))) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}
q(save="no")
```


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

gradle -Pcustom.jar.dir=${HDF5_DIR}/HDFView/lib build.gradle shadowJar
```

### Get ```gradle test``` to work.

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

### Get ```java -jar hellbender-protected.jar``` to work.

If you didn't need to indicate the location of ```libjhdf5.jnilib``` explicitly to get the testing working then you are all set.

Otherwise you will need to tell the VM each time as well like so:

```
java -Djava.library.path=wherever-i-downloaded-the-jnilib -jar hellbender-protected.jar ...
```

The CNV case and PoN workflows (description and examples)
---------------------------------------------------------

#### Requirements
1. Java 1.8
2. A functioning GATK4-protected jar (hellbender-protected.jar)
3. HDF5 1.8.13 
4. The location of the HDF5-Java JNI Libraries Release 2.9 (2.11 for Macs).  (See platform instructions above for typical locations)  
5. Reference genome (fasta files) with fai and dict files.  This can be downloaded as part of the GATK resource bundle: http://www.broadinstitute.org/gatk/guide/article?id=1213
6. Target BED file that was used to create the PoN file.  Format details can be found `here <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ .  **NOTE:**  For the CNV tools, you will need a fourth column for target name, which must be unique across rows.
```
1       12200   12275   target1
1       13505   13600   target2
1       31000   31500   target3
1       35138   35174   target4
....snip....
```
7. PoN file (when running case samples only).  This file should be created using the Create PoN workflow (see below).
<a name="SampleName"/>
8. Sample name (when running case samples only).  This can be extracted from the input bam file by looking at the ``SM`` tag.  Use ``samtools view -H my_sample.bam | egrep SM``.

#### Case sample workflow

This workflow requires a PoN file generated by the Create PoN workflow.

If you do not have a PoN, please skip to the [Create PoN workflow](#create-pon-workflow), below ....

##### Overview of steps

- Step 1. Collect proportional coverage
- Step 2. Create coverage profile
- Step 3. Segment coverage profile
- Step 4. Plot coverage profile
- Step 5. Call segments

##### Step 1. Collect proportional coverage

###### Inputs
- bam file
- target bed file -- must be the same that was used for the PoN
- reference_sequence (required by GATK) -- fasta file with b37 reference.  

###### Outputs
- Proportional coverage tsv file -- Mx5 matrix of proportional coverage, where M is the number of targets.  The fifth column will be named for the sample in the bam file (found in the bam file ``SM`` tag).  If the file exists, it will be overwritten.

```
##fileFormat  = tsv
##commandLine = org.broadinstitute.hellbender.tools.exome.ExomeReadCounts  ...snip...
##title       = Read counts per target and sample
CONTIG  START   END     NAME    SAMPLE1
1       12200   12275   target1    1.150e-05
1       13505   13600   target2    1.500e-05
1       31000   31500   target3    7.000e-05
....snip....
```

###### Invocation

```
 java -Xmx8g -jar <path_to_hellbender_protected_jar> ExomeReadCounts -I <input_bam_file> -O <pcov_output_file_path>  -exome <target_BED> -R <ref_genome> \ 
       -transform PCOV -exonInfo FULL -groupBy SAMPLE -keepdups
```

##### Step 2. Create coverage profile

###### Inputs
- proportional coverage file from Step 1
- target BED file -- must be the same that was used for the PoN
- PoN file
- directory containing the HDF5 JNI native libraries

###### Outputs
- normalized coverage file (tsv) -- details each target with chromosome, start, end, and log copy ratio estimate
```
#fileFormat = tsv
#commandLine = ....snip....
#title = ....snip....
name    contig  start   stop    SAMPLE1
target1    1       12200   12275   -0.5958351605220968
target2    1       13505   13600   -0.2855054918109098
target3    1       31000   31500   -0.11450116047248263
....snip....
```
- pre-tangent-normalization coverage file (tsv) -- same as normalized coverage file (tsv) above, but copy ratio estimates are before the noise reduction step.  The file format is the same as the normalized coverage file (tsv).
- fnt file (tsv) -- proportional coverage divided by the target factors contained in the PoN.  The file format is the same as the proportional coverage in step 1.
- betaHats (tsv) -- used by developers and evaluators, typically, but output location must be specified.  These are the 
 coefficients used in the projection of the case sample into the (reducued) PoN.  This will be a Mx1 matrix where M is the number of targets.

###### Invocation
```
java -Djava.library.path=<hdf_jni_native_dir> -Xmx8g -jar <path_to_hellbender_protected_jar> NormalizeSomaticReadCounts -I <pcov_input_file_path> -T <target_BED> -pon <pon_file> \
 -O <output_target_cr_file> -FNO <output_target_fnt_file> -BHO <output_beta_hats_file> -PTNO <output_pre_tangent_normalization_cr_file>
```


##### Step 3. Segment coverage profile

###### Inputs
- normalized coverage file (tsv) -- from step 2.
- sample name

###### Outputs
- seg file (tsv) -- segment file (tsv) detailing contig, start, end, and copy ratio (segment_mean) for each detected segment.  Note that this is a different format than python recapseg, since the segment mean no longer has log2 applied.
```
Sample  Chromosome      Start   End     Num_Probes      Segment_Mean
SAMPLE1        1       12200   70000   18       0.841235
SAMPLE1        1       300600  1630000 337     1.23232323
....snip....
```

###### Invocation

```
java -Xmx8g -jar <path_to_hellbender_protected_jar>  PerformSegmentation  -S <sample_name> -T <normalized_coverage_file> -O <output_seg_file> -log
```


##### Step 4. Plot coverage profile

###### Inputs
- normalized coverage file (tsv) -- from step 2.
- pre-normalized coverage file (tsv) -- from step 2.
- segmented coverage file (seg) -- from step 3.
- sample name, [see above](#SampleName)

###### Outputs
- beforeAfterTangentLimPlot (png) -- Output before/after tangent normalization plot up to copy-ratio 4
- beforeAfterTangentPlot (png) -- Output before/after tangent normalization plot
- fullGenomePlot (png) -- Full genome plot after tangent normalization
- preQc (txt) -- Median absolute differences of targets before normalization
- postQc (txt) -- Median absolute differences of targets after normalization
- dQc (txt) -- Difference in median absolute differences of targets before and after normalization

###### Invocation

```
java -Xmx8g -jar <path_to_hellbender_protected_jar>  PlotSegmentedCopyRatio  -S <sample_name> -T <normalized_coverage_file> -P <pre_normalized_coverage_file> -seg <segmented_coverage_file> -O <output_seg_file> -log
```


##### Step 5. Call segments

###### Inputs
- normalized coverage file (tsv) -- from step 2.
- seg file (tsv) -- from step 3.
- sample name

###### Outputs
- called file (tsv) -- output is exactly the same as in seg file (step 3), except Segment_Call column is added.  Calls are either "+", "0", or "-" (no quotes). 
```
Sample  Chromosome      Start   End     Num_Probes      Segment_Mean      Segment_Call
SAMPLE1        1       12200   70000   18       0.841235      -
SAMPLE1        1       300600  1630000 337     1.23232323     0 
....snip....
```

###### Invocation
```
java -Xmx8g -jar <path_to_hellbender_protected_jar> CallSegments -T <normalized_coverage_file> -S <seg_file> -O <output_called_seg_file> -sample <sample_name> 
```

#### Create PoN workflow

This workflow can take some time to run depending on how many samples are going into your PoN and the number of targets you are covering.  Basic time estimates are found in the [Overview of Steps](#overview-of-steps).

##### Additional requirements
- Normal sample bam files to be used in the PoN.  The index files (.bai) must be local to all of the associated bam files.

##### Overview of steps

- Step 1. Collect proportional coverage.  (~20 minutes for mean 150x coverage and 150k targets, per sample)
- Step 2. Combine proportional coverage files  (< 5 minutes for 150k targets and 300 samples)
- Step 3. Create the PoN file (~1.75 hours for 150k targets and 300 samples) 

All time estimates are using the internal Broad infrastructure.

##### Step 1. Collect proportional coverage on each bam file

This is exactly the same as the case sample workflow, except that this needs to be run once for each input bam file, each with a different output file name.  Otherwise, the inputs should be the same for each bam file.

Please see documentation [above](#step-1-collect-proportional-coverage).

IMPORTANT NOTE: You must create a list of the proportional coverage files (i.e. output files) that you create in this step.  One output file per line in a text file (see step 2)

##### Step 2. Merge proportional coverage files

This step merges the proportional coverage files into one large file with a separate column for each samples.

###### Inputs
- list of proportional coverage files generated (possibly manually) in step 1.  This is a text file.
```
/path/to/pcov_file1.txt
/path/to/pcov_file2.txt
/path/to/pcov_file3.txt
....snip....
```

###### Outputs
- merged tsv of proportional coverage 
```
CONTIG  START   END     NAME    SAMPLE1    SAMPLE2 SAMPLE3 ....snip....
1       12191   12227   target1    8.835E-6  1.451E-5     1.221E-5    ....snip....
1       12596   12721   target2    1.602E-5  1.534E-5     1.318E-5   ....snip....
....snip....
```

###### Invocation
```
java -Xmx8g -jar  <path_to_hellbender_protected_jar> CombineReadCounts --inputList <text_file_list_of_proportional_coverage_files> \
    -O <output_merged_file> -MOF 200 
```

##### Step 3. Create the PoN file

###### Inputs
- merged tsv of proportional coverage -- generated in step 2.

###### Outputs
- PoN file -- HDF5 format.  This file can be used for running case samples sequenced with the same process.

###### Invocation
```
java -Xmx16g -Djava.library.path=<hdf_jni_native_dir> -jar <path_to_hellbender_protected_jar> CreatePanelOfNormals -I <merged_pcov_file> \
       -O <output_pon_file_full_path>
```


Running the CNV case and PoN creation Workflows with premade Queue scripts
--------------------------------------------------------------------------

*These workflows are not officially supported, but are used extensively, internally, for our evaluations.*

#### Requirements
1. Java 1.8
2. A functioning GATK4-protected jar (hellbender-protected.jar)
3. HDF5 1.8.13 
4. The location of the HDF5-Java JNI Libraries Release 2.9 (2.11 for Macs).  (See platform instructions above for typical locations)  
5. Queue scala scripts (see instructions below)
6. Reference genome (fasta files) with fai and dict files. This can be downloaded as part of the GATK resource bundle: http://www.broadinstitute.org/gatk/guide/article?id=1213
7. Target BED file that was used to create the PoN file.
8. PoN file (when running case samples only).  This file should be created using the CreatePoNPipeline workflow (example below).


#### Download Queue scala scripts
Navigate to [https://www.broadinstitute.org/gatk/download/gatk4].

Download the following files to the same directory.  This step only needs to be done once.
- CaseSampleHBExomePipeline.scala
- CreatePoNPipeline.scala
- recapseg-hb-eval.jar

#### Test that invoking the help documentation works
```
java -jar recapseg-hb-eval.jar -S CreatePoNPipeline.scala --help
```

```
java -jar recapseg-hb-eval.jar -S CaseSampleHBExomePipeline.scala --help
```

#### Example case sample run

*Parameters to the Queue jar (recapseg-hb-eval.jar) may need to be adjusted depending on your execution environment*

```
### Modify parameters below this line

# A simple text file listing each bam file on a separate line
INPUT_BAMS=case_bam_list.txt

# Output location (absolute paths recommended).  
#  Seg file will appear in $OUTDIR/tumor_pcov
#  Calls will appear in $OUTDIR/caller
OUTDIR=/home/lichtens/evals/out_case/

# Requirement #2 above
GATK4PJAR=/home/lichtens/hellbender-protected.jar

# Requirement #4 above.  Directory of HDF5 JNI shared libraries (.so/.dynlib)
HDFLOC=/usr/lib/jni/

# Requirement #6 above (This is b37 at the Broad Institute)
REF=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta

# Requirement #7 above.  Must be the same as was used to create the PoN.
TARGETS=/home/lichtens/my_target_list.bed

# Requirement #8 above.  PoN file must be created with the same target file (above) and settings.
PON=/home/lichtens/evals/out_pon/create_pon/my_blood_normals.pon

# How much memory to allocate to each job, in GB
MEM=8

# See --help for descriptions of these parameters.  
#  If -keepDups was specified for the PoN, it must be specified for case samples as well.
# -rawcov will generate a separate file for the raw counts as a separate process.
OTHER_OPTS=" -keepDups -rawcov "

#### Do not modify below this line

# Run the sample(s)
java -jar recapseg-hb-eval.jar -S CaseSampleHBExomePipeline.scala -mem ${MEM} -pon ${PON} -i ${INPUT_BAMS} -o ${OUTDIR} -hbJar ${GATK4PJAR} -r ${REF} -L ${TARGETS} -qsub -run -logDir ${OUTDIR} -hvl ${HDFLOC} ${OTHER_OPTS}

```

#### Example create PoN run

*Parameters to the Queue jar (recapseg-hb-eval.jar) may need to be adjusted depending on your execution environment*

```
### Modify parameters below this line

# A simple text file listing each bam file on a separate line
INPUT_BAMS=blood_normals_bam_list.txt

# Output location (absolute paths recommended).  
#  PoN file will appear as $OUTDIR/create_pon/$PON_FILENAME
OUTDIR=/home/lichtens/evals/out_pon/

# Requirement #2 above
GATK4PJAR=/home/lichtens/hellbender-protected.jar

# Requirement #4 above
HDFLOC=/usr/lib/jni/

# Requirement #6 above (This is b37 at the Broad Institute)
REF=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta

# Requirement #7 above.  PoN file must be created with the same target file (below) and settings.
#  This parameter is only the base filename.  Do not specify a directory.  This file will appear as
#  $OUTDIR/create_pon/$PON_FILENAME
PON_FILENAME=my_blood_normals.pon

# Requirement #8 above
TARGETS=/home/lichtens/my_target_list.bed

# How much memory to allocate to each job, in GB.  
#  Larger PoNs can require a lot of RAM (208k targets x 300 samples needs approx. 14GB RAM)
MEM=14

# See --help for descriptions of these parameters.  
#  If -keepDups was specified for the PoN, it must be specified for case samples as well.
# -rawcov will generate a separate file for the raw counts as a separate process.
OTHER_OPTS=" -keepDups -rawcov "

#### Do not modify below this line

# Run the sample(s)
java -jar recapseg-hb-eval.jar -S CreatePoNPipeline.scala -pon ${PON_FILENAME} -i ${INPUT_BAMS} -o ${OUTDIR} -hbJar ${GATK4PJAR} -r ${REF} -L ${TARGETS} -qsub -run -logDir ${OUTDIR} -mem ${MEM} -hvl ${HDFLOC} ${OTHER_OPTS}

```

