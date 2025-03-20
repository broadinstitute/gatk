# GATK Developer Guide
This document is intended to serve as a basic introduction to developing for the GATK.

TODO

## A basic tool walkthrough
Let's start with a walkthrough of one of the most basic GATK tools: PrintReads.  Here is all the code for PrintReads (minus imports and doc comments):

```java
@CommandLineProgramProperties(
        summary = "Prints reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file.",
        oneLineSummary = "Print reads in the SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
@WorkflowProperties
public final class PrintReads extends ReadWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;
    private SAMFileGATKReadWriter outputWriter;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        outputWriter.addRead(read);
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
```

PrintReads accepts a file in SAM/BAM/CRAM format, and iterates through the reads in the file, printing each to a specified output file.

TODO: finish the description, or get rid of this, I don't know

## Classes you might want to extend
This section describes classes you might want to extend when making new tools in GATK.

### CommandLineProgram
`CommandLineProgram` is the highest level class you should extend when making a new tool.  You should really only extend `CommandLineProgram` directly if your tool does not use any of the infrastructure in the GATK engine for interfacing with genomic data files and performing operations on genomic data types.

#### Methods to override when extending CommandLineProgram
```java
protected void onStartup() {}
```
`onStartup` is your setup method.  It runs after the command-line arguments are parsed.  If you need to open a file reader/writer, validate command-line arguments, or initialize any other variables you need to use, this is probably the place to do it.  By default, `onStartup` does nothing, and it is not actually a requirement to override it.

```java
protected abstract Object doWork();
```
`doWork` is where most of your logic will probably go.  It's the only method you are required to implement if you directly extend `CommandLineProgram`.  If any of your code can throw a `RuntimeException`, it's also helpful to do that here, because then it will be logged appropriately.  

```java
protected void onShutdown() {}
```
`onShutdown` pairs well with `onStartup` because it is your cleanup method.  Closing readers and writers, logging results, and reporting failures and their implications are all operations that fit well into this method.  `onShutdown` will execute even if `doWork` throws an exception.  By default, `onShutdown` does nothing, and it is not a requirement to override it.

### GATKTool
`GATKTool` is the base class for (almost) all tools in GATK.  It extends `CommandLineProgram`.  Extending `GATKTool` is the right choice if extending one of the Walker classes isn't the right fit for what you want your tool to do.  It provides base functionality from the GATK engine for performing operations with different types of genomic data files.  This includes a number of default arguments for things like read files, references, and interval lists, along with quality of life features like a progress meter that can be invoked in subclasses to log progress to the user.

#### Methods to override when extending GATKTool
```java
public void onTraversalStart() {}
```
`onTraversalStart` is your initialization method for tools that extend GATKTool.  It is called after command-line arguments have been parsed and after `GATKTool` performs some initial setup such as initializing the progress meter.  To examine all of the initialization performed by `GATKTool`, see `GATKTool`'s implementation of `onStartup`.
Opening readers and writers, and processing command-line arguments specific to your tool are tasks you will likely want to implement in `onTraversalStart`.

```java
public abstract void traverse();
```
`traverse` is where most of your tool logic is likely to go.  This is the only method you are required to implement if you directly extend `GATKTool`.  If code within this method throws a `RuntimeException`, it will be properly logged.

```java
public Object onTraversalSuccess() { return null; }
```
`onTraversalSuccess` is called after `traverse` finishes successfully (i.e with no uncaught exceptions).  This can be used to close resources and can also return a value to be printed by the engine.

```java
public void closeTool() {}
```
`closeTool` is called whether `traverse` has succeeded or not.  If `traverse` has succeeded, it is called after `onTraversalSuccess`.  It can be used for cleanup activities that have to happen regardless of the outcome of `traverse`.


### The Walker classes
GATK provides a number of Walker base classes that can be extended for building tools that iterate over entities to perform some operation.  For most cases, when deciding what class to extend for your tool, you should start with checking to see if there is a walker that serves your purposes.

The walker classes you are mostly likely to use are:

#### VariantWalker
A `VariantWalker` iterates over a source of variants (e.g. one or more VCF files) and processes one variant at a time.  This processing is done in the `apply` method, which is the primary method you will override if extending `VariantWalker`.  For an example demonstrating in more depth how to extend `VariantWalker`, see `ExampleVariantWalker`.

#### ReadWalker
A `ReadWalker` iterates over a source of reads (e.g. one or more BAM files) and processes one read at a time.  This processing is done in the `apply` method, which is the primary method you will override if extending `ReadWalker`.  For an example demonstrating in more depth how to extend `ReadWalker`, see `ExampleReadWalkerWithVariants` or `ExampleReadWAlkerWithReference`.

#### FeatureWalker
A `FeatureWalker` iterates over a source of features and processes one at a time.  A `Feature`, in the context of GATK, is basically a locus on a reference sequence and some other information tied to that locus.  For example, a `BEDFeature` is an entry in a BED file.  For an example demonstrating how to extend `FeatureWalker`, see `ExampleFeatureWalker`.

