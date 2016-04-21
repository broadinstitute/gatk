package org.broadinstitute.hellbender.utils.test;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Builder for command line argument lists
 * It will convert old style "Argument=Value" into new style "--Argument value" strings
 * Use this only in test code.
 */
public final class ArgumentsBuilder {
    private final List<String> args= new ArrayList<>();

    public ArgumentsBuilder(){}

    public ArgumentsBuilder(Object[] args){
        for (Object arg: args){
            if (arg instanceof String){
                this.add((String) arg);
            } else {
                this.add(arg);
            }
        }
    }

    /**
     * Add a string to the arguments list
     * Strings are processed specially, they are reformatted to match the new unix style arguments
     * @param arg A string representing one or more arguments
     * @return the ArgumentsBuilder
     */
    public ArgumentsBuilder add(String arg){
        List<String> chunks = Arrays.asList(StringUtils.split(arg.trim()));
        for (String chunk : chunks){
            if(chunk.contains("=")){
                String tmp = "--"+chunk;
                args.addAll(Arrays.asList(tmp.split("=")));
            }
            else{
                args.add(chunk);
            }
        }
        return this;
    }

    /**
     * add an input file argument {@link StandardArgumentDefinitions#INPUT_LONG_NAME}
     */
    public ArgumentsBuilder addInput(File input) {
        addFileArgument(StandardArgumentDefinitions.INPUT_LONG_NAME, input);
        return this;
    }

    /**
     * add an output file argument using {@link StandardArgumentDefinitions#OUTPUT_LONG_NAME}
     */
    public ArgumentsBuilder addOutput(File output) {
        addFileArgument(StandardArgumentDefinitions.OUTPUT_LONG_NAME, output);
        return this;
    }

    /**
     * add a reference file argument using {@link StandardArgumentDefinitions#REFERENCE_LONG_NAME}
     */
    public ArgumentsBuilder addReference(File reference){
        addFileArgument(StandardArgumentDefinitions.REFERENCE_LONG_NAME, reference);
        return this;
    }

    /**
     * add a vcf file argument using {@link StandardArgumentDefinitions#VARIANT_LONG_NAME}
     */
    public ArgumentsBuilder addVCF(File fileIn) {
        addFileArgument(StandardArgumentDefinitions.VARIANT_LONG_NAME, fileIn);
        return this;
    }

    /**
     * add an argument with a file as its parameter
     */
    public ArgumentsBuilder addFileArgument(String argumentName, File file){
        Utils.nonNull(file);
        Utils.nonNull(argumentName);
        add("--" + argumentName);
        add(file.getAbsolutePath());
        return this;
    }

    /**
     * add an argument with a boolean as its parameter
     */
    public ArgumentsBuilder addBooleanArgument(String argumentName, boolean yes){
        Utils.nonNull(argumentName);
        add("--" + argumentName);
        add(yes);
        return this;
    }

    /**
     * add an argument with a given value to this builder
     */
    public ArgumentsBuilder addArgument(final String argumentName, final String argumentValue) {
        Utils.nonNull(argumentValue);
        Utils.nonNull(argumentName);
        add("--" + argumentName);
        add(argumentValue);
        return this;
    }

    /**
     * add a positional argument to this builder
     */
    public ArgumentsBuilder addPositionalArgument(final String argumentName) {
        Utils.nonNull(argumentName);
        add(argumentName);
        return this;
    }

    /**
     * Add any object's string representation to the arguments list
     */
    public ArgumentsBuilder add(Object arg) {
        this.args.add(arg.toString());
        return this;
    }

    /**
     * @return the arguments as List<String>
     */
    public List<String> getArgsList(){
        return this.args;
    }

    /**
     * @return the arguments as String[]
     */
    public String[] getArgsArray(){
        return this.args.toArray(new String[this.args.size()]);
    }


    /**
     * @return the arguments as a single String
     */
    public String getString() {
        return String.join(" ", args);
    }

}
