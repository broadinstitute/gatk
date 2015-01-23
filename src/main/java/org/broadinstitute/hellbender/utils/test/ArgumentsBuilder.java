package org.broadinstitute.hellbender.utils.test;

import org.apache.commons.lang3.StringUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Builder for command line argument lists
 * It will convert old style "Argument=Value" into new style "--Argument value" strings
 * Use this only in test code.
 */
public class ArgumentsBuilder {
    final private List<String> args= new ArrayList<>();

    public ArgumentsBuilder(){
    }

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
     * Strings are processed specially, they are reformated to match the new unix style arguments
     * @param arg A string representing one or more arguments
     * @return
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
     * Add any object's string representation to the arguments list
     * @param arg
     * @return
     */
    public ArgumentsBuilder add(Object arg) {
        this.args.add(arg.toString());
        return this;
    }

    /**
     * get the arguments as List<String>
     * @return
     */
    public List<String> getArgsList(){
        return this.args;
    }

    /**
     * get the arguments as String[]
     * @return
     */
    public String[] getArgsArray(){
        return this.args.toArray(new String[this.args.size()]);
    }



}
