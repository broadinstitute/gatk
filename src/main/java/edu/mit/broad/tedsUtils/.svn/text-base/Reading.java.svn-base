/*
 * $Id: Reading.java 44052 2007-09-04 21:12:29Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * A Sanger reading.
 * Has methods for reading and writing SCF-formatted files.
 * Note that the collections returned are the "live" data that are part of the Reading's state,
 * not independent copies.  Therefore, you can modify the Reading by modifying the collections.
 * Helpful, if what you want to do is to read, modify, and write a Reading.  Dangerous, otherwise.
 * Want to change a base call?  Now you can.
 *
 * @author tsharpe
 * @version $Revision: 44052 $
 */
public class Reading
    implements Cloneable
{
    /**
     * Make one.  Any of the args can be null.
     */
    public Reading( List<Sample> samples, List<Call> calls, List<Comment> comments, byte[] privateData )
    {
        mSamples = samples;
        mCalls = calls;
        mComments = comments;
        mPrivateData = privateData;
    }

    /**
     * Get the called sequence.
     */
    public String getSequence()
    {
        String result = null;
        if ( mCalls != null )
        {
            StringBuilder sb = new StringBuilder(mCalls.size());
            for ( Call call : mCalls )
            {
                sb.append(call.getCall());
            }
            result = sb.toString();
        }
        return result;
    }

    /**
     * Get the quality scores.
     * You own the byte array that's returned -- it's not part of the Reading's state.
     */
    public byte[] getQuality()
    {
        byte[] result = new byte[getNCalls()];
        for ( int idx = 0; idx < result.length; ++idx )
        {
            result[idx] = (byte)mCalls.get(idx).getQuality();
        }
        return result;
    }

    /**
     * The number of samples available.
     */
    public int getNSamples()
    {
        return mSamples == null ? 0 : mSamples.size();
    }

    /**
     * The intensity data (or null, if there is none).
     */
    public List<Sample> getSamples()
    {
        return mSamples;
    }

    /**
     * Get the samples as an array.
     */
    public Sample[] getSampleArray()
    {
        Sample[] result = EMPTY_SAMPLE_ARRAY;
        int nnn = getNSamples();
        if ( nnn > 0 )
        {
            result = mSamples.toArray(new Sample[nnn]);
        }
        return result;
    }

    /**
     * Set the intensity data.
     */
    public void setSamples( List<Sample> samples )
    {
        mSamples = samples;
    }

    /**
     * The number of calls available.
     */
    public int getNCalls()
    {
        return mCalls == null ? 0 : mCalls.size();
    }

    /**
     * The base-call data (or null, if there is none).
     */
    public List<Call> getCalls()
    {
        return mCalls;
    }

    /**
     * Set the calls.
     */
    public void setCalls( List<Call> calls )
    {
        mCalls = calls;
    }

    /**
     * Get the index of the first call that has a sample index greater than or equal to the
     * specified sample index.  Returns the total number of calls if there is no such call.
     */
    public int getCallIndex( int sampleIdx )
    {
        int nCalls = mCalls.size();
        for ( int idx = 0; idx < nCalls; ++idx )
        {
            if ( mCalls.get(idx).getSampleIndex() >= sampleIdx )
            {
                nCalls = idx;
            }
        }
        return nCalls;
    }

    /**
     * The number of comments available.
     */
    public int getNComments()
    {
        return mComments == null ? 0 : mComments.size();
    }

    /**
     * The <field-ID,value> pairs meant to comment on the data (or null, if there are none).
     * There doesn't seem to be any requirement in the SCF spec that the comment field-IDs be
     * unique, so this method doesn't return a map.
     */
    public List<Comment> getComments()
    {
        return mComments;
    }

    /**
     * Add a comment.
     */
    public void addComment( String key, String value )
    {
        if ( mComments == null )
        {
            mComments = new ArrayList<Comment>();
        }
        mComments.add(new Comment(key,value));
    }

    /**
     * Set the comments.
     */
    public void setComments( List<Comment> comments )
    {
        mComments = comments;
    }

    /**
     * The length of the private data.
     */
    public int getPrivateDataLength()
    {
        return mPrivateData == null ? 0 : mPrivateData.length;
    }

    /**
     * Some extra info that the trace writer stuck in the file (or null, if there is none).
     */
    public byte[] getPrivateData()
    {
        return mPrivateData;
    }

    /**
     * Set the extra info.
     * @param privateData The extra info, or null.
     */
    public void setPrivateData( byte[] privateData )
    {
        mPrivateData = privateData;
    }

    /**
     * Copy a Reading.
     * The lists are cloned, but the objects they list are supposed to be immutable,
     * and they are not cloned.
     */
    @Override
    public Reading clone()
    {
        try
        {
            Reading that = (Reading)super.clone();
            if ( mSamples != null )
                that.mSamples = new ArrayList<Sample>(mSamples);
            if ( mCalls != null )
                that.mCalls = new ArrayList<Call>(mCalls);
            if ( mComments != null )
                that.mComments = new ArrayList<Comment>(mComments);
            if ( mPrivateData != null )
                that.mPrivateData = mPrivateData.clone();
            return that;
        }
        catch ( CloneNotSupportedException e )
        {
            throw new IllegalStateException(e);
        }
    }

    /**
     * Write a trace in SCF format.
     * The output stream is closed by this method.
     * @param raw Where to write it.
     */
    public void writeSCF( OutputStream raw )
        throws IOException
    {
        DataOutputStream os = new DataOutputStream(new BufferedOutputStream(raw));
        try
        {
            // write header

            int offset = 128;
            os.writeInt(SCF_MAGIC_NUMBER);

            int nSamples = getNSamples();
            os.writeInt(nSamples);
            os.writeInt(offset);

            int bytesPerSample = 2;
            if ( nSamples > 0 )
            {
                if ( maxSampleIntensity() < 255.5 )
                {
                    bytesPerSample = 1;
                }
                offset += 4*bytesPerSample*nSamples;
            }

            int nCalls = getNCalls();
            os.writeInt(nCalls);
            os.writeInt(0);
            os.writeInt(0);
            os.writeInt(offset);
            offset += 12*nCalls;

            byte[] commentsBuf = null;
            if ( getNComments() == 0 )
            {
                os.writeInt(0);
                os.writeInt(offset);
            }
            else
            {
                commentsBuf = encodeComments();
                os.writeInt(commentsBuf.length);
                os.writeInt(offset);
                offset += commentsBuf.length;
            }

            os.writeInt(SCF_VERSION_3_1);
            os.writeInt(bytesPerSample);
            os.writeInt(SCF_CODE_SET); // code set

            int nPrivateBytes = getPrivateDataLength();
            os.writeInt(nPrivateBytes);
            os.writeInt(offset);

            for ( int idx = 0; idx < 18; ++idx )
                os.writeInt(0);

            // write the samples

            if ( nSamples > 0 )
            {
                int[] samples = doubleDelta();
                if ( bytesPerSample == 1 )
                    for ( int sss : samples )
                        os.writeByte(sss);
                else
                    for ( int sss : samples )
                        os.writeShort(sss);
            }

            // write the calls

            if ( nCalls > 0 )
            {
                for ( Call call : mCalls )
                    os.writeInt(call.getSampleIndex());
                for ( Call call : mCalls )
                    os.writeByte(call.getQuality(Base.A));
                for ( Call call : mCalls )
                    os.writeByte(call.getQuality(Base.C));
                for ( Call call : mCalls )
                    os.writeByte(call.getQuality(Base.G));
                for ( Call call : mCalls )
                    os.writeByte(call.getQuality(Base.T));
                for ( Call call : mCalls )
                    os.writeByte(call.getCall());
                for ( Call call : mCalls )
                    os.writeByte(call.getSubstitutionScore());
                for ( Call call : mCalls )
                    os.writeByte(call.getInsertionScore());
                for ( Call call : mCalls )
                    os.writeByte(call.getDeletionScore());
            }

            // write the comments

            if ( commentsBuf != null )
                os.write(commentsBuf);

            // write the private data

            if ( nPrivateBytes > 0 )
                os.write(mPrivateData);
        }
        finally
        {
            os.close();
        }
    }

    /**
     * Create a trace from a file name.
     */
    public static Reading readSCF( String fileName )
        throws IOException
    {
        return readSCF(new File(fileName));
    }

    /**
     * Create a trace from an SCF file.
     */
    public static Reading readSCF( File file )
        throws IOException
    {
        return readSCF(getInputStream(file));
    }

    /**
     * Create a trace from an SCF-formatted input stream.
     * The stream is closed by this method.
     * @param raw The SCF-formatted input stream.
     * @return A Reading.
     */
    public static Reading readSCF( InputStream raw )
        throws IOException
    {
        DataInputStream is = new DataInputStream(new BufferedInputStream(raw));
        try
        {
            int magic = is.readInt();
            if ( magic != SCF_MAGIC_NUMBER )
            {
                throw new IOException("Not an SCF file: wrong magic number");
            }

            int nSamples = is.readInt();
            int sampleAreaIdx = is.readInt();

            int nBases = is.readInt();
            is.skipBytes(8); // skip over obsolete left and right clips
            int baseAreaIdx = is.readInt();
            int baseAreaEnd = baseAreaIdx + 12*nBases;

            int nCommentBytes = is.readInt();
            int commentAreaIdx = is.readInt();
            int commentAreaEnd = commentAreaIdx + nCommentBytes;

            int version = is.readInt();
            int nBytesPerSample = is.readInt();
            is.skipBytes(4); // skip over codeset

            if ( version < SCF_VERSION_2 )
            {
                nBytesPerSample = 1;
            }
            int sampleAreaEnd = sampleAreaIdx + 4 * nBytesPerSample * nSamples;

            int nPrivateBytes = is.readInt();
            int privateAreaIdx = is.readInt();
            int privateAreaEnd = privateAreaIdx + nPrivateBytes;

            is.skipBytes(4*18); // header spare bytes

            int len = Math.max(privateAreaEnd,Math.max(commentAreaEnd,Math.max(baseAreaEnd,sampleAreaEnd)));
            byte[] buffer = new byte[len];
            is.readFully(buffer,128,len-128);

            List<Sample> samples = null;
            if ( nSamples > 0 )
            {
                if ( version < SCF_VERSION_3 )
                {
                    samples = readSamplesV2(nSamples,nBytesPerSample,buffer,sampleAreaIdx);
                }
                else
                {
                    samples = readSamplesV3(nSamples,nBytesPerSample,buffer,sampleAreaIdx);
                }
            }

            List<Call> calls = null;
            if ( nBases > 0 )
            {
                if ( version < SCF_VERSION_3 )
                {
                    calls = readCallsV2(nBases,buffer,baseAreaIdx);
                }
                else
                {
                    calls = readCallsV3(nBases,buffer,baseAreaIdx);
                }
            }

            List<Comment> comments = null;
            if ( nCommentBytes > 0 )
            {
                comments = readComments(nCommentBytes,buffer,commentAreaIdx);
            }

            byte[] privateData = null;
            if ( nPrivateBytes > 0 )
            {
                privateData = new byte[nPrivateBytes];
                System.arraycopy(buffer,privateAreaIdx,privateData,0,nPrivateBytes);
            }

            return new Reading(samples,calls,comments,privateData);
        }
        finally
        {
            is.close();
        }
    }


    public static InputStream getInputStream( File scfFile )
        throws IOException
    {
        InputStream result = new FileInputStream(scfFile);
        if ( scfFile.getName().toLowerCase().endsWith(".gz") )
        {
            result = new GZIPInputStream(result);
        }
        return result;
    }

    public static String getTraceName( File file )
    {
        return getTraceName(file.getName());
    }

    public static String getTraceName( String fileName )
    {
        int idx = fileName.indexOf(".scf");
        if ( idx == -1 )
        {
            idx = fileName.length();
        }
        return fileName.substring(0,idx);
    }

    private int[] doubleDelta()
    {
        int[] encodedSamples = new int[4*mSamples.size()];
        int offset = 0;
        int z = 0;
        int p1 = 0;
        int p2 = 0;
        for ( Sample sample : mSamples )
        {
            z = 2 * p1 - p2;
            p2 = p1;
            p1 = Math.min(MAX_INTENSITY_VALUE,Math.max(0,Math.round((float)sample.getIntensity(Base.A))));
            encodedSamples[offset++] = p1 - z;
        }
        z = p1 = p2 = 0;
        for ( Sample sample : mSamples )
        {
            z = 2 * p1 - p2;
            p2 = p1;
            p1 = Math.min(MAX_INTENSITY_VALUE,Math.max(0,Math.round((float)sample.getIntensity(Base.C))));
            encodedSamples[offset++] = p1 - z;
        }
        z = p1 = p2 = 0;
        for ( Sample sample : mSamples )
        {
            z = 2 * p1 - p2;
            p2 = p1;
            p1 = Math.min(MAX_INTENSITY_VALUE,Math.max(0,Math.round((float)sample.getIntensity(Base.G))));
            encodedSamples[offset++] = p1 - z;
        }
        z = p1 = p2 = 0;
        for ( Sample sample : mSamples )
        {
            z = 2 * p1 - p2;
            p2 = p1;
            p1 = Math.min(MAX_INTENSITY_VALUE,Math.max(0,Math.round((float)sample.getIntensity(Base.T))));
            encodedSamples[offset++] = p1 - z;
        }
        return encodedSamples;
    }

    private byte[] encodeComments()
        throws IOException
    {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        Iterator<Comment> itr = mComments.iterator();
        while ( itr.hasNext() )
        {
            Comment entry = itr.next();
            if ( entry.getKey() != null )
                os.write(entry.getKey().getBytes());
            os.write('=');
            if ( entry.getValue() != null )
                os.write(entry.getValue().getBytes());
            os.write(itr.hasNext()?'\n':0);
        }
        return os.toByteArray();
    }

    private double maxSampleIntensity()
    {
        double result = 0;
        for ( Base channel : Base.values() )
        {
            for ( Sample sample : mSamples )
            {
                double vvv = sample.getIntensity(channel);
                if ( vvv > result )
                {
                    result = vvv;
                }
            }
        }
        return result;
    }

    private static List<Sample> readSamplesV2( int nSamples, int nBytesPerSample, byte[] buffer, int offset )
        throws IOException
    {
        List<Sample> samples = new ArrayList<Sample>(nSamples);
        if ( nBytesPerSample == 1 )
        {
            while ( nSamples-- > 0 )
            {
                int aVal = buffer[offset++] & 0xFF;
                int cVal = buffer[offset++] & 0xFF;
                int gVal = buffer[offset++] & 0xFF;
                int tVal = buffer[offset++] & 0xFF;
                samples.add(new Sample(aVal,cVal,gVal,tVal));
            }
        }
        else if ( nBytesPerSample == 2 )
        {
            while ( nSamples-- > 0 )
            {
                int aVal = readShort(buffer,offset);
                int cVal = readShort(buffer,offset += 2);
                int gVal = readShort(buffer,offset += 2);
                int tVal = readShort(buffer,offset += 2);
                offset += 2;
                samples.add(new Sample(aVal,cVal,gVal,tVal));
            }
        }
        else
        {
            throw new IOException("Bytes per sample is " + nBytesPerSample);
        }
        return samples;
    }

    private static List<Sample> readSamplesV3( int nSamples, int nBytesPerSample, byte[] buffer, int offset )
        throws IOException
    {
        List<Sample> samples = new ArrayList<Sample>(nSamples);
        int a0 = 0;
        int a1 = 0;
        int c0 = 0;
        int c1 = 0;
        int g0 = 0;
        int g1 = 0;
        int t0 = 0;
        int t1 = 0;
        if ( nBytesPerSample == 1 )
        {
            int end = offset + nSamples;
            while ( offset < end )
            {
                a0 += buffer[offset];
                a1 += a0;
                c0 += buffer[offset+nSamples];
                c1 += c0;
                g0 += buffer[offset+2*nSamples];
                g1 += g0;
                t0 += buffer[offset+3*nSamples];
                t1 += t0;
                samples.add(new Sample(a1&0xFF,c1&0xFF,g1&0xFF,t1&0xFF));
                offset += 1;
            }
        }
        else if ( nBytesPerSample == 2 )
        {
            int end = offset + 2*nSamples;
            while ( offset < end )
            {
                a0 += readShort(buffer,offset);
                a1 += a0;
                c0 += readShort(buffer,offset+2*nSamples);
                c1 += c0;
                g0 += readShort(buffer,offset+4*nSamples);
                g1 += g0;
                t0 += readShort(buffer,offset+6*nSamples);
                t1 += t0;
                samples.add(new Sample(a1&0xFFFF,c1&0xFFFF,g1&0xFFFF,t1&0xFFFF));
                offset += 2;
            }
        }
        else
        {
            throw new IOException("Bytes per sample is " + nBytesPerSample);
        }
        return samples;
    }

    private static List<Call> readCallsV2( int nCalls, byte[] buffer, int offset )
    {
        List<Call> calls = new ArrayList<Call>(nCalls);
        while ( nCalls-- > 0 )
        {
            int peakIndex = readInt(buffer,offset);
            int qA = buffer[offset += 4] & 0xFF;
            int qC = buffer[++offset] & 0xFF;
            int qG = buffer[++offset] & 0xFF;
            int qT = buffer[++offset] & 0xFF;
            char call = (char)(buffer[++offset] & 0xFF);
            offset += 4;
            calls.add( new Call(call,peakIndex,qA,qC,qG,qT,0,0,0) );
        }
        return calls;
    }

    private static List<Call> readCallsV3( int nCalls, byte[] buffer, int offset )
    {
        List<Call> calls = new ArrayList<Call>(nCalls);
        for ( int iii = 0; iii < nCalls; ++iii )
        {
            int peakIndex = readInt(buffer,offset+4*iii);
            int idx = offset + 4*nCalls + iii;
            int qA = buffer[idx] & 0xFF;
            int qC = buffer[idx += nCalls] & 0xFF;
            int qG = buffer[idx += nCalls] & 0xFF;
            int qT = buffer[idx += nCalls] & 0xFF;
            char call = (char)(buffer[idx += nCalls] & 0xFF);
            int qSub = buffer[idx += nCalls] & 0xFF;
            int qIns = buffer[idx += nCalls] & 0xFF;
            int qDel = buffer[idx += nCalls] & 0xFF;
            calls.add( new Call(call,peakIndex,qA,qC,qG,qT,qSub,qIns,qDel) );
        }
        return calls;
    }

    private static List<Comment> readComments( int nCommentBytes, byte[] buffer, int offset )
        throws IOException
    {
        List<Comment> comments = new ArrayList<Comment>();

        BufferedReader rdr = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(buffer,offset,nCommentBytes)));
        String line;
        while ( (line = rdr.readLine()) != null )
        {
            int eqIdx = line.indexOf('=');
            if ( eqIdx == -1 )
            {
                if ( line.length() != 0 && !line.equals("\000") )
                {
                    System.err.println("Badly formatted comment ignored: '" + line + "'.");
                }
                continue;
            }
            int endIdx = line.indexOf(0,eqIdx);
            if ( endIdx != -1 )
            {
                comments.add(new Comment(line.substring(0,eqIdx),line.substring(eqIdx+1,endIdx)));
                break;
            }
            comments.add(new Comment(line.substring(0,eqIdx),line.substring(eqIdx+1)));
        }
        rdr.close();

        return comments;
    }

    private static int readShort( byte[] buffer, int idx )
    {
        return ((buffer[idx] & 0xFF) << 8) | (buffer[idx+1] & 0xFF);
    }

    private static int readInt( byte[] buffer, int idx )
    {
        return (buffer[idx] << 24) | ((buffer[idx+1] & 0xFF) << 16) | ((buffer[idx+2] & 0xFF) << 8) | (buffer[idx+3] & 0xFF);
    }

    private List<Sample> mSamples;
    private List<Call> mCalls;
    private List<Comment> mComments;
    private byte[] mPrivateData;

    private static final int SCF_MAGIC_NUMBER = 0x2e736366;
    private static final int SCF_VERSION_3_1 = 0x332e3130;
    private static final int SCF_VERSION_3 = 0x332e3030;
    private static final int SCF_VERSION_2 = 0x322e3030;
    private static final int SCF_CODE_SET = 2;
    private static final int MAX_INTENSITY_VALUE = 0x7FFF;
    private static final Sample[] EMPTY_SAMPLE_ARRAY = new Sample[0];

    /**
     * Four channels of intensity data.
     * This class is, and must remain immutable, since it's shared between clones of Readings.
     *
     * @author tsharpe
     * @version $Revision: 44052 $
     */
    public static class Sample
    {
        /**
         * Make one.
         */
        public Sample( double aVal, double cVal, double gVal, double tVal )
        {
            mVals = new double[5];
            mVals[0] = aVal;
            mVals[1] = cVal;
            mVals[2] = gVal;
            mVals[3] = tVal;
            mVals[4] = aVal + cVal + gVal + tVal;
            test();
        }

        /**
         * Make one by scaling another.
         */
        public Sample( double scale, Sample sample )
        {
            mVals = new double[5];
            for ( int iii = 0; iii < 5; ++iii )
            {
                mVals[iii] = scale * sample.mVals[iii];
            }
            test();
        }

        /**
         * Returns a new sample, scaled by the argument. 
         */
        public Sample scale( double scale )
        {
            Sample result = this;
            if ( this != ZERO )
            {
                result = new Sample(scale,this);
            }
            return result;
        }

        /**
         * Return a new sample with a constant added to each intensity.
         */
        public Sample bias( double bias )
        {
            Sample result = new Sample();

            double[] vals = result.mVals;
            for ( int iii = 0; iii < 4; ++iii )
            {
                vals[iii] = mVals[iii] + bias;
            }
            vals[4] = mVals[4] + 4*bias;

            result.test();
            return result;
        }

        /**
         * Get the intensity data for one of the channels.
         */
        public double getIntensity( Base channel )
        {
            return mVals[channel.ordinal()];
        }

        /**
         * Get the sum of the intensities over all channels. 
         */
        public double getTotalIntensity()
        {
            return mVals[4];
        }

        /**
         * Get a number that represents the degree of heterozygosity of the sample.
         * Ranges from 0 (which means that all of the signal is in one channel) to 1
         * (which means that all of the signal is equally distributed between two channels).
         * Special case: an all-zero sample has a heterozygosity of 0.5.
         */
        public double heterozygosity()
        {
            double result = 0.5;
            if ( mVals[4] != 0 )
            {
                double best = mVals[0];
                double nextBest = mVals[1];
                if ( nextBest > best )
                {
                    best = nextBest;
                    nextBest = mVals[0];
                }

                if ( mVals[2] > best )
                {
                    nextBest = best;
                    best = mVals[2];
                }
                else if ( mVals[2] > nextBest )
                {
                    nextBest = mVals[2];
                }

                if ( mVals[3] > best )
                {
                    nextBest = best;
                    best = mVals[3];
                }
                else if ( mVals[3] > nextBest )
                {
                    nextBest = mVals[3];
                }

                result = 4. * best * nextBest / mVals[4] / mVals[4];
            }
            return result;
        }

        /**
         * Maximum intensity in any channel.
         */
        public double getMaxIntensity()
        {
            return Math.max(Math.max(mVals[0],mVals[1]),Math.max(mVals[2],mVals[3]));
        }

        /**
         * Returns the intensities of the four channels delimited by tabs.
         */
        @Override
        public String toString()
        {
            return mVals[0] + "\t" + mVals[1] + "\t" + mVals[2] + "\t" + mVals[3];
        }

        /**
         * Make a new sample by summing two others.
         */
        public static Sample sum( Sample sample1, Sample sample2 )
        {
            Sample result;
            if ( sample1 == ZERO )
            {
                result = sample2;
            }
            else if ( sample2 == ZERO )
            {
                result = sample1;
            }
            else
            {
                result = new Sample();
                for ( int iii = 0; iii < 5; ++iii )
                {
                    result.mVals[iii] = sample1.mVals[iii] + sample2.mVals[iii];
                }
            }
            result.test();
            return result;
        }

        /**
         * Make a new sample from the difference of two others.
         */
        public static Sample diff( Sample sample1, Sample sample2 )
        {
            Sample result;
            if( sample2 == ZERO )
            {
                result = sample1;
            }
            else
            {
                result = new Sample();
                for ( int iii = 0; iii < 5; ++iii )
                {
                    result.mVals[iii] = sample1.mVals[iii] - sample2.mVals[iii];
                }
            }
            result.test();
            return result;
        }

        /**
         * Make a new sample from the product of two others.
         */
        public static Sample product( Sample sample1, Sample sample2 )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                tot += (result.mVals[iii] = sample1.mVals[iii] * sample2.mVals[iii]);
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        /**
         * Make a new sample from the ratio of two others.
         */
        public static Sample ratio( Sample sample1, Sample sample2 )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                tot += (result.mVals[iii] = sample1.mVals[iii] / sample2.mVals[iii]);
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        /**
         * Make a new sample from the minimum of two others.
         */
        public static Sample min( Sample sample1, Sample sample2 )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                tot += (result.mVals[iii] = Math.min(sample1.mVals[iii],sample2.mVals[iii]));
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        /**
         * Make a new sample from the maximum of two others.
         */
        public static Sample max( Sample sample1, Sample sample2 )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                tot += (result.mVals[iii] = Math.max(sample1.mVals[iii],sample2.mVals[iii]));
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        /**
         * Make a new sample from another, capping the intensities at a specified value.
         */
        public static Sample clip( Sample sample, double maxVal )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                tot += (result.mVals[iii] = Math.min(sample.mVals[iii],maxVal));
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        /**
         * Make a new sample from another, but set the minimum intensities to a specified value.
         */
        public static Sample rectify( Sample sample, double minVal )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                tot += (result.mVals[iii] = Math.max(sample.mVals[iii],minVal));
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        public static Sample noNaN( Sample replacement, Sample toTest )
        {
            Sample result = new Sample();
            double tot = 0.;
            for ( int iii = 0; iii < 4; ++iii )
            {
                double val = toTest.mVals[iii];
                if ( Double.isNaN(val) || Double.isInfinite(val) )
                {
                    val = replacement.mVals[iii];
                }
                tot += (result.mVals[iii] = val);
            }
            result.mVals[4] = tot;
            result.test();
            return result;
        }

        public static List<Sample> biasToZero( List<Sample> samples )
        {
            int nnn = samples.size();
            List<Sample> result = new ArrayList<Sample>(nnn);

            double minVal = Double.MAX_VALUE;
            for ( Sample sample : samples )
            {
                minVal = Math.min(minVal,Math.min(Math.min(sample.mVals[0],sample.mVals[1]),Math.min(sample.mVals[2],sample.mVals[3])));
            }
            minVal = -minVal;
            for ( Sample sample : samples )
            {
                result.add(sample.bias(minVal));
            }
            return result;
        }

        private void test()
        {
            assert(mVals[4] == mVals[0] + mVals[1] + mVals[2] + mVals[3]);
        }

        private Sample()
        {
            mVals = new double[5];
        }

        private double[] mVals;
        public static final Sample ZERO = new Sample(0,0,0,0);
    }

    /**
     * A base call.
     * This class is, and must remain immutable, since it's shared between clones of Readings.
     *
     * @author tsharpe
     * @version $Revision: 44052 $
     */
    public static class Call
    {
        public Call( char call, int sampleIndex, int qA, int qC, int qG, int qT, int qSub, int qIns, int qDel )
        {
            mCall = call;
            mSampleIndex = sampleIndex;
            mQuals = new int[7];
            mQuals[0] = qA;
            mQuals[1] = qC;
            mQuals[2] = qG;
            mQuals[3] = qT;
            mQuals[4] = qSub;
            mQuals[5] = qIns;
            mQuals[6] = qDel;
        }

        /**
         * The code called.
         * @return The call.
         */
        public char getCall()
        {
            return mCall;
        }

        /**
         * The index of the sample where the call was made.
         */
        public int getSampleIndex()
        {
            return mSampleIndex;
        }

        /**
         * The quality score associated with the call.
         */
        public int getQuality()
        {
            return getQuality(Base.valueOf(mCall));
        }

        /**
         * Makes a clone with an updated quality value.
         */
        public Call setQuality( int quality )
        {
            Call clone = new Call(mCall,mSampleIndex,mQuals[0],mQuals[1],mQuals[2],mQuals[3],mQuals[4],mQuals[5],mQuals[6]);
            Base base = Base.valueOf(mCall);
            if ( base == null )
            {
                base = Base.A;
            }
            clone.mQuals[base.ordinal()] = quality;
            return clone;
        }

        /**
         * The quality score associated with a specified channel.
         * Returns the max quality for a null channel.
         */
        public int getQuality( Base channel )
        {
            int result;
            if ( channel == null )
            {
                result = Math.max(Math.max(mQuals[0],mQuals[1]),Math.max(mQuals[2],mQuals[3]));
            }
            else
            {
                result = mQuals[channel.ordinal()];
            }
            return result;
        }

        /**
         * New in version 3.10 -- probability that this call is a substitution for another base (whatever that means).
         */
        public int getSubstitutionScore()
        {
            return mQuals[4];
        }

        /**
         * New in version 3.10 -- probability that this call is an overcall.
         */
        public int getInsertionScore()
        {
            return mQuals[5];
        }

        /**
         * New in version 3.10 -- probability that this call is an undercall.
         */
        public int getDeletionScore()
        {
            return mQuals[6];
        }

        /**
         * Returns the base called, and its quality.
         */
        @Override
        public String toString()
        {
            return mCall + "(" + getQuality() + ")";
        }

        private char mCall;
        private int mSampleIndex;
        private int[] mQuals;
    }

    /**
     * A comment on the data or on how it was processed.
     * Almost a java.util.Map.Entry, but you can't change the value.
     * This class is, and must remain immutable, since it's shared between clones of Readings.
     *
     * @author tsharpe
     * @version $Revision: 44052 $
     */
    public static class Comment
    {
        /**
         * Make one.
         */
        public Comment( String key, String value )
        {
            mKey = key;
            mValue = value;
        }

        /**
         * Get the field ID.
         */
        public String getKey()
        {
            return mKey;
        }

        /**
         * Get the commentary for the associated field ID.
         */
        public String getValue()
        {
            return mValue;
        }

        /**
         * Same definition as for Map.Entry, q.v.
         */
        @Override
        public boolean equals( Object obj )
        {
            boolean result = false;
            if ( this == obj )
            {
                result = true;
            }
            else if ( obj instanceof Comment )
            {
                Comment that = (Comment)obj;
                if ( (mKey == that.getKey() || mKey != null && mKey.equals(that.getKey())) &&
                     (mValue == that.getValue() || mValue != null && mValue.equals(that.mValue)) )
                {
                        result = true;
                }
            }
            return result;
        }

        /**
         * Same definition as for Map.Entry, q.v.
         */
        @Override
        public int hashCode()
        {
            return (mKey == null ? 0 : mKey.hashCode()) ^ (mValue == null ? 0 : mValue.hashCode());
        }

        /**
         * Returns fieldID = value.
         */
        @Override
        public String toString()
        {
            return mKey + " = " + mValue;
        }

        private String mKey;
        private String mValue;
    }
}
