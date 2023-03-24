/**
 * $Id: BaseCall.java 72515 2008-09-12 17:43:22Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import edu.mit.broad.tedsUtils.Base;

/**
 * A called base (including ambiguity codes).
 * I've invented the symbol X to represent no base (as opposed to N which means any base).
 * This class enumerates all the possible sets of the enum Base.
 *
 * @author tsharpe
 * @version $Revision$
 */
public enum BaseCall
{
    X(EnumSet.noneOf(Base.class)),
    A(EnumSet.of(Base.A)),
    C(EnumSet.of(Base.C)),
    M(EnumSet.of(Base.A,Base.C)),
    G(EnumSet.of(Base.G)),
    R(EnumSet.of(Base.A,Base.G)),
    S(EnumSet.of(Base.C,Base.G)),
    V(EnumSet.of(Base.A,Base.C,Base.G)),
    T(EnumSet.of(Base.T)),
    W(EnumSet.of(Base.A,Base.T)),
    Y(EnumSet.of(Base.C,Base.T)),
    H(EnumSet.of(Base.A,Base.C,Base.T)),
    K(EnumSet.of(Base.G,Base.T)),
    D(EnumSet.of(Base.A,Base.G,Base.T)),
    B(EnumSet.of(Base.C,Base.G,Base.T)),
    N(EnumSet.allOf(Base.class));

    /**
     * Returns the character code for the BaseCall.
     * (I.e., the first, and only, letter of its name.)
     */
    public char code()
    {
        return name().charAt(0);
    }

    /**
     * Returns the complement of the call.
     * E.g., T is returned for A, and W is returned for W (since the complement of A or T is T or A).
     */
    public BaseCall complement()
    {
        if ( mComplement == null )
        {
            EnumSet<Base> bases = EnumSet.noneOf(Base.class);
            for ( Base base : mBases )
            {
                bases.add(base.complement());
            }
            mComplement = valueOf(bases);
        }
        return mComplement;
    }

    /**
     * Returns a regular expression for the call.
     * E.g., "[AT]" is returned for W.
     */
    public String regExp()
    {
        String result = name();
        if ( mNBases == 0 )
        {
            result = "";
        }
        else if ( mNBases > 1 )
        {
            StringBuilder sb = new StringBuilder(mNBases+2);
            sb.append('[');
            for ( Base base : mBases )
            {
                sb.append(base.name());
            }
            sb.append(']');
            result = sb.toString();
        }
        return result;
    }

    /**
     * The size of the bases() EnumSet.
     */
    public int getNBases()
    {
        return mNBases;
    }

    /**
     * The set of Bases represented by this (possibly ambiguous) BaseCall.
     * You can monkey with it -- it's a copy of the one true set.
     */
    public EnumSet<Base> bases()
    {
        return EnumSet.copyOf(mBases);
    }

    /**
     * Returns the corresponding Base for unambiguous BaseCalls.
     * Returns null unless this BaseCall unambiguously represents a single Base.
     */
    public Base toBase()
    {
        return mBase;
    }
    
    /**
     * Convert a character into the ordinal of one of the enumerated BaseCalls.
     * Upper and lower case characters in the set XACMGRSVTWYHKDBN are converted
     * to the corresponding BaseCall's ordinal value.
     * U is converted to the ordinal for T.
     * All other characters are converted to the value NO_BASE.
     */
    public static int ordinalOf( char code )
    {
        return CHAR_TO_VAL[code];
    }

    /**
     * Translates an unambiguous Base into a BaseCall.
     * The null Base gets transformed into an X.
     */
    public static BaseCall valueOf( Base base )
    {
        switch ( base )
        {
        case A:
            return A;
        case C:
            return C;
        case G:
            return G;
        case T:
            return T;
        default:
            return X;
        }
    }

    /**
     * Translates a character code into a BaseCall.
     * Null is returned for unrecognized codes.
     */
    public static BaseCall valueOf( char code )
    {
        BaseCall result = null;
        int idx = ordinalOf(code);
        if ( idx >= 0 )
        {
            result = values()[idx];
        }
        return result;
    }

    /**
     * Translates a set of Bases into the BaseCall representing that set.
     */
    public static BaseCall valueOf( Set<Base> bases )
    {
        return BASE_SET_TO_BASECALL.get(bases);
    }

    /**
     * Returns the BaseCall that has, as its set of bases, the union of the sets of bases associated with each argument.
     */
    public static BaseCall unionOf( BaseCall bc1, BaseCall bc2 )
    {
        EnumSet<Base> bases = EnumSet.copyOf(bc1.mBases);
        bases.addAll(bc2.mBases);
        return valueOf(bases);
    }

    /**
     * Returns the BaseCall that has, as its set of bases, the intersection of the sets of bases associated with each argument.
     */
    public static BaseCall intersectionOf( BaseCall bc1, BaseCall bc2 )
    {
        EnumSet<Base> bases = EnumSet.copyOf(bc1.mBases);
        bases.retainAll(bc2.mBases);
        return valueOf(bases);
    }

    BaseCall( EnumSet<Base> bases )
    {
        mBases = bases;
        mNBases = mBases.size();
        if ( mNBases == 1 )
        {
            mBase = mBases.iterator().next();
        }
    }

    private EnumSet<Base> mBases;
    private int mNBases; // mBases.size(), cached for performance.
    private Base mBase; // the corresponding Base, if there is one that's unique.
    private BaseCall mComplement;

    /**
     * Shorthand for BaseCall.values().length
     */
    public static final int N_BASECALLS = BaseCall.values().length;
    public static final byte NO_BASE = -1;

    private static final byte[] CHAR_TO_VAL;
    private static final Map<EnumSet<Base>,BaseCall> BASE_SET_TO_BASECALL;

    static
    {
        CHAR_TO_VAL = new byte[Character.MAX_VALUE+1];
        Arrays.fill(CHAR_TO_VAL,NO_BASE);
        BASE_SET_TO_BASECALL = new HashMap<EnumSet<Base>,BaseCall>(N_BASECALLS);
        for ( BaseCall call : values() )
        {
            char code = call.code();
            CHAR_TO_VAL[Character.toLowerCase(code)] = CHAR_TO_VAL[Character.toUpperCase(code)] = (byte)call.ordinal();
            BASE_SET_TO_BASECALL.put(call.mBases, call);
        }
        CHAR_TO_VAL['U'] = CHAR_TO_VAL['u'] = (byte)T.ordinal();
        assert(BASE_SET_TO_BASECALL.size() == N_BASECALLS);
    }
}
