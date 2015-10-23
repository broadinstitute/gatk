package org.broadinstitute.hellbender.utils.mcmc;

import java.util.List;

/**
 * Basic implementation of {@link AbstractParameterizedState} with no convenience fields or getters.  See
 * example of use in GibbsSamplerSingleGaussianUnitTest.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ParameterizedState extends AbstractParameterizedState {
    /**
     * Subclasses of {@link AbstractParameterizedState} must implement the
     * {@link AbstractParameterizedState#copy(Class)} method.  This boilerplate implementation simply calls the
     * copy constructor {@link ParameterizedState#ParameterizedState(ParameterizedState)}.
     */
    @Override
    protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
        return stateClass.cast(new ParameterizedState(this));
    }

    /**
     * Copy constructor.
     */
    public ParameterizedState(final ParameterizedState state) {
        super(state);
    }

    /**
     * Constructor for creating a ParameterizedState given a List of {@link Parameter} objects with values of mixed type.
     * @param parameters    List of Parameters
     */
    public ParameterizedState(final List<Parameter<?>> parameters) {
        super(parameters);
    }
}
