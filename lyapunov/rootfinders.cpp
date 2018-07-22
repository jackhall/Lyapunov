#ifndef ROOTFINDERS_CPP
#define ROOTFINDERS_CPP

#include "common.cpp"
#include "steppers.cpp"


namespace lyapunov {
    namespace bp = boost::python;

    struct Interval {
        /* A utility to handle the contracting interval of a bisection rootfinder.
         */
        double lower, upper;
        double length() const { return upper - lower; }
        double midpoint() const { return (upper + lower) / 2.0; }
    };

    class RootFinder : public ODEIterator {
    /* Abstract base class for any rootfinders.
     */
    protected:
        Stepper _stepper;

        virtual void _step_across(bp::object new_state) = 0;

    public:
        RootFinder() = delete;
        RootFinder(Stepper stepper, bp::object events);
        RootFinder(bp::object stepper, bp::object events) {
            // Wrap python stepper with a subclass of Stepper.
        }

        // RootFinders do not directly manage a system.
        // All system accesses go through Steppers.
        bp::object get_system() const {
            return stepper.get_system();
        }
        void set_system(bp::object new_system) {
            stepper.set_system(new_system);
        }

        double get_step_size() const {
            return stepper.get_step_size();
        }

        void step_across(bp::object new_state = bp::object()) {
            _step_across(new_state);
        }

        virtual bp::tuple next() = 0;
    };

}

#endif
