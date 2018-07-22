/*
	Lyapunov: a library for integrating nonlinear dynamical systems
	Copyright (C) 2013  John Wendell Hall

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

	The author may be reached at jackhall@utexas.edu.
*/

#include <boost/numeric/odeint.hpp>

namespace lyapunov {
    namespace bp = boost::python;

    class Stepper : public ODEIterator {
    /* Abstract base class for any ODE solver object.
     */
    public:
        Stepper() = delete;
        Stepper(bp::object sys) = 0;

        virtual bool step(double step_size) = 0;
    };


    template<typename system_type, typename stepper_type>
    class multistep_iterator {
        /* Not to be used with rootfinding!
         */
    protected:
        stepper_type stepper;
        typedef stepper_type::state_type state_type;

        typedef std::pair<double, state_type> pair_type;
        pair_type current;

        system_type system;
		std::function<void(const state_type&, state_type&, double)> eval_system;

    public:
        multistep_iterator(system_type system)
            : stepper(), current(system.get()), system(system),
              eval_system([this](const state_type& x, state_type& dxdt, const double t) {
                  system.set(t, x);
                  dxdt = system.dxdt();
              }) {}

        void synchronize() {
            /* Moves the managed system to the current state.
             */
            system.set(current.first, current.second);
        }
        pair_type state() const { return current; }
        bool step(double step_size) {
            stepper.do_step(eval_system, current.second, current.first, step_size);
            current.first += step_size;
            return true;
        }
        void reset(state_type new_state) {
            current = std::make_pair(current.first, new_state);
            stepper.reset()
            synchronize();
        }
    };
    template<typename system_type, typename stepper_type>
    class fixed_step_iterator
      : public multistep_iterator<system_type, stepper_type> {
        /* An iterator to control the flow of a basic simulation. It simulates
         * systems defined as objects, such that state is set separately from
         * from derivative evaluation. Allowing the user to split up that
         * computation means that setting the state is potentially expensive,
         * so it shouldn't be done more often than necessary.
         */
        typedef multistep_iterator<system_type, stepper_type> base_type;
        typedef stepper_type::state_type state_type;
        typedef std::pair<double, state_type> pair_type;
        pair_type saved;

    public:
        fixed_step_iterator(system_type system)
            : base_type(system), saved(current) {}

        bool step(double step_size) {
            /* Take a step forward, and return whether that step was successful. Does
             * not synchronize after, because sometimes step is called repeatedly before
             * returning control to the user.
             */
            swap();  // current state saved; do_step overwrites previous saved state
            stepper.do_step(eval_system, saved.second, saved.first,
                            current.second, step_size);
            current.first = saved.first + step_size;
            return true;
        }
        void reset(state_type new_state) {
            /* Erases any saved state and synchronizes.
             */
            saved = current;
            synchronize();
        }
        void swap() {
            /* Does not synchronize, because swap is most commonly called just before
             * stepper.do_step, which also synchronizes.
             */
            std::swap(current, saved);
        }
    };
    template<typename system_type, typename stepper_type>
    class free_step_iterator
      : public fixed_step_iterator<system_type, stepper_type> {
        /*
         */
        typedef fixed_step_iterator<system_type, stepper_type> base_type;
        typedef stepper_type::state_type state_type;
        typedef std::pair<double, state_type> pair_type;

        double step_size;

    public:
        free_step_iterator(system_type system) : base_type(system) {}

        bool step(double step_size) {
            swap();
            auto result = stepper.try_step(eval_system, saved.second, saved.first,
                                           current.second, step_size);
            if(result == ode::success) {
                current.first = saved.first + step_size;
            } else return false;

        }
    };

}
