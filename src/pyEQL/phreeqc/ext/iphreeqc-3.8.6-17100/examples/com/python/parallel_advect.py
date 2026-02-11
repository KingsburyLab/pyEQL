# -*- coding: utf-8 -*-

"""A coupled advection-reaction model using the COM server.

The sketch below. shows the setup. The coupled model contains a advection and a
reaction model. The reaction model can have one or more Phreeqc calculators.
The calculators can work in parallel using the module `multiprocessing`.

+--------------------------------------------------------------------------+
|                            CoupledModel                                  |
|  +--------------------------------------------------------------------+  |
|  |                         AdvectionModel                             |  |
|  |                                                                    |  |
|  +--------------------------------------------------------------------+  |
|  |                         ReactionModel                              |  |
|  |  +--------------------+--------------------+--------------------+  |  |
|  |  | PhreeqcCalculator1 | PhreeqcCalculator2 | PhreeqcCalculator3 |  |  |
+--------------------------------------------------------------------------+

Author: Mike Müller, mmueller@hydrocomputing.com
"""

import multiprocessing
import os
import time

import matplotlib.pyplot as plt
from win32com.client import Dispatch


class CoupledModel(object):
    """This is a coupled advection model.

    Since it is just a simple example, we use a 1D model.
    The same approach can be applied to a 2d or 3D model
    as long as the advection model supports it.
    The PHREEQC part is the same.
    Furthermore, instead of a simple advection model,
    we can have a more sophisticated transport model.
    """

    def __init__(self, ncells, nshifts, initial_conditions, processes):
        self.nshifts = nshifts
        self.reaction_model = ReactionModel(ncells, initial_conditions,
                                            processes)
        self.reaction_model.make_initial_state()
        init_conc = dict([(name, [value] * ncells) for name, value in
                          list(self.reaction_model.init_conc.items())])
        self.advection_model = AdvectionModel(init_conc,
                                              self.reaction_model.inflow_conc)
        self.component_names = self.reaction_model.component_names
        self.results = {}
        for name in self.component_names:
            self.results[name] = []

    def run(self):
        """Go over all time steps (shifts).
        """
        for shift in range(self.nshifts):
            self.advection_model.advect()
            self.advection_model.save_results(self.results)
            self.reaction_model.modify(self.advection_model.conc)
            self.advection_model.update_conc(self.reaction_model.conc)
        self.reaction_model.finish()


class AdvectionModel(object):
    """Very simple 1D advection model.

    This model can be replaced by a more sophisticated transport
    model with two or three dimensions.
    This could be an external code such as MT3D or anything
    else that is moving concentrations through the subsurface.
    """

    def __init__(self, init_conc, inflow_conc):
        """Set the initial and inflow concentrations.

        Both concentrations are dictionaries with the specie names as keys.
        Values are 1D arrays for `init_conc` and scalars for `inflow_conc`.
        """
        self.conc = init_conc
        self.inflow_conc = inflow_conc
        self.outflow = {}

    def update_conc(self, new_conc):
        """Update the concentrations after the reactions step.

        This is very simple but could be more involved
        if the transport model is more complex.
        """
        self.conc = new_conc

    def advect(self):
        """Shift one cell.
        """
        for name in self.conc:
            self.outflow[name] = self.conc[name][-1]
            self.conc[name][1:] = self.conc[name][:-1]
            self.conc[name][0] = self.inflow_conc[name]

    def save_results(self, results):
        """Save the calculation results.

        Typically, we would write our results into a file.
        For simplicity we just add the current outflow that
        we stored in `self.outflow` and add it to `results`,
        which is a dictionary with all specie names as keys
        and lists as values.
        """
        for name in self.conc:
            results[name].append(self.outflow[name])


class ReactionModel(object):
    """Calculate reactions using PHREEQC as computational engine.

    We have no direct contact with IPhreeqc here.
    We make one or more instances of `PhreeqcCalculator`
    that are actually using IPhreeqc.
    We can use more than one processor with `multiprocessing`.
    """

    def __init__(self, ncells, initial_conditions, processes):
        if processes > ncells:
            raise ValueError('Number of processes needs to be less or equal '
                             'than number of cells. %d processes %d cells.'
                             % (processes, ncells))
        if processes < 1:
            raise ValueError('Need at least one process got %d' % processes)
        self.parallel = False
        if processes > 1:
            self.parallel = True
        self.ncells = ncells
        self.initial_conditions = initial_conditions
        self.processes = processes
        self.inflow_conc = {}
        self.init_conc = {}
        self.conc = {}
        self.component_names = []
        self.calculators = []
        self.cell_ranges = []
        self._init_calculators()
        self.make_initial_state()

    def _init_calculators(self):
        """If we are going parallel we need several calculators.
        """
        if self.parallel:
            # Domain decomposition.
            slave_ncells, reminder = divmod(self.ncells, self.processes)
            root_ncells = slave_ncells + reminder
            current_cell = root_ncells
            root_calculator = PhreeqcCalculator(root_ncells,
                                                self.initial_conditions)
            self.calculators = [root_calculator]
            self.cell_ranges = [(0, root_ncells)]
            for process in range(self.processes - 1):
                self.calculators.append(PhreeqcCalculatorProxy(slave_ncells,
                                                    self.initial_conditions))
                self.cell_ranges.append((current_cell,
                                         current_cell + slave_ncells))
                current_cell += slave_ncells
            assert current_cell == self.ncells
            self.calculators.reverse()
            self.cell_ranges.reverse()
        else:
            root_calculator = PhreeqcCalculator(self.ncells,
                                                self.initial_conditions)
            # Just one calculator and the entire range but still use a list
            # to provide the same interface as the parallel case.
            self.calculators = [root_calculator]
            self.cell_ranges = [(0, self.ncells)]

    def make_initial_state(self):
        """Get the initial values from the calculator(s).
        """
        self.inflow_conc = self.calculators[0].inflow_conc
        self.init_conc = self.calculators[0].init_conc
        self.component_names = self.calculators[0].component_names
        if self.parallel:
            # Make sure all calculators are initialized the same.
            for calculator in self.calculators[1:]:
                assert self.inflow_conc == calculator.inflow_conc
                assert self.init_conc == calculator.init_conc
                assert self.component_names == calculator.component_names

    def modify(self, new_conc):
        """Pass new conc after advection to the calculator.
        """
        self.conc = {}
        for name in self.component_names:
            self.conc[name] = []
        for cell_range, calculator in zip(self.cell_ranges, self.calculators):
            current_conc = dict([(name, value[cell_range[0]:cell_range[1]]) for
                                  name, value in list(new_conc.items())])
            calculator.modify(current_conc)
        for calculator in self.calculators:
            conc = calculator.get_modified()
            for name in self.component_names:
                self.conc[name].extend(conc[name])

    def finish(self):
        """This is necessary for multiprocessing.

        Multiprocessing uses external processes. These need to be
        explicitly closed to avoid hanging of the program at
        the end.
        """
        for calculator in self.calculators:
            calculator.finish()


class PhreeqcCalculator(object):
    """All PHREEQC calculations happen here.

    This is the only place where we interact with IPhreeqc.
    Each instance of this class might run in a different
    process using `multiprocessing`.
    """

    def __init__(self, ncells, initial_conditions):
        """
        ncells - number of cells
        initial_conditions - string containing PHREEQC input for
                             solution and exchange, see example below
        """
        self.ncells = ncells
        self.initial_conditions = initial_conditions
        self.inflow_conc = {}
        self.init_conc = {}
        self.conc = {}
        self.phreeqc = Dispatch('IPhreeqcCOM.Object')
        self.phreeqc.LoadDatabase(r"phreeqc.dat")
        self.components = []
        self.component_names = []
        self._make_initial_state()

    def _make_initial_state(self):
        """Copy solution to all cells and calculate initial conditions.
        """
        self.phreeqc.RunString(self.initial_conditions)
        self.components = self.phreeqc.GetComponentList()
        start = 1
        end = self.ncells
        code = ''
        code += "COPY solution 1 %d-%d\n" % (start, end)
        code += "COPY exchange 1 %d-%d\n" % (start, end)
        code += "END\n"
        code += "RUN_CELLS; -cells %d-%d\n" % (start, end)
        code += self.make_selected_output(self.components)
        self.phreeqc.RunString(code)
        self.conc = self.get_selected_output()
        all_names = list(self.conc.keys())
        self.component_names = [name for name in all_names if name not in
                                ('cb', 'H', 'O')]
        code = ''
        code += self.make_selected_output(self.components)
        code += "RUN_CELLS; -cells 0-1\n"
        self.phreeqc.RunString(code)
        start_conc = self.get_selected_output()
        for name in self.component_names:
            self.inflow_conc[name] = start_conc[name][0]
            self.init_conc[name] = start_conc[name][1]

    def modify(self, new_conc):
        """Set new concentration after advection and re-calculate.
        """
        conc = self.conc
        end = self.ncells + 1
        conc.update(new_conc)
        modify = []
        for index, cell in enumerate(range(1, end)):
            modify.append("SOLUTION_MODIFY %d" % cell)
            modify.append("\t-cb      %e" % conc['cb'][index])
            modify.append("\t-total_h %s" % conc['H'][index])
            modify.append("\t-total_o %s" % conc['O'][index])
            modify.append("\t-totals")
            for name in self.component_names:
                modify.append("\t\t%s\t%s" % (name, conc[name][index]))
        modify.append("RUN_CELLS; -cells %d-%d\n" % (1, self.ncells))
        code = '\n'.join(modify)
        self.phreeqc.RunString(code)
        self.conc = self.get_selected_output()

    def get_modified(self):
        """Return calculated conc.
        """
        return self.conc

    @ staticmethod # this is just a function but belongs here
    def make_selected_output(components):
        """
        Build SELECTED_OUTPUT data block.
        """
        headings = "-headings    cb    H    O    "
        headings += '\t'.join(components)
        selected_output = """
        SELECTED_OUTPUT
            -reset false
        USER_PUNCH
        """
        selected_output += headings + "\n"
        # charge balance, H, and O
        code = '10 w = TOT("water")\n'
        code += '20 PUNCH CHARGE_BALANCE, TOTMOLE("H"), TOTMOLE("O")\n'
        # All other elements
        lino = 30
        for component in components:
            code += '%d PUNCH w*TOT(\"%s\")\n' % (lino, component)
            lino += 10
        selected_output += code
        return selected_output

    def get_selected_output(self):
        """Return calculation result as dict.

        Header entries are the keys and the columns
        are the values as lists of numbers.
        """
        output = self.phreeqc.GetSelectedOutputArray()
        header = output[0]
        conc = {}
        for head in header:
            conc[head] = []
        for row in output[1:]:
            for col, head in enumerate(header):
                conc[head].append(row[col])
        return conc

    def finish(self):
        """Placeholder to give same interface as the multiprocessing version.
        """
        pass


class PhreeqcCalculatorProxy(object):
    """Proxy that communicates with other processes.

    We uses this proxy for parallel computations.
    All code that is specific for parallel computing is located
    in here.
    """

    def __init__(self, ncells, initial_conditions):
        """Go parallel.
        """
        self.in_queue = multiprocessing.JoinableQueue()
        self.out_queue = multiprocessing.JoinableQueue()
        self.process = multiprocessing.Process(
            target=process_worker,
            args=(ncells, initial_conditions, self.in_queue, self.out_queue))
        self.process.start()
        (self.inflow_conc,
         self.init_conc,
         self.component_names) = self.out_queue.get()

    def modify(self, new_conc):
        """Run PHREEQC in another process.
        """
        self.in_queue.put(new_conc)

    def get_modified(self):
        """Return calculated conc.
        """
        return self.out_queue.get()

    def finish(self):
        """Terminate the process.
        """
        self.in_queue.put(None)
        self.process.join()


def process_worker(ncells, initial_conditions, in_queue, out_queue):
    """This runs in another process.
    """
    print('Started process with ID', os.getpid())
    calculator = PhreeqcCalculator(ncells, initial_conditions)
    out_queue.put((calculator.inflow_conc, calculator.init_conc,
                  calculator.component_names))
    while True:
        new_conc = in_queue.get()
        # None is the sentinel. We are done
        if new_conc is None:
            break
        calculator.modify(new_conc)
        out_queue.put(calculator.conc)


def plot(ncells, outflow, specie_names):
    """Plot the results.
    """
    colors = {'Ca': 'r', 'Cl': 'b', 'K': 'g', 'N': 'y', 'Na': 'm'}
    x = [i / float(ncells) for i in
         range(1, len(outflow[specie_names[0]]) + 1)]
    args = []
    for name in specie_names:
        args.extend([x, outflow[name], colors[name]])
    # pylint: disable-msg=W0142
    # we do want *
    plt.plot(*args)
    plt.legend(specie_names, loc=(0.8, 0.5))
    plt.ylabel('MILLIMOLES PER KILOGRAM WATER')
    plt.xlabel('PORE VOLUME')
    plt.show()


def measure_time(func, *args, **kwargs):
    """Convenience function to measure run times.
    """
    import sys
    start = time.perf_counter()
    result = func(*args, **kwargs)
    return result, time.perf_counter() - start

if __name__ == '__main__':

    def main(ncells, nshifts, processes=2):
        """
        Specify initial conditions data blocks.

        Uniform initial conditions are assumed.
        """
        initial_conditions = """
        TITLE Example 11.--Transport and ion exchange.
        SOLUTION 0  CaCl2
            units            mmol/kgw
            temp             25.0
            pH               7.0     charge
            pe               12.5    O2(g)   -0.68
            Ca               0.6
            Cl               1.2
        SOLUTION 1  Initial solution for column
            units            mmol/kgw
            temp             25.0
            pH               7.0     charge
            pe               12.5    O2(g)   -0.68
            Na               1.0
            K                0.2
            N(5)             1.2
            END
        EXCHANGE 1
            equilibrate 1
            X                0.0011
        END
            """

        def run():
            """Do the work.
            """
            model = CoupledModel(ncells, nshifts, initial_conditions,
                                 processes)
            model.run()
            return model, model.results
        (model, outflow), run_time = measure_time(run)
        print('Statistics')
        print('==========')
        print('number of cells:    ', ncells)
        print('number of shifts:   ', nshifts)
        print('number of processes:', processes)
        print('run_time:           ', run_time)
        plot(ncells, outflow, model.component_names)

    main(ncells=400, nshifts=1200, processes=2)
