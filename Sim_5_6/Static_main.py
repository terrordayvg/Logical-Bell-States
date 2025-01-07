#### Data for bell states for Surface 3,5 and BS 3,5 codes.
#### Static noise --------
#### Author: Vladlen Galetsky, last update: 6/1/2025 
###########################################################
"""
Requirements: 

Define: 
1-j: amount of QEC cycles
2-pp: vector of physical error rate
3-max_shots: 1 000 000


Codes:
circuit3up: d=3 Bacon Shor code
circ5d: d=5 rotated Surface code
circsurf_notr: d=3 unrotated Surface code
circsurf_notr2: d=3 rotated Surface code
circbs5: d=5 Bacon Shor code
circuit2: Unencoded Bell states

"""

import itertools
import stim
import pytest
import sys 
import sinter 
import pymatching
import random

sys.path.insert(1, '...') #input path here ----------------------------------------- 
from typing import List
import matplotlib.pyplot as plt

import sinter
import stim
from typing import List
import matplotlib.pyplot as plt
import numpy as np

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info import state_fidelity
import math
import sys
from typing import Callable, TypeVar, List, Any, Iterable, Optional, TYPE_CHECKING, Dict, Union, Literal, Tuple
from typing import cast

import numpy as np

from sinter._probability_util import fit_binomial, shot_error_rate_to_piece_error_rate, Fit
from qiskit import QuantumCircuit, transpile
#from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
from qiskit_aer.noise import NoiseModel
from qiskit import transpile, assemble
#from qiskit_ibm_provider import IBMProvider

TCurveId = TypeVar('TCurveId')

MARKERS: str = "ov*sp^<>8PhH+xXDd|" * 100
T = TypeVar('T')
TVal = TypeVar('TVal')
TKey = TypeVar('TKey')
import gen

def better_sorted_str_terms(val: Any) -> Any:

    if isinstance(val, tuple):
        return tuple(better_sorted_str_terms(e) for e in val)
    if not isinstance(val, str):
        return val
    terms = split_by(val, lambda c: c in '.0123456789')
    result = []
    for term in terms:
        term = ''.join(term)
        if '.' in term:
            try:
                term = float(term)
            except ValueError:
                try:
                    term = tuple(int(e) for e in term.split('.'))
                except ValueError:
                    pass
        else:
            try:
                term = int(term)
            except ValueError:
                pass
        result.append(term)
    return tuple(LooseCompare(e) for e in result)

def group_by(items: Iterable[TVal],
             *,
             key: Callable[[TVal], TKey],
             ) -> Dict[TKey, List[TVal]]:
    

    result: Dict[TKey, List[TVal]] = {}

    for item in items:
        curve_id = key(item)
        result.setdefault(curve_id, []).append(item)

    return result

def plot_custom(
        *,
        stats: 'Iterable[sinter.TaskStats]',
        x_func: Callable[['sinter.TaskStats'], Any],
        y_func: Callable[['sinter.TaskStats'], Union['sinter.Fit', float, int]],
        group_func: Callable[['sinter.TaskStats'], TCurveId] = lambda _: None,
        filter_func: Callable[['sinter.TaskStats'], Any] = lambda _: True,
        plot_args_func: Callable[[int, TCurveId, List['sinter.TaskStats']], Dict[str, Any]] = lambda index, group_key, group_stats: dict(),
        line_fits: Optional[Tuple[Literal['linear', 'log', 'sqrt'], Literal['linear', 'log', 'sqrt']]] = None,
) -> None:


    # Backwards compatibility to when the group stats argument wasn't present.
    import inspect
    if len(inspect.signature(plot_args_func).parameters) == 2:
        old_plot_args_func = cast(Callable[[int, TCurveId], Any], plot_args_func)
        plot_args_func = lambda a, b, _: old_plot_args_func(a, b)

    filtered_stats: List['sinter.TaskStats'] = [
        stat
        for stat in stats
        if filter_func(stat)
    ]

    curve_groups = group_by(filtered_stats, key=group_func)
    for k, curve_id in enumerate(sorted(curve_groups.keys(), key=better_sorted_str_terms)):
        this_group_stats = sorted(curve_groups[curve_id], key=x_func)

        xs = []
        ys = []
        xs_range = []
        ys_low = []
        ys_high = []
        saw_fit = False
        for stat in this_group_stats:
            num_kept = stat.shots - stat.discards
            if num_kept == 0:
                continue
            x = float(x_func(stat))
            y = y_func(stat)
            if isinstance(y, Fit):
                xs_range.append(x)
                ys_low.append(y.low)
                ys_high.append(y.high)
                saw_fit = True
                y = y.best
            if not math.isnan(y):
                xs.append(x)
                ys.append(y)
        curve_id="BC code d=4"
        kwargs: Dict[str, Any] = dict(plot_args_func(k, curve_id, this_group_stats))
        kwargs.setdefault('marker', MARKERS[k])
        if curve_id is not None:
            kwargs.setdefault('label', "BC code d=4")
            kwargs.setdefault('color', f'C{k}')
        kwargs.setdefault('color', 'black')
        #ax.plot(xs, ys, **kwargs)

        print("-------------------------------------")
        print("Real output:")
        print(xs,ys)

        if line_fits is not None and len(set(xs)) >= 2:
            x_scale, y_scale = line_fits
            fit_xs = _rescale(xs, x_scale, False)
            fit_ys = _rescale(ys, y_scale, False)

            from scipy.stats import linregress
            line_fit = linregress(fit_xs, fit_ys)

            x0 = fit_xs[0]
            x1 = fit_xs[-1]
            dx = x1 - x0
            x0 -= dx*10
            x1 += dx*10
            if x0 < 0 <= fit_xs[0] > x0 and x_scale == 'sqrt':
                x0 = 0

            out_xs = np.linspace(x0, x1, 1000)
            out_ys = out_xs * line_fit.slope + line_fit.intercept
            out_xs = _rescale(out_xs, x_scale, True)
            out_ys = _rescale(out_ys, y_scale, True)

            line_kwargs = kwargs.copy()
            line_kwargs.pop('marker', None)
            line_kwargs.pop('label', "BC code d=4")
            line_kwargs['linestyle'] = '--'
            line_kwargs.setdefault('linewidth', 1)
            line_kwargs['linewidth'] /= 2
            ax.plot(out_xs, out_ys, **line_kwargs)

        if saw_fit:
            fit_kwargs = kwargs.copy()
            fit_kwargs.setdefault('zorder', 0)
            fit_kwargs.setdefault('alpha', 1)
            fit_kwargs['zorder'] -= 100
            fit_kwargs['alpha'] *= 0.25
            fit_kwargs.pop('marker', None)
            fit_kwargs.pop('linestyle', None)
            fit_kwargs.pop('label', "BC code d=4")
            #ax.fill_between(xs_range, ys_low, ys_high, **fit_kwargs)

            print("-------------------------------------")
            print("Area ranges")
            print(xs_range,ys_low,ys_high)
            return xs, ys, ys_low, ys_high

def plot_error_rate2(
        *,
        stats: 'Iterable[sinter.TaskStats]',
        x_func: Callable[['sinter.TaskStats'], Any],
        failure_units_per_shot_func: Callable[['sinter.TaskStats'], Any] = lambda _: 1,
        failure_values_func: Callable[['sinter.TaskStats'], Any] = lambda _: 1,
        group_func: Callable[['sinter.TaskStats'], TCurveId] = lambda _: None,
        filter_func: Callable[['sinter.TaskStats'], Any] = lambda _: True,
        plot_args_func: Callable[[int, TCurveId, List['sinter.TaskStats']], Dict[str, Any]] = lambda index, group_key, group_stats: dict(),
        highlight_max_likelihood_factor: Optional[float] = 1e3,
        line_fits: Optional[Tuple[Literal['linear', 'log', 'sqrt'], Literal['linear', 'log', 'sqrt']]] = None,
) -> None:
    """Plots error rates in curves with uncertainty highlights.

    Args:
        ax: The plt.Axes to plot onto. For example, the `ax` value from `fig, ax = plt.subplots(1, 1)`.
        stats: The collected statistics to plot.
        x_func: The X coordinate to use for each stat's data point. For example, this could be
            `x_func=lambda stat: stat.json_metadata['physical_error_rate']`.
        failure_units_per_shot_func: How many error chances there are per shot. This rescales what the
            logical error rate means. By default, it is the logical error rate per shot, but this allows
            you to instead make it the logical error rate per round. For example, if the metadata
            associated with a shot has a field 'r' which is the number of rounds, then this can be
            achieved with `failure_units_per_shot_func=lambda stats: stats.metadata['r']`.
        failure_values_func: How many independent ways there are for a shot to fail, such as
            the number of independent observables in a memory experiment. This affects how the failure
            units rescaling plays out (e.g. with 1 independent failure the "center" of the conversion
            is at 50% whereas for 2 independent failures the "center" is at 75%).
        group_func: Optional. When specified, multiple curves will be plotted instead of one curve.
            The statistics are grouped into curves based on whether or not they get the same result
            out of this function. For example, this could be `group_func=lambda stat: stat.decoder`.
        filter_func: Optional. When specified, some curves will not be plotted.
            The statistics are filtered and only plotted if filter_func(stat) returns True.
            For example, `filter_func=lambda s: s.json_metadata['basis'] == 'x'` would plot only stats
            where the saved metadata indicates the basis was 'x'.
        plot_args_func: Optional. Specifies additional arguments to give the the underlying calls to
            `plot` and `fill_between` used to do the actual plotting. For example, this can be used
            to specify markers and colors. Takes the index of the curve in sorted order and also a
            curve_id (these will be 0 and None respectively if group_func is not specified). For example,
            this could be:

                plot_args_func=lambda index, curve_id: {'color': 'red'
                                                        if curve_id == 'pymatching'
                                                        else 'blue'}

        highlight_max_likelihood_factor: Controls how wide the uncertainty highlight region around curves is.
            Must be 1 or larger. Hypothesis probabilities at most that many times as unlikely as the max likelihood
            hypothesis will be highlighted.
        line_fits: Defaults to None. Set this to a tuple (x_scale, y_scale) to include a dashed line
            fit to every curve. The scales determine how to transform the coordinates before
            performing the fit, and can be set to 'linear', 'sqrt', or 'log'.
    """
    if highlight_max_likelihood_factor is None:
        highlight_max_likelihood_factor = 1
    if not (highlight_max_likelihood_factor >= 1):
        raise ValueError(f"not (highlight_max_likelihood_factor={highlight_max_likelihood_factor} >= 1)")

    def y_func(stat: 'sinter.TaskStats') -> Union[float, 'sinter.Fit']:
        result = fit_binomial(
            num_shots=stat.shots - stat.discards,
            num_hits=stat.errors,
            max_likelihood_factor=highlight_max_likelihood_factor,
        )

        pieces = failure_units_per_shot_func(stat)
        values = failure_values_func(stat)
        result = Fit(
            low=shot_error_rate_to_piece_error_rate(result.low, pieces=pieces, values=values),
            best=shot_error_rate_to_piece_error_rate(result.best, pieces=pieces, values=values),
            high=shot_error_rate_to_piece_error_rate(result.high, pieces=pieces, values=values),
        )

        if stat.errors == 0:
            result = Fit(low=result.low, high=result.high, best=float('nan'))

        if highlight_max_likelihood_factor == 1:
            return result.best
        return result

    print("Dataset:============================================")
    print(x_func)
    print(y_func)

    X,Y,Y_l,Y_h=plot_custom(
        stats=stats,
        x_func=x_func,
        y_func=y_func,
        group_func=group_func,
        filter_func=filter_func,
        plot_args_func=plot_args_func,
        line_fits=line_fits,
    )

    #print(X,Y,Y_l,Y_h)
    return X,Y,Y_l,Y_h


@pytest.mark.parametrize('width,height,basis,rounds', itertools.product(
    [4, 6],
    [6, 12],
    ['X', 'Z'],
    [1, 5, 6],
))
def test_make_bacon_shor_xx_lattice_surgery_circuit(width: int, height: int, basis: str, rounds: int):
    circuit = make_bacon_shor_xx_lattice_surgery_circuit(
        width=width,
        height=height,
        basis=basis,
        rounds=rounds,
    )
    circuit = gen.NoiseModel.uniform_depolarizing(1e-3).noisy_circuit(circuit)
    circuit.detector_error_model()

    expected_determined = circuit.num_detectors + circuit.num_observables
    assert gen.count_determined_measurements_in_circuit(circuit) == expected_determined

    assert circuit.num_ticks == (rounds + 2) * 4 + 1

    expected_distance = min(width // 2, rounds) if basis == 'X' else height // 2
    actual_distance = len(circuit.shortest_graphlike_error())
    assert actual_distance == expected_distance

    assert gen.gates_used_by_circuit(circuit) <= {
        'R',
        'M',
        'RX',
        'MX',
        'MXX',
        'MZZ',

        'TICK',
        'DEPOLARIZE1',
        'DEPOLARIZE2',
        'DETECTOR',
        'X_ERROR',
        'Z_ERROR',
        'OBSERVABLE_INCLUDE',
        'QUBIT_COORDS',
        'SHIFT_COORDS',
    }


def test_gate_counts():
    circuit = make_bacon_shor_xx_lattice_surgery_circuit(
        width=5,
        height=5,
        basis='X',
        rounds=3,
    )
    no_noise = gen.NoiseRule(after={})
    circuit = gen.NoiseModel(
        idle_depolarization=1e-3,
        gate_rules={
            'RX': 0.01,
            'MX': 0.01,
            'MPP': 0.01,
        }
    ).noisy_circuit(circuit)

    assert gen.gate_counts_for_circuit(circuit) == {
        'MXX': 215,
        'MZZ': 200,
        'DEPOLARIZE1': 170,  # Counts idles.

        'RX': 50,
        'MX': 50,

        'DETECTOR': 102,
        'QUBIT_COORDS': 50,
        'SHIFT_COORDS': 5,
        'OBSERVABLE_INCLUDE': 3,
        'TICK': 21,
    }



def count_determined_measurements_in_circuit(circuit: stim.Circuit) -> int:
    """Simulates the circuit, counting how many measurements were determined.

    In most cases, for a quantum error correcting code, the result should be
    related to the number of detectors plus the number of observables declared
    in the circuit.
    """

    num_determined_measurements = 0
    sim = stim.TableauSimulator()
    n = circuit.num_qubits

    def run_block(block: stim.Circuit, reps: int):
        nonlocal num_determined_measurements
        for _ in range(reps):
            for inst in block:
                if isinstance(inst, stim.CircuitRepeatBlock):
                    run_block(inst.body_copy(), inst.repeat_count)
                elif inst.name == 'M' or inst.name == 'MR':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_z(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MX' or inst.name == 'MRX':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_x(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MY' or inst.name == 'MRY':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_y(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MPP':
                    args = inst.gate_args_copy()
                    targets = inst.targets_copy()
                    start = 0
                    while start < len(targets):
                        end = start + 1
                        while end < len(targets) and targets[end].is_combiner:
                            end += 2

                        p = stim.PauliString(n)
                        for t in targets[start:end:2]:
                            if t.is_x_target:
                                p[t.value] = 'X'
                            elif t.is_y_target:
                                p[t.value] = 'Y'
                            elif t.is_z_target:
                                p[t.value] = 'Z'
                            else:
                                raise NotImplementedError(f'{t=} {inst=}')

                        known = sim.peek_observable_expectation(p) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, targets[start:end], args))

                        start = end
                else:
                    sim.do(inst)

    run_block(circuit, 1)
    return num_determined_measurements

def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)

    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors

if __name__ == '__main__':
  
    j=1 #for the first QEC cycle in memory
    pp=np.arange(0.001,0.01,0.001) #Physical error
    pp=[0.001,0.01]

    Vec_Y=[]
    Vec_X=pp       #Iteration vector
    Std_p=[]
    Std_m=[]

    Vec_Y2=[]
    Vec_X2=pp      #Iteration vector
    Std_p2=[]
    Std_m2=[]

    Vec_Y3=[]
    Vec_X3=pp      #Iteration vector
    Std_p3=[]
    Std_m3=[]


    Vec_Y4=[]
    Vec_X4=pp      #Iteration vector
    Std_p4=[]
    Std_m4=[]


    Vec_Y5=[]
    Vec_X5=pp      #Iteration vector
    Std_p5=[]
    Std_m5=[]

    Vec_Y6=[]
    Vec_X6=pp      #Iteration vector
    Std_p6=[]
    Std_m6=[]

    S_Bell=[]

    for i in range(len(pp)):

        #All physical errors are the same:

        p_f=pp[i] # Depolarize 2
        p_f1=p_f  # Depolarize 1
        p_f2=p_f  # Measurement error


        ################### d=3 Bacon-Shor code (BS)
        circuit3up=stim.Circuit("""
QUBIT_COORDS(0, 0) 0
QUBIT_COORDS(0, 1) 1
QUBIT_COORDS(0, 2) 2
QUBIT_COORDS(0, 3) 3
QUBIT_COORDS(0, 4) 4
QUBIT_COORDS(1, 0) 5
QUBIT_COORDS(1, 1) 6
QUBIT_COORDS(1, 2) 7
QUBIT_COORDS(1, 3) 8
QUBIT_COORDS(1, 4) 9
QUBIT_COORDS(2, 0) 10
QUBIT_COORDS(2, 1) 11
QUBIT_COORDS(2, 2) 12
QUBIT_COORDS(2, 3) 13
QUBIT_COORDS(2, 4) 14
QUBIT_COORDS(3, 0) 15
QUBIT_COORDS(3, 1) 16
QUBIT_COORDS(3, 2) 17
QUBIT_COORDS(3, 3) 18
QUBIT_COORDS(3, 4) 19
QUBIT_COORDS(4, 0) 20
QUBIT_COORDS(4, 1) 21
QUBIT_COORDS(4, 2) 22
QUBIT_COORDS(4, 3) 23
QUBIT_COORDS(4, 4) 24
QUBIT_COORDS(5, 0) 25
QUBIT_COORDS(5, 1) 26
QUBIT_COORDS(5, 2) 27
QUBIT_COORDS(5, 3) 28
QUBIT_COORDS(5, 4) 29
QUBIT_COORDS(6, 0) 30
QUBIT_COORDS(6, 1) 31
QUBIT_COORDS(6, 2) 32
QUBIT_COORDS(6, 3) 33
QUBIT_COORDS(6, 4) 34
QUBIT_COORDS(7, 0) 35
QUBIT_COORDS(7, 1) 36
QUBIT_COORDS(7, 2) 37
QUBIT_COORDS(7, 3) 38
QUBIT_COORDS(7, 4) 39
QUBIT_COORDS(8, 0) 40
QUBIT_COORDS(8, 1) 41
QUBIT_COORDS(8, 2) 42
QUBIT_COORDS(8, 3) 43
QUBIT_COORDS(8, 4) 44
QUBIT_COORDS(9, 0) 45
QUBIT_COORDS(9, 1) 46
QUBIT_COORDS(9, 2) 47
QUBIT_COORDS(9, 3) 48
QUBIT_COORDS(9, 4) 49
#ignore creation in first phase of merging
RX 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
TICK
MPP X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49

TICK
MPP X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44

TICK
MPP Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48

TICK
MPP Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
DETECTOR(0.5, 0, 0) rec[-80]
DETECTOR(0.5, 1, 0) rec[-79]
DETECTOR(0.5, 2, 0) rec[-78]
DETECTOR(0.5, 3, 0) rec[-77]
DETECTOR(0.5, 4, 0) rec[-76]
DETECTOR(1.5, 0, 0) rec[-60]
DETECTOR(1.5, 1, 0) rec[-59]
DETECTOR(1.5, 2, 0) rec[-58]
DETECTOR(1.5, 3, 0) rec[-57]
DETECTOR(1.5, 4, 0) rec[-56]
DETECTOR(2.5, 0, 0) rec[-75]
DETECTOR(2.5, 1, 0) rec[-74]
DETECTOR(2.5, 2, 0) rec[-73]
DETECTOR(2.5, 3, 0) rec[-72]
DETECTOR(2.5, 4, 0) rec[-71]
DETECTOR(3.5, 0, 0) rec[-55]
DETECTOR(3.5, 1, 0) rec[-54]
DETECTOR(3.5, 2, 0) rec[-53]
DETECTOR(3.5, 3, 0) rec[-52]
DETECTOR(3.5, 4, 0) rec[-51]
DETECTOR(5.5, 0, 0) rec[-50]
DETECTOR(5.5, 1, 0) rec[-49]
DETECTOR(5.5, 2, 0) rec[-48]
DETECTOR(5.5, 3, 0) rec[-47]
DETECTOR(5.5, 4, 0) rec[-46]
DETECTOR(6.5, 0, 0) rec[-70]
DETECTOR(6.5, 1, 0) rec[-69]
DETECTOR(6.5, 2, 0) rec[-68]
DETECTOR(6.5, 3, 0) rec[-67]
DETECTOR(6.5, 4, 0) rec[-66]
DETECTOR(7.5, 0, 0) rec[-45]
DETECTOR(7.5, 1, 0) rec[-44]
DETECTOR(7.5, 2, 0) rec[-43]
DETECTOR(7.5, 3, 0) rec[-42]
DETECTOR(7.5, 4, 0) rec[-41]
DETECTOR(8.5, 0, 0) rec[-65]
DETECTOR(8.5, 1, 0) rec[-64]
DETECTOR(8.5, 2, 0) rec[-63]
DETECTOR(8.5, 3, 0) rec[-62]
DETECTOR(8.5, 4, 0) rec[-61]
SHIFT_COORDS(0, 0, 1)

TICK #complete

MPP X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X20*X25 X21*X26 X22*X27 X23*X28 X24*X29 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49

TICK #complete

MPP X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44


TICK #complete

MPP Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48

TICK #complete



MPP Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49


DETECTOR(0.5, 0, 0) rec[-165] rec[-164] rec[-163] rec[-162] rec[-161] rec[-85] rec[-84] rec[-83] rec[-82] rec[-81]
DETECTOR(1.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]
DETECTOR(2.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]
DETECTOR(3.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]
DETECTOR(5.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
DETECTOR(6.5, 0, 0) rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]
DETECTOR(7.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
DETECTOR(8.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]
DETECTOR(0, 0.5, 0) rec[-125] rec[-123] rec[-121] rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-109] rec[-107] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]
DETECTOR(0, 1.5, 0) rec[-105] rec[-103] rec[-101] rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-89] rec[-87] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]
DETECTOR(0, 2.5, 0) rec[-124] rec[-122] rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-110] rec[-108] rec[-106] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]
DETECTOR(0, 3.5, 0) rec[-104] rec[-102] rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-90] rec[-88] rec[-86] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]
OBSERVABLE_INCLUDE(2) rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]
SHIFT_COORDS(0, 0, 1)

TICK
REPEAT 3 {

    MPP("""+str(p_f)+""") X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X20*X25 X21*X26 X22*X27 X23*X28 X24*X29 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
    
    DEPOLARIZE2("""+str(p_f)+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 20 25 21 26 22 27 23 28 24 29 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49

    TICK

    MPP("""+str(p_f)+""") X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
    
    DEPOLARIZE2("""+str(p_f)+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1 
    DEPOLARIZE1("""+str(p_f)+""") 0 1 2 3 4
    TICK

    MPP("""+str(p_f)+""") Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48

    DEPOLARIZE2("""+str(p_f)+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
    DEPOLARIZE1("""+str(p_f)+""") 4 9 14 19 24 29 34 39 44






    TICK

    MPP("""+str(p_f)+""") Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49

    DEPOLARIZE2("""+str(p_f)+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
    DEPOLARIZE1("""+str(p_f)+""") 5 10 15 20 25 30 35 40








    DETECTOR(0.5, 0, 0) rec[-170] rec[-169] rec[-168] rec[-167] rec[-166] rec[-85] rec[-84] rec[-83] rec[-82] rec[-81]
    DETECTOR(1.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]
    DETECTOR(2.5, 0, 0) rec[-165] rec[-164] rec[-163] rec[-162] rec[-161] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]
    DETECTOR(3.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]
    DETECTOR(4.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]
    DETECTOR(5.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
    DETECTOR(6.5, 0, 0) rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]
    DETECTOR(7.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
    DETECTOR(8.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]
    DETECTOR(0, 0.5, 0) rec[-125] rec[-123] rec[-121] rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-109] rec[-107] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]
    DETECTOR(0, 1.5, 0) rec[-105] rec[-103] rec[-101] rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-89] rec[-87] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]
    DETECTOR(0, 2.5, 0) rec[-124] rec[-122] rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-110] rec[-108] rec[-106] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]
    DETECTOR(0, 3.5, 0) rec[-104] rec[-102] rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-90] rec[-88] rec[-86] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]
    SHIFT_COORDS(0, 0, 1)
    TICK
}




MPP("""+str(p_f)+""") X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49

DEPOLARIZE2("""+str(p_f)+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
DEPOLARIZE1("""+str(p_f)+""") 20 21 22 23 24 25 26 27 28 29


TICK

MPP("""+str(p_f)+""") X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44

DEPOLARIZE2("""+str(p_f)+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1 
DEPOLARIZE1("""+str(p_f)+""") 0 1 2 3 4


TICK

MPP("""+str(p_f)+""") Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48

DEPOLARIZE2("""+str(p_f)+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
DEPOLARIZE1("""+str(p_f)+""") 4 9 14 19 24 29 34 39 44


TICK

MPP("""+str(p_f)+""") Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
DEPOLARIZE2("""+str(p_f)+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
DEPOLARIZE1("""+str(p_f)+""") 5 10 15 20 25 30 35 40

DETECTOR(0.5, 0, 0) rec[-165] rec[-164] rec[-163] rec[-162] rec[-161] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]
DETECTOR(1.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]
DETECTOR(2.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]
DETECTOR(3.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]
DETECTOR(5.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
DETECTOR(6.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]
DETECTOR(7.5, 0, 0) rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
DETECTOR(8.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]
DETECTOR(0, 0.5, 0) rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32]
DETECTOR(0, 1.5, 0) rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12]
DETECTOR(0, 2.5, 0) rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31]
DETECTOR(0, 3.5, 0) rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11]
DETECTOR(5, 0.5, 0) rec[-110] rec[-108] rec[-106] rec[-104] rec[-102] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]
DETECTOR(5, 1.5, 0) rec[-90] rec[-88] rec[-86] rec[-84] rec[-82] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]
DETECTOR(5, 2.5, 0) rec[-109] rec[-107] rec[-105] rec[-103] rec[-101] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]
DETECTOR(5, 3.5, 0) rec[-89] rec[-87] rec[-85] rec[-83] rec[-81] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]
SHIFT_COORDS(0, 0, 1)


REPEAT """+str(j)+"""{
    

    MPP("""+str(p_f)+""") X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
    DEPOLARIZE2("""+str(p_f)+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
    DEPOLARIZE1("""+str(p_f)+""") 20 21 22 23 24 25 26 27 28 29







    TICK
    MPP("""+str(p_f)+""") X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
    DEPOLARIZE2("""+str(p_f)+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1 
    DEPOLARIZE1("""+str(p_f)+""") 0 1 2 3 4 45 46 47 48 49


    TICK
    MPP("""+str(p_f)+""") Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48
    DEPOLARIZE2("""+str(p_f)+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
    DEPOLARIZE1("""+str(p_f)+""") 4 9 14 19 24 29 34 39







    TICK
   
    MPP("""+str(p_f)+""") Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
    DEPOLARIZE2("""+str(p_f)+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
    DEPOLARIZE1("""+str(p_f)+""") 5 10 15 20 25 30 35 40


    
    DETECTOR(0.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]#
    DETECTOR(1.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]#
    DETECTOR(2.5, 0, 0) rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]#
    DETECTOR(3.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]#
    DETECTOR(5.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]#
    DETECTOR(6.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]#
    DETECTOR(7.5, 0, 0) rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]#
    DETECTOR(8.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]#
    DETECTOR(0, 0.5, 0) rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32]#
    DETECTOR(0, 1.5, 0) rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12]#
    DETECTOR(0, 2.5, 0) rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31]#
    DETECTOR(0, 3.5, 0) rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11]#
    DETECTOR(5, 0.5, 0) rec[-110] rec[-108] rec[-106] rec[-104] rec[-102] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]#
    DETECTOR(5, 1.5, 0) rec[-90] rec[-88] rec[-86] rec[-84] rec[-82] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]#
    DETECTOR(5, 2.5, 0) rec[-109] rec[-107] rec[-105] rec[-103] rec[-101] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]#
    DETECTOR(5, 3.5, 0) rec[-89] rec[-87] rec[-85] rec[-83] rec[-81] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]#

}

TICK
MX("""+str(p_f)+""") 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
DETECTOR(0.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
DETECTOR(1.5, 0, 0) rec[-110] rec[-109] rec[-108] rec[-107] rec[-106] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41] rec[-40] rec[-39] rec[-38] rec[-37] rec[-36]
DETECTOR(2.5, 0, 0) rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-40] rec[-39] rec[-38] rec[-37] rec[-36] rec[-35] rec[-34] rec[-33] rec[-32] rec[-31]
DETECTOR(3.5, 0, 0) rec[-105] rec[-104] rec[-103] rec[-102] rec[-101] rec[-35] rec[-34] rec[-33] rec[-32] rec[-31] rec[-30] rec[-29] rec[-28] rec[-27] rec[-26]
DETECTOR(5.5, 0, 0) rec[-100] rec[-99] rec[-98] rec[-97] rec[-96] rec[-25] rec[-24] rec[-23] rec[-22] rec[-21] rec[-20] rec[-19] rec[-18] rec[-17] rec[-16]
DETECTOR(6.5, 0, 0) rec[-120] rec[-119] rec[-118] rec[-117] rec[-116] rec[-20] rec[-19] rec[-18] rec[-17] rec[-16] rec[-15] rec[-14] rec[-13] rec[-12] rec[-11]
DETECTOR(7.5, 0, 0) rec[-95] rec[-94] rec[-93] rec[-92] rec[-91] rec[-15] rec[-14] rec[-13] rec[-12] rec[-11] rec[-10] rec[-9] rec[-8] rec[-7] rec[-6]
DETECTOR(8.5, 0, 0) rec[-115] rec[-114] rec[-113] rec[-112] rec[-111] rec[-10] rec[-9] rec[-8] rec[-7] rec[-6] rec[-5] rec[-4] rec[-3] rec[-2] rec[-1]
OBSERVABLE_INCLUDE(0) rec[-30] rec[-29] rec[-28] rec[-27] rec[-26]
OBSERVABLE_INCLUDE(1) rec[-25] rec[-24] rec[-23] rec[-22] rec[-21]
               
""")
        ################### d=5 rotated Surface code (S)

        circ5d=stim.Circuit("""
#system 1
QUBIT_COORDS(1, 1) 1 #65
QUBIT_COORDS(2, 0) 2 #66
QUBIT_COORDS(3, 1) 3 #67
QUBIT_COORDS(5, 1) 5 #69
QUBIT_COORDS(6, 0) 6 #70
QUBIT_COORDS(7, 1) 7 #71
QUBIT_COORDS(9, 1) 9 #73
QUBIT_COORDS(1, 3) 12 #76
QUBIT_COORDS(2, 2) 13 #77
QUBIT_COORDS(3, 3) 14 #78
QUBIT_COORDS(4, 2) 15 #79
QUBIT_COORDS(5, 3) 16 #80
QUBIT_COORDS(6, 2) 17 #81
QUBIT_COORDS(7, 3) 18 #82
QUBIT_COORDS(8, 2) 19 #83
QUBIT_COORDS(9, 3) 20 #84
QUBIT_COORDS(10, 2) 21 #85
QUBIT_COORDS(0, 4) 22 #86
QUBIT_COORDS(1, 5) 23 #87
QUBIT_COORDS(2, 4) 24 #88
QUBIT_COORDS(3, 5) 25 #89
QUBIT_COORDS(4, 4) 26 #90
QUBIT_COORDS(5, 5) 27 #91
QUBIT_COORDS(6, 4) 28 #92
QUBIT_COORDS(7, 5) 29 #93
QUBIT_COORDS(8, 4) 30 #94
QUBIT_COORDS(9, 5) 31 #95
QUBIT_COORDS(1, 7) 34 #98
QUBIT_COORDS(2, 6) 35 #99
QUBIT_COORDS(3, 7) 36 #100
QUBIT_COORDS(4, 6) 37 #101
QUBIT_COORDS(5, 7) 38 #102
QUBIT_COORDS(6, 6) 39 #103
QUBIT_COORDS(7, 7) 40 #104
QUBIT_COORDS(8, 6) 41 #105
QUBIT_COORDS(9, 7) 42 #106
QUBIT_COORDS(10, 6) 43 #107
QUBIT_COORDS(0, 8) 44 #108
QUBIT_COORDS(1, 9) 45 #109
QUBIT_COORDS(2, 8) 46 #110
QUBIT_COORDS(3, 9) 47 #111
QUBIT_COORDS(4, 8) 48 #112
QUBIT_COORDS(5, 9) 49 #113
QUBIT_COORDS(6, 8) 50 #114
QUBIT_COORDS(7, 9) 51 #115
QUBIT_COORDS(8, 8) 52 #116
QUBIT_COORDS(9, 9) 53 #117
QUBIT_COORDS(4, 10) 59 #123
QUBIT_COORDS(8, 10) 63 #127


#11 intermediate qubits, 5 data, 6 X syndrome
#data
QUBIT_COORDS(11, 1) 64 
QUBIT_COORDS(11, 3) 65
QUBIT_COORDS(11, 5) 66
QUBIT_COORDS(11, 7) 67 
QUBIT_COORDS(11, 9) 68

#syndrome
QUBIT_COORDS(10, 0) 69 
QUBIT_COORDS(10, 4) 70
QUBIT_COORDS(10, 8) 71
QUBIT_COORDS(12, 2) 72
QUBIT_COORDS(12, 6) 73
QUBIT_COORDS(12, 10) 74

#system 2
QUBIT_COORDS(13, 1) 75
QUBIT_COORDS(14, 0) 76
QUBIT_COORDS(15, 1) 77
QUBIT_COORDS(17, 1) 79
QUBIT_COORDS(18, 0) 80
QUBIT_COORDS(19, 1) 81
QUBIT_COORDS(21, 1) 83
QUBIT_COORDS(13, 3) 86
QUBIT_COORDS(14, 2) 87
QUBIT_COORDS(15, 3) 88
QUBIT_COORDS(16, 2) 89
QUBIT_COORDS(17, 3) 90
QUBIT_COORDS(18, 2) 91
QUBIT_COORDS(19, 3) 92
QUBIT_COORDS(20, 2) 93
QUBIT_COORDS(21, 3) 94
QUBIT_COORDS(22, 2) 95
QUBIT_COORDS(12, 4) 96
QUBIT_COORDS(13, 5) 97
QUBIT_COORDS(14, 4) 98
QUBIT_COORDS(15, 5) 99
QUBIT_COORDS(16, 4) 100
QUBIT_COORDS(17, 5) 101
QUBIT_COORDS(18, 4) 102
QUBIT_COORDS(19, 5) 103
QUBIT_COORDS(20, 4) 104
QUBIT_COORDS(21, 5) 105
QUBIT_COORDS(13, 7) 108
QUBIT_COORDS(14, 6) 109
QUBIT_COORDS(15, 7) 110
QUBIT_COORDS(16, 6) 111
QUBIT_COORDS(17, 7) 112
QUBIT_COORDS(18, 6) 113
QUBIT_COORDS(19, 7) 114
QUBIT_COORDS(20, 6) 115
QUBIT_COORDS(21, 7) 116
QUBIT_COORDS(22, 6) 117
QUBIT_COORDS(12, 8) 118
QUBIT_COORDS(13, 9) 119
QUBIT_COORDS(14, 8) 120
QUBIT_COORDS(15, 9) 121
QUBIT_COORDS(16, 8) 122
QUBIT_COORDS(17, 9) 123
QUBIT_COORDS(18, 8) 124
QUBIT_COORDS(19, 9) 125
QUBIT_COORDS(20, 8) 126
QUBIT_COORDS(21, 9) 127
QUBIT_COORDS(16, 10) 133
QUBIT_COORDS(20, 10) 137


TICK
RX 1 3 5 7 9 12 14 16 18 20 23 25 27 29 31 34 36 38 40 42 45 47 49 51 53 75 77 79 81 83 86 88 90 92 94 97 99 101 103 105 108 110 112 114 116 119 121 123 125 127 64 65 66 67 68
R 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63 76 80 87 89 91 93 95 96 98 100 102 104 109 111 113 115 117 118 120 122 124 126 133 137 69 70 71 72 73 74
TICK
H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74

#DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
#DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136


TICK
CX 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 69 64 72 86 70 66 73 108 71 68 65 21 67 43 
#DEPOLARIZE2("""+str(p_f)+""") 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 69 64 72 86 70 66 73 108 71 68 65 21 67 43 
#DEPOLARIZE1("""+str(p_f1)+""") 0 1 4 5 8 9 10 11 12 32 33 34 40 54 55 56 57 58 59 60 61 62 63 74 75 78 79 82 83 84 85 95 106 107 117 128 129 130 131 132 133 134 135 136 137




TICK
CX 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 69 9 72 65 70 31 73 67 71 53 64 21 66 43
#DEPOLARIZE2("""+str(p_f)+""") 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 69 9 72 65 70 31 73 67 71 53 64 21 66 43
#DEPOLARIZE1("""+str(p_f1)+""") 0 4 8 10 11 26 32 33 47 51 54 55 56 57 58 59 60 61 62 63 68 74 78 82 83 84 85 95 105 106 107 117 127 128 129 130 131 132 133 134 135 136 137


TICK
CX 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 72 75 70 65 73 97 71 67 74 119 66 96 68 118
#DEPOLARIZE2("""+str(p_f)+""") 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 72 75 70 65 73 97 71 67 74 119 66 96 68 118
#DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 4 6 7 8 10 11 22 23 32 33 44 45 54 55 56 57 58 60 61 62 64 69 76 77 78 80 81 82 84 85 106 107 128 129 130 131 132 134 135 136


TICK
CX 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 72 64 70 20 73 66 71 42 74 68 65 96 67 118
#DEPOLARIZE2("""+str(p_f)+""") 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 72 64 70 20 73 66 71 42 74 68 65 96 67 118
#DEPOLARIZE1("""+str(p_f1)+""") 0 2 4 6 8 10 11 22 32 33 44 45 53 54 55 56 57 58 60 61 62 76 78 80 82 84 85 94 106 107 116 127 128 129 130 131 132 134 135 136



TICK
H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
#DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
#DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

TICK
MR 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63 76 80 87 89 91 93 95 96 98 100 102 104 109 111 113 115 117 118 120 122 124 126 133 137 69 70 71 72 73 74
#DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 14 16 18 20 23 25 27 29 31 32 33 34 36 38 40 42 45 47 49 51 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 88 90 92 94 97 99 101 103 105 106 107 108 110 112 114 116 119 121 123 125 127 128 129 130 131 132 134 135 136


DETECTOR(2, 0, 0) rec[-54]
DETECTOR(2, 4, 0) rec[-46]
DETECTOR(2, 8, 0) rec[-36]
DETECTOR(4, 2, 0) rec[-51]
DETECTOR(4, 6, 0) rec[-41]
DETECTOR(4, 10, 0) rec[-32]
DETECTOR(6, 0, 0) rec[-53]
DETECTOR(6, 4, 0) rec[-44]
DETECTOR(6, 8, 0) rec[-34]
DETECTOR(8, 2, 0) rec[-49]
DETECTOR(8, 6, 0) rec[-39]
DETECTOR(8, 10, 0) rec[-31]

DETECTOR(14, 0, 0) rec[-30]
DETECTOR(14, 4, 0) rec[-22]
DETECTOR(14, 8, 0) rec[-12]
DETECTOR(16, 2, 0) rec[-27]
DETECTOR(16, 6, 0) rec[-17]
DETECTOR(16, 10, 0) rec[-8]
DETECTOR(18, 0, 0) rec[-29]
DETECTOR(18, 4, 0) rec[-20]
DETECTOR(18, 8, 0) rec[-10]
DETECTOR(20, 2, 0) rec[-25]
DETECTOR(20, 6, 0) rec[-15]
DETECTOR(20, 10, 0) rec[-7]

DETECTOR(10, 0, 0) rec[-6]
DETECTOR(10, 4, 0) rec[-5]
DETECTOR(10, 8, 0) rec[-4]
DETECTOR(12, 2, 0) rec[-3]
DETECTOR(12, 6, 0) rec[-2]
DETECTOR(12, 10, 0) rec[-1]


REPEAT 3 {
    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
    DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

    TICK
    CX 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 69 64 72 86 70 66 73 108 71 68 65 21 67 43 
    DEPOLARIZE2("""+str(p_f)+""") 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 69 64 72 86 70 66 73 108 71 68 65 21 67 43 
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 4 5 8 9 10 11 12 32 33 34 40 54 55 56 57 58 59 60 61 62 63 74 75 78 79 82 83 84 85 95 106 107 117 128 129 130 131 132 133 134 135 136 137



    TICK
    CX 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 69 9 72 65 70 31 73 67 71 53 64 21 66 43
    DEPOLARIZE2("""+str(p_f)+""") 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 69 9 72 65 70 31 73 67 71 53 64 21 66 43
    DEPOLARIZE1("""+str(p_f1)+""") 0 4 8 10 11 26 32 33 47 51 54 55 56 57 58 59 60 61 62 63 68 74 78 82 83 84 85 95 105 106 107 117 127 128 129 130 131 132 133 134 135 136 137


    TICK
    CX 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 72 75 70 65 73 97 71 67 74 119 66 96 68 118
    DEPOLARIZE2("""+str(p_f)+""") 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 72 75 70 65 73 97 71 67 74 119 66 96 68 118
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 4 6 7 8 10 11 22 23 32 33 44 45 54 55 56 57 58 60 61 62 64 69 76 77 78 80 81 82 84 85 106 107 128 129 130 131 132 134 135 136


    TICK
    CX 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 72 64 70 20 73 66 71 42 74 68 65 96 67 118
    DEPOLARIZE2("""+str(p_f)+""") 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 72 64 70 20 73 66 71 42 74 68 65 96 67 118
    DEPOLARIZE1("""+str(p_f1)+""") 0 2 4 6 8 10 11 22 32 33 44 45 53 54 55 56 57 58 60 61 62 76 78 80 82 84 85 94 106 107 116 127 128 129 130 131 132 134 135 136



    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
    DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 69 70 71 72 73 74
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

    TICK
    MR("""+str(p_f2)+""") 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63 76 80 87 89 91 93 95 96 98 100 102 104 109 111 113 115 117 118 120 122 124 126 133 137 69 70 71 72 73 74
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 14 16 18 20 23 25 27 29 31 32 33 34 36 38 40 42 45 47 49 51 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 88 90 92 94 97 99 101 103 105 106 107 108 110 112 114 116 119 121 123 125 127 128 129 130 131 132 134 135 136


    SHIFT_COORDS(0, 0, 1)
    DETECTOR(2, 0, 0) rec[-54] rec[-108]
    DETECTOR(6, 0, 0) rec[-53] rec[-107]
    DETECTOR(2, 2, 0) rec[-52] rec[-106]
    DETECTOR(4, 2, 0) rec[-51] rec[-105]
    DETECTOR(6, 2, 0) rec[-50] rec[-104]
    DETECTOR(8, 2, 0) rec[-49] rec[-103]
    DETECTOR(10, 2, 0) rec[-48] rec[-102]
    DETECTOR(0, 4, 0) rec[-47] rec[-101]
    DETECTOR(2, 4, 0) rec[-46] rec[-100]
    DETECTOR(4, 4, 0) rec[-45] rec[-99]
    DETECTOR(6, 4, 0) rec[-44] rec[-98]
    DETECTOR(8, 4, 0) rec[-43] rec[-97]
    DETECTOR(2, 6, 0) rec[-42] rec[-96]
    DETECTOR(4, 6, 0) rec[-41] rec[-95]
    DETECTOR(6, 6, 0) rec[-40] rec[-94]
    DETECTOR(8, 6, 0) rec[-39] rec[-93]
    DETECTOR(10, 6, 0) rec[-38] rec[-92]
    DETECTOR(0, 8, 0) rec[-37] rec[-91]
    DETECTOR(2, 8, 0) rec[-36] rec[-90]
    DETECTOR(4, 8, 0) rec[-35] rec[-89]
    DETECTOR(6, 8, 0) rec[-34] rec[-88]
    DETECTOR(8, 8, 0) rec[-33] rec[-87]
    DETECTOR(4, 10, 0) rec[-32] rec[-86]
    DETECTOR(8, 10, 0) rec[-31] rec[-85]
    
    DETECTOR(14, 0, 0) rec[-30] rec[-84]
    DETECTOR(18, 0, 0) rec[-29] rec[-83]
    DETECTOR(14, 2, 0) rec[-28] rec[-82]
    DETECTOR(16, 2, 0) rec[-27] rec[-81]
    DETECTOR(18, 2, 0) rec[-26] rec[-80]
    DETECTOR(20, 2, 0) rec[-25] rec[-79]
    DETECTOR(22, 2, 0) rec[-24] rec[-78]
    DETECTOR(12, 4, 0) rec[-23] rec[-77]
    DETECTOR(14, 4, 0) rec[-22] rec[-76]
    DETECTOR(16, 4, 0) rec[-21] rec[-75]
    DETECTOR(18, 4, 0) rec[-20] rec[-74]
    DETECTOR(20, 4, 0) rec[-19] rec[-73]
    DETECTOR(14, 6, 0) rec[-18] rec[-72]
    DETECTOR(16, 6, 0) rec[-17] rec[-71]
    DETECTOR(18, 6, 0) rec[-16] rec[-70]
    DETECTOR(20, 6, 0) rec[-15] rec[-69]
    DETECTOR(22, 6, 0) rec[-14] rec[-68]
    DETECTOR(12, 8, 0) rec[-13] rec[-67]
    DETECTOR(14, 8, 0) rec[-12] rec[-66]
    DETECTOR(16, 8, 0) rec[-11] rec[-65]
    DETECTOR(18, 8, 0) rec[-10] rec[-64]
    DETECTOR(20, 8, 0) rec[-9] rec[-63]
    DETECTOR(16, 10, 0) rec[-8] rec[-62]
    DETECTOR(20, 10, 0) rec[-7] rec[-61]
    
    DETECTOR(10, 0, 0) rec[-6] rec[-60]
    DETECTOR(10, 4, 0) rec[-5] rec[-59]
    DETECTOR(10, 8, 0) rec[-4] rec[-58]
    DETECTOR(12, 2, 0) rec[-3] rec[-57]
    DETECTOR(12, 6, 0) rec[-2] rec[-56]
    DETECTOR(12, 10, 0) rec[-1] rec[-55]
}
#[53-49]
MZ("""+str(p_f2)+""") 64 65 66 67 68
DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137

R 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63 76 80 87 89 91 93 95 96 98 100 102 104 109 111 113 115 117 118 120 122 124 126 133 137

TICK
H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

TICK
CX 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 
DEPOLARIZE2("""+str(p_f)+""") 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 
DEPOLARIZE1("""+str(p_f1)+""") 0 1 4 5 8 9 10 11 12 17 20 21 33 34 43 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 78 79 82 83 84 85 86 95 106 107 108 120 128 129 130 131 132 133 134 135 136 137

TICK
CX 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 
DEPOLARIZE2("""+str(p_f)+""") 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 
DEPOLARIZE1("""+str(p_f1)+""") 0 4 8 10 11 26 32 33 47 51 54 55 56 57 58 59 60 61 62 63 68 74 78 82 83 84 85 95 105 106 107 117 127 128 129 130 131 132 133 134 135 136 137


TICK
CX 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 
DEPOLARIZE2("""+str(p_f)+""") 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 
DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 4 6 7 8 10 11 22 23 32 33 44 45 54 55 56 57 58 60 61 62 64 69 76 77 78 80 81 82 84 85 106 107 128 129 130 131 132 134 135 136

TICK
CX 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 
DEPOLARIZE2("""+str(p_f)+""") 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 
DEPOLARIZE1("""+str(p_f1)+""") 0 2 4 6 8 10 11 22 32 33 44 45 53 54 55 56 57 58 60 61 62 76 78 80 82 84 85 94 106 107 116 127 128 129 130 131 132 134 135 136

TICK
H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

TICK
MR("""+str(p_f2)+""") 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63 76 80 87 89 91 93 95 96 98 100 102 104 109 111 113 115 117 118 120 122 124 126 133 137 
DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 14 16 18 20 23 25 27 29 31 32 33 34 36 38 40 42 45 47 49 51 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 88 90 92 94 97 99 101 103 105 106 107 108 110 112 114 116 119 121 123 125 127 128 129 130 131 132 134 135 136

SHIFT_COORDS(0, 0, 1)
DETECTOR(2, 0, 0) rec[-48] rec[-107]#
DETECTOR(6, 0, 0) rec[-47] rec[-106]
DETECTOR(2, 2, 0) rec[-46] rec[-105]
DETECTOR(4, 2, 0) rec[-45] rec[-104]
DETECTOR(6, 2, 0) rec[-44] rec[-103]
DETECTOR(8, 2, 0) rec[-43] rec[-102]
DETECTOR(10, 2, 0) rec[-42] rec[-101] rec[-53] rec[-52]
DETECTOR(0, 4, 0) rec[-41] rec[-100]
DETECTOR(2, 4, 0) rec[-40] rec[-99]
DETECTOR(4, 4, 0) rec[-39] rec[-98]
DETECTOR(6, 4, 0) rec[-38] rec[-97]
DETECTOR(8, 4, 0) rec[-37] rec[-96]
DETECTOR(2, 6, 0) rec[-36] rec[-95]
DETECTOR(4, 6, 0) rec[-35] rec[-94]
DETECTOR(6, 6, 0) rec[-34] rec[-93]
DETECTOR(8, 6, 0) rec[-33] rec[-92]
DETECTOR(10, 6, 0) rec[-32] rec[-91] rec[-51] rec[-50]
DETECTOR(0, 8, 0) rec[-31] rec[-90]
DETECTOR(2, 8, 0) rec[-30] rec[-89]
DETECTOR(4, 8, 0) rec[-29] rec[-88]
DETECTOR(6, 8, 0) rec[-28] rec[-87]
DETECTOR(8, 8, 0) rec[-27] rec[-86]
DETECTOR(4, 10, 0) rec[-26] rec[-85]
DETECTOR(8, 10, 0) rec[-25] rec[-84]
    
DETECTOR(14, 0, 0) rec[-24] rec[-83]
DETECTOR(18, 0, 0) rec[-23] rec[-82]
DETECTOR(14, 2, 0) rec[-22] rec[-81]
DETECTOR(16, 2, 0) rec[-21] rec[-80]
DETECTOR(18, 2, 0) rec[-20] rec[-79]
DETECTOR(20, 2, 0) rec[-19] rec[-78]
DETECTOR(22, 2, 0) rec[-18] rec[-77]
DETECTOR(12, 4, 0) rec[-17] rec[-76] rec[-51] rec[-52]
DETECTOR(14, 4, 0) rec[-16] rec[-75]
DETECTOR(16, 4, 0) rec[-15] rec[-74]
DETECTOR(18, 4, 0) rec[-14] rec[-73]
DETECTOR(20, 4, 0) rec[-13] rec[-72]
DETECTOR(14, 6, 0) rec[-12] rec[-71]
DETECTOR(16, 6, 0) rec[-11] rec[-70]
DETECTOR(18, 6, 0) rec[-10] rec[-69]
DETECTOR(20, 6, 0) rec[-9] rec[-68]
DETECTOR(22, 6, 0) rec[-8] rec[-67]
DETECTOR(12, 8, 0) rec[-7] rec[-66] rec[-49] rec[-50]
DETECTOR(14, 8, 0) rec[-6] rec[-65]
DETECTOR(16, 8, 0) rec[-5] rec[-64]
DETECTOR(18, 8, 0) rec[-4] rec[-63]
DETECTOR(20, 8, 0) rec[-3] rec[-62]
DETECTOR(16, 10, 0) rec[-2] rec[-61]
DETECTOR(20, 10, 0) rec[-1] rec[-60]
    
    
REPEAT """+str(j)+"""{

    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
    DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

    TICK
    CX 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 
    DEPOLARIZE2("""+str(p_f)+""") 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52 76 77 98 99 120 121 89 90 111 112 80 81 102 103 124 125 93 94 115 116 97 96 119 118 88 87 110 109 101 100 123 122 92 91 114 113 105 104 127 126 
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 4 5 8 9 10 11 12 17 20 21 33 34 43 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 78 79 82 83 84 85 86 95 106 107 108 120 128 129 130 131 132 133 134 135 136 137

    TICK
    CX 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 
    DEPOLARIZE2("""+str(p_f)+""") 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52 76 75 98 97 120 119 89 88 111 110 80 79 102 101 124 123 93 92 115 114 86 96 108 118 77 87 99 109 90 100 112 122 81 91 103 113 94 104 116 126 
    DEPOLARIZE1("""+str(p_f1)+""") 0 4 8 10 11 26 32 33 47 51 54 55 56 57 58 59 60 61 62 63 68 74 78 82 83 84 85 95 105 106 107 117 127 128 129 130 131 132 133 134 135 136 137


    TICK
    CX 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 
    DEPOLARIZE2("""+str(p_f)+""") 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43 98 88 120 110 89 79 111 101 133 123 102 92 124 114 93 83 115 105 137 127 86 87 108 109 99 100 121 122 90 91 112 113 103 104 125 126 94 95 116 117 
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 4 6 7 8 10 11 22 23 32 33 44 45 54 55 56 57 58 60 61 62 64 69 76 77 78 80 81 82 84 85 106 107 128 129 130 131 132 134 135 136

    TICK
    CX 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 
    DEPOLARIZE2("""+str(p_f)+""") 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43 98 86 120 108 89 77 111 99 133 121 102 90 124 112 93 81 115 103 137 125 75 87 97 109 88 100 110 122 79 91 101 113 92 104 114 126 83 95 105 117 
    DEPOLARIZE1("""+str(p_f1)+""") 0 2 4 6 8 10 11 22 32 33 44 45 53 54 55 56 57 58 60 61 62 76 78 80 82 84 85 94 106 107 116 127 128 129 130 131 132 134 135 136

    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
    DEPOLARIZE1("""+str(p_f1)+""") 2 6 15 19 24 28 37 41 46 50 59 63 76 80 89 93 98 102 111 115 120 124 133 137 
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 13 14 16 17 18 20 21 22 23 25 26 27 29 30 31 32 33 34 35 36 38 39 40 42 43 44 45 47 48 49 51 52 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 87 88 90 91 92 94 95 96 97 99 100 101 103 104 105 106 107 108 109 110 112 113 114 116 117 118 119 121 122 123 125 126 127 128 129 130 131 132 134 135 136

    TICK
    MR("""+str(p_f2)+""") 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63 76 80 87 89 91 93 95 96 98 100 102 104 109 111 113 115 117 118 120 122 124 126 133 137  
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 7 8 9 10 11 12 14 16 18 20 23 25 27 29 31 32 33 34 36 38 40 42 45 47 49 51 53 54 55 56 57 58 60 61 62 64 65 66 67 68 75 77 78 79 81 82 83 84 85 86 88 90 92 94 97 99 101 103 105 106 107 108 110 112 114 116 119 121 123 125 127 128 129 130 131 132 134 135 136

    SHIFT_COORDS(0, 0, 1)
    DETECTOR(2, 0, 0) rec[-48] rec[-96]#
    DETECTOR(6, 0, 0) rec[-47] rec[-95]
    DETECTOR(2, 2, 0) rec[-46] rec[-94]
    DETECTOR(4, 2, 0) rec[-45] rec[-93]
    DETECTOR(6, 2, 0) rec[-44] rec[-92]
    DETECTOR(8, 2, 0) rec[-43] rec[-91]
    DETECTOR(10, 2, 0) rec[-42] rec[-90] 
    DETECTOR(0, 4, 0) rec[-41] rec[-89]
    DETECTOR(2, 4, 0) rec[-40] rec[-88]
    DETECTOR(4, 4, 0) rec[-39] rec[-87]
    DETECTOR(6, 4, 0) rec[-38] rec[-86]
    DETECTOR(8, 4, 0) rec[-37] rec[-85]
    DETECTOR(2, 6, 0) rec[-36] rec[-84]
    DETECTOR(4, 6, 0) rec[-35] rec[-83]
    DETECTOR(6, 6, 0) rec[-34] rec[-82]
    DETECTOR(8, 6, 0) rec[-33] rec[-81]
    DETECTOR(10, 6, 0) rec[-32] rec[-80] 
    DETECTOR(0, 8, 0) rec[-31] rec[-79]
    DETECTOR(2, 8, 0) rec[-30] rec[-78]
    DETECTOR(4, 8, 0) rec[-29] rec[-77]
    DETECTOR(6, 8, 0) rec[-28] rec[-76]
    DETECTOR(8, 8, 0) rec[-27] rec[-75]
    DETECTOR(4, 10, 0) rec[-26] rec[-74]
    DETECTOR(8, 10, 0) rec[-25] rec[-73]
    
    DETECTOR(14, 0, 0) rec[-24] rec[-72]
    DETECTOR(18, 0, 0) rec[-23] rec[-71]
    DETECTOR(14, 2, 0) rec[-22] rec[-70]
    DETECTOR(16, 2, 0) rec[-21] rec[-69]
    DETECTOR(18, 2, 0) rec[-20] rec[-68]
    DETECTOR(20, 2, 0) rec[-19] rec[-67]
    DETECTOR(22, 2, 0) rec[-18] rec[-66]
    DETECTOR(12, 4, 0) rec[-17] rec[-65] 
    DETECTOR(14, 4, 0) rec[-16] rec[-64]
    DETECTOR(16, 4, 0) rec[-15] rec[-63]
    DETECTOR(18, 4, 0) rec[-14] rec[-62]
    DETECTOR(20, 4, 0) rec[-13] rec[-61]
    DETECTOR(14, 6, 0) rec[-12] rec[-60]
    DETECTOR(16, 6, 0) rec[-11] rec[-59]
    DETECTOR(18, 6, 0) rec[-10] rec[-58]
    DETECTOR(20, 6, 0) rec[-9] rec[-57]
    DETECTOR(22, 6, 0) rec[-8] rec[-56]
    DETECTOR(12, 8, 0) rec[-7] rec[-55] 
    DETECTOR(14, 8, 0) rec[-6] rec[-54]
    DETECTOR(16, 8, 0) rec[-5] rec[-53]
    DETECTOR(18, 8, 0) rec[-4] rec[-52]
    DETECTOR(20, 8, 0) rec[-3] rec[-51]
    DETECTOR(16, 10, 0) rec[-2] rec[-50]
    DETECTOR(20, 10, 0) rec[-1] rec[-49]

}

MX("""+str(p_f2)+""") 1 3 5 7 9 12 14 16 18 20 23 25 27 29 31 34 36 38 40 42 45 47 49 51 53 75 77 79 81 83 86 88 90 92 94 97 99 101 103 105 108 110 112 114 116 119 121 123 125 127

DETECTOR(2, 0, 1) rec[-49] rec[-50] rec[-98]#
DETECTOR(2, 4, 1) rec[-39] rec[-40] rec[-44] rec[-45] rec[-90]#
DETECTOR(2, 8, 1) rec[-29] rec[-30] rec[-34] rec[-35] rec[-80]#
DETECTOR(4, 2, 1) rec[-43] rec[-44] rec[-48] rec[-49] rec[-95]##
DETECTOR(4, 6, 1) rec[-33] rec[-34] rec[-38] rec[-39] rec[-85]#
DETECTOR(4, 10, 1) rec[-28] rec[-29] rec[-76]#
DETECTOR(6, 0, 1) rec[-47] rec[-48] rec[-97]#
DETECTOR(6, 4, 1) rec[-37] rec[-38] rec[-42] rec[-43] rec[-88]#
DETECTOR(6, 8, 1) rec[-27] rec[-28] rec[-32] rec[-33] rec[-78]#
DETECTOR(8, 2, 1) rec[-41] rec[-42] rec[-46] rec[-47] rec[-93]#
DETECTOR(8, 6, 1) rec[-31] rec[-32] rec[-36] rec[-37] rec[-83]
DETECTOR(8, 10, 1) rec[-26] rec[-27] rec[-75] # 

OBSERVABLE_INCLUDE(0) rec[-30] rec[-35] rec[-40] rec[-45] rec[-50]#

DETECTOR(14, 0, 1) rec[-24] rec[-25] rec[-74]#
DETECTOR(14, 4, 1) rec[-14] rec[-15] rec[-19] rec[-20] rec[-66]#
DETECTOR(14, 8, 1) rec[-4] rec[-5] rec[-9] rec[-10] rec[-56]#
DETECTOR(16, 2, 1) rec[-18] rec[-19] rec[-23] rec[-24] rec[-71]#
DETECTOR(16, 6, 1) rec[-8] rec[-9] rec[-13] rec[-14] rec[-61]#
DETECTOR(16, 10, 1) rec[-3] rec[-4] rec[-52]#
DETECTOR(18, 0, 1) rec[-22] rec[-23] rec[-73]#
DETECTOR(18, 4, 1) rec[-12] rec[-13] rec[-17] rec[-18] rec[-64]#
DETECTOR(18, 8, 1) rec[-2] rec[-3] rec[-7] rec[-8] rec[-54]#
DETECTOR(20, 2, 1) rec[-16] rec[-17] rec[-21] rec[-22] rec[-69]#
DETECTOR(20, 6, 1) rec[-6] rec[-7] rec[-11] rec[-12] rec[-59]#
DETECTOR(20, 10, 1) rec[-1] rec[-2] rec[-51]#

OBSERVABLE_INCLUDE(1) rec[-5] rec[-10] rec[-15] rec[-20] rec[-25]

""")

   


        #d=3 unoratated Surface code
        circsurf_notr=stim.Circuit('''
    QUBIT_COORDS(0, 0) 0
    QUBIT_COORDS(1, 0) 1
    QUBIT_COORDS(2, 0) 2
    QUBIT_COORDS(3, 0) 3
    QUBIT_COORDS(4, 0) 4
    QUBIT_COORDS(0, 1) 5
    QUBIT_COORDS(1, 1) 6
    QUBIT_COORDS(2, 1) 7
    QUBIT_COORDS(3, 1) 8
    QUBIT_COORDS(4, 1) 9
    QUBIT_COORDS(0, 2) 10
    QUBIT_COORDS(1, 2) 11
    QUBIT_COORDS(2, 2) 12
    QUBIT_COORDS(3, 2) 13
    QUBIT_COORDS(4, 2) 14
    QUBIT_COORDS(0, 3) 15
    QUBIT_COORDS(1, 3) 16
    QUBIT_COORDS(2, 3) 17
    QUBIT_COORDS(3, 3) 18
    QUBIT_COORDS(4, 3) 19
    QUBIT_COORDS(0, 4) 20
    QUBIT_COORDS(1, 4) 21
    QUBIT_COORDS(2, 4) 22
    QUBIT_COORDS(3, 4) 23
    QUBIT_COORDS(4, 4) 24
    QUBIT_COORDS(5, 0) 25
    QUBIT_COORDS(5, 1) 26
    QUBIT_COORDS(5, 2) 27
    QUBIT_COORDS(5, 3) 28
    QUBIT_COORDS(5, 4) 29
    QUBIT_COORDS(6, 0) 30
    QUBIT_COORDS(7, 0) 31
    QUBIT_COORDS(8, 0) 32
    QUBIT_COORDS(9, 0) 33
    QUBIT_COORDS(10, 0) 34
    QUBIT_COORDS(6, 1) 35
    QUBIT_COORDS(7, 1) 36
    QUBIT_COORDS(8, 1) 37
    QUBIT_COORDS(9, 1) 38
    QUBIT_COORDS(10, 1) 39
    QUBIT_COORDS(6, 2) 40
    QUBIT_COORDS(7, 2) 41
    QUBIT_COORDS(8, 2) 42
    QUBIT_COORDS(9, 2) 43
    QUBIT_COORDS(10, 2) 44
    QUBIT_COORDS(6, 3) 45
    QUBIT_COORDS(7, 3) 46
    QUBIT_COORDS(8, 3) 47
    QUBIT_COORDS(9, 3) 48
    QUBIT_COORDS(10, 3) 49
    QUBIT_COORDS(6, 4) 50
    QUBIT_COORDS(7, 4) 51
    QUBIT_COORDS(8, 4) 52
    QUBIT_COORDS(9, 4) 53
    QUBIT_COORDS(10, 4) 54

    #data qubits into |+> state
    RX 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54
    TICK
    
    #x,z stabilizer qubits
    R 1 3 5 7 9 11 13 15 17 19 21 23
    R 31 33 35 37 39 41 43 45 47 49 51 53
    R 25 27 29
    TICK
    #Merging------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
    TICK
    CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 29 50 27 40 25 30 28 19 26 9
    TICK
    CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 27 28 25 26
    TICK
    CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 29 28 27 26
    TICK
    CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 29 24 27 14 25 4 26 35 28 45

    TICK
    H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53


    TICK
    MR  1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53

    DETECTOR(1, 0, 0) rec[-27]
    DETECTOR(1, 2, 0) rec[-22]
    DETECTOR(1, 4, 0) rec[-17]

    DETECTOR(3, 0, 0) rec[-26]
    DETECTOR(3, 2, 0) rec[-21]
    DETECTOR(3, 4, 0) rec[-16]

    DETECTOR(5, 0, 0) rec[-15]
    DETECTOR(5, 2, 0) rec[-14]
    DETECTOR(5, 4, 0) rec[-13]

    DETECTOR(7, 0, 0) rec[-12]
    DETECTOR(7, 2, 0) rec[-7]
    DETECTOR(7, 4, 0) rec[-2]

    DETECTOR(9, 0, 0) rec[-11]
    DETECTOR(9, 2, 0) rec[-6]
    DETECTOR(9, 4, 0) rec[-1]


    REPEAT 3 {
        
        H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''')  1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54
        TICK

        CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 29 50 27 40 25 30 28 19 26 9
        DEPOLARIZE2('''+str(p_f)+''')  1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 29 50 27 40 25 30 28 19 26 9
        TICK

        CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 27 28 25 26
        DEPOLARIZE2('''+str(p_f)+''')  1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 27 28 25 26
        DEPOLARIZE1('''+str(p_f)+''')  0 2 4 21 23 29 30 32 34 51 53


        TICK
        CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 29 28 27 26
        DEPOLARIZE2('''+str(p_f)+''')  11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 29 28 27 26
        DEPOLARIZE1('''+str(p_f)+''')  1 3 20 22 24 25 31 33 50 52 54


        TICK
        CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 29 24 27 14 25 4 26 35 28 45
        DEPOLARIZE2('''+str(p_f)+''')  1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 29 24 27 14 25 4 26 35 28 45
        DEPOLARIZE1('''+str(p_f)+''')  5 15 34 44

        TICK
        H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''')  1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''')  0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54


        TICK
        MR('''+str(p_f)+''')  1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53
        DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54

        SHIFT_COORDS(0, 0, 1)
        DETECTOR(1, 0, 0) rec[-27] rec[-54]
        DETECTOR(3, 0, 0) rec[-26] rec[-53]
        DETECTOR(0, 1, 0) rec[-25] rec[-52]
        DETECTOR(2, 1, 0) rec[-24] rec[-51]
        DETECTOR(4, 1, 0) rec[-23] rec[-50]
        DETECTOR(1, 2, 0) rec[-22] rec[-49]
        DETECTOR(3, 2, 0) rec[-21] rec[-48]
        DETECTOR(0, 3, 0) rec[-20] rec[-47]
        DETECTOR(2, 3, 0) rec[-19] rec[-46]
        DETECTOR(4, 3, 0) rec[-18] rec[-45] 
        DETECTOR(1, 4, 0) rec[-17] rec[-44]
        DETECTOR(3, 4, 0) rec[-16] rec[-43]

        DETECTOR(5, 0, 0) rec[-15] rec[-42] 
        DETECTOR(5, 2, 0) rec[-14] rec[-41]
        DETECTOR(5, 4, 0) rec[-13] rec[-40]

        DETECTOR(7, 0, 0) rec[-12] rec[-39]
        DETECTOR(9, 0, 0) rec[-11] rec[-38]
        DETECTOR(6, 1, 0) rec[-10] rec[-37]


        DETECTOR(8, 1, 0) rec[-9] rec[-36]
        DETECTOR(10, 1, 0) rec[-8] rec[-35]
        DETECTOR(7, 2, 0) rec[-7] rec[-34]
        DETECTOR(9, 2, 0) rec[-6] rec[-33]
        DETECTOR(6, 3, 0) rec[-5] rec[-32]
        DETECTOR(8, 3, 0) rec[-4] rec[-31]
        DETECTOR(10, 3, 0) rec[-3] rec[-30] 
        DETECTOR(7, 4, 0) rec[-2] rec[-29]
        DETECTOR(9, 4, 0) rec[-1] rec[-28]


    }
    ##################################################################################
    #Measurement of the data qubits in the intermediate section
    MZ('''+str(p_f)+''')  26 28    
    DEPOLARIZE1('''+str(p_f1)+''') 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54

    ##################################################################################

    R 1 3 5 7 9 11 13 15 17 19 21 23
    R 31 33 35 37 39 41 43 45 47 49 51 53

    TICK
    H 1 3 11 13 21 23 31 33 41 43 51 53
    DEPOLARIZE1('''+str(p_f1)+''')  1 3 11 13 21 23 31 33 41 43 51 53
    DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54


    TICK
    CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 
    DEPOLARIZE2('''+str(p_f)+''')  1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47
    DEPOLARIZE1('''+str(p_f1)+''') 0 9 10 19 20 25 26 29 30 39 40 49 50

    TICK
    CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 
    DEPOLARIZE2('''+str(p_f)+''')  1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49  
    DEPOLARIZE1('''+str(p_f)+''')  0 2 4 21 23 29 30 32 34 51 53
    TICK

    CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49
    DEPOLARIZE2('''+str(p_f)+''')  11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 
    DEPOLARIZE1('''+str(p_f)+''') 1 3 20 22 24 25 31 33 50 52 54
    TICK
    CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49
    DEPOLARIZE2('''+str(p_f)+''')  1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 
    DEPOLARIZE1('''+str(p_f)+''')  4 5 14 15 24 25 26 29 34 44 45



    TICK
    H 1 3 11 13 21 23 31 33 41 43 51 53
    DEPOLARIZE1('''+str(p_f1)+''')  1 3 11 13 21 23 31 33 41 43 51 53
    DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54

    TICK
    MR('''+str(p_f2)+''')  1 3 5 7 9 11 13 15 17 19 21 23 31 33 35 37 39 41 43 45 47 49 51 53
    DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 6 8 10 12 14 16 18 20 22 24 25 26 27 28 29 30 32 34 36 38 40 42 44 46 48 50 52 54



    DETECTOR(1, 0, 0) rec[-24] rec[-53]
    DETECTOR(3, 0, 0) rec[-23] rec[-52]
    DETECTOR(0, 1, 0) rec[-22] rec[-51]
    DETECTOR(2, 1, 0) rec[-21] rec[-50]
    DETECTOR(4, 1, 0) rec[-20] rec[-49] rec[-26]
    DETECTOR(1, 2, 0) rec[-19] rec[-48]
    DETECTOR(3, 2, 0) rec[-18] rec[-47]
    DETECTOR(0, 3, 0) rec[-17] rec[-46]
    DETECTOR(2, 3, 0) rec[-16] rec[-45]
    DETECTOR(4, 3, 0) rec[-15] rec[-44] rec[-25] 
    DETECTOR(1, 4, 0) rec[-14] rec[-43]
    DETECTOR(3, 4, 0) rec[-13] rec[-42]

    DETECTOR(7, 0, 0) rec[-12] rec[-38]
    DETECTOR(9, 0, 0) rec[-11] rec[-37]
    DETECTOR(6, 1, 0) rec[-10] rec[-36] rec[-26]
    DETECTOR(8, 1, 0) rec[-9] rec[-35]
    DETECTOR(10, 1, 0) rec[-8] rec[-34]
    DETECTOR(7, 2, 0) rec[-7] rec[-33]
    DETECTOR(9, 2, 0) rec[-6] rec[-32]
    DETECTOR(6, 3, 0) rec[-5] rec[-31] rec[-25]
    DETECTOR(8, 3, 0) rec[-4] rec[-30]
    DETECTOR(10, 3, 0) rec[-3] rec[-29] 
    DETECTOR(7, 4, 0) rec[-2] rec[-28]
    DETECTOR(9, 4, 0) rec[-1] rec[-27]



    #Delay for information to send----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    R 1 3 5 7 9 11 13 15 17 19 21 23
    R 31 33 35 37 39 41 43 45 47 49 51 53

    REPEAT '''+str(j)+''' {
        
        TICK
        H 1 3 11 13 21 23 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''')  1 3 11 13 21 23 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54

        
        TICK
        CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 
        DEPOLARIZE2('''+str(p_f)+''')  1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47  
        DEPOLARIZE1('''+str(p_f1)+''') 0 9 10 19 20 25 26 29 30 39 40 49 50

        TICK
        CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 
        DEPOLARIZE2('''+str(p_f)+''')  1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49  
        DEPOLARIZE1('''+str(p_f)+''')  0 2 4 21 23 29 30 32 34 51 53

        TICK
        CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 
        DEPOLARIZE2('''+str(p_f)+''')  11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 
        DEPOLARIZE1('''+str(p_f)+''') 1 3 20 22 24 25 31 33 50 52 54

        TICK
        CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 
        DEPOLARIZE2('''+str(p_f)+''')  1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49  
        DEPOLARIZE1('''+str(p_f)+''')  4 5 14 15 24 25 26 29 34 44 45

        TICK
        H 1 3 11 13 21 23 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''')  1 3 11 13 21 23 31 33 41 43 51 53
        DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54

        TICK
        MR('''+str(p_f2)+''')  1 3  5  7  9 11 13 15 17 19 21 23 31 33 35 37 39 41 43 45 47 49 51 53
        DEPOLARIZE1('''+str(p_f1)+''') 0 2 4 6 8 10 12 14 16 18 20 22 24 25 26 28 29 30 32 34 36 38 40 42 44 46 48 50 52 54


        SHIFT_COORDS(0, 0, 1)
        DETECTOR(1, 0, 0) rec[-24] rec[-48]
        DETECTOR(3, 0, 0) rec[-23] rec[-47]
        DETECTOR(0, 1, 0) rec[-22] rec[-46]
        DETECTOR(2, 1, 0) rec[-21] rec[-45]
        DETECTOR(4, 1, 0) rec[-20] rec[-44]
        DETECTOR(1, 2, 0) rec[-19] rec[-43]
        DETECTOR(3, 2, 0) rec[-18] rec[-42]
        DETECTOR(0, 3, 0) rec[-17] rec[-41]
        DETECTOR(2, 3, 0) rec[-16] rec[-40]
        DETECTOR(4, 3, 0) rec[-15] rec[-39] 
        DETECTOR(1, 4, 0) rec[-14] rec[-38]
        DETECTOR(3, 4, 0) rec[-13] rec[-37]

      

        DETECTOR(7, 0, 0) rec[-12] rec[-36]
        DETECTOR(9, 0, 0) rec[-11] rec[-35]
        DETECTOR(6, 1, 0) rec[-10] rec[-34]
        DETECTOR(8, 1, 0) rec[-9] rec[-33]
        DETECTOR(10, 1, 0) rec[-8] rec[-32]
        DETECTOR(7, 2, 0) rec[-7] rec[-31]
        DETECTOR(9, 2, 0) rec[-6] rec[-30]
        DETECTOR(6, 3, 0) rec[-5] rec[-29]
        DETECTOR(8, 3, 0) rec[-4] rec[-28]
        DETECTOR(10, 3, 0) rec[-3] rec[-27] 
        DETECTOR(7, 4, 0) rec[-2] rec[-26]
        DETECTOR(9, 4, 0) rec[-1] rec[-25]


    }

    MX('''+str(p_f2)+''')  0 2 4 6 8  10 12 14 16 18 20 22 24 30 32 34 36 38 40 42 44 46 48 50 52 54

    DETECTOR(1, 0, 1) rec[-23] rec[-25] rec[-26] rec[-50]                   
    DETECTOR(1, 2, 1) rec[-18] rec[-20] rec[-21] rec[-23] rec[-45]
    DETECTOR(1, 4, 1) rec[-15] rec[-16] rec[-18] rec[-40]
    
    DETECTOR(3, 0, 1) rec[-22] rec[-24] rec[-25] rec[-49]
    DETECTOR(3, 2, 1) rec[-17] rec[-19] rec[-20] rec[-22] rec[-44]
    DETECTOR(3, 4, 1) rec[-14] rec[-15] rec[-17] rec[-39]

    DETECTOR(7, 0, 1) rec[-10] rec[-12] rec[-13] rec[-38]
    DETECTOR(7, 2, 1) rec[-5] rec[-7] rec[-8] rec[-10] rec[-33]
    DETECTOR(7, 4, 1) rec[-2] rec[-3] rec[-5] rec[-28]

    DETECTOR(9, 0, 1) rec[-9] rec[-11] rec[-12] rec[-37]
    DETECTOR(9, 2, 1) rec[-4] rec[-6] rec[-7] rec[-9] rec[-32]
    DETECTOR(9, 4, 1) rec[-1] rec[-2] rec[-4] rec[-27]


    OBSERVABLE_INCLUDE(0) rec[-16] rec[-21] rec[-26] #rec[-3] rec[-8] rec[-13] 
    OBSERVABLE_INCLUDE(1) rec[-3] rec[-8] rec[-13] 


''')    
####################################################################################################################################################################################################################        
        
        #d=3 rotated Surface code (S)
        circsurf_notr2=stim.Circuit("""
QUBIT_COORDS(1, 1) 1
QUBIT_COORDS(2, 0) 2
QUBIT_COORDS(3, 1) 3
QUBIT_COORDS(5, 1) 5
QUBIT_COORDS(1, 3) 8
QUBIT_COORDS(2, 2) 9
QUBIT_COORDS(3, 3) 10
QUBIT_COORDS(4, 2) 11
QUBIT_COORDS(5, 3) 12
QUBIT_COORDS(6, 2) 13
QUBIT_COORDS(0, 4) 14
QUBIT_COORDS(1, 5) 15
QUBIT_COORDS(2, 4) 16
QUBIT_COORDS(3, 5) 17
QUBIT_COORDS(4, 4) 18
QUBIT_COORDS(5, 5) 19
QUBIT_COORDS(4, 6) 25

#data
QUBIT_COORDS(7, 1) 0
QUBIT_COORDS(7, 3) 4
QUBIT_COORDS(7, 5) 7

#stabilizers
QUBIT_COORDS(6, 0) 26
QUBIT_COORDS(8, 2) 30
QUBIT_COORDS(6, 4) 32
QUBIT_COORDS(8, 6) 33


QUBIT_COORDS(9, 1) 27
QUBIT_COORDS(10, 0) 28
QUBIT_COORDS(11, 1) 29
QUBIT_COORDS(13, 1) 31
QUBIT_COORDS(9, 3) 34
QUBIT_COORDS(10, 2) 35
QUBIT_COORDS(11, 3) 36
QUBIT_COORDS(12, 2) 37
QUBIT_COORDS(13, 3) 38
QUBIT_COORDS(14, 2) 39
QUBIT_COORDS(8, 4) 40
QUBIT_COORDS(9, 5) 41
QUBIT_COORDS(10, 4) 42
QUBIT_COORDS(11, 5) 43
QUBIT_COORDS(12, 4) 44
QUBIT_COORDS(13, 5) 45
QUBIT_COORDS(12, 6) 51



RX 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45 0 4 7
R 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33 
TICK
H 2 11 16 25 28 37 42 51 26 30 32 33


TICK
CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44 26 0 32 7 30 34 4 13



TICK
CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44 26 5 32 19 30 4 0 13



TICK
CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39 30 27 32 4 33 41 7 40


TICK
CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39 32 12 30 0 33 7 4 40



TICK
H 2 11 16 25 28 37 42 51 26 30 32 33


TICK
MR 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33


DETECTOR(2, 0, 0) rec[-20]
DETECTOR(2, 4, 0) rec[-15]
DETECTOR(4, 2, 0) rec[-18]
DETECTOR(4, 6, 0) rec[-13]

DETECTOR(6, 0, 0) rec[-4]
DETECTOR(8, 2, 0) rec[-3]
DETECTOR(6, 4, 0) rec[-2]
DETECTOR(8, 6, 0) rec[-1]

DETECTOR(10, 0, 0) rec[-12]  
DETECTOR(10, 4, 0) rec[-7]  
DETECTOR(12, 2, 0) rec[-10] 
DETECTOR(12, 6, 0) rec[-5] 




REPEAT 3 {

    TICK
    H 2 11 16 25 28 37 42 51 26 30 32 33
    DEPOLARIZE1("""+str(p_f1)+""") 2 11 16 25 28 37 42 51 26 30 32 33
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21 

    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44 26 0 32 7 30 34 4 13
    DEPOLARIZE2("""+str(p_f)+""") 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44 26 0 32 7 30 34 4 13
    DEPOLARIZE1("""+str(p_f1)+""") 1 5 6 8 20 21 22 23 24 25 27 31 33 39 46 47 48 49 50 51

    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44 26 5 32 19 30 4 0 13
    DEPOLARIZE2("""+str(p_f)+""") 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44 26 5 32 19 30 4 0 13
    DEPOLARIZE1("""+str(p_f1)+""") 6 7 20 21 22 23 24 25 31 33 43 45 46 47 48 49 50 51


    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39 30 27 32 4 33 41 7 40
    DEPOLARIZE2("""+str(p_f)+""") 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39 30 27 32 4 33 41 7 40
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 6 14 15 20 21 22 23 24 26 28 29 46 47 48 49 50

    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39 32 12 30 0 33 7 4 40
    DEPOLARIZE2("""+str(p_f)+""") 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39 32 12 30 0 33 7 4 40
    DEPOLARIZE1("""+str(p_f1)+""") 2 6 14 15 20 21 22 23 24 26 28 38 41 45 46 47 48 49 50

    TICK
    H 2 11 16 25 28 37 42 51 26 30 32 33
    DEPOLARIZE1("""+str(p_f1)+""") 2 11 16 25 28 37 42 51
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21 

    TICK
    MR("""+str(p_f2)+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21 
    SHIFT_COORDS(0, 0, 1)

    DETECTOR(2, 0, 0) rec[-20] rec[-40]
    DETECTOR(2, 2, 0) rec[-19] rec[-39]
    DETECTOR(4, 2, 0) rec[-18] rec[-38]
    DETECTOR(6, 2, 0) rec[-17] rec[-37]
    DETECTOR(0, 4, 0) rec[-16] rec[-36]
    DETECTOR(2, 4, 0) rec[-15] rec[-35]
    DETECTOR(4, 4, 0) rec[-14] rec[-34]
    DETECTOR(4, 6, 0) rec[-13] rec[-33]

    DETECTOR(10, 0, 0) rec[-12] rec[-32]   
    DETECTOR(10, 2, 0) rec[-11] rec[-31]   
    DETECTOR(12, 2, 0) rec[-10] rec[-30]  
    DETECTOR(14, 2, 0) rec[-9] rec[-29]  
    DETECTOR(8, 4, 0) rec[-8] rec[-28]   
    DETECTOR(10, 4, 0) rec[-7] rec[-27]   
    DETECTOR(12, 4, 0) rec[-6] rec[-26]  
    DETECTOR(12, 6, 0) rec[-5] rec[-25]

    DETECTOR(6, 0, 0) rec[-4] rec[-24]   
    DETECTOR(8, 2, 0) rec[-3] rec[-23]   
    DETECTOR(6, 4, 0) rec[-2] rec[-22]  
    DETECTOR(8, 6, 0) rec[-1] rec[-21]

}

#################### splitting ################################################################
MZ("""+str(p_f)+""") 0 4 7
DEPOLARIZE1("""+str(p_f1)+""") 1 2 3 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51


#Reset auxiliary qubits
R 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 
#

TICK
H 2 11 16 25 28 37 42 51
DEPOLARIZE1("""+str(p_f1)+""") 2 11 16 25 28 37 42 51
DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21  

TICK
CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
DEPOLARIZE2("""+str(p_f)+""") 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
DEPOLARIZE1("""+str(p_f1)+""") 0 1 4 5 6 7 8 13 20 21 22 23 24 25 26 27 30 31 32 33 34 46 47 48 49 50 51
TICK
CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44
DEPOLARIZE2("""+str(p_f)+""") 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44
DEPOLARIZE1("""+str(p_f1)+""") 6 7 20 21 22 23 24 25 31 33 43 45 46 47 48 49 50 51

TICK
CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39
DEPOLARIZE2("""+str(p_f)+""") 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39
DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 6 14 15 20 21 22 23 24 26 28 29 46 47 48 49 50

TICK
CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39
DEPOLARIZE2("""+str(p_f)+""") 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39
DEPOLARIZE1("""+str(p_f1)+""") 2 6 14 15 20 21 22 23 24 26 28 38 41 45 46 47 48 49 50

TICK
H 2 11 16 25 28 37 42 51
DEPOLARIZE1("""+str(p_f1)+""") 2 11 16 25 28 37 42 51
DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21  

TICK
MR("""+str(p_f2)+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51
DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21 
SHIFT_COORDS(0, 0, 1)

DETECTOR(2, 0, 0) rec[-16] rec[-39]
DETECTOR(2, 2, 0) rec[-15] rec[-38]
DETECTOR(4, 2, 0) rec[-14] rec[-37]
DETECTOR(6, 2, 0) rec[-13] rec[-36] rec[-18] rec[-19] 
DETECTOR(0, 4, 0) rec[-12] rec[-35]
DETECTOR(2, 4, 0) rec[-11] rec[-34]
DETECTOR(4, 4, 0) rec[-10] rec[-33]
DETECTOR(4, 6, 0) rec[-9] rec[-32]

DETECTOR(10, 0, 0) rec[-8] rec[-31]   
DETECTOR(10, 2, 0) rec[-7] rec[-30]   
DETECTOR(12, 2, 0) rec[-6] rec[-29]  
DETECTOR(14, 2, 0) rec[-5] rec[-28]  
DETECTOR(8, 4, 0) rec[-4] rec[-27] rec[-17] rec[-18] 
DETECTOR(10, 4, 0) rec[-3] rec[-26]   
DETECTOR(12, 4, 0) rec[-2] rec[-25]  
DETECTOR(12, 6, 0) rec[-1] rec[-24]

#Delay for information to send----------------------------------------------------------------------------------------------------------------------------------------------------------------------

R 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 

REPEAT """+str(j)+""" {
    TICK
    H 2 11 16 25 28 37 42 51
    DEPOLARIZE1("""+str(p_f1)+""") 2 11 16 25 28 37 42 51
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21  
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
    DEPOLARIZE2("""+str(p_f)+""") 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 4 5 6 7 8 13 20 21 22 23 24 25 26 27 30 31 32 33 34 46 47 48 49 50 51

    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44
    DEPOLARIZE2("""+str(p_f)+""") 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44
    DEPOLARIZE1("""+str(p_f1)+""") 6 7 20 21 22 23 24 25 31 33 43 45 46 47 48 49 50 51

    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39
    DEPOLARIZE2("""+str(p_f)+""") 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 2 3 6 14 15 20 21 22 23 24 26 28 29 46 47 48 49 50

    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39
    DEPOLARIZE2("""+str(p_f)+""") 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39
    DEPOLARIZE1("""+str(p_f1)+""") 2 6 14 15 20 21 22 23 24 26 28 38 41 45 46 47 48 49 50

    TICK
    H 2 11 16 25 28 37 42 51
    DEPOLARIZE1("""+str(p_f1)+""") 2 11 16 25 28 37 42 51
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21  
    TICK

    MR("""+str(p_f2)+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51
    DEPOLARIZE1("""+str(p_f1)+""") 0 1 3 4 5 6 7 8 10 12 15 17 19 20 21 22 23 24 26 27 29 30 31 32 33 34 36 38 41 43 45 46 47 48 49 50 52 53 54

    SHIFT_COORDS(0, 0, 1)
    DETECTOR(2, 0, 0) rec[-16] rec[-32]
    DETECTOR(2, 2, 0) rec[-15] rec[-31]
    DETECTOR(4, 2, 0) rec[-14] rec[-30]
    DETECTOR(6, 2, 0) rec[-13] rec[-29]
    DETECTOR(0, 4, 0) rec[-12] rec[-28]
    DETECTOR(2, 4, 0) rec[-11] rec[-27]
    DETECTOR(4, 4, 0) rec[-10] rec[-26]
    DETECTOR(4, 6, 0) rec[-9] rec[-25]

    DETECTOR(10, 0, 0) rec[-8] rec[-24]   
    DETECTOR(10, 2, 0) rec[-7] rec[-23]   
    DETECTOR(12, 2, 0) rec[-6] rec[-22]  
    DETECTOR(14, 2, 0) rec[-5] rec[-21]  
    DETECTOR(8, 4, 0) rec[-4] rec[-20]   
    DETECTOR(10, 4, 0) rec[-3] rec[-19]   
    DETECTOR(12, 4, 0) rec[-2] rec[-18]  
    DETECTOR(12, 6, 0) rec[-1] rec[-17]
}

MX("""+str(p_f2)+""") 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45

DETECTOR(2, 0, 1) rec[-17] rec[-18] rec[-34]#
DETECTOR(2, 4, 1) rec[-11] rec[-12] rec[-14] rec[-15] rec[-29]#

DETECTOR(4, 2, 1) rec[-13] rec[-14] rec[-16] rec[-17] rec[-32]
DETECTOR(4, 6, 1) rec[-10] rec[-11] rec[-27] #


DETECTOR(10, 0, 1) rec[-8] rec[-9] rec[-26] #
DETECTOR(10, 4, 1) rec[-2] rec[-3] rec[-5] rec[-6] rec[-21] #

DETECTOR(12, 2, 1) rec[-4] rec[-5] rec[-7] rec[-8] rec[-24]#
DETECTOR(12, 6, 1) rec[-1] rec[-2] rec[-19] #

OBSERVABLE_INCLUDE(0) rec[-12] rec[-15] rec[-18] #rec[-3] rec[-6] rec[-9] 

OBSERVABLE_INCLUDE(1) rec[-3] rec[-6] rec[-9] 
               
""")
#######################################################################################################################################################################################################################################
        
        #Unencoded Bell states
        circuit2=stim.Circuit("""
            QUBIT_COORDS(0,1) 0
            QUBIT_COORDS(0,2) 1

            H 0
            DEPOLARIZE1("""+str(p_f)+""") 0
            DEPOLARIZE1("""+str(p_f)+""") 1

            CX 0 1
            DEPOLARIZE2("""+str(p_f)+""") 0 1


            #Travel distance

          

            TICK
            M("""+str(p_f2)+""") 0 1
            DETECTOR rec[-1] rec[-2]

            """)

        #d=5 Bacon-Shor code (BS) ########################################################################
        circbs5=stim.Circuit(""" QUBIT_COORDS(0, 0) 0
QUBIT_COORDS(0, 1) 1
QUBIT_COORDS(0, 2) 2
QUBIT_COORDS(0, 3) 3
QUBIT_COORDS(0, 4) 4
QUBIT_COORDS(0, 5) 5
QUBIT_COORDS(0, 6) 6
QUBIT_COORDS(0, 7) 7
QUBIT_COORDS(0, 8) 8
QUBIT_COORDS(1, 0) 9
QUBIT_COORDS(1, 1) 10
QUBIT_COORDS(1, 2) 11
QUBIT_COORDS(1, 3) 12
QUBIT_COORDS(1, 4) 13
QUBIT_COORDS(1, 5) 14
QUBIT_COORDS(1, 6) 15
QUBIT_COORDS(1, 7) 16
QUBIT_COORDS(1, 8) 17
QUBIT_COORDS(2, 0) 18
QUBIT_COORDS(2, 1) 19
QUBIT_COORDS(2, 2) 20
QUBIT_COORDS(2, 3) 21
QUBIT_COORDS(2, 4) 22
QUBIT_COORDS(2, 5) 23
QUBIT_COORDS(2, 6) 24
QUBIT_COORDS(2, 7) 25
QUBIT_COORDS(2, 8) 26
QUBIT_COORDS(3, 0) 27
QUBIT_COORDS(3, 1) 28
QUBIT_COORDS(3, 2) 29
QUBIT_COORDS(3, 3) 30
QUBIT_COORDS(3, 4) 31
QUBIT_COORDS(3, 5) 32
QUBIT_COORDS(3, 6) 33
QUBIT_COORDS(3, 7) 34
QUBIT_COORDS(3, 8) 35
QUBIT_COORDS(4, 0) 36
QUBIT_COORDS(4, 1) 37
QUBIT_COORDS(4, 2) 38
QUBIT_COORDS(4, 3) 39
QUBIT_COORDS(4, 4) 40
QUBIT_COORDS(4, 5) 41
QUBIT_COORDS(4, 6) 42
QUBIT_COORDS(4, 7) 43
QUBIT_COORDS(4, 8) 44
QUBIT_COORDS(5, 0) 45
QUBIT_COORDS(5, 1) 46
QUBIT_COORDS(5, 2) 47
QUBIT_COORDS(5, 3) 48
QUBIT_COORDS(5, 4) 49
QUBIT_COORDS(5, 5) 50
QUBIT_COORDS(5, 6) 51
QUBIT_COORDS(5, 7) 52
QUBIT_COORDS(5, 8) 53
QUBIT_COORDS(6, 0) 54
QUBIT_COORDS(6, 1) 55
QUBIT_COORDS(6, 2) 56
QUBIT_COORDS(6, 3) 57
QUBIT_COORDS(6, 4) 58
QUBIT_COORDS(6, 5) 59
QUBIT_COORDS(6, 6) 60
QUBIT_COORDS(6, 7) 61
QUBIT_COORDS(6, 8) 62
QUBIT_COORDS(7, 0) 63
QUBIT_COORDS(7, 1) 64
QUBIT_COORDS(7, 2) 65
QUBIT_COORDS(7, 3) 66
QUBIT_COORDS(7, 4) 67
QUBIT_COORDS(7, 5) 68
QUBIT_COORDS(7, 6) 69
QUBIT_COORDS(7, 7) 70
QUBIT_COORDS(7, 8) 71
QUBIT_COORDS(8, 0) 72
QUBIT_COORDS(8, 1) 73
QUBIT_COORDS(8, 2) 74
QUBIT_COORDS(8, 3) 75
QUBIT_COORDS(8, 4) 76
QUBIT_COORDS(8, 5) 77
QUBIT_COORDS(8, 6) 78
QUBIT_COORDS(8, 7) 79
QUBIT_COORDS(8, 8) 80
QUBIT_COORDS(9, 0) 81
QUBIT_COORDS(9, 1) 82
QUBIT_COORDS(9, 2) 83
QUBIT_COORDS(9, 3) 84
QUBIT_COORDS(9, 4) 85
QUBIT_COORDS(9, 5) 86
QUBIT_COORDS(9, 6) 87
QUBIT_COORDS(9, 7) 88
QUBIT_COORDS(9, 8) 89
QUBIT_COORDS(10, 0) 90
QUBIT_COORDS(10, 1) 91
QUBIT_COORDS(10, 2) 92
QUBIT_COORDS(10, 3) 93
QUBIT_COORDS(10, 4) 94
QUBIT_COORDS(10, 5) 95
QUBIT_COORDS(10, 6) 96
QUBIT_COORDS(10, 7) 97
QUBIT_COORDS(10, 8) 98
QUBIT_COORDS(11, 0) 99
QUBIT_COORDS(11, 1) 100
QUBIT_COORDS(11, 2) 101
QUBIT_COORDS(11, 3) 102
QUBIT_COORDS(11, 4) 103
QUBIT_COORDS(11, 5) 104
QUBIT_COORDS(11, 6) 105
QUBIT_COORDS(11, 7) 106
QUBIT_COORDS(11, 8) 107
QUBIT_COORDS(12, 0) 108
QUBIT_COORDS(12, 1) 109
QUBIT_COORDS(12, 2) 110
QUBIT_COORDS(12, 3) 111
QUBIT_COORDS(12, 4) 112
QUBIT_COORDS(12, 5) 113
QUBIT_COORDS(12, 6) 114
QUBIT_COORDS(12, 7) 115
QUBIT_COORDS(12, 8) 116
QUBIT_COORDS(13, 0) 117
QUBIT_COORDS(13, 1) 118
QUBIT_COORDS(13, 2) 119
QUBIT_COORDS(13, 3) 120
QUBIT_COORDS(13, 4) 121
QUBIT_COORDS(13, 5) 122
QUBIT_COORDS(13, 6) 123
QUBIT_COORDS(13, 7) 124
QUBIT_COORDS(13, 8) 125
QUBIT_COORDS(14, 0) 126
QUBIT_COORDS(14, 1) 127
QUBIT_COORDS(14, 2) 128
QUBIT_COORDS(14, 3) 129
QUBIT_COORDS(14, 4) 130
QUBIT_COORDS(14, 5) 131
QUBIT_COORDS(14, 6) 132
QUBIT_COORDS(14, 7) 133
QUBIT_COORDS(14, 8) 134
QUBIT_COORDS(15, 0) 135
QUBIT_COORDS(15, 1) 136
QUBIT_COORDS(15, 2) 137
QUBIT_COORDS(15, 3) 138
QUBIT_COORDS(15, 4) 139
QUBIT_COORDS(15, 5) 140
QUBIT_COORDS(15, 6) 141
QUBIT_COORDS(15, 7) 142
QUBIT_COORDS(15, 8) 143
QUBIT_COORDS(16, 0) 144
QUBIT_COORDS(16, 1) 145
QUBIT_COORDS(16, 2) 146
QUBIT_COORDS(16, 3) 147
QUBIT_COORDS(16, 4) 148
QUBIT_COORDS(16, 5) 149
QUBIT_COORDS(16, 6) 150
QUBIT_COORDS(16, 7) 151
QUBIT_COORDS(16, 8) 152
QUBIT_COORDS(17, 0) 153
QUBIT_COORDS(17, 1) 154
QUBIT_COORDS(17, 2) 155
QUBIT_COORDS(17, 3) 156
QUBIT_COORDS(17, 4) 157
QUBIT_COORDS(17, 5) 158
QUBIT_COORDS(17, 6) 159
QUBIT_COORDS(17, 7) 160
QUBIT_COORDS(17, 8) 161
RX 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161
TICK
MPP X0*X9 X1*X10 X2*X11 X3*X12 X4*X13 X5*X14 X6*X15 X7*X16 X8*X17 X18*X27 X19*X28 X20*X29 X21*X30 X22*X31 X23*X32 X24*X33 X25*X34 X26*X35 X36*X45 X37*X46 X38*X47 X39*X48 X40*X49 X41*X50 X42*X51 X43*X52 X44*X53 X54*X63 X55*X64 X56*X65 X57*X66 X58*X67 X59*X68 X60*X69 X61*X70 X62*X71 X90*X99 X91*X100 X92*X101 X93*X102 X94*X103 X95*X104 X96*X105 X97*X106 X98*X107 X108*X117 X109*X118 X110*X119 X111*X120 X112*X121 X113*X122 X114*X123 X115*X124 X116*X125 X126*X135 X127*X136 X128*X137 X129*X138 X130*X139 X131*X140 X132*X141 X133*X142 X134*X143 X144*X153 X145*X154 X146*X155 X147*X156 X148*X157 X149*X158 X150*X159 X151*X160 X152*X161
TICK
MPP X9*X18 X10*X19 X11*X20 X12*X21 X13*X22 X14*X23 X15*X24 X16*X25 X17*X26 X27*X36 X28*X37 X29*X38 X30*X39 X31*X40 X32*X41 X33*X42 X34*X43 X35*X44 X45*X54 X46*X55 X47*X56 X48*X57 X49*X58 X50*X59 X51*X60 X52*X61 X53*X62 X63*X72 X64*X73 X65*X74 X66*X75 X67*X76 X68*X77 X69*X78 X70*X79 X71*X80 X81*X90 X82*X91 X83*X92 X84*X93 X85*X94 X86*X95 X87*X96 X88*X97 X89*X98 X99*X108 X100*X109 X101*X110 X102*X111 X103*X112 X104*X113 X105*X114 X106*X115 X107*X116 X117*X126 X118*X127 X119*X128 X120*X129 X121*X130 X122*X131 X123*X132 X124*X133 X125*X134 X135*X144 X136*X145 X137*X146 X138*X147 X139*X148 X140*X149 X141*X150 X142*X151 X143*X152
TICK
MPP Z0*Z1 Z2*Z3 Z4*Z5 Z6*Z7 Z9*Z10 Z11*Z12 Z13*Z14 Z15*Z16 Z18*Z19 Z20*Z21 Z22*Z23 Z24*Z25 Z27*Z28 Z29*Z30 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48 Z49*Z50 Z51*Z52 Z54*Z55 Z56*Z57 Z58*Z59 Z60*Z61 Z63*Z64 Z65*Z66 Z67*Z68 Z69*Z70 Z72*Z73 Z74*Z75 Z76*Z77 Z78*Z79 Z81*Z82 Z83*Z84 Z85*Z86 Z87*Z88 Z90*Z91 Z92*Z93 Z94*Z95 Z96*Z97 Z99*Z100 Z101*Z102 Z103*Z104 Z105*Z106 Z108*Z109 Z110*Z111 Z112*Z113 Z114*Z115 Z117*Z118 Z119*Z120 Z121*Z122 Z123*Z124 Z126*Z127 Z128*Z129 Z130*Z131 Z132*Z133 Z135*Z136 Z137*Z138 Z139*Z140 Z141*Z142 Z144*Z145 Z146*Z147 Z148*Z149 Z150*Z151 Z153*Z154 Z155*Z156 Z157*Z158 Z159*Z160
TICK
MPP Z1*Z2 Z3*Z4 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z14*Z15 Z16*Z17 Z19*Z20 Z21*Z22 Z23*Z24 Z25*Z26 Z28*Z29 Z30*Z31 Z32*Z33 Z34*Z35 Z37*Z38 Z39*Z40 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49 Z50*Z51 Z52*Z53 Z55*Z56 Z57*Z58 Z59*Z60 Z61*Z62 Z64*Z65 Z66*Z67 Z68*Z69 Z70*Z71 Z73*Z74 Z75*Z76 Z77*Z78 Z79*Z80 Z82*Z83 Z84*Z85 Z86*Z87 Z88*Z89 Z91*Z92 Z93*Z94 Z95*Z96 Z97*Z98 Z100*Z101 Z102*Z103 Z104*Z105 Z106*Z107 Z109*Z110 Z111*Z112 Z113*Z114 Z115*Z116 Z118*Z119 Z120*Z121 Z122*Z123 Z124*Z125 Z127*Z128 Z129*Z130 Z131*Z132 Z133*Z134 Z136*Z137 Z138*Z139 Z140*Z141 Z142*Z143 Z145*Z146 Z147*Z148 Z149*Z150 Z151*Z152 Z154*Z155 Z156*Z157 Z158*Z159 Z160*Z161
DETECTOR(0.5, 0, 0) rec[-288]
DETECTOR(0.5, 1, 0) rec[-287]
DETECTOR(0.5, 2, 0) rec[-286]
DETECTOR(0.5, 3, 0) rec[-285]
DETECTOR(0.5, 4, 0) rec[-284]
DETECTOR(0.5, 5, 0) rec[-283]
DETECTOR(0.5, 6, 0) rec[-282]
DETECTOR(0.5, 7, 0) rec[-281]
DETECTOR(0.5, 8, 0) rec[-280]
DETECTOR(1.5, 0, 0) rec[-216]
DETECTOR(1.5, 1, 0) rec[-215]
DETECTOR(1.5, 2, 0) rec[-214]
DETECTOR(1.5, 3, 0) rec[-213]
DETECTOR(1.5, 4, 0) rec[-212]
DETECTOR(1.5, 5, 0) rec[-211]
DETECTOR(1.5, 6, 0) rec[-210]
DETECTOR(1.5, 7, 0) rec[-209]
DETECTOR(1.5, 8, 0) rec[-208]
DETECTOR(2.5, 0, 0) rec[-279]
DETECTOR(2.5, 1, 0) rec[-278]
DETECTOR(2.5, 2, 0) rec[-277]
DETECTOR(2.5, 3, 0) rec[-276]
DETECTOR(2.5, 4, 0) rec[-275]
DETECTOR(2.5, 5, 0) rec[-274]
DETECTOR(2.5, 6, 0) rec[-273]
DETECTOR(2.5, 7, 0) rec[-272]
DETECTOR(2.5, 8, 0) rec[-271]
DETECTOR(3.5, 0, 0) rec[-207]
DETECTOR(3.5, 1, 0) rec[-206]
DETECTOR(3.5, 2, 0) rec[-205]
DETECTOR(3.5, 3, 0) rec[-204]
DETECTOR(3.5, 4, 0) rec[-203]
DETECTOR(3.5, 5, 0) rec[-202]
DETECTOR(3.5, 6, 0) rec[-201]
DETECTOR(3.5, 7, 0) rec[-200]
DETECTOR(3.5, 8, 0) rec[-199]
DETECTOR(4.5, 0, 0) rec[-270]
DETECTOR(4.5, 1, 0) rec[-269]
DETECTOR(4.5, 2, 0) rec[-268]
DETECTOR(4.5, 3, 0) rec[-267]
DETECTOR(4.5, 4, 0) rec[-266]
DETECTOR(4.5, 5, 0) rec[-265]
DETECTOR(4.5, 6, 0) rec[-264]
DETECTOR(4.5, 7, 0) rec[-263]
DETECTOR(4.5, 8, 0) rec[-262]
DETECTOR(5.5, 0, 0) rec[-198]
DETECTOR(5.5, 1, 0) rec[-197]
DETECTOR(5.5, 2, 0) rec[-196]
DETECTOR(5.5, 3, 0) rec[-195]
DETECTOR(5.5, 4, 0) rec[-194]
DETECTOR(5.5, 5, 0) rec[-193]
DETECTOR(5.5, 6, 0) rec[-192]
DETECTOR(5.5, 7, 0) rec[-191]
DETECTOR(5.5, 8, 0) rec[-190]
DETECTOR(6.5, 0, 0) rec[-261]
DETECTOR(6.5, 1, 0) rec[-260]
DETECTOR(6.5, 2, 0) rec[-259]
DETECTOR(6.5, 3, 0) rec[-258]
DETECTOR(6.5, 4, 0) rec[-257]
DETECTOR(6.5, 5, 0) rec[-256]
DETECTOR(6.5, 6, 0) rec[-255]
DETECTOR(6.5, 7, 0) rec[-254]
DETECTOR(6.5, 8, 0) rec[-253]
DETECTOR(7.5, 0, 0) rec[-189]
DETECTOR(7.5, 1, 0) rec[-188]
DETECTOR(7.5, 2, 0) rec[-187]
DETECTOR(7.5, 3, 0) rec[-186]
DETECTOR(7.5, 4, 0) rec[-185]
DETECTOR(7.5, 5, 0) rec[-184]
DETECTOR(7.5, 6, 0) rec[-183]
DETECTOR(7.5, 7, 0) rec[-182]
DETECTOR(7.5, 8, 0) rec[-181]
DETECTOR(9.5, 0, 0) rec[-180]
DETECTOR(9.5, 1, 0) rec[-179]
DETECTOR(9.5, 2, 0) rec[-178]
DETECTOR(9.5, 3, 0) rec[-177]
DETECTOR(9.5, 4, 0) rec[-176]
DETECTOR(9.5, 5, 0) rec[-175]
DETECTOR(9.5, 6, 0) rec[-174]
DETECTOR(9.5, 7, 0) rec[-173]
DETECTOR(9.5, 8, 0) rec[-172]
DETECTOR(10.5, 0, 0) rec[-252]
DETECTOR(10.5, 1, 0) rec[-251]
DETECTOR(10.5, 2, 0) rec[-250]
DETECTOR(10.5, 3, 0) rec[-249]
DETECTOR(10.5, 4, 0) rec[-248]
DETECTOR(10.5, 5, 0) rec[-247]
DETECTOR(10.5, 6, 0) rec[-246]
DETECTOR(10.5, 7, 0) rec[-245]
DETECTOR(10.5, 8, 0) rec[-244]
DETECTOR(11.5, 0, 0) rec[-171]
DETECTOR(11.5, 1, 0) rec[-170]
DETECTOR(11.5, 2, 0) rec[-169]
DETECTOR(11.5, 3, 0) rec[-168]
DETECTOR(11.5, 4, 0) rec[-167]
DETECTOR(11.5, 5, 0) rec[-166]
DETECTOR(11.5, 6, 0) rec[-165]
DETECTOR(11.5, 7, 0) rec[-164]
DETECTOR(11.5, 8, 0) rec[-163]
DETECTOR(12.5, 0, 0) rec[-243]
DETECTOR(12.5, 1, 0) rec[-242]
DETECTOR(12.5, 2, 0) rec[-241]
DETECTOR(12.5, 3, 0) rec[-240]
DETECTOR(12.5, 4, 0) rec[-239]
DETECTOR(12.5, 5, 0) rec[-238]
DETECTOR(12.5, 6, 0) rec[-237]
DETECTOR(12.5, 7, 0) rec[-236]
DETECTOR(12.5, 8, 0) rec[-235]
DETECTOR(13.5, 0, 0) rec[-162]
DETECTOR(13.5, 1, 0) rec[-161]
DETECTOR(13.5, 2, 0) rec[-160]
DETECTOR(13.5, 3, 0) rec[-159]
DETECTOR(13.5, 4, 0) rec[-158]
DETECTOR(13.5, 5, 0) rec[-157]
DETECTOR(13.5, 6, 0) rec[-156]
DETECTOR(13.5, 7, 0) rec[-155]
DETECTOR(13.5, 8, 0) rec[-154]
DETECTOR(14.5, 0, 0) rec[-234]
DETECTOR(14.5, 1, 0) rec[-233]
DETECTOR(14.5, 2, 0) rec[-232]
DETECTOR(14.5, 3, 0) rec[-231]
DETECTOR(14.5, 4, 0) rec[-230]
DETECTOR(14.5, 5, 0) rec[-229]
DETECTOR(14.5, 6, 0) rec[-228]
DETECTOR(14.5, 7, 0) rec[-227]
DETECTOR(14.5, 8, 0) rec[-226]
DETECTOR(15.5, 0, 0) rec[-153]
DETECTOR(15.5, 1, 0) rec[-152]
DETECTOR(15.5, 2, 0) rec[-151]
DETECTOR(15.5, 3, 0) rec[-150]
DETECTOR(15.5, 4, 0) rec[-149]
DETECTOR(15.5, 5, 0) rec[-148]
DETECTOR(15.5, 6, 0) rec[-147]
DETECTOR(15.5, 7, 0) rec[-146]
DETECTOR(15.5, 8, 0) rec[-145]
DETECTOR(16.5, 0, 0) rec[-225]
DETECTOR(16.5, 1, 0) rec[-224]
DETECTOR(16.5, 2, 0) rec[-223]
DETECTOR(16.5, 3, 0) rec[-222]
DETECTOR(16.5, 4, 0) rec[-221]
DETECTOR(16.5, 5, 0) rec[-220]
DETECTOR(16.5, 6, 0) rec[-219]
DETECTOR(16.5, 7, 0) rec[-218]
DETECTOR(16.5, 8, 0) rec[-217]
SHIFT_COORDS(0, 0, 1)
TICK
MPP X0*X9 X1*X10 X2*X11 X3*X12 X4*X13 X5*X14 X6*X15 X7*X16 X8*X17 X18*X27 X19*X28 X20*X29 X21*X30 X22*X31 X23*X32 X24*X33 X25*X34 X26*X35 X36*X45 X37*X46 X38*X47 X39*X48 X40*X49 X41*X50 X42*X51 X43*X52 X44*X53 X54*X63 X55*X64 X56*X65 X57*X66 X58*X67 X59*X68 X60*X69 X61*X70 X62*X71 X72*X81 X73*X82 X74*X83 X75*X84 X76*X85 X77*X86 X78*X87 X79*X88 X80*X89 X90*X99 X91*X100 X92*X101 X93*X102 X94*X103 X95*X104 X96*X105 X97*X106 X98*X107 X108*X117 X109*X118 X110*X119 X111*X120 X112*X121 X113*X122 X114*X123 X115*X124 X116*X125 X126*X135 X127*X136 X128*X137 X129*X138 X130*X139 X131*X140 X132*X141 X133*X142 X134*X143 X144*X153 X145*X154 X146*X155 X147*X156 X148*X157 X149*X158 X150*X159 X151*X160 X152*X161
TICK
MPP X9*X18 X10*X19 X11*X20 X12*X21 X13*X22 X14*X23 X15*X24 X16*X25 X17*X26 X27*X36 X28*X37 X29*X38 X30*X39 X31*X40 X32*X41 X33*X42 X34*X43 X35*X44 X45*X54 X46*X55 X47*X56 X48*X57 X49*X58 X50*X59 X51*X60 X52*X61 X53*X62 X63*X72 X64*X73 X65*X74 X66*X75 X67*X76 X68*X77 X69*X78 X70*X79 X71*X80 X81*X90 X82*X91 X83*X92 X84*X93 X85*X94 X86*X95 X87*X96 X88*X97 X89*X98 X99*X108 X100*X109 X101*X110 X102*X111 X103*X112 X104*X113 X105*X114 X106*X115 X107*X116 X117*X126 X118*X127 X119*X128 X120*X129 X121*X130 X122*X131 X123*X132 X124*X133 X125*X134 X135*X144 X136*X145 X137*X146 X138*X147 X139*X148 X140*X149 X141*X150 X142*X151 X143*X152
TICK
MPP Z0*Z1 Z2*Z3 Z4*Z5 Z6*Z7 Z9*Z10 Z11*Z12 Z13*Z14 Z15*Z16 Z18*Z19 Z20*Z21 Z22*Z23 Z24*Z25 Z27*Z28 Z29*Z30 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48 Z49*Z50 Z51*Z52 Z54*Z55 Z56*Z57 Z58*Z59 Z60*Z61 Z63*Z64 Z65*Z66 Z67*Z68 Z69*Z70 Z72*Z73 Z74*Z75 Z76*Z77 Z78*Z79 Z81*Z82 Z83*Z84 Z85*Z86 Z87*Z88 Z90*Z91 Z92*Z93 Z94*Z95 Z96*Z97 Z99*Z100 Z101*Z102 Z103*Z104 Z105*Z106 Z108*Z109 Z110*Z111 Z112*Z113 Z114*Z115 Z117*Z118 Z119*Z120 Z121*Z122 Z123*Z124 Z126*Z127 Z128*Z129 Z130*Z131 Z132*Z133 Z135*Z136 Z137*Z138 Z139*Z140 Z141*Z142 Z144*Z145 Z146*Z147 Z148*Z149 Z150*Z151 Z153*Z154 Z155*Z156 Z157*Z158 Z159*Z160
TICK
MPP Z1*Z2 Z3*Z4 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z14*Z15 Z16*Z17 Z19*Z20 Z21*Z22 Z23*Z24 Z25*Z26 Z28*Z29 Z30*Z31 Z32*Z33 Z34*Z35 Z37*Z38 Z39*Z40 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49 Z50*Z51 Z52*Z53 Z55*Z56 Z57*Z58 Z59*Z60 Z61*Z62 Z64*Z65 Z66*Z67 Z68*Z69 Z70*Z71 Z73*Z74 Z75*Z76 Z77*Z78 Z79*Z80 Z82*Z83 Z84*Z85 Z86*Z87 Z88*Z89 Z91*Z92 Z93*Z94 Z95*Z96 Z97*Z98 Z100*Z101 Z102*Z103 Z104*Z105 Z106*Z107 Z109*Z110 Z111*Z112 Z113*Z114 Z115*Z116 Z118*Z119 Z120*Z121 Z122*Z123 Z124*Z125 Z127*Z128 Z129*Z130 Z131*Z132 Z133*Z134 Z136*Z137 Z138*Z139 Z140*Z141 Z142*Z143 Z145*Z146 Z147*Z148 Z149*Z150 Z151*Z152 Z154*Z155 Z156*Z157 Z158*Z159 Z160*Z161

SHIFT_COORDS(0, 0, 1)

DETECTOR(0.5, 0, 0) rec[-585] rec[-584] rec[-583] rec[-582] rec[-581] rec[-580] rec[-579] rec[-578] rec[-577] rec[-297] rec[-296] rec[-295] rec[-294] rec[-293] rec[-292] rec[-291] rec[-290] rec[-289]
DETECTOR(1.5, 0, 0) rec[-513] rec[-512] rec[-511] rec[-510] rec[-509] rec[-508] rec[-507] rec[-506] rec[-505] rec[-216] rec[-215] rec[-214] rec[-213] rec[-212] rec[-211] rec[-210] rec[-209] rec[-208]
DETECTOR(2.5, 0, 0) rec[-576] rec[-575] rec[-574] rec[-573] rec[-572] rec[-571] rec[-570] rec[-569] rec[-568] rec[-288] rec[-287] rec[-286] rec[-285] rec[-284] rec[-283] rec[-282] rec[-281] rec[-280]
DETECTOR(3.5, 0, 0) rec[-504] rec[-503] rec[-502] rec[-501] rec[-500] rec[-499] rec[-498] rec[-497] rec[-496] rec[-207] rec[-206] rec[-205] rec[-204] rec[-203] rec[-202] rec[-201] rec[-200] rec[-199]
DETECTOR(4.5, 0, 0) rec[-567] rec[-566] rec[-565] rec[-564] rec[-563] rec[-562] rec[-561] rec[-560] rec[-559] rec[-279] rec[-278] rec[-277] rec[-276] rec[-275] rec[-274] rec[-273] rec[-272] rec[-271]
DETECTOR(5.5, 0, 0) rec[-495] rec[-494] rec[-493] rec[-492] rec[-491] rec[-490] rec[-489] rec[-488] rec[-487] rec[-198] rec[-197] rec[-196] rec[-195] rec[-194] rec[-193] rec[-192] rec[-191] rec[-190]
DETECTOR(6.5, 0, 0) rec[-558] rec[-557] rec[-556] rec[-555] rec[-554] rec[-553] rec[-552] rec[-551] rec[-550] rec[-270] rec[-269] rec[-268] rec[-267] rec[-266] rec[-265] rec[-264] rec[-263] rec[-262]
DETECTOR(7.5, 0, 0) rec[-486] rec[-485] rec[-484] rec[-483] rec[-482] rec[-481] rec[-480] rec[-479] rec[-478] rec[-189] rec[-188] rec[-187] rec[-186] rec[-185] rec[-184] rec[-183] rec[-182] rec[-181]
DETECTOR(9.5, 0, 0) rec[-477] rec[-476] rec[-475] rec[-474] rec[-473] rec[-472] rec[-471] rec[-470] rec[-469] rec[-180] rec[-179] rec[-178] rec[-177] rec[-176] rec[-175] rec[-174] rec[-173] rec[-172]
DETECTOR(10.5, 0, 0) rec[-549] rec[-548] rec[-547] rec[-546] rec[-545] rec[-544] rec[-543] rec[-542] rec[-541] rec[-252] rec[-251] rec[-250] rec[-249] rec[-248] rec[-247] rec[-246] rec[-245] rec[-244]
DETECTOR(11.5, 0, 0) rec[-468] rec[-467] rec[-466] rec[-465] rec[-464] rec[-463] rec[-462] rec[-461] rec[-460] rec[-171] rec[-170] rec[-169] rec[-168] rec[-167] rec[-166] rec[-165] rec[-164] rec[-163]
DETECTOR(12.5, 0, 0) rec[-540] rec[-539] rec[-538] rec[-537] rec[-536] rec[-535] rec[-534] rec[-533] rec[-532] rec[-243] rec[-242] rec[-241] rec[-240] rec[-239] rec[-238] rec[-237] rec[-236] rec[-235]
DETECTOR(13.5, 0, 0) rec[-459] rec[-458] rec[-457] rec[-456] rec[-455] rec[-454] rec[-453] rec[-452] rec[-451] rec[-162] rec[-161] rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-155] rec[-154]
DETECTOR(14.5, 0, 0) rec[-531] rec[-530] rec[-529] rec[-528] rec[-527] rec[-526] rec[-525] rec[-524] rec[-523] rec[-234] rec[-233] rec[-232] rec[-231] rec[-230] rec[-229] rec[-228] rec[-227] rec[-226]
DETECTOR(15.5, 0, 0) rec[-450] rec[-449] rec[-448] rec[-447] rec[-446] rec[-445] rec[-444] rec[-443] rec[-442] rec[-153] rec[-152] rec[-151] rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-145]
DETECTOR(16.5, 0, 0) rec[-522] rec[-521] rec[-520] rec[-519] rec[-518] rec[-517] rec[-516] rec[-515] rec[-514] rec[-225] rec[-224] rec[-223] rec[-222] rec[-221] rec[-220] rec[-219] rec[-218] rec[-217]
DETECTOR(0, 0.5, 0) rec[-441] rec[-437] rec[-433] rec[-429] rec[-425] rec[-421] rec[-417] rec[-413] rec[-409] rec[-405] rec[-401] rec[-397] rec[-393] rec[-389] rec[-385] rec[-381] rec[-377] rec[-373] rec[-144] rec[-140] rec[-136] rec[-132] rec[-128] rec[-124] rec[-120] rec[-116] rec[-112] rec[-108] rec[-104] rec[-100] rec[-96] rec[-92] rec[-88] rec[-84] rec[-80] rec[-76]
DETECTOR(0, 1.5, 0) rec[-369] rec[-365] rec[-361] rec[-357] rec[-353] rec[-349] rec[-345] rec[-341] rec[-337] rec[-333] rec[-329] rec[-325] rec[-321] rec[-317] rec[-313] rec[-309] rec[-305] rec[-301] rec[-72] rec[-68] rec[-64] rec[-60] rec[-56] rec[-52] rec[-48] rec[-44] rec[-40] rec[-36] rec[-32] rec[-28] rec[-24] rec[-20] rec[-16] rec[-12] rec[-8] rec[-4]
DETECTOR(0, 2.5, 0) rec[-440] rec[-436] rec[-432] rec[-428] rec[-424] rec[-420] rec[-416] rec[-412] rec[-408] rec[-404] rec[-400] rec[-396] rec[-392] rec[-388] rec[-384] rec[-380] rec[-376] rec[-372] rec[-143] rec[-139] rec[-135] rec[-131] rec[-127] rec[-123] rec[-119] rec[-115] rec[-111] rec[-107] rec[-103] rec[-99] rec[-95] rec[-91] rec[-87] rec[-83] rec[-79] rec[-75]
DETECTOR(0, 3.5, 0) rec[-368] rec[-364] rec[-360] rec[-356] rec[-352] rec[-348] rec[-344] rec[-340] rec[-336] rec[-332] rec[-328] rec[-324] rec[-320] rec[-316] rec[-312] rec[-308] rec[-304] rec[-300] rec[-71] rec[-67] rec[-63] rec[-59] rec[-55] rec[-51] rec[-47] rec[-43] rec[-39] rec[-35] rec[-31] rec[-27] rec[-23] rec[-19] rec[-15] rec[-11] rec[-7] rec[-3]
DETECTOR(0, 4.5, 0) rec[-439] rec[-435] rec[-431] rec[-427] rec[-423] rec[-419] rec[-415] rec[-411] rec[-407] rec[-403] rec[-399] rec[-395] rec[-391] rec[-387] rec[-383] rec[-379] rec[-375] rec[-371] rec[-142] rec[-138] rec[-134] rec[-130] rec[-126] rec[-122] rec[-118] rec[-114] rec[-110] rec[-106] rec[-102] rec[-98] rec[-94] rec[-90] rec[-86] rec[-82] rec[-78] rec[-74]
DETECTOR(0, 5.5, 0) rec[-367] rec[-363] rec[-359] rec[-355] rec[-351] rec[-347] rec[-343] rec[-339] rec[-335] rec[-331] rec[-327] rec[-323] rec[-319] rec[-315] rec[-311] rec[-307] rec[-303] rec[-299] rec[-70] rec[-66] rec[-62] rec[-58] rec[-54] rec[-50] rec[-46] rec[-42] rec[-38] rec[-34] rec[-30] rec[-26] rec[-22] rec[-18] rec[-14] rec[-10] rec[-6] rec[-2]
DETECTOR(0, 6.5, 0) rec[-438] rec[-434] rec[-430] rec[-426] rec[-422] rec[-418] rec[-414] rec[-410] rec[-406] rec[-402] rec[-398] rec[-394] rec[-390] rec[-386] rec[-382] rec[-378] rec[-374] rec[-370] rec[-141] rec[-137] rec[-133] rec[-129] rec[-125] rec[-121] rec[-117] rec[-113] rec[-109] rec[-105] rec[-101] rec[-97] rec[-93] rec[-89] rec[-85] rec[-81] rec[-77] rec[-73]
DETECTOR(0, 7.5, 0) rec[-366] rec[-362] rec[-358] rec[-354] rec[-350] rec[-346] rec[-342] rec[-338] rec[-334] rec[-330] rec[-326] rec[-322] rec[-318] rec[-314] rec[-310] rec[-306] rec[-302] rec[-298] rec[-69] rec[-65] rec[-61] rec[-57] rec[-53] rec[-49] rec[-45] rec[-41] rec[-37] rec[-33] rec[-29] rec[-25] rec[-21] rec[-17] rec[-13] rec[-9] rec[-5] rec[-1]
OBSERVABLE_INCLUDE(2) rec[-261] rec[-260] rec[-259] rec[-258] rec[-257] rec[-256] rec[-255] rec[-254] rec[-253]
TICK
REPEAT 3 {
    MPP("""+str(p_f2)+""")  X0*X9 X1*X10 X2*X11 X3*X12 X4*X13 X5*X14 X6*X15 X7*X16 X8*X17 X18*X27 X19*X28 X20*X29 X21*X30 X22*X31 X23*X32 X24*X33 X25*X34 X26*X35 X36*X45 X37*X46 X38*X47 X39*X48 X40*X49 X41*X50 X42*X51 X43*X52 X44*X53 X54*X63 X55*X64 X56*X65 X57*X66 X58*X67 X59*X68 X60*X69 X61*X70 X62*X71 X72*X81 X73*X82 X74*X83 X75*X84 X76*X85 X77*X86 X78*X87 X79*X88 X80*X89 X90*X99 X91*X100 X92*X101 X93*X102 X94*X103 X95*X104 X96*X105 X97*X106 X98*X107 X108*X117 X109*X118 X110*X119 X111*X120 X112*X121 X113*X122 X114*X123 X115*X124 X116*X125 X126*X135 X127*X136 X128*X137 X129*X138 X130*X139 X131*X140 X132*X141 X133*X142 X134*X143 X144*X153 X145*X154 X146*X155 X147*X156 X148*X157 X149*X158 X150*X159 X151*X160 X152*X161
    DEPOLARIZE2("""+str(p_f)+""") 0 9 1 10 2 11 3 12 4 13 5 14 6 15 7 16 8 17 18 27 19 28 20 29 21 30 22 31 23 32 24 33 25 34 26 35 36 45 37 46 38 47 39 48 40 49 41 50 42 51 43 52 44 53 54 63 55 64 56 65 57 66 58 67 59 68 60 69 61 70 62 71 72 81 73 82 74 83 75 84 76 85 77 86 78 87 79 88 80 89 90 99 91 100 92 101 93 102 94 103 95 104 96 105 97 106 98 107 108 117 109 118 110 119 111 120 112 121 113 122 114 123 115 124 116 125 126 135 127 136 128 137 129 138 130 139 131 140 132 141 133 142 134 143 144 153 145 154 146 155 147 156 148 157 149 158 150 159 151 160 152 161
    DEPOLARIZE1("""+str(p_f)+""") 53 55 56 57 58 59 60 61 62 73 74 75 76 78 79 80 91 92 93 94 95 96 97 98 108 109 110 111 112 113 114 115 116
    TICK
    MPP("""+str(p_f2)+""")  X9*X18 X10*X19 X11*X20 X12*X21 X13*X22 X14*X23 X15*X24 X16*X25 X17*X26 X27*X36 X28*X37 X29*X38 X30*X39 X31*X40 X32*X41 X33*X42 X34*X43 X35*X44 X45*X54 X46*X55 X47*X56 X48*X57 X49*X58 X50*X59 X51*X60 X52*X61 X53*X62 X63*X72 X64*X73 X65*X74 X66*X75 X67*X76 X68*X77 X69*X78 X70*X79 X71*X80 X81*X90 X82*X91 X83*X92 X84*X93 X85*X94 X86*X95 X87*X96 X88*X97 X89*X98 X99*X108 X100*X109 X101*X110 X102*X111 X103*X112 X104*X113 X105*X114 X106*X115 X107*X116 X117*X126 X118*X127 X119*X128 X120*X129 X121*X130 X122*X131 X123*X132 X124*X133 X125*X134 X135*X144 X136*X145 X137*X146 X138*X147 X139*X148 X140*X149 X141*X150 X142*X151 X143*X152
    DEPOLARIZE2("""+str(p_f)+""") 9 18 10 19 11 20 12 21 13 22 14 23 15 24 16 25 17 26 27 36 28 37 29 38 30 39 31 40 32 41 33 42 34 43 35 44 45 54 46 55 47 56 48 57 49 58 50 59 51 60 52 61 53 62 63 72 64 73 65 74 66 75 67 76 68 77 69 78 70 79 71 80 81 90 82 91 83 92 84 93 85 94 86 95 87 96 88 97 89 98 99 108 100 109 101 110 102 111 103 112 104 113 105 114 106 115 107 116 117 126 118 127 119 128 120 129 121 130 122 131 123 132 124 133 125 134 135 144 136 145 137 146 138 147 139 148 140 149 141 150 142 151 143 152
    DEPOLARIZE1("""+str(p_f)+""") 0 2 3 4 5 6 7 8 29 30 32 34 41 43 50 52 55 60 62 64 66 68 70 72 73 75 77 79 82 84 86 88 90 92 94 96 98 100 102 104 106 108 110 112 114 116 118 120 122 124 126 128 130 132 134 136 138 140 142 144 146 148 150 152

    TICK
    MPP("""+str(p_f2)+""")  Z0*Z1 Z2*Z3 Z4*Z5 Z6*Z7 Z9*Z10 Z11*Z12 Z13*Z14 Z15*Z16 Z18*Z19 Z20*Z21 Z22*Z23 Z24*Z25 Z27*Z28 Z29*Z30 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48 Z49*Z50 Z51*Z52 Z54*Z55 Z56*Z57 Z58*Z59 Z60*Z61 Z63*Z64 Z65*Z66 Z67*Z68 Z69*Z70 Z72*Z73 Z74*Z75 Z76*Z77 Z78*Z79 Z81*Z82 Z83*Z84 Z85*Z86 Z87*Z88 Z90*Z91 Z92*Z93 Z94*Z95 Z96*Z97 Z99*Z100 Z101*Z102 Z103*Z104 Z105*Z106 Z108*Z109 Z110*Z111 Z112*Z113 Z114*Z115 Z117*Z118 Z119*Z120 Z121*Z122 Z123*Z124 Z126*Z127 Z128*Z129 Z130*Z131 Z132*Z133 Z135*Z136 Z137*Z138 Z139*Z140 Z141*Z142 Z144*Z145 Z146*Z147 Z148*Z149 Z150*Z151 Z153*Z154 Z155*Z156 Z157*Z158 Z159*Z160
    DEPOLARIZE2("""+str(p_f)+""") 0 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 25 27 28 29 30 31 32 33 34 36 37 38 39 40 41 42 43 45 46 47 48 49 50 51 52 54 55 56 57 58 59 60 61 63 64 65 66 67 68 69 70 72 73 74 75 76 77 78 79 81 82 83 84 85 86 87 88 90 91 92 93 94 95 96 97 99 100 101 102 103 104 105 106 108 109 110 111 112 113 114 115 117 118 119 120 121 122 123 124 126 127 128 129 130 131 132 133 135 136 137 138 139 140 141 142 144 145 146 147 148 149 150 151 153 154 155 156 157 158 159 160
    DEPOLARIZE1("""+str(p_f)+""") 8 17 26 35 44 53 62 71 80 89 98 107 116 125 134 143 152

    TICK
    MPP("""+str(p_f2)+""")  Z1*Z2 Z3*Z4 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z14*Z15 Z16*Z17 Z19*Z20 Z21*Z22 Z23*Z24 Z25*Z26 Z28*Z29 Z30*Z31 Z32*Z33 Z34*Z35 Z37*Z38 Z39*Z40 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49 Z50*Z51 Z52*Z53 Z55*Z56 Z57*Z58 Z59*Z60 Z61*Z62 Z64*Z65 Z66*Z67 Z68*Z69 Z70*Z71 Z73*Z74 Z75*Z76 Z77*Z78 Z79*Z80 Z82*Z83 Z84*Z85 Z86*Z87 Z88*Z89 Z91*Z92 Z93*Z94 Z95*Z96 Z97*Z98 Z100*Z101 Z102*Z103 Z104*Z105 Z106*Z107 Z109*Z110 Z111*Z112 Z113*Z114 Z115*Z116 Z118*Z119 Z120*Z121 Z122*Z123 Z124*Z125 Z127*Z128 Z129*Z130 Z131*Z132 Z133*Z134 Z136*Z137 Z138*Z139 Z140*Z141 Z142*Z143 Z145*Z146 Z147*Z148 Z149*Z150 Z151*Z152 Z154*Z155 Z156*Z157 Z158*Z159 Z160*Z161
    DEPOLARIZE2("""+str(p_f)+""") 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 35 37 38 39 40 41 42 43 44 46 47 48 49 50 51 52 53 55 56 57 58 59 60 61 62 64 65 66 67 68 69 70 71 73 74 75 76 77 78 79 80 82 83 84 85 86 87 88 89 91 92 93 94 95 96 97 98 100 101 102 103 104 105 106 107 109 110 111 112 113 114 115 116 118 119 120 121 122 123 124 125 127 128 129 130 131 132 133 134 136 137 138 139 140 141 142 143 145 146 147 148 149 150 151 152 154 155 156 157 158 159 160 161
    DEPOLARIZE1("""+str(p_f)+""") 0 9 18 27 36 45 54 63 72 81 90 99 108 117 126

    SHIFT_COORDS(0, 0, 1)
    
    DETECTOR(0.5, 0, 0) rec[-594] rec[-593] rec[-592] rec[-591] rec[-590] rec[-589] rec[-588] rec[-587] rec[-586] rec[-297] rec[-296] rec[-295] rec[-294] rec[-293] rec[-292] rec[-291] rec[-290] rec[-289]
    DETECTOR(1.5, 0, 0) rec[-513] rec[-512] rec[-511] rec[-510] rec[-509] rec[-508] rec[-507] rec[-506] rec[-505] rec[-216] rec[-215] rec[-214] rec[-213] rec[-212] rec[-211] rec[-210] rec[-209] rec[-208]
    DETECTOR(2.5, 0, 0) rec[-585] rec[-584] rec[-583] rec[-582] rec[-581] rec[-580] rec[-579] rec[-578] rec[-577] rec[-288] rec[-287] rec[-286] rec[-285] rec[-284] rec[-283] rec[-282] rec[-281] rec[-280]
    DETECTOR(3.5, 0, 0) rec[-504] rec[-503] rec[-502] rec[-501] rec[-500] rec[-499] rec[-498] rec[-497] rec[-496] rec[-207] rec[-206] rec[-205] rec[-204] rec[-203] rec[-202] rec[-201] rec[-200] rec[-199]
    DETECTOR(4.5, 0, 0) rec[-576] rec[-575] rec[-574] rec[-573] rec[-572] rec[-571] rec[-570] rec[-569] rec[-568] rec[-279] rec[-278] rec[-277] rec[-276] rec[-275] rec[-274] rec[-273] rec[-272] rec[-271]
    DETECTOR(5.5, 0, 0) rec[-495] rec[-494] rec[-493] rec[-492] rec[-491] rec[-490] rec[-489] rec[-488] rec[-487] rec[-198] rec[-197] rec[-196] rec[-195] rec[-194] rec[-193] rec[-192] rec[-191] rec[-190]
    DETECTOR(6.5, 0, 0) rec[-567] rec[-566] rec[-565] rec[-564] rec[-563] rec[-562] rec[-561] rec[-560] rec[-559] rec[-270] rec[-269] rec[-268] rec[-267] rec[-266] rec[-265] rec[-264] rec[-263] rec[-262]
    DETECTOR(7.5, 0, 0) rec[-486] rec[-485] rec[-484] rec[-483] rec[-482] rec[-481] rec[-480] rec[-479] rec[-478] rec[-189] rec[-188] rec[-187] rec[-186] rec[-185] rec[-184] rec[-183] rec[-182] rec[-181]
    DETECTOR(8.5, 0, 0) rec[-558] rec[-557] rec[-556] rec[-555] rec[-554] rec[-553] rec[-552] rec[-551] rec[-550] rec[-261] rec[-260] rec[-259] rec[-258] rec[-257] rec[-256] rec[-255] rec[-254] rec[-253]
    DETECTOR(9.5, 0, 0) rec[-477] rec[-476] rec[-475] rec[-474] rec[-473] rec[-472] rec[-471] rec[-470] rec[-469] rec[-180] rec[-179] rec[-178] rec[-177] rec[-176] rec[-175] rec[-174] rec[-173] rec[-172]
    DETECTOR(10.5, 0, 0) rec[-549] rec[-548] rec[-547] rec[-546] rec[-545] rec[-544] rec[-543] rec[-542] rec[-541] rec[-252] rec[-251] rec[-250] rec[-249] rec[-248] rec[-247] rec[-246] rec[-245] rec[-244]
    DETECTOR(11.5, 0, 0) rec[-468] rec[-467] rec[-466] rec[-465] rec[-464] rec[-463] rec[-462] rec[-461] rec[-460] rec[-171] rec[-170] rec[-169] rec[-168] rec[-167] rec[-166] rec[-165] rec[-164] rec[-163]
    DETECTOR(12.5, 0, 0) rec[-540] rec[-539] rec[-538] rec[-537] rec[-536] rec[-535] rec[-534] rec[-533] rec[-532] rec[-243] rec[-242] rec[-241] rec[-240] rec[-239] rec[-238] rec[-237] rec[-236] rec[-235]
    DETECTOR(13.5, 0, 0) rec[-459] rec[-458] rec[-457] rec[-456] rec[-455] rec[-454] rec[-453] rec[-452] rec[-451] rec[-162] rec[-161] rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-155] rec[-154]
    DETECTOR(14.5, 0, 0) rec[-531] rec[-530] rec[-529] rec[-528] rec[-527] rec[-526] rec[-525] rec[-524] rec[-523] rec[-234] rec[-233] rec[-232] rec[-231] rec[-230] rec[-229] rec[-228] rec[-227] rec[-226]
    DETECTOR(15.5, 0, 0) rec[-450] rec[-449] rec[-448] rec[-447] rec[-446] rec[-445] rec[-444] rec[-443] rec[-442] rec[-153] rec[-152] rec[-151] rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-145]
    DETECTOR(16.5, 0, 0) rec[-522] rec[-521] rec[-520] rec[-519] rec[-518] rec[-517] rec[-516] rec[-515] rec[-514] rec[-225] rec[-224] rec[-223] rec[-222] rec[-221] rec[-220] rec[-219] rec[-218] rec[-217]
    DETECTOR(0, 0.5, 0) rec[-441] rec[-437] rec[-433] rec[-429] rec[-425] rec[-421] rec[-417] rec[-413] rec[-409] rec[-405] rec[-401] rec[-397] rec[-393] rec[-389] rec[-385] rec[-381] rec[-377] rec[-373] rec[-144] rec[-140] rec[-136] rec[-132] rec[-128] rec[-124] rec[-120] rec[-116] rec[-112] rec[-108] rec[-104] rec[-100] rec[-96] rec[-92] rec[-88] rec[-84] rec[-80] rec[-76]
    DETECTOR(0, 1.5, 0) rec[-369] rec[-365] rec[-361] rec[-357] rec[-353] rec[-349] rec[-345] rec[-341] rec[-337] rec[-333] rec[-329] rec[-325] rec[-321] rec[-317] rec[-313] rec[-309] rec[-305] rec[-301] rec[-72] rec[-68] rec[-64] rec[-60] rec[-56] rec[-52] rec[-48] rec[-44] rec[-40] rec[-36] rec[-32] rec[-28] rec[-24] rec[-20] rec[-16] rec[-12] rec[-8] rec[-4]
    DETECTOR(0, 2.5, 0) rec[-440] rec[-436] rec[-432] rec[-428] rec[-424] rec[-420] rec[-416] rec[-412] rec[-408] rec[-404] rec[-400] rec[-396] rec[-392] rec[-388] rec[-384] rec[-380] rec[-376] rec[-372] rec[-143] rec[-139] rec[-135] rec[-131] rec[-127] rec[-123] rec[-119] rec[-115] rec[-111] rec[-107] rec[-103] rec[-99] rec[-95] rec[-91] rec[-87] rec[-83] rec[-79] rec[-75]
    DETECTOR(0, 3.5, 0) rec[-368] rec[-364] rec[-360] rec[-356] rec[-352] rec[-348] rec[-344] rec[-340] rec[-336] rec[-332] rec[-328] rec[-324] rec[-320] rec[-316] rec[-312] rec[-308] rec[-304] rec[-300] rec[-71] rec[-67] rec[-63] rec[-59] rec[-55] rec[-51] rec[-47] rec[-43] rec[-39] rec[-35] rec[-31] rec[-27] rec[-23] rec[-19] rec[-15] rec[-11] rec[-7] rec[-3]
    DETECTOR(0, 4.5, 0) rec[-439] rec[-435] rec[-431] rec[-427] rec[-423] rec[-419] rec[-415] rec[-411] rec[-407] rec[-403] rec[-399] rec[-395] rec[-391] rec[-387] rec[-383] rec[-379] rec[-375] rec[-371] rec[-142] rec[-138] rec[-134] rec[-130] rec[-126] rec[-122] rec[-118] rec[-114] rec[-110] rec[-106] rec[-102] rec[-98] rec[-94] rec[-90] rec[-86] rec[-82] rec[-78] rec[-74]
    DETECTOR(0, 5.5, 0) rec[-367] rec[-363] rec[-359] rec[-355] rec[-351] rec[-347] rec[-343] rec[-339] rec[-335] rec[-331] rec[-327] rec[-323] rec[-319] rec[-315] rec[-311] rec[-307] rec[-303] rec[-299] rec[-70] rec[-66] rec[-62] rec[-58] rec[-54] rec[-50] rec[-46] rec[-42] rec[-38] rec[-34] rec[-30] rec[-26] rec[-22] rec[-18] rec[-14] rec[-10] rec[-6] rec[-2]
    DETECTOR(0, 6.5, 0) rec[-438] rec[-434] rec[-430] rec[-426] rec[-422] rec[-418] rec[-414] rec[-410] rec[-406] rec[-402] rec[-398] rec[-394] rec[-390] rec[-386] rec[-382] rec[-378] rec[-374] rec[-370] rec[-141] rec[-137] rec[-133] rec[-129] rec[-125] rec[-121] rec[-117] rec[-113] rec[-109] rec[-105] rec[-101] rec[-97] rec[-93] rec[-89] rec[-85] rec[-81] rec[-77] rec[-73]
    DETECTOR(0, 7.5, 0) rec[-366] rec[-362] rec[-358] rec[-354] rec[-350] rec[-346] rec[-342] rec[-338] rec[-334] rec[-330] rec[-326] rec[-322] rec[-318] rec[-314] rec[-310] rec[-306] rec[-302] rec[-298] rec[-69] rec[-65] rec[-61] rec[-57] rec[-53] rec[-49] rec[-45] rec[-41] rec[-37] rec[-33] rec[-29] rec[-25] rec[-21] rec[-17] rec[-13] rec[-9] rec[-5] rec[-1]
}
#   288     287   286   285    284      283  282     281    280     279    278     277     276      275      274   273      272     271   270      269     268     267     266    265      264     263     262     261     260      259     258    257     256     255     254     253     252     251     250      249       248     247       246      245      244      243      242        241       240        239      238      237        236      235        234       233      232       231        230      229      228        227       226       225       224       223        222       221       220       219      218        217      
TICK
MPP("""+str(p_f2)+""")  X0*X9 X1*X10 X2*X11 X3*X12 X4*X13 X5*X14 X6*X15 X7*X16 X8*X17 X18*X27 X19*X28 X20*X29 X21*X30 X22*X31 X23*X32 X24*X33 X25*X34 X26*X35 X36*X45 X37*X46 X38*X47 X39*X48 X40*X49 X41*X50 X42*X51 X43*X52 X44*X53 X54*X63 X55*X64 X56*X65 X57*X66 X58*X67 X59*X68 X60*X69 X61*X70 X62*X71 X90*X99 X91*X100 X92*X101 X93*X102 X94*X103 X95*X104 X96*X105 X97*X106 X98*X107 X108*X117 X109*X118 X110*X119 X111*X120 X112*X121 X113*X122 X114*X123 X115*X124 X116*X125 X126*X135 X127*X136 X128*X137 X129*X138 X130*X139 X131*X140 X132*X141 X133*X142 X134*X143 X144*X153 X145*X154 X146*X155 X147*X156 X148*X157 X149*X158 X150*X159 X151*X160 X152*X161
DEPOLARIZE2("""+str(p_f)+""") 0 9 1 10 2 11 3 12 4 13 5 14 6 15 7 16 8 17 18 27 19 28 20 29 21 30 22 31 23 32 24 33 25 34 26 35 36 45 37 46 38 47 39 48 40 49 41 50 42 51 43 52 44 53 54 63 55 64 56 65 57 66 58 67 59 68 60 69 61 70 62 71 90 99 91 100 92 101 93 102 94 103 95 104 96 105 97 106 98 107 108 117 109 118 110 119 111 120 112 121 113 122 114 123 115 124 116 125 126 135 127 136 128 137 129 138 130 139 131 140 132 141 133 142 134 143 144 153 145 154 146 155 147 156 148 157 149 158 150 159 151 160 152 161
DEPOLARIZE1("""+str(p_f)+""") 5 6 7 8 9 10 18 27 29 32 33 35 36 40
TICK
MPP("""+str(p_f2)+""")  X9*X18 X10*X19 X11*X20 X12*X21 X13*X22 X14*X23 X15*X24 X16*X25 X17*X26 X27*X36 X28*X37 X29*X38 X30*X39 X31*X40 X32*X41 X33*X42 X34*X43 X35*X44 X45*X54 X46*X55 X47*X56 X48*X57 X49*X58 X50*X59 X51*X60 X52*X61 X53*X62 X63*X72 X64*X73 X65*X74 X66*X75 X67*X76 X68*X77 X69*X78 X70*X79 X71*X80 X81*X90 X82*X91 X83*X92 X84*X93 X85*X94 X86*X95 X87*X96 X88*X97 X89*X98 X99*X108 X100*X109 X101*X110 X102*X111 X103*X112 X104*X113 X105*X114 X106*X115 X107*X116 X117*X126 X118*X127 X119*X128 X120*X129 X121*X130 X122*X131 X123*X132 X124*X133 X125*X134 X135*X144 X136*X145 X137*X146 X138*X147 X139*X148 X140*X149 X141*X150 X142*X151 X143*X152
DEPOLARIZE2("""+str(p_f)+""") 9 18 10 19 11 20 12 21 13 22 14 23 15 24 16 25 17 26 27 36 28 37 29 38 30 39 31 40 32 41 33 42 34 43 35 44 45 54 46 55 47 56 48 57 49 58 50 59 51 60 52 61 53 62 63 72 64 73 65 74 66 75 67 76 68 77 69 78 70 79 71 80 81 90 82 91 83 92 84 93 85 94 86 95 87 96 88 97 89 98 99 108 100 109 101 110 102 111 103 112 104 113 105 114 106 115 107 116 117 126 118 127 119 128 120 129 121 130 122 131 123 132 124 133 125 134 135 144 136 145 137 146 138 147 139 148 140 149 141 150 142 151 143 152
DEPOLARIZE1("""+str(p_f)+""") 0 2 3 4 5 6 7 8 29 30 32 34 41 43 50 52 55 60 62 64 66 68 70 72 73 75 77 79 82 84 86 88 90 92 94 96 98 100 102 104 106 108 110 112 114 116 118 120 122 124 126 128 130 132 134 136 138 140 142 144 146 148 150 152


TICK
MPP("""+str(p_f2)+""")  Z0*Z1 Z2*Z3 Z4*Z5 Z6*Z7 Z9*Z10 Z11*Z12 Z13*Z14 Z15*Z16 Z18*Z19 Z20*Z21 Z22*Z23 Z24*Z25 Z27*Z28 Z29*Z30 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48 Z49*Z50 Z51*Z52 Z54*Z55 Z56*Z57 Z58*Z59 Z60*Z61 Z63*Z64 Z65*Z66 Z67*Z68 Z69*Z70 Z72*Z73 Z74*Z75 Z76*Z77 Z78*Z79 Z81*Z82 Z83*Z84 Z85*Z86 Z87*Z88 Z90*Z91 Z92*Z93 Z94*Z95 Z96*Z97 Z99*Z100 Z101*Z102 Z103*Z104 Z105*Z106 Z108*Z109 Z110*Z111 Z112*Z113 Z114*Z115 Z117*Z118 Z119*Z120 Z121*Z122 Z123*Z124 Z126*Z127 Z128*Z129 Z130*Z131 Z132*Z133 Z135*Z136 Z137*Z138 Z139*Z140 Z141*Z142 Z144*Z145 Z146*Z147 Z148*Z149 Z150*Z151 Z153*Z154 Z155*Z156 Z157*Z158 Z159*Z160
DEPOLARIZE2("""+str(p_f)+""") 0 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 25 27 28 29 30 31 32 33 34 36 37 38 39 40 41 42 43 45 46 47 48 49 50 51 52 54 55 56 57 58 59 60 61 63 64 65 66 67 68 69 70 72 73 74 75 76 77 78 79 81 82 83 84 85 86 87 88 90 91 92 93 94 95 96 97 99 100 101 102 103 104 105 106 108 109 110 111 112 113 114 115 117 118 119 120 121 122 123 124 126 127 128 129 130 131 132 133 135 136 137 138 139 140 141 142 144 145 146 147 148 149 150 151 153 154 155 156 157 158 159 160
DEPOLARIZE1("""+str(p_f)+""") 8 17 26 35 44 53 62 71 80 89 98 107 116 125 134 143 152




TICK
MPP("""+str(p_f2)+""")  Z1*Z2 Z3*Z4 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z14*Z15 Z16*Z17 Z19*Z20 Z21*Z22 Z23*Z24 Z25*Z26 Z28*Z29 Z30*Z31 Z32*Z33 Z34*Z35 Z37*Z38 Z39*Z40 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49 Z50*Z51 Z52*Z53 Z55*Z56 Z57*Z58 Z59*Z60 Z61*Z62 Z64*Z65 Z66*Z67 Z68*Z69 Z70*Z71 Z73*Z74 Z75*Z76 Z77*Z78 Z79*Z80 Z82*Z83 Z84*Z85 Z86*Z87 Z88*Z89 Z91*Z92 Z93*Z94 Z95*Z96 Z97*Z98 Z100*Z101 Z102*Z103 Z104*Z105 Z106*Z107 Z109*Z110 Z111*Z112 Z113*Z114 Z115*Z116 Z118*Z119 Z120*Z121 Z122*Z123 Z124*Z125 Z127*Z128 Z129*Z130 Z131*Z132 Z133*Z134 Z136*Z137 Z138*Z139 Z140*Z141 Z142*Z143 Z145*Z146 Z147*Z148 Z149*Z150 Z151*Z152 Z154*Z155 Z156*Z157 Z158*Z159 Z160*Z161
DEPOLARIZE2("""+str(p_f)+""") 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 35 37 38 39 40 41 42 43 44 46 47 48 49 50 51 52 53 55 56 57 58 59 60 61 62 64 65 66 67 68 69 70 71 73 74 75 76 77 78 79 80 82 83 84 85 86 87 88 89 91 92 93 94 95 96 97 98 100 101 102 103 104 105 106 107 109 110 111 112 113 114 115 116 118 119 120 121 122 123 124 125 127 128 129 130 131 132 133 134 136 137 138 139 140 141 142 143 145 146 147 148 149 150 151 152 154 155 156 157 158 159 160 161
DEPOLARIZE1("""+str(p_f)+""") 0 9 18 27 36 45 54 63 72 81 90 99 108 117 126


SHIFT_COORDS(0, 0, 1)

DETECTOR(0.5, 0, 0) rec[-585] rec[-584] rec[-583] rec[-582] rec[-581] rec[-580] rec[-579] rec[-578] rec[-577] rec[-288] rec[-287] rec[-286] rec[-285] rec[-284] rec[-283] rec[-282] rec[-281] rec[-280]
DETECTOR(1.5, 0, 0) rec[-504] rec[-503] rec[-502] rec[-501] rec[-500] rec[-499] rec[-498] rec[-497] rec[-496] rec[-216] rec[-215] rec[-214] rec[-213] rec[-212] rec[-211] rec[-210] rec[-209] rec[-208]
DETECTOR(2.5, 0, 0) rec[-576] rec[-575] rec[-574] rec[-573] rec[-572] rec[-571] rec[-570] rec[-569] rec[-568] rec[-279] rec[-278] rec[-277] rec[-276] rec[-275] rec[-274] rec[-273] rec[-272] rec[-271]
DETECTOR(3.5, 0, 0) rec[-495] rec[-494] rec[-493] rec[-492] rec[-491] rec[-490] rec[-489] rec[-488] rec[-487] rec[-207] rec[-206] rec[-205] rec[-204] rec[-203] rec[-202] rec[-201] rec[-200] rec[-199]
DETECTOR(4.5, 0, 0) rec[-567] rec[-566] rec[-565] rec[-564] rec[-563] rec[-562] rec[-561] rec[-560] rec[-559] rec[-270] rec[-269] rec[-268] rec[-267] rec[-266] rec[-265] rec[-264] rec[-263] rec[-262]
DETECTOR(5.5, 0, 0) rec[-486] rec[-485] rec[-484] rec[-483] rec[-482] rec[-481] rec[-480] rec[-479] rec[-478] rec[-198] rec[-197] rec[-196] rec[-195] rec[-194] rec[-193] rec[-192] rec[-191] rec[-190]
DETECTOR(6.5, 0, 0) rec[-558] rec[-557] rec[-556] rec[-555] rec[-554] rec[-553] rec[-552] rec[-551] rec[-550] rec[-261] rec[-260] rec[-259] rec[-258] rec[-257] rec[-256] rec[-255] rec[-254] rec[-253]
DETECTOR(7.5, 0, 0) rec[-477] rec[-476] rec[-475] rec[-474] rec[-473] rec[-472] rec[-471] rec[-470] rec[-469] rec[-189] rec[-188] rec[-187] rec[-186] rec[-185] rec[-184] rec[-183] rec[-182] rec[-181]
DETECTOR(9.5, 0, 0) rec[-468] rec[-467] rec[-466] rec[-465] rec[-464] rec[-463] rec[-462] rec[-461] rec[-460] rec[-180] rec[-179] rec[-178] rec[-177] rec[-176] rec[-175] rec[-174] rec[-173] rec[-172]
DETECTOR(10.5, 0, 0) rec[-540] rec[-539] rec[-538] rec[-537] rec[-536] rec[-535] rec[-534] rec[-533] rec[-532] rec[-252] rec[-251] rec[-250] rec[-249] rec[-248] rec[-247] rec[-246] rec[-245] rec[-244]
DETECTOR(11.5, 0, 0) rec[-459] rec[-458] rec[-457] rec[-456] rec[-455] rec[-454] rec[-453] rec[-452] rec[-451] rec[-171] rec[-170] rec[-169] rec[-168] rec[-167] rec[-166] rec[-165] rec[-164] rec[-163]
DETECTOR(12.5, 0, 0) rec[-531] rec[-530] rec[-529] rec[-528] rec[-527] rec[-526] rec[-525] rec[-524] rec[-523] rec[-243] rec[-242] rec[-241] rec[-240] rec[-239] rec[-238] rec[-237] rec[-236] rec[-235]
DETECTOR(13.5, 0, 0) rec[-450] rec[-449] rec[-448] rec[-447] rec[-446] rec[-445] rec[-444] rec[-443] rec[-442] rec[-162] rec[-161] rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-155] rec[-154]
DETECTOR(14.5, 0, 0) rec[-522] rec[-521] rec[-520] rec[-519] rec[-518] rec[-517] rec[-516] rec[-515] rec[-514] rec[-234] rec[-233] rec[-232] rec[-231] rec[-230] rec[-229] rec[-228] rec[-227] rec[-226]
DETECTOR(15.5, 0, 0) rec[-441] rec[-440] rec[-439] rec[-438] rec[-437] rec[-436] rec[-435] rec[-434] rec[-433] rec[-153] rec[-152] rec[-151] rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-145]
DETECTOR(16.5, 0, 0) rec[-513] rec[-512] rec[-511] rec[-510] rec[-509] rec[-508] rec[-507] rec[-506] rec[-505] rec[-225] rec[-224] rec[-223] rec[-222] rec[-221] rec[-220] rec[-219] rec[-218] rec[-217]
DETECTOR(0, 0.5, 0) rec[-432] rec[-428] rec[-424] rec[-420] rec[-416] rec[-412] rec[-408] rec[-404] rec[-400] rec[-144] rec[-140] rec[-136] rec[-132] rec[-128] rec[-124] rec[-120] rec[-116] rec[-112]
DETECTOR(0, 1.5, 0) rec[-360] rec[-356] rec[-352] rec[-348] rec[-344] rec[-340] rec[-336] rec[-332] rec[-328] rec[-72] rec[-68] rec[-64] rec[-60] rec[-56] rec[-52] rec[-48] rec[-44] rec[-40]
DETECTOR(0, 2.5, 0) rec[-431] rec[-427] rec[-423] rec[-419] rec[-415] rec[-411] rec[-407] rec[-403] rec[-399] rec[-143] rec[-139] rec[-135] rec[-131] rec[-127] rec[-123] rec[-119] rec[-115] rec[-111]
DETECTOR(0, 3.5, 0) rec[-359] rec[-355] rec[-351] rec[-347] rec[-343] rec[-339] rec[-335] rec[-331] rec[-327] rec[-71] rec[-67] rec[-63] rec[-59] rec[-55] rec[-51] rec[-47] rec[-43] rec[-39]
DETECTOR(0, 4.5, 0) rec[-430] rec[-426] rec[-422] rec[-418] rec[-414] rec[-410] rec[-406] rec[-402] rec[-398] rec[-142] rec[-138] rec[-134] rec[-130] rec[-126] rec[-122] rec[-118] rec[-114] rec[-110]
DETECTOR(0, 5.5, 0) rec[-358] rec[-354] rec[-350] rec[-346] rec[-342] rec[-338] rec[-334] rec[-330] rec[-326] rec[-70] rec[-66] rec[-62] rec[-58] rec[-54] rec[-50] rec[-46] rec[-42] rec[-38]
DETECTOR(0, 6.5, 0) rec[-429] rec[-425] rec[-421] rec[-417] rec[-413] rec[-409] rec[-405] rec[-401] rec[-397] rec[-141] rec[-137] rec[-133] rec[-129] rec[-125] rec[-121] rec[-117] rec[-113] rec[-109]
DETECTOR(0, 7.5, 0) rec[-357] rec[-353] rec[-349] rec[-345] rec[-341] rec[-337] rec[-333] rec[-329] rec[-325] rec[-69] rec[-65] rec[-61] rec[-57] rec[-53] rec[-49] rec[-45] rec[-41] rec[-37]
DETECTOR(9, 0.5, 0) rec[-396] rec[-392] rec[-388] rec[-384] rec[-380] rec[-376] rec[-372] rec[-368] rec[-364] rec[-108] rec[-104] rec[-100] rec[-96] rec[-92] rec[-88] rec[-84] rec[-80] rec[-76]
DETECTOR(9, 1.5, 0) rec[-324] rec[-320] rec[-316] rec[-312] rec[-308] rec[-304] rec[-300] rec[-296] rec[-292] rec[-36] rec[-32] rec[-28] rec[-24] rec[-20] rec[-16] rec[-12] rec[-8] rec[-4]
DETECTOR(9, 2.5, 0) rec[-395] rec[-391] rec[-387] rec[-383] rec[-379] rec[-375] rec[-371] rec[-367] rec[-363] rec[-107] rec[-103] rec[-99] rec[-95] rec[-91] rec[-87] rec[-83] rec[-79] rec[-75]
DETECTOR(9, 3.5, 0) rec[-323] rec[-319] rec[-315] rec[-311] rec[-307] rec[-303] rec[-299] rec[-295] rec[-291] rec[-35] rec[-31] rec[-27] rec[-23] rec[-19] rec[-15] rec[-11] rec[-7] rec[-3]
DETECTOR(9, 4.5, 0) rec[-394] rec[-390] rec[-386] rec[-382] rec[-378] rec[-374] rec[-370] rec[-366] rec[-362] rec[-106] rec[-102] rec[-98] rec[-94] rec[-90] rec[-86] rec[-82] rec[-78] rec[-74]
DETECTOR(9, 5.5, 0) rec[-322] rec[-318] rec[-314] rec[-310] rec[-306] rec[-302] rec[-298] rec[-294] rec[-290] rec[-34] rec[-30] rec[-26] rec[-22] rec[-18] rec[-14] rec[-10] rec[-6] rec[-2]
DETECTOR(9, 6.5, 0) rec[-393] rec[-389] rec[-385] rec[-381] rec[-377] rec[-373] rec[-369] rec[-365] rec[-361] rec[-105] rec[-101] rec[-97] rec[-93] rec[-89] rec[-85] rec[-81] rec[-77] rec[-73]
DETECTOR(9, 7.5, 0) rec[-321] rec[-317] rec[-313] rec[-309] rec[-305] rec[-301] rec[-297] rec[-293] rec[-289] rec[-33] rec[-29] rec[-25] rec[-21] rec[-17] rec[-13] rec[-9] rec[-5] rec[-1]

REPEAT """+str(j)+"""{
   
    TICK
    MPP("""+str(p_f2)+""")  X0*X9 X1*X10 X2*X11 X3*X12 X4*X13 X5*X14 X6*X15 X7*X16 X8*X17 X18*X27 X19*X28 X20*X29 X21*X30 X22*X31 X23*X32 X24*X33 X25*X34 X26*X35 X36*X45 X37*X46 X38*X47 X39*X48 X40*X49 X41*X50 X42*X51 X43*X52 X44*X53 X54*X63 X55*X64 X56*X65 X57*X66 X58*X67 X59*X68 X60*X69 X61*X70 X62*X71 X90*X99 X91*X100 X92*X101 X93*X102 X94*X103 X95*X104 X96*X105 X97*X106 X98*X107 X108*X117 X109*X118 X110*X119 X111*X120 X112*X121 X113*X122 X114*X123 X115*X124 X116*X125 X126*X135 X127*X136 X128*X137 X129*X138 X130*X139 X131*X140 X132*X141 X133*X142 X134*X143 X144*X153 X145*X154 X146*X155 X147*X156 X148*X157 X149*X158 X150*X159 X151*X160 X152*X161
    DEPOLARIZE2("""+str(p_f)+""") 0 9 1 10 2 11 3 12 4 13 5 14 6 15 7 16 8 17 18 27 19 28 20 29 21 30 22 31 23 32 24 33 25 34 26 35 36 45 37 46 38 47 39 48 40 49 41 50 42 51 43 52 44 53 54 63 55 64 56 65 57 66 58 67 59 68 60 69 61 70 62 71 90 99 91 100 92 101 93 102 94 103 95 104 96 105 97 106 98 107 108 117 109 118 110 119 111 120 112 121 113 122 114 123 115 124 116 125 126 135 127 136 128 137 129 138 130 139 131 140 132 141 133 142 134 143 144 153 145 154 146 155 147 156 148 157 149 158 150 159 151 160 152 161
    DEPOLARIZE1("""+str(p_f)+""") 55 56 57 58 59 60 61 62 63 72 74 75 76 77 78 79 80 82 83 84 85 86 87 88 89 93 94 95 96 97 98 109 110 111 112 113 114 115 116 127 128 129 130 131 132 133 134 145 146 147 148 149 150 151 152

    TICK
    MPP("""+str(p_f2)+""")  X9*X18 X10*X19 X11*X20 X12*X21 X13*X22 X14*X23 X15*X24 X16*X25 X17*X26 X27*X36 X28*X37 X29*X38 X30*X39 X31*X40 X32*X41 X33*X42 X34*X43 X35*X44 X45*X54 X46*X55 X47*X56 X48*X57 X49*X58 X50*X59 X51*X60 X52*X61 X53*X62 X63*X72 X64*X73 X65*X74 X66*X75 X67*X76 X68*X77 X69*X78 X70*X79 X71*X80 X81*X90 X82*X91 X83*X92 X84*X93 X85*X94 X86*X95 X87*X96 X88*X97 X89*X98 X99*X108 X100*X109 X101*X110 X102*X111 X103*X112 X104*X113 X105*X114 X106*X115 X107*X116 X117*X126 X118*X127 X119*X128 X120*X129 X121*X130 X122*X131 X123*X132 X124*X133 X125*X134 X135*X144 X136*X145 X137*X146 X138*X147 X139*X148 X140*X149 X141*X150 X142*X151 X143*X152
    DEPOLARIZE2("""+str(p_f)+""") 9 18 10 19 11 20 12 21 13 22 14 23 15 24 16 25 17 26 27 36 28 37 29 38 30 39 31 40 32 41 33 42 34 43 35 44 45 54 46 55 47 56 48 57 49 58 50 59 51 60 52 61 53 62 63 72 64 73 65 74 66 75 67 76 68 77 69 78 70 79 71 80 81 90 82 91 83 92 84 93 85 94 86 95 87 96 88 97 89 98 99 108 100 109 101 110 102 111 103 112 104 113 105 114 106 115 107 116 117 126 118 127 119 128 120 129 121 130 122 131 123 132 124 133 125 134 135 144 136 145 137 146 138 147 139 148 140 149 141 150 142 151 143 152
    DEPOLARIZE1("""+str(p_f)+""") 0 2 3 4 5 6 7 8 29 30 32 34 41 43 50 52 55 60 62 64 66 68 70 72 73 75 77 79 82 84 86 88 90 92 94 96 98 100 102 104 106 108 110 112 114 116 118 120 122 124 126 128 130 132 134 136 138 140 142 144 146 148 150 152

    TICK
    MPP("""+str(p_f2)+""")  Z0*Z1 Z2*Z3 Z4*Z5 Z6*Z7 Z9*Z10 Z11*Z12 Z13*Z14 Z15*Z16 Z18*Z19 Z20*Z21 Z22*Z23 Z24*Z25 Z27*Z28 Z29*Z30 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48 Z49*Z50 Z51*Z52 Z54*Z55 Z56*Z57 Z58*Z59 Z60*Z61 Z63*Z64 Z65*Z66 Z67*Z68 Z69*Z70 Z72*Z73 Z74*Z75 Z76*Z77 Z78*Z79 Z81*Z82 Z83*Z84 Z85*Z86 Z87*Z88 Z90*Z91 Z92*Z93 Z94*Z95 Z96*Z97 Z99*Z100 Z101*Z102 Z103*Z104 Z105*Z106 Z108*Z109 Z110*Z111 Z112*Z113 Z114*Z115 Z117*Z118 Z119*Z120 Z121*Z122 Z123*Z124 Z126*Z127 Z128*Z129 Z130*Z131 Z132*Z133 Z135*Z136 Z137*Z138 Z139*Z140 Z141*Z142 Z144*Z145 Z146*Z147 Z148*Z149 Z150*Z151 Z153*Z154 Z155*Z156 Z157*Z158 Z159*Z160
    DEPOLARIZE2("""+str(p_f)+""") 0 1 2 3 4 5 6 7 9 10 11 12 13 14 15 16 18 19 20 21 22 23 24 25 27 28 29 30 31 32 33 34 36 37 38 39 40 41 42 43 45 46 47 48 49 50 51 52 54 55 56 57 58 59 60 61 63 64 65 66 67 68 69 70 72 73 74 75 76 77 78 79 81 82 83 84 85 86 87 88 90 91 92 93 94 95 96 97 99 100 101 102 103 104 105 106 108 109 110 111 112 113 114 115 117 118 119 120 121 122 123 124 126 127 128 129 130 131 132 133 135 136 137 138 139 140 141 142 144 145 146 147 148 149 150 151 153 154 155 156 157 158 159 160
    DEPOLARIZE1("""+str(p_f)+""") 8 17 26 35 44 53 62 71 80 89 98 107 116 125 134 143 152


    TICK
    MPP("""+str(p_f2)+""")  Z1*Z2 Z3*Z4 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z14*Z15 Z16*Z17 Z19*Z20 Z21*Z22 Z23*Z24 Z25*Z26 Z28*Z29 Z30*Z31 Z32*Z33 Z34*Z35 Z37*Z38 Z39*Z40 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49 Z50*Z51 Z52*Z53 Z55*Z56 Z57*Z58 Z59*Z60 Z61*Z62 Z64*Z65 Z66*Z67 Z68*Z69 Z70*Z71 Z73*Z74 Z75*Z76 Z77*Z78 Z79*Z80 Z82*Z83 Z84*Z85 Z86*Z87 Z88*Z89 Z91*Z92 Z93*Z94 Z95*Z96 Z97*Z98 Z100*Z101 Z102*Z103 Z104*Z105 Z106*Z107 Z109*Z110 Z111*Z112 Z113*Z114 Z115*Z116 Z118*Z119 Z120*Z121 Z122*Z123 Z124*Z125 Z127*Z128 Z129*Z130 Z131*Z132 Z133*Z134 Z136*Z137 Z138*Z139 Z140*Z141 Z142*Z143 Z145*Z146 Z147*Z148 Z149*Z150 Z151*Z152 Z154*Z155 Z156*Z157 Z158*Z159 Z160*Z161
    DEPOLARIZE2("""+str(p_f)+""") 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 35 37 38 39 40 41 42 43 44 46 47 48 49 50 51 52 53 55 56 57 58 59 60 61 62 64 65 66 67 68 69 70 71 73 74 75 76 77 78 79 80 82 83 84 85 86 87 88 89 91 92 93 94 95 96 97 98 100 101 102 103 104 105 106 107 109 110 111 112 113 114 115 116 118 119 120 121 122 123 124 125 127 128 129 130 131 132 133 134 136 137 138 139 140 141 142 143 145 146 147 148 149 150 151 152 154 155 156 157 158 159 160 161
    DEPOLARIZE1("""+str(p_f)+""") 0 9 18 27 36 45 54 63 72 81 90 99 108 117 126

    SHIFT_COORDS(0, 0, 1)

    DETECTOR(0.5, 0, 0) rec[-576] rec[-575] rec[-574] rec[-573] rec[-572] rec[-571] rec[-570] rec[-569] rec[-568] rec[-288] rec[-287] rec[-286] rec[-285] rec[-284] rec[-283] rec[-282] rec[-281] rec[-280]
    DETECTOR(1.5, 0, 0) rec[-504] rec[-503] rec[-502] rec[-501] rec[-500] rec[-499] rec[-498] rec[-497] rec[-496] rec[-216] rec[-215] rec[-214] rec[-213] rec[-212] rec[-211] rec[-210] rec[-209] rec[-208]
    DETECTOR(2.5, 0, 0) rec[-567] rec[-566] rec[-565] rec[-564] rec[-563] rec[-562] rec[-561] rec[-560] rec[-559] rec[-279] rec[-278] rec[-277] rec[-276] rec[-275] rec[-274] rec[-273] rec[-272] rec[-271]
    DETECTOR(3.5, 0, 0) rec[-495] rec[-494] rec[-493] rec[-492] rec[-491] rec[-490] rec[-489] rec[-488] rec[-487] rec[-207] rec[-206] rec[-205] rec[-204] rec[-203] rec[-202] rec[-201] rec[-200] rec[-199]
    DETECTOR(4.5, 0, 0) rec[-558] rec[-557] rec[-556] rec[-555] rec[-554] rec[-553] rec[-552] rec[-551] rec[-550] rec[-270] rec[-269] rec[-268] rec[-267] rec[-266] rec[-265] rec[-264] rec[-263] rec[-262]
    DETECTOR(5.5, 0, 0) rec[-486] rec[-485] rec[-484] rec[-483] rec[-482] rec[-481] rec[-480] rec[-479] rec[-478] rec[-198] rec[-197] rec[-196] rec[-195] rec[-194] rec[-193] rec[-192] rec[-191] rec[-190]
    DETECTOR(6.5, 0, 0) rec[-549] rec[-548] rec[-547] rec[-546] rec[-545] rec[-544] rec[-543] rec[-542] rec[-541] rec[-261] rec[-260] rec[-259] rec[-258] rec[-257] rec[-256] rec[-255] rec[-254] rec[-253]
    DETECTOR(7.5, 0, 0) rec[-477] rec[-476] rec[-475] rec[-474] rec[-473] rec[-472] rec[-471] rec[-470] rec[-469] rec[-189] rec[-188] rec[-187] rec[-186] rec[-185] rec[-184] rec[-183] rec[-182] rec[-181]
    DETECTOR(9.5, 0, 0) rec[-468] rec[-467] rec[-466] rec[-465] rec[-464] rec[-463] rec[-462] rec[-461] rec[-460] rec[-180] rec[-179] rec[-178] rec[-177] rec[-176] rec[-175] rec[-174] rec[-173] rec[-172] 
    DETECTOR(10.5, 0, 0) rec[-540] rec[-539] rec[-538] rec[-537] rec[-536] rec[-535] rec[-534] rec[-533] rec[-532] rec[-252] rec[-251] rec[-250] rec[-249] rec[-248] rec[-247] rec[-246] rec[-245] rec[-244]
    DETECTOR(11.5, 0, 0) rec[-459] rec[-458] rec[-457] rec[-456] rec[-455] rec[-454] rec[-453] rec[-452] rec[-451] rec[-171] rec[-170] rec[-169] rec[-168] rec[-167] rec[-166] rec[-165] rec[-164] rec[-163]
    DETECTOR(12.5, 0, 0) rec[-531] rec[-530] rec[-529] rec[-528] rec[-527] rec[-526] rec[-525] rec[-524] rec[-523] rec[-243] rec[-242] rec[-241] rec[-240] rec[-239] rec[-238] rec[-237] rec[-236] rec[-235]
    DETECTOR(13.5, 0, 0) rec[-450] rec[-449] rec[-448] rec[-447] rec[-446] rec[-445] rec[-444] rec[-443] rec[-442] rec[-162] rec[-161] rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-155] rec[-154]
    DETECTOR(14.5, 0, 0) rec[-522] rec[-521] rec[-520] rec[-519] rec[-518] rec[-517] rec[-516] rec[-515] rec[-514] rec[-234] rec[-233] rec[-232] rec[-231] rec[-230] rec[-229] rec[-228] rec[-227] rec[-226]
    DETECTOR(15.5, 0, 0) rec[-441] rec[-440] rec[-439] rec[-438] rec[-437] rec[-436] rec[-435] rec[-434] rec[-433] rec[-153] rec[-152] rec[-151] rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-145]
    DETECTOR(16.5, 0, 0) rec[-513] rec[-512] rec[-511] rec[-510] rec[-509] rec[-508] rec[-507] rec[-506] rec[-505] rec[-225] rec[-224] rec[-223] rec[-222] rec[-221] rec[-220] rec[-219] rec[-218] rec[-217]
    DETECTOR(0, 0.5, 0) rec[-432] rec[-428] rec[-424] rec[-420] rec[-416] rec[-412] rec[-408] rec[-404] rec[-400] rec[-144] rec[-140] rec[-136] rec[-132] rec[-128] rec[-124] rec[-120] rec[-116] rec[-112]
    DETECTOR(0, 1.5, 0) rec[-360] rec[-356] rec[-352] rec[-348] rec[-344] rec[-340] rec[-336] rec[-332] rec[-328] rec[-72] rec[-68] rec[-64] rec[-60] rec[-56] rec[-52] rec[-48] rec[-44] rec[-40]
    DETECTOR(0, 2.5, 0) rec[-431] rec[-427] rec[-423] rec[-419] rec[-415] rec[-411] rec[-407] rec[-403] rec[-399] rec[-143] rec[-139] rec[-135] rec[-131] rec[-127] rec[-123] rec[-119] rec[-115] rec[-111]
    DETECTOR(0, 3.5, 0) rec[-359] rec[-355] rec[-351] rec[-347] rec[-343] rec[-339] rec[-335] rec[-331] rec[-327] rec[-71] rec[-67] rec[-63] rec[-59] rec[-55] rec[-51] rec[-47] rec[-43] rec[-39]
    DETECTOR(0, 4.5, 0) rec[-430] rec[-426] rec[-422] rec[-418] rec[-414] rec[-410] rec[-406] rec[-402] rec[-398] rec[-142] rec[-138] rec[-134] rec[-130] rec[-126] rec[-122] rec[-118] rec[-114] rec[-110]
    DETECTOR(0, 5.5, 0) rec[-358] rec[-354] rec[-350] rec[-346] rec[-342] rec[-338] rec[-334] rec[-330] rec[-326] rec[-70] rec[-66] rec[-62] rec[-58] rec[-54] rec[-50] rec[-46] rec[-42] rec[-38]
    DETECTOR(0, 6.5, 0) rec[-429] rec[-425] rec[-421] rec[-417] rec[-413] rec[-409] rec[-405] rec[-401] rec[-397] rec[-141] rec[-137] rec[-133] rec[-129] rec[-125] rec[-121] rec[-117] rec[-113] rec[-109]
    DETECTOR(0, 7.5, 0) rec[-357] rec[-353] rec[-349] rec[-345] rec[-341] rec[-337] rec[-333] rec[-329] rec[-325] rec[-69] rec[-65] rec[-61] rec[-57] rec[-53] rec[-49] rec[-45] rec[-41] rec[-37]
    DETECTOR(9, 0.5, 0) rec[-396] rec[-392] rec[-388] rec[-384] rec[-380] rec[-376] rec[-372] rec[-368] rec[-364] rec[-108] rec[-104] rec[-100] rec[-96] rec[-92] rec[-88] rec[-84] rec[-80] rec[-76]
    DETECTOR(9, 1.5, 0) rec[-324] rec[-320] rec[-316] rec[-312] rec[-308] rec[-304] rec[-300] rec[-296] rec[-292] rec[-36] rec[-32] rec[-28] rec[-24] rec[-20] rec[-16] rec[-12] rec[-8] rec[-4]
    DETECTOR(9, 2.5, 0) rec[-395] rec[-391] rec[-387] rec[-383] rec[-379] rec[-375] rec[-371] rec[-367] rec[-363] rec[-107] rec[-103] rec[-99] rec[-95] rec[-91] rec[-87] rec[-83] rec[-79] rec[-75]
    DETECTOR(9, 3.5, 0) rec[-323] rec[-319] rec[-315] rec[-311] rec[-307] rec[-303] rec[-299] rec[-295] rec[-291] rec[-35] rec[-31] rec[-27] rec[-23] rec[-19] rec[-15] rec[-11] rec[-7] rec[-3]
    DETECTOR(9, 4.5, 0) rec[-394] rec[-390] rec[-386] rec[-382] rec[-378] rec[-374] rec[-370] rec[-366] rec[-362] rec[-106] rec[-102] rec[-98] rec[-94] rec[-90] rec[-86] rec[-82] rec[-78] rec[-74]
    DETECTOR(9, 5.5, 0) rec[-322] rec[-318] rec[-314] rec[-310] rec[-306] rec[-302] rec[-298] rec[-294] rec[-290] rec[-34] rec[-30] rec[-26] rec[-22] rec[-18] rec[-14] rec[-10] rec[-6] rec[-2]
    DETECTOR(9, 6.5, 0) rec[-393] rec[-389] rec[-385] rec[-381] rec[-377] rec[-373] rec[-369] rec[-365] rec[-361] rec[-105] rec[-101] rec[-97] rec[-93] rec[-89] rec[-85] rec[-81] rec[-77] rec[-73]
    DETECTOR(9, 7.5, 0) rec[-321] rec[-317] rec[-313] rec[-309] rec[-305] rec[-301] rec[-297] rec[-293] rec[-289] rec[-33] rec[-29] rec[-25] rec[-21] rec[-17] rec[-13] rec[-9] rec[-5] rec[-1]

}

TICK
MX("""+str(p_f2)+""")  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161

DETECTOR(0.5, 0, 0) rec[-450] rec[-449] rec[-448] rec[-447] rec[-446] rec[-445] rec[-444] rec[-443] rec[-442] rec[-162] rec[-161] rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-145]
DETECTOR(1.5, 0, 0) rec[-378] rec[-377] rec[-376] rec[-375] rec[-374] rec[-373] rec[-372] rec[-371] rec[-370] rec[-153] rec[-152] rec[-151] rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-140] rec[-139] rec[-138] rec[-137] rec[-136]
DETECTOR(2.5, 0, 0) rec[-441] rec[-440] rec[-439] rec[-438] rec[-437] rec[-436] rec[-435] rec[-434] rec[-433] rec[-144] rec[-143] rec[-142] rec[-141] rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-130] rec[-129] rec[-128] rec[-127]
DETECTOR(3.5, 0, 0) rec[-369] rec[-368] rec[-367] rec[-366] rec[-365] rec[-364] rec[-363] rec[-362] rec[-361] rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-120] rec[-119] rec[-118]
DETECTOR(4.5, 0, 0) rec[-432] rec[-431] rec[-430] rec[-429] rec[-428] rec[-427] rec[-426] rec[-425] rec[-424] rec[-126] rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-120] rec[-119] rec[-118] rec[-117] rec[-116] rec[-115] rec[-114] rec[-113] rec[-112] rec[-111] rec[-110] rec[-109]
DETECTOR(5.5, 0, 0) rec[-360] rec[-359] rec[-358] rec[-357] rec[-356] rec[-355] rec[-354] rec[-353] rec[-352] rec[-117] rec[-116] rec[-115] rec[-114] rec[-113] rec[-112] rec[-111] rec[-110] rec[-109] rec[-108] rec[-107] rec[-106] rec[-105] rec[-104] rec[-103] rec[-102] rec[-101] rec[-100]
DETECTOR(6.5, 0, 0) rec[-423] rec[-422] rec[-421] rec[-420] rec[-419] rec[-418] rec[-417] rec[-416] rec[-415] rec[-108] rec[-107] rec[-106] rec[-105] rec[-104] rec[-103] rec[-102] rec[-101] rec[-100] rec[-99] rec[-98] rec[-97] rec[-96] rec[-95] rec[-94] rec[-93] rec[-92] rec[-91]
DETECTOR(7.5, 0, 0) rec[-351] rec[-350] rec[-349] rec[-348] rec[-347] rec[-346] rec[-345] rec[-344] rec[-343] rec[-99] rec[-98] rec[-97] rec[-96] rec[-95] rec[-94] rec[-93] rec[-92] rec[-91] rec[-90] rec[-89] rec[-88] rec[-87] rec[-86] rec[-85] rec[-84] rec[-83] rec[-82]
DETECTOR(9.5, 0, 0) rec[-342] rec[-341] rec[-340] rec[-339] rec[-338] rec[-337] rec[-336] rec[-335] rec[-334] rec[-81] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66] rec[-65] rec[-64]
DETECTOR(10.5, 0, 0) rec[-414] rec[-413] rec[-412] rec[-411] rec[-410] rec[-409] rec[-408] rec[-407] rec[-406] rec[-72] rec[-71] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56] rec[-55]
DETECTOR(11.5, 0, 0) rec[-333] rec[-332] rec[-331] rec[-330] rec[-329] rec[-328] rec[-327] rec[-326] rec[-325] rec[-63] rec[-62] rec[-61] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
DETECTOR(12.5, 0, 0) rec[-405] rec[-404] rec[-403] rec[-402] rec[-401] rec[-400] rec[-399] rec[-398] rec[-397] rec[-54] rec[-53] rec[-52] rec[-51] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41] rec[-40] rec[-39] rec[-38] rec[-37]
DETECTOR(13.5, 0, 0) rec[-324] rec[-323] rec[-322] rec[-321] rec[-320] rec[-319] rec[-318] rec[-317] rec[-316] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41] rec[-40] rec[-39] rec[-38] rec[-37] rec[-36] rec[-35] rec[-34] rec[-33] rec[-32] rec[-31] rec[-30] rec[-29] rec[-28]
DETECTOR(14.5, 0, 0) rec[-396] rec[-395] rec[-394] rec[-393] rec[-392] rec[-391] rec[-390] rec[-389] rec[-388] rec[-36] rec[-35] rec[-34] rec[-33] rec[-32] rec[-31] rec[-30] rec[-29] rec[-28] rec[-27] rec[-26] rec[-25] rec[-24] rec[-23] rec[-22] rec[-21] rec[-20] rec[-19]
DETECTOR(15.5, 0, 0) rec[-315] rec[-314] rec[-313] rec[-312] rec[-311] rec[-310] rec[-309] rec[-308] rec[-307] rec[-27] rec[-26] rec[-25] rec[-24] rec[-23] rec[-22] rec[-21] rec[-20] rec[-19] rec[-18] rec[-17] rec[-16] rec[-15] rec[-14] rec[-13] rec[-12] rec[-11] rec[-10]
DETECTOR(16.5, 0, 0) rec[-387] rec[-386] rec[-385] rec[-384] rec[-383] rec[-382] rec[-381] rec[-380] rec[-379] rec[-18] rec[-17] rec[-16] rec[-15] rec[-14] rec[-13] rec[-12] rec[-11] rec[-10] rec[-9] rec[-8] rec[-7] rec[-6] rec[-5] rec[-4] rec[-3] rec[-2] rec[-1]

OBSERVABLE_INCLUDE(0) rec[-90] rec[-89] rec[-88] rec[-87] rec[-86] rec[-85] rec[-84] rec[-83] rec[-82]
OBSERVABLE_INCLUDE(1) rec[-81] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76] rec[-75] rec[-74] rec[-73]""")
        


        # Verification if there are detectors missing or something wrong going on:--------------------------
        # Check conditions for each code:

        print("Count detector sizes, expected and true :=====================================================")
        a=count_determined_measurements_in_circuit(circ5d)
        expected_determined = circ5d.num_detectors + circ5d.num_observables
        print(expected_determined,a) #expected, true

        a=count_determined_measurements_in_circuit(circuit3up)
        expected_determined = circuit3up.num_detectors + circuit3up.num_observables
        print(expected_determined,a) #expected, true

        a=count_determined_measurements_in_circuit(circsurf_notr)
        expected_determined = circsurf_notr.num_detectors + circsurf_notr.num_observables
        print(expected_determined,a) #expected, true

        a=count_determined_measurements_in_circuit(circsurf_notr2)
        expected_determined = circsurf_notr2.num_detectors + circsurf_notr2.num_observables
        print(expected_determined,a) #expected, true


        a=count_determined_measurements_in_circuit(circbs5)
        expected_determined = circbs5.num_detectors + circbs5.num_observables
        print(expected_determined,a) #expected, true

        sampler2 = circuit2.compile_detector_sampler()
        S_Bell.append(np.sum(sampler2.sample(shots=10**6)) / 10**6)


        print("Distances of error :=====================================================")
        actual_distance5 = len(circ5d.shortest_graphlike_error())
        actual_distance1 = len(circsurf_notr.shortest_graphlike_error())
        actual_distance2 = len(circsurf_notr2.shortest_graphlike_error())
        actual_distance4= len(circbs5.shortest_graphlike_error())
        actual_distance3 = len(circuit3up.shortest_graphlike_error())


        print("distance 3?")
        print(actual_distance1,actual_distance2,actual_distance3)
        if (actual_distance1 and actual_distance2 and actual_distance3 !=3):
            raise ValueError('Distance !=3, increase/decrease lattice dimensions')


        print("distance 5?")
        if (actual_distance4 and actual_distance5 !=5):
            raise ValueError('Distance !=5, increase/decrease lattice dimensions')
        print(actual_distance4,actual_distance5)


        print(" ------------------------------------------------------------------------ ")
        print("End of verification...")


        ############################################################

    
        surface_code_tasks6 = [
            sinter.Task(
            circuit = circuit3up,
            json_metadata={'d': 3, 'p': p_f},
            )#,



        ]

        surface_code_tasks4 = [
            sinter.Task(
            circuit = circ5d,
            json_metadata={'d': 5, 'p': p_f},
            )#,


    
        ]



       

        surface_code_tasks2 = [
            sinter.Task(
            circuit = circsurf_notr,
            json_metadata={'d': 3, 'p': p_f},
            )#,


    
        ]

        surface_code_tasks3 = [
            sinter.Task(
            circuit = circsurf_notr2,
            json_metadata={'d': 3, 'p': p_f},
            )#,
            

    
        ]


        surface_code_tasks5 = [
            sinter.Task(
            circuit = circbs5,
            json_metadata={'d': 5, 'p': p_f},
            )#,
            

    
        ]


        collected_stats6: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks6,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True



            )
   

        collected_stats2: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks2,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True


        
            )

        collected_stats3: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks3,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True
            )


        collected_stats4: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks4,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True
            )



        collected_stats5: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks5,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True
            )

        


        X2,Y2,Y_l2,Y_h2=plot_error_rate2(
    stats=collected_stats2,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y2.append(Y2[0])
        Std_p2.append(Y_h2[0])
        Std_m2.append(Y_l2[0])



        X3,Y3,Y_l3,Y_h3=plot_error_rate2(
    stats=collected_stats3,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y3.append(Y3[0])
        Std_p3.append(Y_h3[0])
        Std_m3.append(Y_l3[0])


        X4,Y4,Y_l4,Y_h4=plot_error_rate2(
    stats=collected_stats4,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y4.append(Y4[0])
        Std_p4.append(Y_h4[0])
        Std_m4.append(Y_l4[0])


        X5,Y5,Y_l5,Y_h5=plot_error_rate2(
    stats=collected_stats5,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y5.append(Y5[0])
        Std_p5.append(Y_h5[0])
        Std_m5.append(Y_l5[0])



        X6,Y6,Y_l6,Y_h6=plot_error_rate2(
    stats=collected_stats6,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y6.append(Y6[0])
        Std_p6.append(Y_h6[0])
        Std_m6.append(Y_l6[0])


    
    #Data output for all the codes for the 1QEC static noise case
    #Physical Error, Logical error, Standard deviation lower bound, Standard deviation upper bound
    #Vec_Xx,         Vec_Yx,                        Std_mx,                       Std_px
    print("S3")
    print(Vec_X2,Vec_Y2, Std_m2, Std_p2)
    print("rS3")
    print(Vec_X3,Vec_Y3, Std_m3, Std_p3)
    print("rS5")
    print(Vec_X4,Vec_Y4, Std_m4, Std_p4)
    print("BS5")
    print(Vec_X5,Vec_Y5, Std_m5, Std_p5)
    print("BS3")
    print(Vec_X6,Vec_Y6, Std_m6, Std_p6)
    print("Bell pairs")
    print(S_Bell)



    fig, ax = plt.subplots(1, 1)

    ax.fill_between(Vec_X, Std_m6, Std_p6,alpha=0.2,color="black")
    ax.plot(Vec_X,Vec_Y6,"v",linestyle='-',color="black",label="BS [[18,2,3]]")

    ax.fill_between(Vec_X, Std_m5, Std_p5,alpha=0.2,color="orange")
    ax.plot(Vec_X,Vec_Y5,"v",linestyle='-',color="orange",label="BS [[50,2,5]]")

    ax.fill_between(Vec_X3, Std_m3, Std_p3,alpha=0.2,color="green")
    ax.plot(Vec_X3,Vec_Y3,"v",linestyle='-',color="green",label="R-S [[18,2,3]]")

    ax.fill_between(Vec_X, Std_m4, Std_p4,alpha=0.2,color="purple")
    ax.plot(Vec_X,Vec_Y4,"v",linestyle='-',color="purple",label="R-S [[50,2,5]]")

    ax.fill_between(Vec_X2, Std_m2, Std_p2,alpha=0.2,color="blue")
    ax.plot(Vec_X2,Vec_Y2,"v",linestyle='-',color="blue",label="S [[18,2,3]]")

    ax.plot(Vec_X,S_Bell,"v",linestyle='-',color="red",label="Unencoded Bell-State")

    ax.grid(linestyle='--')
    plt.legend(loc="lower right")
    ax.set_xlabel("Physical Error Rate",fontsize=14)
    ax.set_ylabel("Logical Error Rate",fontsize=14)
    plt.yscale("log")   
    plt.xscale("log")   

    plt.show()