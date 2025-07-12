#!/usr/bin/python3.7
# -*- coding: utf-8 -*-
""" Linix TP - wave_generator class file

This class file provide a way to have nice waves generated.

Author :                Munier Louis
Work :                  LiniX TP
Last Modification :     lmr210331 - 1618
"""
import math
import matplotlib.pyplot as plt

from debug_plot import save_figure


class SignalGenerator:
    """ Signal generator base class

    It provides all the needed tools to have nice control signals.
    """
    PERIOD = None
    STEP = None
    STEP_MIN = 0.007  # minimal timestep : 0.0016307s
    AMPLITUDE = None
    MAX_VELOCITY = None
    MAX_ACCELERATION = None

    SIGNAL = []
    SIGNAL_FILTERED = []

    def __init__(self, m_max_velocity: float = None, m_max_acceleration: float = None):
        """ Function to initialize signal generator with basic constants.

        Args :
            m_max_velocity : maximum motor velocity
            m_max_acceleration : maximum motor acceleration
        """
        self.MAX_VELOCITY = m_max_velocity
        self.MAX_ACCELERATION = m_max_acceleration

    def _generate(self, **params: {}):
        """ Increase possibilities of daughter's functions with some transmitted params.

        Args :
            **params : dictionnary with all the needed parameters to set class, save
                       figure for example
        """
        if params['nb_periods'] > 1:
            self._repeat_period(nb_periods=params['nb_periods'])

        if params['filter']:
            self._filter()

        if params['plot']:
            self._draw(**params)

    def _repeat_period(self, nb_periods: int = 1):
        """ Function to repeat basic signal as many as needed.

        Args :
            nb_periods : number of needed total periods
        """
        sgn = self.SIGNAL.copy()

        for _ in range(nb_periods - 1):
            self.SIGNAL.extend(sgn)

        self.SIGNAL.append(sgn[0])

    def _filter(self):
        """ Function to generate first order signal and plot it if needed. """
        alpha = self.STEP / (self.PERIOD + self.STEP)

        for i in range(len(self.SIGNAL)):
            if i == 0:
                self.SIGNAL_FILTERED.append(alpha * self.SIGNAL[i])
            else:
                self.SIGNAL_FILTERED.append(
                    (1 - alpha) * self.SIGNAL_FILTERED[-1] + alpha * self.SIGNAL[i])

    def _draw(self, **params: {}):
        """ Function to draw generated signal.

        Args :
            **params : dictionnary with all the needed parameters to set class, save
                       figure for example
        """
        # Create and show plot
        fig = plt.figure(figsize=(8, 6), dpi=150)

        times = [i * self.STEP for i in range(len(self.SIGNAL))]
        plt.plot(times, self.SIGNAL, label='Signal')

        if params['filter']:
            times = [i * self.STEP for i in range(len(self.SIGNAL))]
            plt.plot(times, self.SIGNAL_FILTERED, label='1st order filtering')

        plt.xlabel("Time (s)")
        plt.ylabel("Signal (-)")
        plt.title("Signal over time")
        plt.grid(True)
        plt.legend()
        plt.show()

        if params['save']:
            save_figure(fig, **params['fig_options']['signal'])

    def set(self, period: float, amplitude: float):
        """ Set values inside classe's constants.

        Args :
            period : period of the generated signal
            amplitude : amplitude of the generated signal
        """
        self.PERIOD = period
        self.AMPLITUDE = amplitude

        self.STEP = self.PERIOD
        while self.STEP >= 2*self.STEP_MIN:
            self.STEP /= 2

    def get_signal(self):
        """ Return generated signal command. """
        return self.SIGNAL

    def get_signal_filtered(self):
        """ Return filtered generated signal command. """
        if len(self.SIGNAL_FILTERED) == 0:
            self._filter()

        return self.SIGNAL_FILTERED

    def get_timestep(self):
        """ Return generated signal timestep. """
        return self.STEP


class RectangularGenerator(SignalGenerator):
    """ Rectangular generator class

    This derived class implements rectangular generator function.
    """

    def generate(self, **params: {}):
        """ Function to generate rectangular signal and plot it if needed.

        Args :
            **params : dictionnary with all the needed parameters to set class, save
                       figure for example
        """
        times = [i * self.STEP for i in range(int(self.PERIOD / self.STEP))]

        self.SIGNAL.clear()
        self.SIGNAL = [self.AMPLITUDE for i in times[0: int(len(times) / 2)]]
        self.SIGNAL.extend([-x for x in self.SIGNAL])

        super()._generate(**params)


class TriangularGenerator(SignalGenerator):
    """ Triangle generator class

    This derived class implements Triangle generator function.
    """
    N_SECTION = 4

    def generate(self, **params: {}):
        """ Function to generate triangle signal and plot it if needed.

        Args :
            **params : dictionnary with all the needed parameters to set class, save
                       figure for example
        """
        times = [i * self.STEP for i in range(int(self.PERIOD / self.STEP))]
        mid_point = int(len(times) / 2)
        len_section = int(len(times) / self.N_SECTION)

        self.SIGNAL.clear()
        self.SIGNAL = [self.N_SECTION * t * self.AMPLITUDE /
                       self.PERIOD for t in times[0 : len_section]]
        self.SIGNAL.extend([self.AMPLITUDE * self.N_SECTION * (0.5 - t / self.PERIOD)
                           for t in times[len_section : mid_point]])
        self.SIGNAL.extend([-x for x in self.SIGNAL])

        super()._generate(**params)


class TrapezoidalGenerator(SignalGenerator):
    """ Trapezoidal generator class

    This derived class implements Trapezoidal generator function.
    """
    N_SECTION = 6

    def generate(self, **params: {}):
        """ Function to generate trapezoidal signal and plot it if needed.

        Args :
            **params : dictionnary with all the needed parameters to set class, save
                       figure for example
        """
        times = [i * self.STEP for i in range(int(self.PERIOD / self.STEP))]
        mid_point = int(len(times) / 2)
        len_section = int(len(times) / self.N_SECTION)

        self.SIGNAL.clear()
        self.SIGNAL = [self.N_SECTION * t * self.AMPLITUDE /
                       self.PERIOD for t in times[0 : len_section]]
        self.SIGNAL.extend(
            [self.AMPLITUDE for t in times[len_section : -(mid_point + len_section)]])
        self.SIGNAL.extend([self.AMPLITUDE * self.N_SECTION * (0.5 - t / self.PERIOD)
                           for t in times[-(mid_point + len_section) : mid_point]])
        self.SIGNAL.extend([-x for x in self.SIGNAL])

        super()._generate(**params)


class SineGenerator(SignalGenerator):
    """ Sine generator class

    This derived class implements sine generator function.
    """

    def generate(self, **params: {}):
        """ Function to generate sine signal and plot it if needed.

        Args :
            **params : dictionnary with all the needed parameters to set class, save
                       figure for example
        """
        angles = [2 * i * math.pi * self.STEP /
                  self.PERIOD for i in range(int(self.PERIOD / self.STEP))]

        self.SIGNAL.clear()
        self.SIGNAL = [self.AMPLITUDE * math.sin(a) for a in angles]

        super()._generate(**params)
