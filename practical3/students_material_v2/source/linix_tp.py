#!/usr/bin/python3.7
# -*- coding: utf-8 -*-
""" Linix TP - linix_tp main file

This file is made to work with the LiniX TP and provides tools to let students
experiments with the moteus MJBots controller and a belt driven linear axis.

Author :                Munier Louis
Work :                  LiniX TP
Last Modification :     lmr210331 - 1815
"""
import asyncio
import math
import time
import moteus

import debug_plot as dp
import signal_generator as sg


# Provide some constants to avoid magic numbers
POS_EPSILON = 5e-2  # [turn] precision reached pos before starting moving
TIME_STOP = 0.5  # [s] time to smoothely stop the motor after sended signal
DEFAULT_POS_BOUNDARY = 1.0  # [turn] max allowed motor turn in position mode
MAX_TORQUE = 1.7  # [Nm] max allowed torque
MAX_FREQUENCY = 4  # [Hz] max allowed frequency
MAX_FREQUENCY_POS = 1  # [Hz] max allowed rectangular freq in position mode
SIM_TIMEOUT = 1 # [s] After sim_time + SIM_TIMEOUT, simulation is stopped to avoid running for ever


async def main(**params):
    """ Main function to experiment during LiniX TP hours.

    Can be fully modified to obtain needed behavior and shows different plots.

    Args :
        reset : set it to True if original parameters has to be set back
    """

    # Create needed variables
    sim_time = params['period_time'] * params['nb_periods']

    # Create controller and stream, objects then do a stop command in case the
    # controller had faulted previously, send the stop command to clear it.
    c = moteus.Controller()
    s = moteus.Stream(c)
    await c.set_stop()

    # Verify some bounded parameters
    check_params(**params)

    # Check if the parameters has to be resetted to factory defaults
    if params['mode'] == 'reset':
        await reset_params(s)
        await c.set_stop()

        print('Parameters resetted !')
        return

    # Create specified signal generator and generate a new signal
    sgn = generate_signal(**params)

    if params['mode'] == 'stop':
        pass
    elif params['mode'] == 'run':
        # Set controller params
        await set_params(s)

        # Reset motor position
        state = await c.set_position(query=True)
        if abs(state.values[moteus.Register.POSITION]) > DEFAULT_POS_BOUNDARY:
            await c.set_stop()
            raise Exception('Error : Position out of boundaries.')

        while abs(state.values[moteus.Register.POSITION]) > POS_EPSILON:
            state = await c.set_position(position=0, velocity=0.5, maximum_torque=0.04, query=True)

        # Send motor command
        states = await send_command(
            c,
            params['command_type'],
            sgn.get_signal(),
            sim_time,
            sgn.get_timestep()
        )

        # Retrieve positions and velocities :
        kwargs = {}
        kwargs['command'] = sgn.get_signal()
        kwargs['position'] = [state.values[moteus.Register.POSITION]
                              for state in states]
        if params['translate']:
            d_pulley = 20.27 / 1000

            # Convert from turn to m
            if params['command_type'] == 'position':
                kwargs['command'] = [com * d_pulley *
                                     math.pi for com in kwargs['command']]
            kwargs['position'] = [pos * d_pulley *
                                  math.pi for pos in kwargs['position']]

        kwargs['velocity'] = [state.values[moteus.Register.VELOCITY]
                              for state in states]
        kwargs['torque'] = [state.values[moteus.Register.TORQUE]
                            for state in states]

        # Plot and save it
        if params['save']:
            kwargs['fig_options'] = params['fig_options']

        dp.subplot_1d(
            params['command_type'],
            ['position', 'velocity', 'torque'],
            sgn.get_timestep(),
            axis_label=['Time [s]', 'Position [m]',
                        'Velocity [m/s]', 'Torque [Nm]'],
            save=params['save'],
            **kwargs
        )
    else:
        raise Exception('Error : Unexpected mode.')


async def reset_params(stream: moteus.Stream, save: bool = True):
    """ The goal of this function is to reset the configuration values.

    Args :
        stream : moteus stream connected to a specific motor
        save : if needed, save parameters in permanent memory
    """
    params = {}
    params['servopos.position_min'] = -1.0
    params['servopos.position_max'] = 1.0

    params['servo.pid_position.kp'] = 5.0
    params['servo.pid_position.ki'] = 0.0
    params['servo.pid_position.kd'] = 0.07

    await modify_params(stream, save=save, **params)


def check_params(**params):
    """ The goal of this function is to check some basics parameters.

    Args :
        **params : dictionnary of parameters
    """
    command_types = ['position', 'velocity']
    if params['command_type'] not in command_types:
        raise Exception('Error : Unexpected command type.')

    if params['command_type'] == 'position' and params['amplitude'] > DEFAULT_POS_BOUNDARY:
        raise Exception(
            'Error : Amplitude too hight in position control mode.')

    if params['command_type'] == 'position' and 1 / params['period_time'] > MAX_FREQUENCY_POS and params['signal_type'] == 'rectangular':
        raise Exception(
            'Error : Frequency too hight in rectangular position control mode. Max period : {}'.format(1 / MAX_FREQUENCY_POS))

    if 1 / params['period_time'] > MAX_FREQUENCY:
        raise Exception(
            'Error : Frequency too high, increase period_time parameter.')


def generate_signal(**params):
    """ The goal of this function is to generate the command signal.

    Args :
        **params : dictionnary of parameters
    """
    if params['signal_type'] == 'rectangular':
        sgn = sg.RectangularGenerator()
    elif params['signal_type'] == 'triangular':
        sgn = sg.TriangularGenerator()
    elif params['signal_type'] == 'trapezoidal':
        sgn = sg.TrapezoidalGenerator()
    elif params['signal_type'] == 'sine':
        sgn = sg.SineGenerator()
    else:
        raise Exception('Error : Unexpected signal type.')

    sgn.set(params['period_time'], params['amplitude'])
    sgn.generate(**params)
    sgn.get_signal_filtered()

    return sgn


async def send_command(controller: moteus.Controller, command_type: str, commands: [float],
                       sim_time: float, timestep: float):
    """ Function to send generated commands to the controller.

    Args :
        controller : moteus controller instance
        command_type : type of commands position / velocity
        commands : list of positions / velocities commands
        sim_time : simulation time
        timestep : time between each commands
    """
    # Create variable to store all the states of the controller
    states = []

    # Retrieve time to run controller during a certain period
    start_time = time.time()
    print("In send_command(), timestep : {}".format(timestep))

    for c in commands:
        t_before_pos = time.time()

        if command_type == 'position':
            states.append(await controller.set_position(
                position=c,
                maximum_torque=MAX_TORQUE,
                query=True
            ))
        elif command_type == 'velocity':
            states.append(await controller.set_position(
                position=math.nan,
                velocity=c,
                query=True
            ))

        while time.time() - t_before_pos < timestep:
            await asyncio.sleep(0)

        if time.time() - start_time > timestep * len(commands) + SIM_TIMEOUT:
            break

    # Smoother end for the controller
    end_time = time.time()
    while time.time() - end_time < TIME_STOP:
        if command_type == 'position':
            await controller.set_position(
                position=commands[-1],
                maximum_torque=MAX_TORQUE,
                query=True
            )
        elif command_type == 'velocity':
            await controller.set_position(
                position=math.nan,
                velocity=0.0,
                query=True
            )

    await controller.set_stop()

    return states


async def modify_params(stream: moteus.Stream, save: bool = False, **params: {}):
    """ The goal of this function is to modify configuration values.

    Args :
        stream : moteus stream connected to a specific motor
        save : specify if new parameters has to be saved in constant memory
        **params : dictionnary of new parameters
    """
    for k, v in params.items():
        msg = 'conf set {} {}'.format(k, v).encode('latin1')
        await stream.command(msg)

    if save:
        await stream.command(b'conf write')


async def set_params(stream: moteus.Stream, save: bool = True):
    """ The goal of this function is to reset the configuration values.

    Args :
        stream : moteus stream connected to a specific motor
        save : if needed, save parameters in permanent memory
    """
    params = {}
    params['servopos.position_min'] = -1.0  # default : -1.0
    params['servopos.position_max'] = 1.0 # default : 1.0

    params['servo.pid_position.ilimit'] = 10.0  # default : 5.0
    params['servo.pid_position.kp'] = 30.0  # default : 5.0
    params['servo.pid_position.ki'] = 5.0  # default : 0.0
    params['servo.pid_position.kd'] = 0.035  # default : 0.07

    await modify_params(stream, save=save, **params)


if __name__ == "__main__":
    # Set the different parameters to send to the main
    settings = {}
    
    # min : 0.25 (frequency = 1 / period_time => freq_max = 4Hz)
    settings['period_time'] = 1.5
    settings['nb_periods'] = 3
    settings['amplitude'] = 0.4
    settings['command_type'] = 'position'  # 'position' or 'velocity'
    settings['mode'] = 'run'  # 'stop', 'run', 'reset'

    # to plot the position in degree instead of turn
    settings['translate'] = False

    # Type of signal : 'rectangular', 'triangular', 'trapezoidal', 'sine'
    settings['signal_type'] = 'sine'

    # To show and save plots
    settings['plot'] = True
    settings['filter'] = True  # show filter plot
    settings['save'] = True
    settings['fig_options'] = {
        'signal': {
            'name': 'default_name_signal_plot',
            'extension': 'png'
        },
        'result': {
            'name': 'default_name_result_plot',
            'extension': 'png'
        }
    }

    asyncio.run(main(**settings))
