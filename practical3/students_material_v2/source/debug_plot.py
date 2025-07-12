#!/usr/bin/python3.7
# -*- coding: utf-8 -*-
""" Linix TP - debug_plot file

This file is made to work with the LiniX TP and provides tools to let students
create useful plots to populate report that will serve to evaluate their knowledges.

Author :                Munier Louis
Work :                  LiniX TP
Last Modification :     lmr210331 - 1640
"""
import matplotlib.pyplot as plt


def plot_1d_1scale(command_type: str, func_labels: [str], timestep: float, plot_title: str = None,
                   axis_label: [str] = None, save: bool = False, **params: {}):
    """ This function is used to easily plot 1d functions.

    Args :
        command_type : type of desired command (position or velocity)
        func_labels : legend of each plot, except position / velocity commands
        timestep : if needed, the timestep between each position values
        plot_title : title of the plot, if Non, default one is used
        axis_label : label for each axis, if None, default ones are used
        save : if needed the plot is saved in the current folder
        **params : values to plot (desired command, position, velocity ...)
    """

    # Fill necessary variables with default value if not specify by inputs
    if axis_label is None:
        axis_label = ["x axis [s]", "y axis [m]"]

    if plot_title is None:
        plot_title = "Default title."

    times = [i * timestep for i in range(len(params['commands']))]

    # Create and show plot
    fig = plt.figure(figsize=(8, 6), dpi=150)
    plt.tight_layout()
    plt.grid(True)

    l = 0
    for k, v in params.items():
        if k == 'commands':
            plt.plot(times, v, label="Desired {}".format(command_type))
        else:
            plt.plot(times, v, label=func_labels[l])
            l += 1

    plt.title(plot_title)
    plt.xlabel(axis_label[0])
    plt.ylabel(axis_label[1])
    plt.legend()

    # Show and save plot in file
    plt.show()

    if save:
        save_figure(fig, **params['fig_options']['result'])


def plot_1d_2scales(command_type: str, func_labels: [str], timestep: float, plot_title: str = None,
                    axis_label: [str] = None, save: bool = False, **params: {}):
    """ This function is used to easily plot 1d functions with two different scales.

    Args :
        command_type : type of desired command (position or velocity)
        func_labels : legend of each plot, except position / velocity commands
        timestep : if needed, the timestep between each position values
        plot_title : title of the plot, if Non, default one is used
        axis_label : label for each axis, if None, default ones are used
        save : if needed the plot is saved in the current folder
        **params : values to plot (desired command, position, velocity ...)
    """

    # Fill necessary variables with default value if not specify by inputs
    if axis_label is None:
        axis_label = ["x axis [s]", "y1 axis [m]", "y2 axis [m/s]"]

    if plot_title is None:
        plot_title = "Default title."

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    fig.tight_layout()
    plt.grid(True)

    plots = []
    color = ['tab:blue', 'tab:green', 'tab:red', 'tab:orange']
    for i, fl in enumerate(func_labels):
        if i == 0:
            ax.set_xlabel(axis_label[i])
        elif i == 1:
            ax = ax.twinx()

        ax.set_ylabel(axis_label[i + 1])

        if command_type == fl:
            time = [i * timestep for i in range(len(params['command']))]
            p, = ax.plot(
                time,
                params['command'],
                label='Desired ' + fl,
                color=color[-1]
            )

            plots.append(p)

        time = [i * timestep for i in range(len(params[fl]))]
        p, = ax.plot(time, params[fl], label='Achieved ' + fl, color=color[i])
        ax.yaxis.label.set_color(p.get_color())
        plots.append(p)

    ax.legend(handles=plots)

    # Show and save plot in file
    plt.show()

    if save:
        save_figure(fig, **params['fig_options']['result'])


def subplot_1d(command_type: str, func_labels: [str], timestep: float,
               axis_label: [str] = None, save: bool = False, **params: {}):
    """ This function is used to easily plot 1d functions with two different scales.

    Args :
        command_type : type of desired command (position or velocity)
        func_labels : legend of each plot, except position / velocity commands
        timestep : if needed, the timestep between each position values
        axis_label : label for each axis, if None, default ones are used
        save : if needed the plot is saved in the current folder
        **params : values to plot (desired command, position, velocity ...)
    """

    # Fill necessary variables with default value if not specify by inputs
    if axis_label is None:
        axis_label = ["x axis [s]", "y1 axis [m]",
                      "y2 axis [m/s]"]

    # Create plot
    fig, axs = plt.subplots(len(func_labels), figsize=(8, 6), dpi=150)
    fig.tight_layout()

    color = ['tab:blue', 'tab:green', 'tab:orange']
    for i, fl in enumerate(func_labels):
        axs[i].grid(True)
        axs[i].set_title(fl[0].upper() + fl[1:] + ' over time')
        axs[i].set_ylabel(axis_label[i + 1])

        if fl == func_labels[-1]:
            axs[i].set_xlabel(axis_label[0])

        if command_type == fl:
            time = [i * timestep for i in range(len(params['command']))]
            axs[i].plot(
                time,
                params['command'],
                label='Desired ' + fl,
                color=color[-1]
            )

        time = [i * timestep for i in range(len(params[fl]))]
        axs[i].plot(time, params[fl], label='Achieved ' + fl, color=color[i])
        axs[i].legend()

    # Show and save plot in file
    plt.show()

    if save:
        save_figure(fig, **params['fig_options']['result'])


def save_figure(figure: object, **fig_options: {}):
    """
    This function is used to easily save figure on computer.

    Args :
        figure : the figure plot to save
        **fig_options : dictionnary with different options to save the plotted figure
    """
    name = fig_options['name']
    extension = fig_options['extension']

    name = name.replace(" ", "_").replace(".", "dot")
    name = "{}.{}".format(name, extension)
    figure.savefig(name, bbox_inches='tight')
