 # `ideal_traj_gen` Module

 ## Overview

 The `ideal_traj_gen` module provides functions to generate and plot ideal lift motion trajectories based on jerk, acceleration, speed, and distance constraints. It includes methods to handle single-segment and multi-segment trajectories and plot the results using Matplotlib.

 ## Features

 - Generate single-segment trajectories with specified initial conditions, jerk input, duration, and sampling rate.
 - Generate multi-segment trajectories with specified segments, initial conditions, and sampling rate.
 - Generate ideal lift motion trajectories based on jerk, maximum acceleration, maximum speed, and distance constraints.
 - Plot and save trajectory graphs, including displacement, velocity, acceleration, and jerk.

 ## Installation

 To use the `ideal_traj_gen` module, ensure you have the necessary dependencies installed. You can install them using pip:

 ```bash
 pip install numpy matplotlib scipy
 ```

 ## Usage

 ### Importing the Module

 ```python
 from ideal_traj_gen import traj_gen_ideal_lift_motion, plot_trajectories
 ```

 ### Generating and Plotting Trajectories

 #### Case 1: Distance is Too Short to Reach Max Acceleration

 ```python
 # Define parameters
 jerk_const = 1.0  # constant jerk
 max_acc = 1.5
 max_speed = 4.0
 distance = 4

 # Generate the trajectory
 t_sec, disp, vel, acc, jerk = traj_gen_ideal_lift_motion(jerk_const, max_acc, max_speed, distance)

 # Plot and save the trajectory
 fig = plot_trajectories(t_sec, disp, vel, acc, jerk)
 fig.savefig("ideal_traj_gen_case1.png")
 ```

 #### Case 2: Distance is Long Enough to Reach Max Speed

 ```python
 # Define parameters
 distance = 12.0

 # Generate the trajectory
 t_sec, disp, vel, acc, jerk = traj_gen_ideal_lift_motion(jerk_const, max_acc, max_speed, distance)

 # Plot and save the trajectory
 fig = plot_trajectories(t_sec, disp, vel, acc, jerk)
 fig.savefig("ideal_traj_gen_case2.png")
 ```

 #### Case 3: Distance is Long Enough to Reach Max Speed and Maintain It

 ```python
 # Define parameters
 distance = 24

 # Generate the trajectory
 t_sec, disp, vel, acc, jerk = traj_gen_ideal_lift_motion(jerk_const, max_acc, max_speed, distance)

 # Plot and save the trajectory
 fig = plot_trajectories(t_sec, disp, vel, acc, jerk)
 fig.savefig("ideal_traj_gen_case3.png")
 ```

 ### Function Documentation

 #### `traj_gen_single(init_con, jerk_input, duration, sampling_rate=1000)`

 Generates a single-segment trajectory.

 - **Parameters:**
   - `init_con` (list): Initial conditions [displacement, velocity, acceleration]
   - `jerk_input` (float): Jerk input (m/s^3)
   - `duration` (float): Duration of the trajectory (s)
   - `sampling_rate` (int, optional): Sampling rate of the trajectory (Hz)

 - **Returns:**
   - `t` (numpy.ndarray): Time vector (s)
   - `x` (numpy.ndarray): Displacement vector (m)
   - `dx` (numpy.ndarray): Velocity vector (m/s)
   - `ddx` (numpy.ndarray): Acceleration vector (m/s^2)
   - `dddx` (numpy.ndarray): Jerk vector (m/s^3)

 #### `traj_gen_multi(segments, init_con=[0, 0, 0], sampling_rate=1000)`

 Generates a multi-segment trajectory.

 - **Parameters:**
   - `segments` (list of tuples): Each tuple contains (jerk_input, duration) for a segment
   - `init_con` (list, optional): Initial conditions [displacement, velocity, acceleration]
   - `sampling_rate` (int, optional): Sampling rate of the trajectory (Hz)

 - **Returns:**
   - `t_out` (numpy.ndarray): Time vector (s)
   - `disp_out` (numpy.ndarray): Displacement vector (m)
   - `vel_out` (numpy.ndarray): Velocity vector (m/s)
   - `acc_out` (numpy.ndarray): Acceleration vector (m/s^2)
   - `jerk_out` (numpy.ndarray): Jerk vector (m/s^3)

 #### `plot_trajectories(t, displacement, velocity, acceleration, jerk)`

 Plots the trajectories of displacement, velocity, acceleration, and jerk.

 - **Parameters:**
   - `t` (numpy.ndarray): Time vector (s)
   - `displacement` (numpy.ndarray): Displacement vector (m)
   - `velocity` (numpy.ndarray): Velocity vector (m/s)
   - `acceleration` (numpy.ndarray): Acceleration vector (m/s^2)
   - `jerk` (numpy.ndarray): Jerk vector (m/s^3)

 - **Returns:**
   - `fig` (matplotlib.figure.Figure): The figure object containing the plot

 #### `traj_gen_ideal_lift_motion(jerk, max_acc, max_speed, distance, sampling_rate=1000)`

 Generates an ideal lift motion trajectory.

 - **Parameters:**
   - `jerk` (float): Jerk input (m/s^3)
   - `max_acc` (float): Maximum acceleration (m/s^2)
   - `max_speed` (float): Maximum speed (m/s)
   - `distance` (float): Total distance of the trajectory (m)
   - `sampling_rate` (int, optional): Sampling rate of the trajectory (Hz)

 - **Returns:**
   - `t_out` (numpy.ndarray): Time vector (s)
   - `disp_out` (numpy.ndarray): Displacement vector (m)
   - `vel_out` (numpy.ndarray): Velocity vector (m/s)
   - `acc_out` (numpy.ndarray): Acceleration vector (m/s^2)
   - `jerk_out` (numpy.ndarray): Jerk vector (m/s^3)

 ## License

 This module is released under the MIT License. See the LICENSE file for details.
