import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def traj_gen_single(init_con, jerk_input, duration, sampling_rate=1000):
    """
    Generate a trajectory given initial conditions, jerk input, duration, and sampling rate.

    Parameters
    ----------
    init_con : list
        Initial conditions [displacement, velocity, acceleration]
    jerk_input : float
        Jerk input (m/s^3)
    duration : float
        Duration of the trajectory (s)
    sampling_rate : int, optional
        Sampling rate of the trajectory (Hz)

    Returns
    -------
    t : numpy.ndarray
        Time vector (s)
    x : numpy.ndarray
        Displacement vector (m)
    dx : numpy.ndarray
        Velocity vector (m/s)
    ddx : numpy.ndarray
        Acceleration vector (m/s^2)
    dddx : numpy.ndarray
        Jerk vector (m/s^3)
    """
    assert len(
        init_con) == 3, "Initial conditions must be a list of 3 elements"
    assert duration > 0, "Duration must be greater than 0"
    assert sampling_rate > 0, "Sampling rate must be greater than 0"

    def state_space_model(t, y, jerk):
        return [y[1], y[2], jerk]

    t = np.linspace(0, duration, int(sampling_rate * duration))
    sol = solve_ivp(lambda t, y: state_space_model(t, y, jerk_input),
                    [0, duration],
                    init_con,
                    t_eval=t)

    x = sol.y[0]
    dx = sol.y[1]
    ddx = sol.y[2]
    dddx = np.full_like(t, jerk_input)

    return t, x, dx, ddx, dddx


def traj_gen_multi(segments, init_con=[0, 0, 0], sampling_rate=1000):
    """
    Generate a multi-segment trajectory.

    Parameters
    ----------
    segments : list of tuples
        Each tuple contains (jerk_input, duration) for a segment
    init_con : list, optional
        Initial conditions [displacement, velocity, acceleration]
    sampling_rate : int, optional
        Sampling rate of the trajectory (Hz)

    Returns
    -------
    t_out : numpy.ndarray
        Time vector (s)
    disp_out : numpy.ndarray
        Displacement vector (m)
    vel_out : numpy.ndarray
        Velocity vector (m/s)
    acc_out : numpy.ndarray
        Acceleration vector (m/s^2)
    jerk_out : numpy.ndarray
        Jerk vector (m/s^3)
    """
    t_out = np.array([])
    disp_out = np.array([])
    vel_out = np.array([])
    acc_out = np.array([])
    jerk_out = np.array([])

    for jerk_in, duration in segments:
        if duration == 0:
            continue

        t, disp_temp, vel_temp, acc_temp, jerk_temp = traj_gen_single(
            init_con, jerk_in, duration, sampling_rate)

        if t_out.size > 0:
            t += t_out[-1]

        t_out = np.concatenate((t_out, t))
        disp_out = np.concatenate((disp_out, disp_temp))
        vel_out = np.concatenate((vel_out, vel_temp))
        acc_out = np.concatenate((acc_out, acc_temp))
        jerk_out = np.concatenate((jerk_out, jerk_temp))

        init_con = [disp_temp[-1], vel_temp[-1], acc_temp[-1]]

    return t_out, disp_out, vel_out, acc_out, jerk_out


def plot_trajectories(t, displacement, velocity, acceleration, jerk):
    """
    Plot the trajectories of displacement, velocity, acceleration, and jerk.

    Parameters
    ----------
    t : numpy.ndarray
        Time vector (s)
    displacement : numpy.ndarray
        Displacement vector (m)
    velocity : numpy.ndarray
        Velocity vector (m/s)
    acceleration : numpy.ndarray
        Acceleration vector (m/s^2)
    jerk : numpy.ndarray
        Jerk vector (m/s^3)

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot
    """

    fig, axs = plt.subplots(4, 1, figsize=(6, 6), sharex=True)

    axs[0].plot(t, jerk, '-', label='Jerk')
    axs[0].set_ylabel('Jerk (m/s^3)')
    axs[0].legend(loc='upper right')
    axs[0].grid(True)

    axs[1].plot(t, acceleration, '-', label='Acceleration', color='orange')
    axs[1].set_ylabel('Acceleration (m/s^2)')
    axs[1].legend(loc='upper right')
    axs[1].grid(True)

    axs[2].plot(t, velocity, '-', label='Velocity', color='green')
    axs[2].set_ylabel('Velocity (m/s)')
    axs[2].legend(loc='upper right')
    axs[2].grid(True)

    axs[3].plot(t, displacement, '-', label='Displacement', color='red')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Displacement (m)')
    axs[3].legend(loc='upper right')
    axs[3].grid(True)

    # fig.suptitle('Ideal Lift Kinematics', fontsize=16)
    plt.tight_layout()

    plt.show()
    return fig


def traj_gen_ideal_lift_motion(jerk,
                               max_acc,
                               max_speed,
                               distance,
                               sampling_rate=1000):
    """
    Generate a connected trajectory based on jerk, maximum acceleration, maximum speed, and distance constraints.

    Parameters
    ----------
    jerk : float
        Jerk input (m/s^3)
    max_acc : float
        Maximum acceleration (m/s^2)
    max_speed : float
        Maximum speed (m/s)
    distance : float
        Total distance of the trajectory (m)
    sampling_rate : int, optional
        Sampling rate of the trajectory (Hz)

    Returns
    -------
    t_out : numpy.ndarray
        Time vector (s)
    disp_out : numpy.ndarray
        Displacement vector (m)
    vel_out : numpy.ndarray
        Velocity vector (m/s)
    acc_out : numpy.ndarray
        Acceleration vector (m/s^2)
    jerk_out : numpy.ndarray
        Jerk vector (m/s^3)
    """
    init_con = [0, 0, 0]

    t1 = max_acc / jerk  # Time to reach max acceleration
    v1 = 0.5 * jerk * t1**2  # Velocity reached at max acceleration
    d1 = v1 * t1 * 4  # Minimum total distance to reach max acceleration

    assert max_speed > 2 * v1, "Max acceleration must be reached before max speed."

    t2 = (max_speed -
          2 * v1) / max_acc  # Time to reach max speed from max acceleration
    d2 = max_speed * (2 * t1 + t2)  # Minimum total distance to reach max speed

    if distance < d1:
        print("case 1: distance is too short to reach max acceleration")
        t_j = (distance / 2 / jerk)**(1 / 3)
        segments = [
            (jerk, t_j),
            (-jerk, 2 * t_j),
            (jerk, t_j),
        ]
    elif distance < d2:
        print("case 2: distance is too short to reach max speed")
        const_a = 1
        const_b = (2 * v1 / max_acc + 2 * t1)
        const_c = (4 * v1 * t1 - distance) / max_acc
        discriminant = const_b**2 - 4 * const_a * const_c
        assert discriminant >= 0, "discriminant must be positive"
        t_acc = (-const_b + np.sqrt(discriminant)) / (2 * const_a)

        segments = [
            (jerk, t1),
            (0.0, t_acc),
            (-jerk, 2 * t1),
            (0.0, t_acc),
            (jerk, t1),
        ]
    else:
        print("case 3: distance is long enough to reach max speed")
        t_vel = (distance - d2) / max_speed
        segments = [
            (jerk, t1),
            (0.0, t2),
            (-jerk, t1),
            (0.0, t_vel),
            (-jerk, t1),
            (0.0, t2),
            (jerk, t1),
        ]

    return traj_gen_multi(segments, init_con, sampling_rate)
