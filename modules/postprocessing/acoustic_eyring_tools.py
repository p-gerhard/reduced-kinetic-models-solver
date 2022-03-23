import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


import logging
import json

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


logging.basicConfig(
    format="[%(asctime)s] - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)


SOUND_SPEED = 1.0

# SOUND_SPEED = 343.


class RoomSimulation:
    def __init__(
        self, simu_file, room_volume, room_area, room_alpha, src_toff, t_max=1.0, e0=1, is_2d=False
    ):
        if SOUND_SPEED != 343:
            logger.warning(
                "{:<30}: {:>6.8f} m/s".format(
                    "Sound Speed is NOT 343 m/s but equal to:", SOUND_SPEED
                )
            )

        self.simu_file = simu_file
        self.room_volume = room_volume
        self.room_area = room_area
        self.room_alpha = room_alpha
        self.src_toff = src_toff
        self.t_max = t_max

        self.n_points = 3000
        self.e0 = e0
        self.is_2d = is_2d

        if self.is_2d:
            self.mean_free_path = self.get_2d_mean_free_path()
        else:
            self.mean_free_path = self.get_3d_mean_free_path()

        # Get caracteristic times T30 T60
        self.t_30, self.t_60 = self.get_eyring_t_char()

        self.t_list, self.eyring_energy_decay = self.get_eyring_energy_decay()

        self.write_eyring_data()

    def get_3d_mean_free_path(self):
        return (4.0 * self.room_volume) / self.room_area

    def get_2d_mean_free_path(self):
        return (np.pi * self.room_volume) / self.room_area

    def get_eyring_t_char(self):

        rt_type = 30
        t_30 = (
            -(rt_type / 10.0)
            * np.log(10.0)
            * self.mean_free_path
            / (SOUND_SPEED * np.log(1.0 - self.room_alpha))
        )

        rt_type = 60
        t_60 = (
            -(rt_type / 10.0)
            * np.log(10.0)
            * self.mean_free_path
            / (SOUND_SPEED * np.log(1.0 - self.room_alpha))
        )
        return t_30, t_60

    # Compute total energy decay using ergodic theory using
    # E(t) = E_0 * exp( c * ln(1 - alpha) * t / <l>)
    # or
    # E(t) = E_0 * (1 - alpha) ^ (ct / <l>)

    def get_eyring_energy_decay(self):

        n_points = 3000
        t_list = np.linspace(0.0, self.t_max + 0.01, n_points)
        w_list = np.zeros(n_points)

        for idx, t in enumerate(t_list):
            if t < self.src_toff:
                w_list[idx] = self.e0
            else:
                w_list[idx] = self.e0 * np.exp(
                    SOUND_SPEED
                    * np.log(1.0 - self.room_alpha)
                    * (t - self.src_toff)
                    / self.mean_free_path
                )

        return (
            t_list,
            w_list,
        )

    def write_eyring_data(self, do_plot=True):
        print("Theoritical Eyring T30: {:<6.8} s".format(self.t_30))
        print("Theoritical Eyring T60: {:<6.8} s".format(self.t_60))

        with open("eyring_t30.dat", "w") as f:
            f.write("{:6.12f}".format(self.t_30))

        with open("eyring_t60.dat", "w") as f:
            f.write("{:6.12f}".format(self.t_60))

        eyring_data = np.column_stack((self.t_list, self.eyring_energy_decay))

        np.savetxt("eyring_energy_decay.csv", eyring_data, delimiter=",")

        num_data = np.genfromtxt(self.simu_file)

        if do_plot:
            plt.plot(self.t_list, self.eyring_energy_decay)
            plt.plot(num_data[:,0], num_data[:,1] / np.max(num_data[:,1]))
            plt.xlabel("time")
            plt.ylabel("E_tot")
            plt.title("Energy decay (Eyring)")
            plt.show()
