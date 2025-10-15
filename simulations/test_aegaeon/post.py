#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Mplanet  = 5.68683765495e26   
Rplanet  = 60330e3

def wrap_0_360(x):
    x = np.asarray(x) % 360.0
    return x.item() if np.ndim(x) == 0 else x

cols = ["id", "a", "e", "inc", "wpi", "Omega", "lambda", "mass"]

mimas = pd.read_csv("mimas.aei_out", delim_whitespace=True, names=cols)
aegaeon = pd.read_csv("aegaeon.aei_out", delim_whitespace=True, names=cols)

phi = 7.0 * mimas["lambda"].values - 6.0 * aegaeon["lambda"].values - mimas["wpi"].values
phi = wrap_0_360(phi)

tempo_anos = mimas["id"].values / (365.25 * 24 * 60 * 60)
tempo_anos_a = aegaeon["id"].values / (365.25 * 24 * 60 * 60)

plt.figure(figsize=(8, 4))
plt.plot(tempo_anos, phi, '.', markersize=2)
plt.xlabel("Tempo (anos)")
plt.ylim(0,360)
plt.yticks(np.linspace(0,360,7))
plt.ylabel("Ângulo ressonante φ (graus)")
plt.title("Evolução do ângulo ressonante 7:6")
plt.tight_layout()
plt.savefig("phi.png", dpi=300)
#plt.show()

fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
fig.subplots_adjust(hspace=0.1)
axs[0].plot(tempo_anos, mimas["a"]/Rplanet)
axs[0].set_ylabel("a ($R_S$)")
axs[1].plot(tempo_anos, mimas["e"])
axs[1].set_ylabel("e")
axs[2].plot(tempo_anos, wrap_0_360(np.degrees(mimas["inc"])))
axs[2].set_ylabel("i (°)")
plt.tight_layout()
plt.savefig("ae_mimas.png", dpi=300)
#plt.show()

fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
fig.subplots_adjust(hspace=0.1)
axs[0].plot(tempo_anos_a, aegaeon["a"]/Rplanet)
axs[0].set_ylabel("a ($R_S$)")
axs[1].plot(tempo_anos_a, aegaeon["e"])
axs[1].set_ylabel("e")
axs[2].plot(tempo_anos_a, wrap_0_360(np.degrees(aegaeon["inc"])))
axs[2].set_ylabel("i (°)")
plt.tight_layout()
plt.savefig("ae_aegaeon.png", dpi=300)
#plt.show()