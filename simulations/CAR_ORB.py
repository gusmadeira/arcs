#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import math
import numpy as np

# Constantes (mesmas do código C)
PI = np.pi
G = 6.6743e-11
J2 = 16290.543820e-6
J4 = 0.0
J6 = 0.0
RS=60330e3
MS=5.68683765495e26

def wrap_0_2pi(x):
        twopi = 2.0 * np.pi
        while x > twopi or x < 0.0:
            if x > twopi:
                x -= twopi
            elif x < 0.0:
                x += twopi
        return x

def read_float_line(line: str):
    try:
        parts = line.strip().split()
        if len(parts) < 10:
            return None
        vals = [float(p) for p in parts[:10]]
        return vals
    except Exception:
        return None

def main():
    if len(sys.argv) < 2:
        print("Uso: python CAR_ORB.py <arquivo_entrada>")
        sys.exit(1)

    inpath = sys.argv[1]
    outpath = inpath + "_out"

    try:
        fin = open(inpath, "r")
    except OSError as e:
        print(f"Erro ao abrir arquivo de entrada: {e}")
        sys.exit(1)

    try:
        fout = open(outpath, "w")
    except OSError as e:
        print(f"Erro ao abrir arquivo de saída: {e}")
        fin.close()
        sys.exit(1)

    # Leitura da primeira linha válida
    line = fin.readline()
    vals = read_float_line(line)
    # No C, havia uma versão comentada com R, aqui seguimos a de 8 valores
    while vals is None and line:
        line = fin.readline()
        vals = read_float_line(line)

    while line and vals is not None:
        # Desempacotar a linha
        id_, Pos1, Pos2, Pos3, Vel1, Vel2, Vel3, mass, an, en = vals

        # Variáveis conforme o C
        Erre = np.sqrt(Pos1**2.+Pos2**2.)
        z = Pos3
        zdot = Vel3

        # L: ângulo no plano XY (usaremos atan2 para simplificar e depois normalizar)
        if Pos1 < 0.0:
            L = math.atan(Pos2/Pos1) + 1.0*PI
        elif math.atan(Pos2/Pos1) < 0.0:
            L = math.atan(Pos2/Pos1) + 2.0*PI
        elif (Pos1 == 0.0) and (Pos2 > 0.0):
            L = (PI/2.0)
        elif (Pos1 == 0.0) and (Pos2 < 0.0):
            L = ((3.0*PI)/2.0)
        else:
            L = math.atan(Pos2/Pos1)
        L = wrap_0_2pi(L)

        # projs radiais e angulares no plano XY
        Erredot = Vel1 * np.cos(L) + Vel2 * np.sin(L)
        Ldot = (-Vel1 * np.sin(L) + Vel2 * np.cos(L)) / Erre

        # Inicializações
        itera = 0
        SMA = Erre
        ECC = 0.0
        INC = 0.0

        Rc = 0.0
        Lc = 0.0
        zc = 0.0
        Rcdot = 0.0
        Lcdot = 0.0
        zcdot = 0.0

        MU = G * (MS + mass)
        epson = 1.0e-15

        # Loop iterativo (máx 1000 iterações)
        while True:
            ai = SMA

            # frequências médias e épicíclicas (mantendo fidelidade às potências)
            # NMEDI
            SMA2 = SMA * SMA
            SMA3 = SMA2 * SMA
            RS_SMA = RS / SMA
            RS_SMA2 = RS_SMA * RS_SMA
            RS_SMA4 = RS_SMA2 * RS_SMA2
            RS_SMA6 = RS_SMA4 * RS_SMA2

            base_n = np.sqrt(MU / SMA3)
            NMEDI = base_n * (
                1
                + (3.0/4.0) * J2 * RS_SMA2
                - (15.0/16.0) * J4 * RS_SMA4
                + (35.0/32.0) * J6 * RS_SMA6
                - (9.0/32.0) * (J2**2) * RS_SMA4
                + (45.0/64.0) * J2 * J4 * RS_SMA6
                + (27.0/128.0) * (J2**3) * RS_SMA6
                + 3.0 * J2 * (ECC**2) * RS_SMA2
                - 12.0 * RS_SMA2 * J2 * (INC**2)
            )

            # kappa (epicíclica horizontal)
            kappa = base_n * (
                1
                - (3.0/4.0) * J2 * RS_SMA2
                + (45.0/16.0) * J4 * RS_SMA4
                - (175.0/32.0) * J6 * RS_SMA6
                - (9.0/32.0) * (J2**2) * RS_SMA4
                + (135.0/64.0) * J2 * J4 * RS_SMA6
                - (27.0/128.0) * (J2**3) * RS_SMA6
                - 9.0 * RS_SMA2 * J2 * (INC**2)
            )

            # nu (epicíclica vertical)
            nu = base_n * (
                1
                + (9.0/4.0) * J2 * RS_SMA2
                - (75.0/16.0) * J4 * RS_SMA4
                + (245.0/32.0) * J6 * RS_SMA6
                - (81.0/32.0) * (J2**2) * RS_SMA4
                + (675.0/64.0) * J2 * J4 * RS_SMA6
                + (729.0/128.0) * (J2**3) * RS_SMA6
                + 6.0 * RS_SMA2 * J2 * (ECC**2)
                - (51.0/4.0) * RS_SMA2 * J2 * (INC**2)
            )

            eta2 = (MU / SMA3) * (
                1
                - 2.0 * J2 * RS_SMA2
                + (75.0/8.0) * J4 * RS_SMA4
                - (175.0/8.0) * J6 * RS_SMA6
            )

            chi2 = (MU / SMA3) * (
                1
                + (15.0/2.0) * J2 * RS_SMA2
                - (175.0/8.0) * RS_SMA4 * J4
                + (735.0/16.0) * RS_SMA6 * J6
            )

            # alpha's
            # No C: alphaa[1] = (2*nu + kappa)/3, alphaa[2] = 2*nu - kappa
            alpha1 = (2.0 * nu + kappa) / 3.0
            alpha2_ = 2.0 * nu - kappa
            alpha2_prod = alpha1 * alpha2_
            # alpha = sqrt(alpha2_prod) (não usado diretamente depois)
            # alpha = math.sqrt(alpha2_prod) if alpha2_prod > 0 else 0.0

            # parâmetros intermediários
            # sem, exc, inclina
            denom = 1.0 - (Ldot - Lcdot - NMEDI) / (2.0 * NMEDI)
            sem = (Erre - Rc) / denom

            if sem == 0.0 or kappa == 0.0:
                exc = 0.0
            else:
                exc = math.hypot(
                    (Ldot - Lcdot - NMEDI) / (2.0 * NMEDI),
                    (Erredot - Rcdot) / (sem * kappa)
                )

            inclina = np.sqrt( ((z - zc)/sem)**2 + ((zdot-zcdot)/(sem*nu))**2 )

            # Ajuste iterativo idêntico ao C
            while (L < 0.0) or (inclina > PI):
                if inclina > PI:
                    inclina = inclina - PI
                elif inclina < 0.0:
                    inclina = inclina + PI

            # lambda
            lambda_ = L - Lc - 2.0 * (NMEDI / kappa) * ((Erredot - Rcdot) / (sem * kappa))
            lambda_ = wrap_0_2pi(lambda_)

            denomM = sem*kappa*(1-((Erre-Rc)/sem))
            M = math.atan((Erredot - Rcdot) / denomM)

            if (-Erre+Rc+sem)/(sem*exc)<0.0:
                M=M+1.0*PI
            elif M<0.:
                M=M+2.0*PI
            elif ((-Erre+Rc+sem)/(sem*exc))==0.0 and (Erredot-Rcdot)/(sem*kappa*kappa)>0.0:
                M=PI/2.
            elif ((-Erre+Rc+sem)/(sem*exc))==0.0 and (Erredot-Rcdot)/(sem*kappa*kappa)<0.0:
                M=3./2.*PI
            M = wrap_0_2pi(M)

            # lom (argumento relacionado ao nó/longitude do nó)
            denom_lom = (zdot - zcdot)
            numer_lom = nu * (z - zc)
            if denom_lom == 0.0 and numer_lom == 0.0:
                lom = 0.0
            else:
                if sem != 0.0 and inclina != 0.0 and nu != 0.0:
                    lom = math.atan2(numer_lom, denom_lom)
                else:
                    lom = 0.0

            # Ajustes de sinal conforme o C
            if sem != 0.0 and inclina != 0.0 and nu != 0.0:
                if (zdot - zcdot) / (sem * inclina * nu) < 0.0:
                    lom += PI
            if lom < 0.0:
                lom += 2.0 * PI
            lom = wrap_0_2pi(lom)

            # Atualiza variáveis de iteração
            SMA = sem
            ECC = exc
            INC = inclina

            # Termos de correção Rc, Lc, zc e derivadas (mantendo fórmulas do C)
            # evitar divisões por kappa=0, alpha2_prod=0, nu=0
            kappa2 = kappa * kappa
            NMEDI2 = NMEDI * NMEDI

            cos2M = np.cos(2.0 * M)
            sin2M = np.sin(2.0 * M)
            cos2lom = np.cos(2.0 * lom)
            sin2lom = np.sin(2.0 * lom)

            term1_Rc = SMA * (ECC**2) * ((3.0/2.0) * (eta2 / kappa2) - 1.0 - 0.5 * (eta2 / kappa2) * cos2M)
            term2_Rc = 0.0
            if alpha2_prod != 0.0:
                term2_Rc = SMA * (INC**2) * ((3.0/4.0) * (chi2 / kappa2) - 1.0 + (chi2 / (4.0 * alpha2_prod)) * cos2lom)
            Rc = term1_Rc + term2_Rc

            term1_Lc = 0.0
            if kappa != 0.0:
                term1_Lc = (NMEDI / kappa) * (ECC**2) * ((3.0/4.0) + eta2 / (2.0 * kappa2)) * sin2M
            term2_Lc = 0.0
            if alpha2_prod != 0.0 and nu != 0.0:
                term2_Lc = -(INC**2) * (chi2 / (4.0 * alpha2_prod)) * (NMEDI / nu) * sin2lom
            Lc = term1_Lc + term2_Lc

            zc = 0.0
            if sem != 0.0 and ECC != 0.0:
                t1 = 0.0
                if kappa != 0.0 and alpha1 != 0.0:
                    t1 = (chi2 / (2.0 * kappa * alpha1)) * np.sin(lom + M)
                t2 = 0.0
                if kappa != 0.0 and alpha2_ != 0.0:
                    t2 = (3.0/2.0) * (chi2 / (kappa * alpha2_)) * np.sin(lom - M)
                zc = SMA * INC * ECC * (t1 - t2)

            Rcdot = 0.0
            if kappa != 0.0:
                Rcdot = SMA * (ECC**2) * (eta2 / kappa) * np.sin(2.0 * M) - SMA * (INC**2) * (chi2 / (2.0 * alpha2_prod if alpha2_prod != 0.0 else 1.0)) * nu * sin2lom

            Lcdot = (
                NMEDI * (ECC**2) * (
                    (7.0/2.0)
                    - 3.0 * (eta2 / kappa2 if kappa2 != 0.0 else 0.0)
                    - (kappa2 / (2.0 * NMEDI2) if NMEDI2 != 0.0 else 0.0)
                    + ((3.0/2.0) + (eta2 / kappa2 if kappa2 != 0.0 else 0.0)) * cos2M
                )
                + (INC**2) * NMEDI * (
                    2.0
                    - (kappa2 / (2.0 * NMEDI2) if NMEDI2 != 0.0 else 0.0)
                    - (3.0/2.0) * (chi2 / kappa2 if kappa2 != 0.0 else 0.0)
                    - (chi2 / (2.0 * alpha2_prod) if alpha2_prod != 0.0 else 0.0) * cos2lom
                )
            )

            t1 = (chi2 * (kappa + nu)) / (2.0 * kappa * alpha1) * np.cos(lom + M)
            t2 = (3.0/2.0) * (chi2 * (kappa - nu)) / (kappa * alpha2_) * np.cos(lom - M)
            zcdot = SMA * INC * ECC * (t1 + t2)

            modulo = abs(SMA - ai)
            itera += 1

            if (modulo < epson or itera >= 1000):
                # Cálculos finais e escrita
                longnodo = wrap_0_2pi(lambda_ - lom)
                wpi = wrap_0_2pi(lambda_ - M)

                # Formatação semelhante ao printf do C
                # "% 10.5lf % 18.15lf % 18.15lf % 18.15lf  %18.15lf  %18.15lf  %18.15lf  %10.5lf\n"
                fout.write(f"{id_: 10.5f}"
                           f" {SMA: 18.15f} {ECC: 18.15f} {INC: 18.15f}"
                           f" {wpi * 180.0 / PI: 18.15f} {longnodo * 180.0 / PI: 18.15f} {lambda_ * 180.0 / PI: 18.15f}"
                           f" {mass: 10.5f}\n")
                break

        # próxima linha
        line = fin.readline()
        vals = read_float_line(line)

    fin.close()
    fout.close()
    print(f"Feito. Saída em: {outpath}")

if __name__ == "__main__":
    main()
