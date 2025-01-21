#Developed by Gustavo Gonçalves dos Santos
#g.gustavo.santos@gmail.com

import os
from math import radians as rad
from cmath import exp as exp
from numpy import ndarray, zeros, pi, conj
from scipy.sparse import csr_matrix


class system_FaultData:
    def __init__(
        self, systemName,
        system: str = "",
        nbus: int = 0,
        bus: dict = dict(),
        nger: int = 0,
        ger: dict = dict(),
        nline: int = 0,
        branch: dict = dict(),
        ntrafo: int = 0,
        trafo: dict = dict(),
        y1: ndarray = [],
        y2: ndarray = [],
        y0: ndarray = [],
    ):

        # Inicialization
        self.system = os.path.splitext(systemName)[0]
        self.nbus = nbus
        self.bus = bus
        bus['num'] = list()
        bus['name'] = list()
        bus['voltage_magnitude'] = list()
        bus['voltage_angle'] = list()
        bus['voltage_base'] = list()
        self.nger = nger
        self.ger = ger
        ger['bus'] = list()
        ger['name'] = list()
        ger['R1'] = list()
        ger['X1'] = list()
        ger['R2'] = list()
        ger['X2'] = list()
        ger['R0'] = list()
        ger['X0'] = list()
        self.nline = nline
        self.branch = branch
        branch['from'] = list()
        branch['to'] = list()
        branch['R1'] = list()
        branch['X1'] = list()
        branch['R0'] = list()
        branch['X0'] = list()
        self.ntrafo = ntrafo
        self.trafo = trafo
        trafo['from'] = list()
        trafo['to'] = list()
        trafo['R1'] = list()
        trafo['X1'] = list()
        trafo['R0'] = list()
        trafo['X0'] = list()
        trafo['Vm_tap'] = list()
        trafo['Va_tap'] = list()
        trafo['conx'] = list()
        self.y1 = y1
        self.y2 = y2
        self.y0 = y0

    def add_values(self):

        Sbase = 100

        folder = os.path.join(os.path.dirname(__file__), 'Systems')
        fname = os.path.join(folder, self.system + '.txt')

        f = open(f'{fname}', 'r', encoding='latin-1')

        lines = f.readlines()
        f.close()

        i = 0
        while lines[i].strip() != 'END':
            i += 1
            if 'BUS DATA' in lines[i].strip():
                i += 1
                while 'END BUS' != lines[i][0:7]:
                    if lines[i][0] != '(':
                        self.nbus += 1
                        self.bus['num'].append(int(lines[i][0:8].strip()))
                        self.bus['name'].append(lines[i][8:20].strip())
                        self.bus['voltage_magnitude'].append(
                            float(lines[i][20:28].strip()))
                        self.bus['voltage_angle'].append(
                            float(lines[i][28:36].strip())*pi/180)
                        self.bus['voltage_base'].append(
                            float(lines[i][36:44].strip()))
                    i += 1

            if 'GER DATA' in lines[i].strip():
                i += 1
                while 'END GER DATA' != lines[i][0:12]:
                    if lines[i][0] != '(':
                        self.nger += 1
                        self.ger['bus'].append(int(lines[i][0:8].strip()))
                        self.ger['R1'].append(float(lines[i][8:16].strip()))
                        self.ger['X1'].append(float(lines[i][16:24].strip()))
                        self.ger['R2'].append(float(lines[i][24:32].strip()))
                        self.ger['X2'].append(float(lines[i][32:40].strip()))
                        self.ger['R0'].append(float(lines[i][40:48].strip()))
                        self.ger['X0'].append(float(lines[i][48:56].strip()))
                    i += 1

            if 'BRANCH DATA' in lines[i].strip():
                i += 1
                while 'END BRANCH DATA' != lines[i][0:15]:
                    if lines[i][0] != '(':
                        self.nline += 1
                        self.branch['from'].append(int(lines[i][0:8].strip()))
                        self.branch['to'].append(int(lines[i][8:16].strip()))
                        self.branch['R1'].append(
                            float(lines[i][16:24].strip()))
                        self.branch['X1'].append(
                            float(lines[i][24:32].strip()))
                        self.branch['R0'].append(
                            float(lines[i][32:40].strip()))
                        self.branch['X0'].append(
                            float(lines[i][40:48].strip()))
                    i += 1

            if 'TRAFO DATA' in lines[i].strip():
                i += 1
                while 'END TRAFO DATA' != lines[i][0:14]:
                    if lines[i][0] != '(':
                        self.ntrafo += 1
                        self.trafo['from'].append(int(lines[i][0:8].strip()))
                        self.trafo['to'].append(int(lines[i][8:16].strip()))
                        self.trafo['R1'].append(float(lines[i][16:24].strip()))
                        self.trafo['X1'].append(float(lines[i][24:32].strip()))
                        self.trafo['R0'].append(float(lines[i][32:40].strip()))
                        self.trafo['X0'].append(float(lines[i][40:48].strip()))
                        self.trafo['Vm_tap'].append(
                            float(lines[i][48:56].strip()))
                        self.trafo['Va_tap'].append(
                            float(lines[i][56:64].strip()))
                        self.trafo['conx'].append(lines[i][64:72].strip())
                    i += 1

        self.y1 = zeros(shape=[self.nbus, self.nbus], dtype='complex')
        self.y2 = zeros(shape=[self.nbus, self.nbus], dtype='complex')
        self.y0 = zeros(shape=[self.nbus, self.nbus], dtype='complex')

        for i in range(self.nger):
            j = 0
            while True:
                if self.ger['bus'][i] == self.bus['num'][j]:
                    break
                j += 1
            self.y1[j, j] = 1/(complex(
                real=self.ger['R1'][i], imag=self.ger['X1'][i]))
            self.y2[j, j] = 1/(complex(
                real=self.ger['R2'][i], imag=self.ger['X2'][i]))
            self.y0[j, j] = 1/(complex(
                real=self.ger['R0'][i], imag=self.ger['X0'][i]))

        for i in range(self.nline):
            j = 0
            while True:
                if self.branch['from'][i] == self.bus['num'][j]:
                    from_bus = j
                    break
                j += 1
            j = 0
            while True:
                if self.branch['to'][i] == self.bus['num'][j]:
                    to_bus = j
                    break
                j += 1
            self.y1[from_bus, from_bus] += 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y1[from_bus, to_bus] -= 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y1[to_bus, from_bus] -= 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y1[to_bus, to_bus] += 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y2[from_bus, from_bus] += 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y2[from_bus, to_bus] -= 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y2[to_bus, from_bus] -= 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y2[to_bus, to_bus] += 1 / \
                (complex(real=self.branch['R1'][i], imag=self.branch['X1'][i]))
            self.y0[from_bus, from_bus] += 1 / \
                (complex(real=self.branch['R0'][i], imag=self.branch['X0'][i]))
            self.y0[from_bus, to_bus] -= 1 / \
                (complex(real=self.branch['R0'][i], imag=self.branch['X0'][i]))
            self.y0[to_bus, from_bus] -= 1 / \
                (complex(real=self.branch['R0'][i], imag=self.branch['X0'][i]))
            self.y0[to_bus, to_bus] += 1 / \
                (complex(real=self.branch['R0'][i], imag=self.branch['X0'][i]))

        for i in range(self.ntrafo):
            j = 0
            while True:
                if self.trafo['from'][i] == self.bus['num'][j]:
                    from_bus = j
                    break
                j += 1
            j = 0
            while True:
                if self.trafo['to'][i] == self.bus['num'][j]:
                    to_bus = j
                    break
                j += 1
            if self.trafo['conx'][i] == 'Yt-Yt':
                conx = 1
            elif self.trafo['conx'][i] == 'Y-Y':
                conx = 2
            elif self.trafo['conx'][i] == 'Yt-D':
                conx = 3
            elif self.trafo['conx'][i] == 'Y-D':
                conx = 4
            elif self.trafo['conx'][i] == 'D-Y':
                conx = 5
            elif self.trafo['conx'][i] == 'D-Yt':
                conx = 6
            elif self.trafo['conx'][i] == 'D-D':
                conx = 7
            elif self.trafo['conx'][i] == 'Yt-Y':
                conx = 8
            elif self.trafo['conx'][i] == 'Y-Yt':
                conx = 9
            if conx == 1 or conx == 2 or conx == 7 or conx == 8 or conx == 9:
                # t = self.trafo['Vm_tap'][i]
                t = 1
            elif conx == 3 or conx == 4:
                # t = self.trafo['Vm_tap'][i] * exp(-1j*30*pi/180)
                t = exp(1j*30*pi/180)
            elif conx == 5 or conx == 6:
                # t = self.trafo['Vm_tap'][i] * exp(1j*30*pi/180)
                t = exp(-1j*30*pi/180)

            self.y1[from_bus, from_bus] += (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))
            self.y1[from_bus, to_bus] -= (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))*t
            self.y1[to_bus, from_bus] -= (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))*conj(t)
            self.y1[to_bus, to_bus] += (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))*pow(abs(t), 2)

            self.y2[from_bus, from_bus] += (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))
            self.y2[from_bus, to_bus] -= (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))*conj(t)
            self.y2[to_bus, from_bus] -= (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))*t
            self.y2[to_bus, to_bus] += (
                1 / (complex(real=self.trafo['R1'][i], imag=self.trafo['X1'][i])))*pow(abs(t), 2)

            if conx == 1:
                self.y0[from_bus, from_bus] += 1 / \
                    (complex(real=self.trafo['R0']
                     [i], imag=self.trafo['X0'][i]))
                self.y0[from_bus, to_bus] -= 1 / \
                    (complex(real=self.trafo['R0']
                     [i], imag=self.trafo['X0'][i]))
                self.y0[to_bus, from_bus] -= 1 / \
                    (complex(real=self.trafo['R0']
                     [i], imag=self.trafo['X0'][i]))
                self.y0[to_bus, to_bus] += 1 / \
                    (complex(real=self.trafo['R0']
                     [i], imag=self.trafo['X0'][i]))

            elif conx == 2 or conx == 4 or conx == 5 or conx == 7 or conx == 8 or conx == 9:
                self.y0[from_bus, to_bus] = 0
                self.y0[to_bus, from_bus] = 0
            elif conx == 3:
                self.y0[from_bus, from_bus] += 1 / \
                    (complex(real=self.trafo['R0']
                     [i], imag=self.trafo['X0'][i]))
                self.y0[from_bus, to_bus] = 0
                self.y0[to_bus, from_bus] = 0
            elif conx == 6:
                self.y0[from_bus, to_bus] = 0
                self.y0[to_bus, from_bus] = 0
                self.y0[to_bus, to_bus] += 1 / \
                    (complex(real=self.trafo['R0']
                     [i], imag=self.trafo['X0'][i]))
