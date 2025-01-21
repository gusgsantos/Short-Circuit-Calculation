#Developed by Gustavo Gonçalves dos Santos
#g.gustavo.santos@gmail.com


# Simplified short circuit calculation
# Disregard of:
# - Transformer tap;
# - Generator grounding impedance;
# - Mutual impedance of lines;
# - branch susceptance;
# - Load current.


from numpy import array, zeros, pi, conj, exp, linalg, round, angle, sqrt
import os
from os.path import dirname, exists, realpath
from copy import copy

from system_FaultData import system_FaultData


systemName = 'ieee14bus_fault.txt'
faultBus = 6  # Bus of fault
Zf = 2  # Impedance of fault
ShortCircuit_type = 4

# Type of short circuit
# 1 - ABC
# 2 - AT
# 3 - AB
# 4 - ABT

d = system_FaultData(systemName)
d.add_values()

Va = array(d.bus['voltage_angle'])
Vm = array(d.bus['voltage_magnitude'])

Vpre_fault = Vm * exp(1j * Va)
Vpos_fault = list()

Z0 = linalg.inv(d.y0)
Z1 = linalg.inv(d.y1)
Z2 = linalg.inv(d.y2)

a = exp(1j*120*pi/180)
a2 = exp(1j*240*pi/180)

Ia0 = zeros((d.nbus, 1), dtype='complex')
Ia1 = zeros((d.nbus, 1), dtype='complex')
Ia2 = zeros((d.nbus, 1), dtype='complex')

Ia = zeros((d.nbus, 1), dtype='complex')
Ib = zeros((d.nbus, 1), dtype='complex')
Ic = zeros((d.nbus, 1), dtype='complex')

Va0 = zeros((d.nbus, 1), dtype='complex')
Va1 = zeros((d.nbus, 1), dtype='complex')
Va2 = zeros((d.nbus, 1), dtype='complex')

matrix_A = zeros(shape=[3, 3], dtype='complex')
matrix_A[0, 0] = 1
matrix_A[0, 1] = 1
matrix_A[0, 2] = 1
matrix_A[1, 0] = 1
matrix_A[1, 1] = a2
matrix_A[1, 2] = a
matrix_A[2, 0] = 1
matrix_A[2, 1] = a
matrix_A[2, 2] = a2

i = 0
while True:
    if d.bus['num'][i] == faultBus:
        id_fault = i
        break
    i += 1


if ShortCircuit_type == 1:  # ABC
    for i in range(d.nbus):
        Ia1[i] = Vpre_fault[i] / (Z1[i, i] + Zf)
        I012 = array([Ia0[i], Ia1[i], Ia2[i]])
        Ia[i], Ib[i], Ic[i] = matrix_A @ I012
    for i in range(d.nbus):
        if i == id_fault:
            Va1[i] = Vpre_fault[i] * Zf / (Z1[i, i] + Zf)
        else:
            Va1[i] = Vpre_fault[i] - \
                (Z1[i, id_fault] / (Zf + Z1[id_fault, id_fault])) * \
                Vpre_fault[id_fault]
        V012 = array([Va0[i], Va1[i], Va2[i]])
        Vpos_fault.append(matrix_A @ V012)

elif ShortCircuit_type == 2:
    for i in range(d.nbus):
        Ia1[i] = Vpre_fault[i]/(Z0[i, i]+Z1[i, i]+Z2[i, i]+3*Zf)
        Ia2[i] = Ia1[i]
        Ia0[i] = Ia1[i]
        I012 = array([Ia0[i], Ia1[i], Ia2[i]])
        Ia[i], Ib[i], Ic[i] = matrix_A @ I012
    for i in range(d.nbus):
        if i == id_fault:
            Va1[i] = Vpre_fault[i] * \
                (Z0[i, i]+Z2[i, i]+3*Zf)/(Z0[i, i]+Z1[i, i]+Z2[i, i]+3*Zf)
            Va2[i] = -Vpre_fault[i]*Z2[i, i] / \
                (Z0[i, i]+Z1[i, i]+Z2[i, i]+3*Zf)
            Va0[i] = -Vpre_fault[i]*Z0[i, i] / \
                (Z0[i, i]+Z1[i, i]+Z2[i, i]+3*Zf)
        else:
            Va1[i] = Vpre_fault[i] - Ia1[id_fault]*Z1[i, id_fault]
            Va2[i] = - Ia2[id_fault]*Z2[i, id_fault]
            Va0[i] = - Ia0[id_fault]*Z0[i, id_fault]
        V012 = array([Va0[i], Va1[i], Va2[i]])
        Vpos_fault.append(matrix_A @ V012)

elif ShortCircuit_type == 3:  # BC
    for i in range(d.nbus):
        Ia1[i] = Vpre_fault[i]/(Z1[i, i]+Z2[i, i]+Zf)
        Ia2[i] = -Ia1[i]
        I012 = array([Ia0[i], Ia1[i], Ia2[i]])
        Ia[i], Ib[i], Ic[i] = matrix_A @ I012
    for i in range(d.nbus):
        if i == id_fault:
            Va1[i] = Vpre_fault[i] - Ia1[i]*Z1[i, i]
            Va2[i] = -Ia2[i]*Z2[i, i]
            Va0[i] = 0
        else:
            Va1[i] = Vpre_fault[i] - Ia1[id_fault]*Z1[i, id_fault]
            Va2[i] = - Ia2[id_fault]*Z2[i, id_fault]
            Va0[i] = 0
        V012 = array([Va0[i], Va1[i], Va2[i]])
        Vpos_fault.append(matrix_A @ V012)

elif ShortCircuit_type == 4:  # BCT
    for i in range(d.nbus):
        Ia1[i] = Vpre_fault[i] / \
            (Z1[i, i]+Z2[i, i]*(Z0[i, i]+3*Zf)/(Z2[i, i]+Z0[i, i]+3*Zf))
        Ia2[i] = -(Vpre_fault[i] - Z1[i, i]*Ia1[i])/Z2[i, i]
        Ia0[i] = -(Vpre_fault[i] - Z1[i, i]*Ia1[i])/(Z0[i, i]+3*Zf)
        I012 = array([Ia0[i], Ia1[i], Ia2[i]])
        Ia[i], Ib[i], Ic[i] = matrix_A @ I012
    for i in range(d.nbus):
        if i == id_fault:
            Va1[i] = Ia1[i]*(Z2[i, i]*(Z0[i, i]+3*Zf)/(Z2[i, i]+Z0[i, i]+3*Zf))
            Va2[i] = -Ia2[i]*Z2[i, i]
            Va0[i] = -Ia0[i]*Z0[i, i]
        else:
            Va1[i] = Vpre_fault[i] - Ia1[id_fault]*Z1[i, id_fault]
            Va2[i] = -Ia2[id_fault]*Z2[i, id_fault]
            Va0[i] = -Ia0[id_fault]*Z0[i, id_fault]
        V012 = array([Va0[i], Va1[i], Va2[i]])
        Vpos_fault.append(matrix_A @ V012)

# REPORT
b = []
nChar = 150
for i in range(nChar):
    b += [" "]
b += ['\n']

# TITLE: RESULTS
reportData = list()
reportData.append('---------- Short circuit results - ' +
                  str(d.nbus) + '-bus system ----------\n\n')
reportData.append('Fault bus: ' + str(faultBus) + '\n\n')

if ShortCircuit_type == 1:
    reportData.append('Type and fault impedance: ABC - ' +
                      str(Zf) + ' ohms  \n\n')
elif ShortCircuit_type == 2:
    reportData.append('Type and fault impedance: AT - ' +
                      str(Zf) + ' ohms  \n\n')
elif ShortCircuit_type == 3:
    reportData.append('Type and fault impedance: AB - ' +
                      str(Zf) + ' ohms  \n\n')
elif ShortCircuit_type == 4:
    reportData.append('Type and fault impedance: ABT - ' +
                      str(Zf) + ' ohms  \n\n')

reportData.append('Fault Current: \n')
aux = copy(b)
reportData.append(aux)
reportData[-1][0:12] = 'Ia Magnitude'
reportData[-1][25:33] = 'Ia Angle'
reportData[-1][50:62] = 'Ib Magnitude'
reportData[-1][75:84] = 'Ib Angle'
reportData[-1][100:112] = 'Ic Magnitude'
reportData[-1][125:133] = 'Ic Angle'

Ibase = 100/(sqrt(3)*d.bus['voltage_base'][id_fault])

aux = copy(b)
reportData.append(aux)
# Ia
strIa_m = str(round(abs(Ia[id_fault][0]), 4)) + ' pu ' + \
    '(' + str(round(abs(Ia[id_fault][0])*Ibase, 2)) + ' kA)'
strIa_a = str(round(angle(Ia[id_fault][0])*180/pi, 2)) + ' degrees'
reportData[-1][0:len(str(strIa_m))] = strIa_m
reportData[-1][25:25+len(str(strIa_a))] = strIa_a
# Ib
strIb_m = str(round(abs(Ib[id_fault][0]), 4)) + ' pu' + \
    '(' + str(round(abs(Ib[id_fault][0])*Ibase, 2)) + ' kA)'
strIb_a = str(round(angle(Ib[id_fault][0])*180/pi, 2)) + ' degrees'
reportData[-1][50:50+len(str(strIb_m))] = strIb_m
reportData[-1][75:75+len(str(strIb_m))] = strIb_a

# Ic
strIc_m = str(round(abs(Ic[id_fault][0]), 4)) + ' pu' + \
    '(' + str(round(abs(Ib[id_fault][0])*Ibase, 2)) + ' kA)'
strIc_a = str(round(angle(Ic[id_fault][0])*180/pi, 2)) + ' degrees'
reportData[-1][100:100+len(str(strIc_m))] = strIc_m
reportData[-1][125:125+len(str(strIc_m))] = strIc_a

reportData.append('\n\nPost-Fault Voltage: \n')
aux = copy(b)
reportData.append(aux)
reportData[-1][0:3] = 'Bus'
reportData[-1][10:22] = 'Va Magnitude'
reportData[-1][35:43] = 'Va Angle'
reportData[-1][60:72] = 'Vb Magnitude'
reportData[-1][85:94] = 'Vb Angle'
reportData[-1][110:122] = 'Vc Magnitude'
reportData[-1][135:143] = 'Vc Angle'

for i in range(d.nbus):
    aux = copy(b)
    reportData.append(aux)
    bus = str(d.bus['num'][i])
    reportData[-1][0:len(bus)] = bus

    Vm_a = str(round(abs(Vpos_fault[i][0][0]), 4)) + ' pu (' + str(
        round(abs(Vpos_fault[i][0][0])*d.bus['voltage_base'][i], 2)) + ' kV)'
    Va_a = str(round(angle(Vpos_fault[i][0][0])*180/pi, 2)) + ' degrees'
    reportData[-1][10:10+len(Vm_a)] = Vm_a
    reportData[-1][35:35+len(Va_a)] = Va_a

    Vm_b = str(round(abs(Vpos_fault[i][1][0]), 4)) + ' pu (' + str(
        round(abs(Vpos_fault[i][1][0])*d.bus['voltage_base'][i], 2)) + ' kV)'
    Va_b = str(round(angle(Vpos_fault[i][1][0])*180/pi, 2)) + ' degrees'
    reportData[-1][60:60+len(Vm_b)] = Vm_b
    reportData[-1][85:85+len(Va_b)] = Va_b

    Vm_c = str(round(abs(Vpos_fault[i][2][0]), 4)) + ' pu (' + str(
        round(abs(Vpos_fault[i][2][0])*d.bus['voltage_base'][i], 2)) + ' kV)'
    Va_c = str(round(angle(Vpos_fault[i][2][0])*180/pi, 2)) + ' degrees'
    reportData[-1][110:110+len(Vm_c)] = Vm_c
    reportData[-1][135:135+len(Va_c)] = Va_c


nlines = len(reportData)
fullpathname = os.path.join(os.path.dirname(
    __file__), 'Results', os.path.splitext(systemName)[0])
if not os.path.exists(fullpathname):
    os.makedirs(fullpathname)
if ShortCircuit_type == 1:
    fullpathname2 = os.path.join(os.path.dirname(
        __file__), 'Results', os.path.splitext(systemName)[0], 'SC_faultbus_' + str(faultBus) + '_type_ABC_ZF_'+str(Zf)+'_ohms.txt')
    f = open(f'{fullpathname2}', 'w', encoding='latin-1')
    for k in range(nlines):
        f.write("".join(reportData[k]))
    f.close()
elif ShortCircuit_type == 2:
    fullpathname2 = os.path.join(os.path.dirname(
        __file__), 'Results', os.path.splitext(systemName)[0], 'SC_faultbus_' + str(faultBus) + '_type_AT_ZF_'+str(Zf)+'_ohms.txt')
    f = open(f'{fullpathname2}', 'w', encoding='latin-1')
    for k in range(nlines):
        f.write("".join(reportData[k]))
    f.close()
elif ShortCircuit_type == 3:
    fullpathname2 = os.path.join(os.path.dirname(
        __file__), 'Results', os.path.splitext(systemName)[0], 'SC_faultbus_' + str(faultBus) + '_type_BC_ZF_'+str(Zf)+'_ohms.txt')
    f = open(f'{fullpathname2}', 'w', encoding='latin-1')
    for k in range(nlines):
        f.write("".join(reportData[k]))
    f.close()
elif ShortCircuit_type == 4:
    fullpathname2 = os.path.join(os.path.dirname(
        __file__), 'Results', os.path.splitext(systemName)[0], 'SC_faultbus_' + str(faultBus) + '_type_BCT_ZF_'+str(Zf)+'_ohms.txt')
    f = open(f'{fullpathname2}', 'w', encoding='latin-1')
    for k in range(nlines):
        f.write("".join(reportData[k]))
    f.close()
