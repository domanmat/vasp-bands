#!/usr/bin/python

#python script for: (1) printing BS figures, (2) generating band_data.dat (only of max 2 projections),
#                   and (3) comparison of band plots from wannier90.dat and VASP

import sys
import os
import matplotlib.pyplot as plt
import time
import numpy as np
import math
import copy
import csv
from termcolor import colored
start_time = time.time()
plt.rcParams["font.family"] = "Arial"
cm = 1/2.54

# #small plot settings
# figure_width=8 #cm
# figure_height=6 #cm
# plot_width=0.6 #originally 0.8

#medium plot settings
figure_width=12 #cm
figure_height=8 #cm
plot_width=0.6 #originally 0.8

# #large plot settings
# # figure_width=12 #cm
# figure_width=14 #cm
# figure_height=12 #cm
# plot_width=0.6
# # figure_width=24 #cm

plot_spin_channels = 0 # 1 - yes, figure with different color lines 0 - no, all black

if_manual_fermi = 1 # 1 - yes, manual EF; 0 - no, read
# manual_fermi = 0 #for wannier90
# manual_fermi = 3.33448347 #doping 50%
# manual_fermi = 4.00886838 #doping 10%
# manual_fermi = 2.06895214 #NiH2--5xNaF--1xCaF GGA_metal z DOS
# manual_fermi = 4.2175         #LiNiH2 GGA+U_metal z DOS
# manual_fermi = 9.6008         #LaNiO2 GGA+U_metal z DOS
# manual_fermi = 9.67408515     #LaNiO2 GGA_metal z DOS
manual_fermi = 9.69701759     #LaNiO2 SCAN_metal z DOS
# manual_fermi = 6.4157         #NdNiO2 GGA+U_metal z DOS
# manual_fermi = 4.5673         #CaCuO2 GGA+U_metal z DOS
# manual_fermi = 4.2175         #LiNiH2 GGA+U_metal z DOS
# manual_fermi = 3.153          #NaNiH2 GGA+U_metal z DOS
# manual_fermi = 4.6209         #SrNiH3 GGA+U_metal z DOS
# manual_fermi = 4.35611968       #LiNiH2 SCAN_metal z DOS
# manual_fermi = 4.36394605     #CaCuO2 GGA+U AF-C z DOS
# manual_fermi = 8.95858846     #LaNiO2 GGA+U AF-C z DOS
# manual_fermi = 2.90579828     #NaNiH2 GGA+U AF-C z DOS
# manual_fermi = 2.93334439     #NaNiH2 GGA+U AF-G z DOS
# manual_fermi = 4.5909059      #LiNiH2 GGA+U AF-C z DOS
# manual_fermi = 1.05250437     #GGA_metal NM CaAgSO4 z DOS
# manual_fermi = 0.52172137     #GGA+U AF CaAgSO4 z DOS


y_min=-10
y_max=6
y_step=2 #krok na osi y
scale=5 #skala markerów projekcji
alpha=1
colors = ['b', 'lime', 'r', 'c', 'm', 'y', 'k', 'g', 'w']

# chosen_folder=input("Podaj ścieżkę folderu:\n")
# folder=chosen_folder

folder=r"C:\Users\Mateusz\Documents\MobaXterm\home\For_Jose_01_08\NiLiH2-gga_u6-is3_sq2_sc112_015_af1-prim_nm-BS_rw"

# input_print="O.p Cu.dxy"
# input_print="O.p Cu.dx2-y2"
# input_print="Cu.dx2-y2.dz2"
# input_print="O.p Cu.dxy Cu.dz2"
# input_print="O.p Ni.dxy Ni.dz2"
# input_print="O.p Ni.dx2-y2"
# input_print="O.p Ni.dx2-y2 Ni.dz2"
# input_print="H.t Ni.t"
input_print="H.s Ni.dx2-y2"
# input_print="H.s Ni.dxy Ni.dz2"
# input_print="H.s Ni.dx2-y2 Ni.dz2"
# input_print="O.p Cu.dx2-y2"
# input_print="O.p Cu.dx2-y2 Cu.dz2"
# input_print="O.p Ni.dx2-y2"
# input_print="O.p Ni.dx2-y2 Ni.s"
# input_print="O.p Ni.dx2-y2 Nd.t"
# input_print="O.p Ni.dx2-y2 Ni.dz2"
# input_print="O.p Ni.dz2"
# input_print="O.p Ni.t"
# input_print="O.p Ca.dz2"
# input_print="H.t Ni.dx2-y2 Ni.s"
# input_print="H.t Ni.dx2-y2 Ni.s"
# input_print="H.t Ni.dx2-y2"
# input_print="H.t Ni.s"
# input_print="H.t Ni.dz2"
# input_print="H.t Ni.dxy"
# input_print="H.t Ni.dx2-y2 Ni.dz2"
# input_print="H.s Ni.dx2-y2"
# input_print="S.t Ni.dx2-y2"
# input_print="H.t Ni.t"
# input_print="H.t Ni.dz2"
# input_print="O.p Ag.d Cu.d"
# input_print="O.p Cu.d"
# input_print="O.p Ag.dxy.dyz"
# input_print="O.p Cu.dx2-y2"
# input_print="O.p Cu.dxy"
# input_print="V.t C.t.N.t"
# input_print="V.t C.t N.t"
# input_print="C.t N.t"
# input_print="V.t"
# input_print="V1.dyz"


# input_read=input("Podaj atomy i orbitale projekcji:\n")
# input_print=input_read
#s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2  fy3x2   fxyz   fyz2    fz3   fxz2   fzx2    fx3    tot

# zrobić modyfikacje ładującą atomy danego typu osobno
#

procar_dir=os.path.join(folder,"PROCAR")
poscar_dir=os.path.join(folder,"POSCAR")
outcar_dir=os.path.join(folder,"OUTCAR")
incar_dir=os.path.join(folder,"INCAR")
doscar_dir=os.path.join(folder,"DOSCAR")
kpoints_dir=os.path.join(folder,"KPOINTS")
vasprun_dir=os.path.join(folder,"vasprun.xml")
wannier_dir=os.path.join(folder,"wannier90_band.dat")
band_data_dir=os.path.join(folder,"band_data "+input_print+".dat")
#printing of band_data.dat works only for len(input_print)=2, i.e. only 2 printed projections
try:
    os.remove(band_data_dir)
except OSError:
    pass
if len(input_print.split()) == 2:
    with open(band_data_dir, 'a', newline='') as f_data:
        writer = csv.writer(f_data, delimiter='\t')
        data=['band_index','kpoint_index','kpoint_dist','kpoint_x','kpoint_y','kpoint_z','eigenvalues','mixing(0-1)','proj_sum']
        writer.writerow(data)

f_incar = open(incar_dir, "r")
incar_lines=f_incar.readlines()
ispin=''
for line in incar_lines:
    if "ISPIN" in line:
        ispin=line.split()[-1]
print("ISPIN =",ispin)

f_outcar = open(outcar_dir, "r")
outcar_lines=f_outcar.readlines()
lmax=[]
for line in outcar_lines:
    if " LMAX " in line:
        lmax.append(int(line.split()[-1]))
print("LMAX =", max(lmax))
if max(lmax) == 6:
    # print("1")
    orbital_no=[0,1,2,3,4,5,6,7,8,9]
    projection_symbol=[['s'], ['py'], ['pz'], ['px'], ['dxy'], ['dyz'], ['dz2'], ['dxz'], ['dx2-y2'], ['tot']]
elif max(lmax) == 8:
    # print("2")
    orbital_no = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    projection_symbol=[['s'], ['py'], ['pz'], ['px'], ['dxy'], ['dyz'], ['dz2'], ['dxz'], ['dx2-y2'],
                       ['fy3x2'], ['fxyz'], ['fyz2'], ['fz3'], ['fxz2'], ['fzx2'], ['fx3'], ['tot']]
elif max(lmax) == 4:
    # print("3")
    orbital_no = [0, 1, 2, 3, 4, 5]
    projection_symbol=[['s'], ['py'], ['pz'], ['px'], ['tot']]
else:
    print("ERROR reading lmax")
# print("Możliwe projekcje:",projection_symbol)
for l in range(len(projection_symbol)):
    projection_symbol[l].append(orbital_no[l])
print(projection_symbol)

f_poscar = open(poscar_dir, "r")
poscar_lines=f_poscar.readlines()

poscar_scale=float(poscar_lines[1])
print(poscar_scale)
p=poscar_lines[2].split() #p to tupla z 1-szej linii rozdzielonej spacjami, w kazdej, kolejne ich elementy to kolumny
q=poscar_lines[3].split() #p to tupla z 2-giej linii
r=poscar_lines[4].split() #p to tupla z 3-ciej linii
print (p,"\n",q,"\n",r)
for i in range(3): #zaokragla liczby 1,2,3, z tupli (wierszy) do zmiennoprzecinkowych
    p[i]=float(p[i])
    q[i]=float(q[i])
    r[i]=float(r[i])
v_length = []
for w in [p,q,r]: #dla kazdej wartosci 'w' z wierszy p,q,r liczy 'const' - stala sieci (dlugosc wektora w kier. 0,1,2)
    const=math.sqrt(w[0]**2+w[1]**2+w[2]**2)*poscar_scale
    print(const)
    v_length.append(const)
alfa=(math.acos((q[0] * r[0] + q[1] * r[1] + q[2] * r[2]) * poscar_scale ** 2 / (v_length[1] * v_length[2])) / math.pi * 180) #printuje katy alfa, beta, gamma
beta =(math.acos((p[0] * r[0] + p[1] * r[1] + p[2] * r[2]) * poscar_scale ** 2 / (v_length[0] * v_length[2])) / math.pi * 180)
gamma=(math.acos((p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * poscar_scale ** 2 / (v_length[0] * v_length[1])) / math.pi * 180)
volume=float(v_length[0] * v_length[1] * v_length[2] * math.sqrt(1 - math.pow(math.cos(0.01745 * alfa),2) -
                  math.pow(math.cos(0.01745 * beta),2) -
                  math.pow(math.cos(0.01745 * gamma),2) +
                 2*math.cos(0.01745 * alfa) * math.cos(0.01745 * beta) *
                                             math.cos(0.01745 * gamma)))
print(volume)
recip_a=math.pi * 2 * v_length[1] * v_length[2] * math.sin(alfa * 0.01745329) / volume
recip_b=math.pi * 2 * v_length[0] * v_length[2] * math.sin(beta * 0.01745329) / volume
recip_c=math.pi * 2 * v_length[0] * v_length[1] * math.sin(gamma * 0.01745329) / volume
recip_vectors=[]
recip_vectors+=(recip_a,recip_b,recip_c)
print(recip_vectors)

def split_list_slicing(input_list):
    length = len(input_list)
    middle = length // 2
    first_half = input_list[:middle]
    second_half = input_list[middle:]
    return first_half, second_half

# atom_types=[]
# atom_types.append(poscar_lines[5].split())
atom_types=poscar_lines[5].split()
print(atom_types)
print(len(atom_types))
atom_list=[]
atoms_no=poscar_lines[6].split()
print(atoms_no)
atom_index=0
for i in range(len(atom_types)):
    atom_list.append([])
    atom_list[i].append(atom_types[i])
    atom_list[i].append(int(atoms_no[i]))
    atom_list[i].append([])
    for j in range(atom_list[i][1]):
        atom_list[i][2].append(atom_index)
        atom_index+=1
print(atom_list)
print('hello')
# atom_list=[i][0] - nazwy typów
# atom_list=[i][1] - liczba atomów danego typu
# atom_list=[i][2] - numery atomów
# atom_list=[['Ni',1,[0]],['H',2,[1,2]],['Li',1,[3]]]

atoms_sum=0
for i in range(len(atom_list)):
    atoms_sum+=int(atom_list[i][1])
print("Number of ions: ",atoms_sum)

f_vasprun=open(vasprun_dir, "r")
vasprun_lines=f_vasprun.readlines()
efermi=0
for line in vasprun_lines:
    if "efermi" in line:
        efermi=float(line.split()[2])
print("Energy of Fermi level:", efermi)
if if_manual_fermi == 1:
    efermi=manual_fermi


#### do DOSCAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
f_doscar=open(doscar_dir, "r")
doscar_data=f_doscar.readlines()
header=doscar_data[5]
# print('lol')
# print('lol')
# print(header)
# print(header)
# print(header)
#### do DOSCAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# with open(procar_dir, "r") as f_procar:
#     data = f_procar.read()
f_procar=open(procar_dir, "r")
procar_data=f_procar.read()
all_kpoints_array=[]
# procar_kpoints=procar_text.split("k-point ")
procar_kpoints=procar_data.split("k-point ")
# print(procar_kpoints[1])
# print(procar_kpoints[1].split("band ")[2])
for kpoint in procar_kpoints:
    kpoint_array=[]
    linesofkpoint=kpoint.split("\n")
    for line in linesofkpoint:
        if "energy" in line:
            kpoint_array.append(float(line.split()[4])-efermi)
    all_kpoints_array.append(kpoint_array)

#print(kpoint_array)
#all_kpoints_array to energie bands dla każdego kpuntku
all_kpoints_array.remove(all_kpoints_array[0])
#print(all_kpoints_array)
no_kpoints=len(all_kpoints_array)
no_bands=len(all_kpoints_array[0])
print("Number of kpoints: ", no_kpoints)
print("Number of bands:   ", no_bands)


kpoint_lst=[]
kpoint_band=[]
#kpoint_band zawiera teksty wszystkich bandów w k-punktach
#w kolejności kpoint_band[no_kpoints][no_bands]
for kpoint_no in range(0,(no_kpoints)):
    kpoint_band.append(procar_kpoints[kpoint_no+1].split("band "))
    kpoint_band[kpoint_no].pop(0)
    kpoint_lst.append(kpoint_no)
#print(kpoint_band[79]) # ok
# print(len(kpoint_band))
    #kpoint_band2.append(kpoint_band[i].split("ion "))

#print(kpoint_band[0][1])
ion_band=[]
for kpoint in range(0,no_kpoints): #dłg. = no_kpoints
    ion_band.append([])
    for band in range(0,no_bands):
        ion_band[kpoint].append([])
        #print(kpoint_band[kpoint][j])
        linesofbands = kpoint_band[kpoint][band].split("\n")
        linesofbands.pop(0)  #usuwanie zbędnych linii, od góry 2 puste i linia z projekcjami
        linesofbands.pop(0)
        linesofbands.pop(0)
        linesofbands.pop(-1) #usuwanie zbędnych linii, od dołu 2 puste
        linesofbands.pop(-1)
        #print(i,j)
        #print(linesofbands)
        for ion in range(0,atoms_sum):
            #print(i,j)
            ion_band[kpoint][band].append([])
            ion_band[kpoint][band][ion]=(linesofbands[ion].split())
            ion_band[kpoint][band][ion].pop(0) #usuwa pierwszy element każdej linii, czyli np. nazwe atomu

# ion_band[k_point][band][proj_atoms][proj_orbital]
# print(ion_band[0][0][13][9])
# print(ion_band[79][59][0][8]) #ok, zgadza się

# #wszystko - liczba k-pointów
# print(len(ion_band))
# #kpoint - liczba pasm
# print(len(ion_band[0]))
# #kpoint, band - liczba jonów
# print(len(ion_band[0][0]))
# #kpoint, band, jon - liczba projekcji
# print(len(ion_band[0][0][0]))
# #kpoint, band, jon, projekcja - wartości, liczba znaków tej wartości
# print(len(ion_band[0][0][0][0]))
# print('hello')
# #
# print(ion_band[-1]) #wszystkie jony, wszystkie bandy, z ostatniego k-pointu
# print(ion_band[-1][-1]) #wszystkie jony, ostatniego bandu, z ostatniego k-pointu
# print(ion_band[-1][-1][-1]) #ostatni jon, ostatniego bandu, z ostatniego k-pointu
# print(ion_band[-1][-1][-1][-1]) #ostatnia wartość, ostatniego jonu, ostatniego bandu, ostatniego k-pointu

# print('hello')
for kpoint in range(len(ion_band)):                      # długość = liczba k-pointów
    for band in range(len(ion_band[kpoint])):               # długość = liczba pasm
        # print(ion_band[i][j][0][9])
        for ion in range(len(ion_band[kpoint][band])):        # długość = liczba jonów
            #print(len(ion_band[i][j][k]))          # długość = liczba projekcji
            if len(ion_band[kpoint][band][ion])!=len(projection_symbol):
                print(len(ion_band[kpoint][band][ion]))
                # print("error")

# #ion_band[i][j][k], i=no_kpoints, [j]=no_bands, [k]=atoms_sum, każda komórka ma 10 projekcji
# print('hello')

#możliwe projekcje
#if
# for i in range(len(projection_symbol)):
#     print(projection_symbol[i][0],"=",projection_symbol[i][1]) #symbole i numery projekcji

# to_print=input("Podaj atomy i orbitale do narysowania w formacie Ag.x2-y2.dxy F.px :\n")
# input_print="Ni.pz.py.px" --- not working
mySeparator="."
input_x=mySeparator.join(atom_types)
print(input_x)

input_strings=input_print.split()
print(input_strings)
input_proj=[]
input_atom=[]
# input_atom=[]
for i in range(len(input_strings)):
    input_proj.append([])
    input_all=input_strings[i].split(".")
    input_atom.append(input_all.pop(0))
    input_proj[i].append(input_atom[i]) #dodawanie nazwy atomu
    for j in range(len(atom_list)):
        if input_proj[i][0] in atom_list[j][0]:
            input_proj[i].append(atom_list[j][2]) #dodawanie indeksów atomów do typu

for i in range(len(input_strings)):
    input_all=input_strings[i].split(".")
    input_all.pop(0) #usuwa symbol atomu
    input_proj[i].append([])
    input_proj[i].append([])
    for j in range(len(input_all)):
        for l in range(len(projection_symbol)):
            if input_all[j] in projection_symbol[l][0]: #szuka w symbolach projekcji
                input_proj[i][2].append(input_all[j])
                input_proj[i][3].append(projection_symbol[l][1])
    # input_proj[i].append(input_all)
print("atom_list:", atom_list)
print("input_proj: ", input_proj)
#input proj w zamierzeniu =[['Ni',[0],'tot',[9]],['H',[1,2],'s',[0]]]
#
# to_print=[[[0],[9]],[[1,2],[0]]]
# to_print=[['Ni',[0],[9]],['H',[1,2],[0]]]

### k-path ticks on x-axis part
#################################
f_kpoints=open(kpoints_dir, "r")
kpoints_line1 = ""
kpoints_line2 = ""
kpoints_line1 += f_kpoints.readlines(1)[0] #czyta 1. linię, ale potem idzie dalej jbc
kpoints_line2 += f_kpoints.readlines(1)[0] #czyta 1. linię, ale potem idzie dalej jbc
kpoints_coordinates = f_kpoints.readlines()[2:]
#
#czytanie współrzędnych k-punktów
k_coordinates=[]
k_coordinates_cart=[]
to_delete=[]
for line in kpoints_coordinates:
    coord=line.split()[:3]
    # print(coord)
    k_coordinates.append(coord)
for i in range(len(k_coordinates)):
    if len(k_coordinates[i]) == 0:
        # print("zero")
        to_delete.append(i)
# print(len(to_delete),"\n")
for i in range(len(to_delete),0,-1):
    # print(i)
    # print(to_delete[i-1],"\n")
    k_coordinates.pop(to_delete[i-1])
# odcinki oryginalne k-path
k_coordinates_recip=copy.copy(k_coordinates)
# print()
# print(len(k_coordinates_recip))
# print(k_coordinates_recip)

#usuwa co drugą linię od końca - aby nie było przerwy
for i in range(int(len(k_coordinates)/2)-1):
    # print(i+2)
    k_coordinates.pop(-(i+2))
# print(k_coordinates)
# for i in range(len(k_coordinates)):
#     # print(k_coordinates[i]) #każda linia z tuplą 3 współrzędnych
#     k_coordinates_cart.append([float(x) * y for x, y in zip(k_coordinates[i], recip_vectors)])
for i in range(len(k_coordinates)):
    # print(k_coordinates[i]) #każda linia z tuplą 3 współrzędnych
    k_coordinates_cart.append([float(x) * y for x, y in zip(k_coordinates[i], recip_vectors)])
    #     k_coordinates_cart.append(((x-y) * b)**2 for x, y, b in zip(k_coordinates[i], k_coordinates[i-1], recip_vectors)])
# print(len(k_coordinates_cart))
# print('\n',k_coordinates_cart)
kpath_cart=[]
k_sum=0
k_dist=[]
for i in range(len(k_coordinates_cart)):
    # kpath_cart.append([])
    if i >= 1:
        # print(k_coordinates_cart[i])
        # print(sum(k_coordinates_cart[i])-sum(k_coordinates_cart[i-1]))
        sum_sq=[]
        for j in range(len(k_coordinates_cart[i])):
            sum_sq.append(math.pow((k_coordinates_cart[i][j]-k_coordinates_cart[i-1][j]),2))
        k_dist.append(math.sqrt(sum_sq[0]+sum_sq[1]+sum_sq[2]))
        k_sum+=k_dist[i-1]
    kpath_cart.append(k_sum)
print()
# print(k_dist)c
# print(k_coordinates)
# print(len(k_coordinates))
#czytanie współrzędnych k-punktów
#
#old part
print(kpoints_line1)
print(kpoints_line2)
pts_per_section=int(kpoints_line2.split()[0])-1
kpaths= kpoints_line1.split()[2:]
print('points per line:', pts_per_section+1)
#
all_kpoints_coord=[]
# all_kpoints_coord+=kpath_cart[0]
for i in range(0,len(kpath_cart)-1):
    odcinek=float(k_dist[i] / pts_per_section)
    # print(odcinek)
    all_kpoints_coord.append(float(kpath_cart[i]))
    for k in range(1, pts_per_section + 1):
        # print(k)
        all_kpoints_coord.append(float(kpath_cart[i]+(int(k)*odcinek)))
# print(all_kpoints_coord[239:280])
# print(all_kpoints_coord) #współrzędne na osi x
print("all_kpoints_coord:",len(all_kpoints_coord))

#printowanie k-puntków w sieci odwrotnej = k_coordinates_recip_all
# print(len(k_coordinates_recip))
# print(k_coordinates_recip)
k_coordinates_recip_all=[]
#pętla po liczbie odcinków w przestrzeni odwrotnej
for i in range(int(len(k_coordinates_recip)/2)):
    # print(i)
    spacing=[]
    for k in range(0, pts_per_section + 1):
        k_coordinate=[]
        for j in range(3):
            if k == 0:
                #tylko dla k=0 bo nie ma sensu dla każdego tego liczyć - będzie takie samo
                spacing.append(float(float(k_coordinates_recip[i*2+1][j]) - float(k_coordinates_recip[i*2][j])) / pts_per_section)
                # print(j, k_coordinates_recip[i*2][j], k_coordinates_recip[i*2+1][j], spacing[j])
            #k_coordinate to 3 współrzędne każdego k-punktu w odcinku
            k_coordinate.append(float(k_coordinates_recip[i*2][j])+(k * spacing[j]))
        k_coordinates_recip_all.append(k_coordinate)
# print(len(k_coordinates_recip_all))
# print(k_coordinates_recip_all)


# print(kpaths)
x_paths=[]
# print(len(kpaths)-1)
for i in range(0, (len(kpaths))):
    # print(i)
    x_paths.append(kpaths[i].split("-"))
    if len(kpaths) > 1 and i < (len(kpaths) - 1):
        x_paths[i][-1] += ' | '
    # print(x_paths[i])
# print(x_paths)

for i in range(0, (len(kpaths) - 1)):
    if i < (len(kpaths) - 1):
        x_paths[i][-1] += x_paths[i + 1][0]
        x_paths[i + 1].pop(0) #usuwa 1-szy element
# print(x_paths)
x_labels=[]
for i in range(0, len(x_paths)):
    for j in range(0, len(x_paths[i])):
        x_labels.append(x_paths[i][j])
print(x_labels)

# x_positions=[]
# for k in range(0, len(x_labels)):
#     x_positions.append(float(k) * pts_per_section)
# print(x_positions)
x_positions=[]
for k in range(0, len(x_labels)):
    x_positions.append(kpath_cart[k])
print(x_positions)


#wykres samych pasm (eigenvalues) na czarno
bands_eigenvalues=[]
for band in range(0,no_bands):
    # print(i)
    bands_eigenvalues.append([])
    for kpoint in range(0,int(no_kpoints)):
        energy_to_add = all_kpoints_array[kpoint][band]
        bands_eigenvalues[band].append(energy_to_add)
    # #długość to liczba k-pointów
    # print(len(band_array[i]))
    # # plt.plot(band_array,color=colors[i])
    # plt.plot(kpoint_lst,band_array[i],color="black",linewidth=0.5)
# #długość to liczba pasm
# print(len(band_array))

projection_print = []
# #robimy tablicę projection_print[k][n][l] - wybieram projekcje do wykresu
for proj_no in range(len(input_proj)):     # pętla po liczbie elementów (atomów) do wyrysowania
    projection_print.append([])
    # print('n',n)
    for l in range(len(input_proj[proj_no][3])):      # pętla po liczbie orbitali danego do dodania
        projection_print[proj_no].append([])
        for x in range(len(input_proj[proj_no][1])):  # pętla po indeksach atomów danego typu
            # print(input_proj[n][3]) #to numery orbitali do projekcji z inputu
            # print(input_proj[n][1]) #nr-y atomów danego typu do projekcji z inputu
            # #n,l,x - wartości z inputu
            projection_print[proj_no][l].append([])
            proj_orbital=input_proj[proj_no][3][l]
            proj_atoms=input_proj[proj_no][1][x]
            for band in range(0, no_bands):        # pętla po liczbie pasm
                projection_print[proj_no][l][x].append([])
                for k_point in range(0,int(no_kpoints)):  # pętla po liczbie k-punktów
                    # #dodajemy j wartości dla każdego projection_print[k][n][l][i]
                    projection_print[proj_no][l][x][band].append(scale * float(
                    ion_band[k_point][band][proj_atoms][proj_orbital]))

                    # print(projection_print[n][l])
# #ion_band[i][j][k][l] to wszystkie dane z PROCAR
# #ion_band[i][j][k], i=no_kpoints, [j]=no_bands, [k]=atoms_sum, [l] każda komórka ma 10 projekcji
# #ion_band[k_point][band][proj_at][proj_orbital]
# #projection_print[proj_no][l][x][band] to wszystkie dane do rysowania
# #projection_print[proj_no][l][x][band] [proj_no]=elementy do projekcji, [l] nr-y orbitali do proj,
# # [x] nr-y atomów do projekcji, [band] pasma do projekcji
# print('orbitals in proj. element 0:', len(projection_print[0]))
# print('atoms in proj. element 0, orbital 0:', len(projection_print[0][0]))
# print('bands in proj. element 0, orbital 0, atom 0:', len(projection_print[0][0][0]))
# print('kpoints in proj. element 0, orbital 0, atom 0, band 0:', len(projection_print[0][0][0][0]))
# #wymiary są ok


fig1, ax1 = plt.subplots()
legend=[]
plot_legend=[]
plot=[]
plot_legend_key=[]
# kpoint_plot=[]
proj_no=0
band=0
#rysowanie - nieaktywne
for proj_no in range(len(input_proj)): # pętla po elementach do projekcji - 3
    final_plot=[]
    plot.append([])
    plot_legend_key.append([])
    legend.append(input_strings[proj_no])
    #PRINTING
    print()
    print('projection no.:', proj_no+1)
    print('atom type:', input_proj[proj_no][0])
    print('included orbitals:', len(input_proj[proj_no][3]), input_proj[proj_no][2],
          "orbital index", input_proj[proj_no][3])
    print('included atoms:', len(input_proj[proj_no][1]),
          "atom index", input_proj[proj_no][1])
    for band in range(0, no_bands):
        plot[proj_no].append([])
        # kpoint_plot.append(kpoint_lst)
        # print('\nband number:', band + 1)
        final_plot.append([])
        for k_point in range(0, no_kpoints):
            final_plot[band].append(0)
            # #final_plot to pusta tablica zer no_bands x no_kpoints, dla każdego pasma
        # print(len(final_plot[band])) - liczba no_kpoints
        for orbital in range(len(input_proj[proj_no][3])): #pętla po indeksach orbitali inputu
            for proj_at in range(len(projection_print[proj_no][orbital])):
                # liczba atomów w projekcji, minimum 1
                # proj_at nie jest zgodny z input_proj[proj_no][1] - bo w projekcji inna numeracja
                # #PRINTING
                # print('atom index:', input_proj[proj_no][1][proj_at])
                # print(len(projection_print[proj_no][orbital][proj_at][band]), "points to add")
                first=final_plot[band]
                second=projection_print[proj_no][orbital][proj_at][band]
                final_plot[band]=[x + y for x, y in zip(first,second)]

        # zmienna plot - niepotrzebna
        # #ważna linia!!!!
        # plot[proj_no][band]=plt.scatter(kpoint_lst, bands_eigenvalues[band], s=final_plot[band],
        #                                 color=colors[proj_no], alpha=alpha)
        # opcjonalnie
        # linewidths=0.2, edgecolors='black')

    # plot_legend.append(plot[proj_no][0])
    #dodatkowy niewidoczny wykres wyłącznie do ładnej legendy
    plot_legend_key[proj_no] = plt.scatter([-10], [0], s=scale/5, color=colors[proj_no])
    plot_legend.append(plot_legend_key[proj_no])

# plt.figure(1)
# 3 wymiar (jako colormap) gdzie będzie np. mixing - color
print()
print("Number of band projections:", len(bands_eigenvalues[band]))
print()
proj_sum=[]
proj_sum_sq=[]
band_sum=[]
band_sum_sq=[]
diff=[]
mixing=[]
for band in range(0, no_bands):
    proj_sum.append([])
    proj_sum_sq.append([])
    band_sum.append([])
    band_sum_sq.append([])
    diff.append([])
    mixing.append([])
    for proj_no in range(len(input_proj)): # pętla po elementach do projekcji
        proj_sum[band].append([])
        proj_sum_sq[band].append([])
    for k_point in range(0, no_kpoints):
        if len(input_proj) != 3:
            mixing.append(0)
        for proj_no in range(len(input_proj)):
            proj_sum[band][proj_no].append(0)
            proj_sum_sq[band][proj_no].append(0)
    for proj_no in range(len(input_proj)): # pętla po elementach do projekcji
        for orbital in range(len(input_proj[proj_no][3])):
            for proj_at in range(len(input_proj[proj_no][1])):
                first = proj_sum[band][proj_no]
                second = projection_print[proj_no][orbital][proj_at][band]
                proj_sum[band][proj_no] = [x + y for x, y in zip(first, second)]
    if len(input_proj) == 2:
        band_sum[band] = [x + y for x, y in zip(proj_sum[band][0], proj_sum[band][1])]
        band_sum_sq[band] = [x**2 for x in band_sum[band]]
        diff[band] = [x - y for x, y in zip(proj_sum[band][0], proj_sum[band][1])]
        mixing[band]= [(sum-diff) / (sum+0.0000001) / 2 for sum, diff in zip(band_sum[band], diff[band])]

        # if (len(all_kpoints_coord)) != len(bands_eigenvalues[band]) :
        #     all_kpoints_coord2=all_kpoints_coord
        #     all_kpoints_coord.extend(all_kpoints_coord2)
        #     print(all_kpoints_coord)
        #plt.scatter(x=kpoint_lst, y=plt_bnds, s=plt_sum, c=plt_mix, cmap='brg', vmin=0., vmax=1.,  alpha=1.0)
        if (len(all_kpoints_coord)) != len(bands_eigenvalues[band]) :
            # first_half = input_list[:middle]
            # second_half = input_list[middle:]
            bands_eigenvalues_spin = split_list_slicing(bands_eigenvalues[band])
            mixing_spin = split_list_slicing(mixing[band])
            band_sum_sq_spin = split_list_slicing(band_sum_sq[band])
            im    = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues_spin[0], c=mixing_spin[0], s=band_sum_sq_spin[0],
                             cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)
            im_s2 = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues_spin[1], c=mixing_spin[1], s=band_sum_sq_spin[1],
                         cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)
        else:
            im = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues[band], c=mixing[band], s=band_sum_sq[band],
                    cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)

        # writing output to file for Antonio
        # f_data = open(band_data_dir, "w")
        with open(band_data_dir, 'a', newline='') as f_data:
            writer = csv.writer(f_data, delimiter='\t')
            band_index=[(band+1) for i in range(len(all_kpoints_coord))]
            kpoint_index = [(i + 1) for i in range(len(all_kpoints_coord))]
            kpoint_x = [(k_coordinates_recip_all[i][0]) for i in range(len(all_kpoints_coord))]
            kpoint_y = [(k_coordinates_recip_all[i][1]) for i in range(len(all_kpoints_coord))]
            kpoint_z = [(k_coordinates_recip_all[i][2]) for i in range(len(all_kpoints_coord))]
            band_sum_norm = [i/scale for i in (band_sum[band])]
            writer.writerows(zip(band_index, kpoint_index, all_kpoints_coord, kpoint_x, kpoint_y, kpoint_z,
                                 bands_eigenvalues[band], mixing[band], band_sum_norm))

    elif len(input_proj) == 3:
        # band_sum_tot = band_sum - to tylko zmienia nazwę, a referencja pozostaje ta sama, aktualizacja tych samych zmiennych
        band_sum_tot = copy.copy(band_sum)
        band_sum_12 = copy.copy(band_sum)
        band_sum_13 = copy.copy(band_sum)
        band_sum_23 = copy.copy(band_sum)
        band_sum_tot[band] = [x + y + z for x, y, z in zip(proj_sum[band][0], proj_sum[band][1], proj_sum[band][2])]
        band_sum_12[band] = [x + y for x, y in zip(proj_sum[band][0], proj_sum[band][1])]
        band_sum_13[band] = [x + y for x, y in zip(proj_sum[band][0], proj_sum[band][2])]
        band_sum_23[band] = [x + y for x, y in zip(proj_sum[band][1], proj_sum[band][2])]
        band_sum_sq[band] = [x**2 for x in band_sum_tot[band]]
        diff_12 = copy.copy(diff)
        diff_13 = copy.copy(diff)
        diff_23 = copy.copy(diff)
        diff_12[band] = [x - y for x, y in zip(proj_sum[band][0], proj_sum[band][1])]
        diff_13[band] = [x - y for x, y in zip(proj_sum[band][0], proj_sum[band][2])]
        diff_23[band] = [x - y for x, y in zip(proj_sum[band][1], proj_sum[band][2])]
        frac_1 = copy.copy(mixing)
        frac_2 = copy.copy(mixing)
        frac_3 = copy.copy(mixing)
        frac_1[band] = [proj/(tot+0.0000001)  for proj, tot in zip(proj_sum[band][0], band_sum_tot[band])]
        frac_2[band] = [proj/(tot+0.0000001)  for proj, tot in zip(proj_sum[band][1], band_sum_tot[band])]
        frac_3[band] = [proj/(tot+0.0000001)  for proj, tot in zip(proj_sum[band][2], band_sum_tot[band])]
        # mixing_12[band]= [(sum-diff) / (sum+0.0000001) / 2 for sum, diff in zip(band_sum_12[band], diff_12[band])]
        # mixing_13[band]= [(sum-diff) / (sum+0.0000001) / 2 for sum, diff in zip(band_sum_13[band], diff_13[band])]
        # mixing_23[band]= [(sum-diff) / (sum+0.0000001) / 2 for sum, diff in zip(band_sum_23[band], diff_23[band])]
        # # print(no_kpoints)
        # print(len(mixing_12[band]))
        for k_point in range(0, no_kpoints):
            mixing[band].append([])
            for i in range(3):
                mixing[band][k_point].append([])
            max_frac=max(frac_1[band][k_point],frac_2[band][k_point],frac_3[band][k_point],0.001)
            # mixing[band][k_point][1] = mixing_12[band][k_point] #r - c=b
            # mixing[band][k_point][2] = mixing_13[band][k_point] #g - c=g
            # mixing[band][k_point][0] = mixing_23[band][k_point] #b - c=r
            mixing[band][k_point][2] = frac_1[band][k_point]/max_frac #r - c=b
            mixing[band][k_point][1] = frac_2[band][k_point]/max_frac #g - c=g
            mixing[band][k_point][0] = frac_3[band][k_point]/max_frac #b - c=r
        # print(mixing[band])
        # print(mixing[band][0])
        # print(len(mixing[band]))
        colors = ['b', 'lime', 'r', 'c', 'g', 'm', 'y', 'k', 'w']
        # colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
        # print(len(all_kpoints_coord),len(bands_eigenvalues[band]),len([i for i in mixing[band]]) )
        if (len(all_kpoints_coord)) != len(bands_eigenvalues[band]): #dla ISPIN = 2
            bands_eigenvalues_spin  = split_list_slicing(bands_eigenvalues[band])
            mixing_spin             = split_list_slicing(mixing[band])
            band_sum_sq_spin        = split_list_slicing(band_sum_sq[band])
            im      = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues_spin[0],
                         c=mixing_spin[0],
                         # c=[i for i in mixing_spin[0]],
                         s=band_sum_sq_spin[0], alpha=1.0)
                         # vmin=0.0, vmax=1.0,
            im_s2   = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues_spin[1],
                             c=mixing_spin[1],
                             # c=[i for i in mixing_spin[1]],
                             s=band_sum_sq_spin[1], alpha=1.0)
            # im_s2 = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues_spin[1], c=mixing_spin[1], s=band_sum_sq_spin[1],
            #                     cmap='brg', vmin=0.0, vmax=1.0, alpha=1.0)
        else:
            im = ax1.scatter(x=all_kpoints_coord, y=bands_eigenvalues[band],
                         # c=mixing[band],
                         c=[i for i in mixing[band]],
                         s=band_sum_sq[band],
                         # vmin=0.0, vmax=1.0,
                             alpha=1.0)
    #plot for (input_proj)!=2:
    else:
        for proj_no in range(len(input_proj)):
            proj_sum_sq[band][proj_no] = [x**2 for x in proj_sum[band][proj_no]]
            ax1.scatter(all_kpoints_coord, bands_eigenvalues[band], s=proj_sum_sq[band][proj_no],
                        color=colors[proj_no], alpha=alpha)

# def print_format_table():
#     """
#     prints table of formatted text format options
#     """
#     for style in range(8):
#         for fg in range(30,38):
#             s1 = ''
#             for bg in range(40,48):
#                 format = ';'.join([str(style), str(fg), str(bg)])
#                 s1 += '\x1b[%sm %s \x1b[0m' % (format, format)
#             print(s1)
#         print('\n')
# print_format_table()
if len(input_proj) != 2:
    print('\x1b[0;30;41m' + 'Saving band_data.dat aborted! \nWorks only for 2 projections!' + '\x1b[0m')

ax1.vlines(x = x_positions, ymin = y_min, ymax = y_max, colors = 'black', linewidth=0.25)
ax1.set(xlim=[0, kpath_cart[-1]], ylim=[y_min, y_max], ylabel='Energy (eV)') #, xlabel='k-path')
ax1.set_xticks(x_positions, labels=x_labels)
if 'y_step' in globals():
    list_up = np.arange(0, y_max+0.001, y_step)
    list_down = np.arange(0, y_min-0.001, -y_step)
    list_yticks = np.concatenate((list_down, list_up))
    ax1.set_yticks(list_yticks)
    # ax1.set_yticks(np.arange(y_min, y_max+0.001, y_step))

### działa
if len(input_proj)==2:
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 1.0, box.height * 1.0])
    fig1.colorbar(im, ax=ax1, ticks=[0, 0.5, 1],
                  shrink=0.25, aspect=10).set_ticklabels([input_strings[0], 'mixed', input_strings[1]])
                  # shrink=0.5).set_ticklabels([input_strings[0], 'mixed', input_strings[1]])
    # ax1.set_position([box.x0, box.y0, box.width * 1.0, box.height * 1.1])
    # ax1.set_aspect('auto')
    # ratio=1.0
    # xleft, xright = ax1.get_xlim()
    # ybottom, ytop = ax1.get_ylim()
    # ax1.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
else:
    box = ax1.get_position()
    ax1.legend(plot_legend, input_strings, loc='upper center', bbox_to_anchor=(1.30, 0.7),
               markerscale=4, frameon=False)
    ax1.set_position([box.x0, box.y0, box.width * 1.0, box.height * 1.0])

#rysowanie wykresu z czarnymi pasmami
# #plt.xlim(0,int(no_kpoints/2)) z ISPIN=2
for band in range(0, no_bands):
    if (len(all_kpoints_coord)) != len(bands_eigenvalues[band]) :
        bands_eigenvalues_spin = split_list_slicing(bands_eigenvalues[band])
        ax1.plot(all_kpoints_coord, bands_eigenvalues_spin[0], color="black", linewidth=0.5)
        ax1.plot(all_kpoints_coord, bands_eigenvalues_spin[1], color="black", linewidth=0.5)
    else:
        ax1.plot(all_kpoints_coord, bands_eigenvalues[band], color="black", linewidth=0.5)
#Fermi level line
ax1.axhline(y=0.0, color='black', linestyle='-', linewidth=0.25)
#rect - tuple (left, bottom, right, top), default: (0, 0, 1, 1)
# fig1.tight_layout(rect=[0.2, 0.2, 0.8, 2])
# fig1.tight_layout(rect=[0.2, 0.2, 0.8, 2])

#final saving plot
figure_title=os.path.join(folder, "Figure_1")
basename = figure_title + "-" + input_print + " " +" h" + str(figure_height) + "w" + str(figure_width) + " " + str(y_min) + "_" + str(y_max)
print()
print(basename + ".pdf",'\n')
fig1.set_size_inches(w=figure_width*cm, h=figure_height*cm)
# changes figure line widths on each axis
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(0.8)
ax1.tick_params(width=0.8)
ax1.set_position([box.x0+0.05, box.y0, box.width * plot_width, box.height * 1.0])
fig1.savefig(basename + ".png", format="png", dpi=500)
fig1.savefig(basename + ".pdf", format="pdf")
# fig1.savefig(basename + ".svg", format="svg")

if plot_spin_channels == 1:
    fig3, ax3 = plt.subplots()
    fig3.set_size_inches(w=figure_width*cm, h=figure_height*cm)
    for band in range(0, no_bands):
        if (len(all_kpoints_coord)) != len(bands_eigenvalues[band]):
            bands_eigenvalues_spin = split_list_slicing(bands_eigenvalues[band])
            ax3.plot(all_kpoints_coord, bands_eigenvalues_spin[0], color="blue", linewidth=0.5)
            ax3.plot(all_kpoints_coord, bands_eigenvalues_spin[1], color="red", linewidth=0.5)
        else:
            break
    # Fermi level line
    ax3.axhline(y=0.0, color='black', linestyle='-', linewidth=0.25)
    figure_spin_split = os.path.join(folder, "Figure_3")
    basename_spin = figure_spin_split + "-" + "spin_split" + " " + str(y_min) + " " + str(y_max)
    ax3.set(xlim=[0, kpath_cart[-1]], ylim=[y_min, y_max], ylabel='Energy (eV)')
    ax3.set_xticks(x_positions, labels=x_labels)
    ax3.axhline(y=0.0, color='black', linestyle='-', linewidth=0.25)
    ax3.set_position([box.x0+0.05, box.y0, box.width * 0.6, box.height * 1.0])
    fig3.savefig(basename_spin + ".png", format="png", dpi=500)
    fig3.savefig(basename_spin + ".pdf", format="pdf")

is_file = os.path.isfile(wannier_dir)
if is_file is False:
    print('\x1b[0;30;41m' + 'There is no wannier90.dat in directory. \nPrinting of the comparison plot aborted!' + '\x1b[0m')
    print("Energy of Fermi level:", efermi)
    print("--- %.2f seconds ---" % (time.time() - start_time))
    sys.exit()
else:
    print("Energy of Fermi level:", efermi)
    print("Printing wannierization plot...")

#plot wannier orbitals
f_wannier=open(wannier_dir, "r")
wan_data=f_wannier.read()
# print(wan_data)
wan_band_data=wan_data.split("\n  \n")
wan_band_data.pop(-1)
# print("hello", len(wan_band_data), len(wan_band_data[0]))
no_wan_bands=len(wan_band_data)
wan_band_en=[]
wan_band_kpts=[]
# wan_band_data_lines=[]
for band in range(no_wan_bands):
    wan_band_en.append([])
    wan_band_kpts.append([])
    wan_band_data_lines = wan_band_data[band].split("\n")
    # print(len(wan_band_data_lines)) - liczba k-pointów
    for line in wan_band_data_lines:
        wan_band_kpts[band].append(float(line.split()[0]))
        no_wan_kpoints=len(wan_band_kpts[band])
        wan_band_en[band].append(float(line.split()[1])-efermi)
print('numer of wannier bandsplot kpoints:', no_wan_kpoints, '\n')
wan_kpts=[]
# print(wan_band_kpts)
# wan_kpts_norm=[]
for kpoint in range(no_wan_kpoints):
    wan_kpts.append(kpoint)
    # wan_kpts_norm.append(kpoint)
    # wan_kpts_norm[kpoint]=(wan_kpts[kpoint]*no_kpoints)/no_wan_kpoints
wan_kpts_norm=[float((x * (kpath_cart[-1])) / (no_wan_kpoints-1)) for x in wan_kpts]
# print(wan_kpts_norm)
# print(kpoint_lst)
# print(wan_band_en[0]
# print(wan_band_kpts[0])

# plt.figure(2)
fig2, ax2 = plt.subplots()
ax2.vlines(x = x_positions, ymin = y_min, ymax = y_max, colors = 'black', linewidth=0.25)
for band in range(0, no_wan_bands):
    # print(band)
    ax2.plot(wan_kpts_norm, wan_band_en[band], color="red", linewidth=0.5)
# print(wan_band_en[-1])
for band in range(0, no_bands):
    # print(band)
    ax2.plot(all_kpoints_coord, bands_eigenvalues[band], color="black", linewidth=0.5)
# print(len(all_kpoints_coord))
str_fermi=round(efermi, 3)

figure_wan_title=os.path.join(folder, "Figure_2")
basename_wan = figure_wan_title + "-" + "wann90" + " " + str(y_min) + "_" + str(y_max) + " w" + str(figure_width) + " h" + str(figure_height) + " ef" + str(str_fermi)
figure_wan_pdf= basename_wan + ".pdf"
figure_wan_png= basename_wan + ".png"
ax2.set(xlim=[0, kpath_cart[-1]], ylim=[y_min, y_max], ylabel='Energy (eV)')
ax2.set_xticks(x_positions, labels=x_labels)
ax2.axhline(y=0.0, color='black', linestyle='-', linewidth=0.25)
box = ax2.get_position()
ax2.set_position([box.x0+0.05, box.y0, box.width * plot_width, box.height * 1.0])
if 'y_step' in globals():
    list_up = np.arange(0, y_max+0.001, y_step)
    list_down = np.arange(0, y_min-0.001, -y_step)
    list_yticks = np.concatenate((list_down, list_up))
    ax2.set_yticks(list_yticks)
    # ax1.set_yticks(np.arange(y_min, y_max+0.001, y_step))
fig2.set_size_inches(w=figure_width*cm, h=figure_height*cm)
fig2.savefig(figure_wan_png, format="png", dpi=500)
fig2.savefig(figure_wan_pdf, format="pdf")

print("--- %.2f seconds ---" % (time.time() - start_time))
#25.3 seconds
#0.8 seconds xddd
