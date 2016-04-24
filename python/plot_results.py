import os, sys, itertools
import numpy as np
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt
from matplotlib import cm

###########################################
# Plot charge density and electric fields #
###########################################

def plot_fields(output_filename):
    input_xml = et.parse(output_filename[0:-8] + ".xml").getroot()
    output_xml = et.parse(output_filename).getroot()
    dumps = output_xml.findall("dump")
    
    number_of_cells = int(output_xml.findtext("data/number_of_points"))
    number_of_times = len(dumps)
    cell_center_positions = np.fromstring(output_xml.findtext("data/position"), sep="\t")
    
    time_values = np.array([float(dump.get("time")) for dump in dumps])
    charge_density = np.empty((number_of_cells, number_of_times))
    electric_field = np.empty((number_of_cells, number_of_times))
    
    for i in range(len(dumps)):
        charge_density[:,i] = np.fromstring(dumps[i].findtext("delta_charge_density"), sep="\t")
        electric_field[:,i] = np.fromstring(dumps[i].findtext("electric_field_x"), sep="\t")

    # plot electric field and charge density
        
    x, y = np.meshgrid(time_values, cell_center_positions)
    
    plt.figure()
    plt.subplot(211)
    plt.contourf(x, y, charge_density)
    plt.xlabel("t")
    plt.ylabel("x")
    plt.title(r"$\delta\rho$")
    plt.colorbar()
    
    plt.subplot(212)
    plt.contourf(x, y, electric_field)
    plt.xlabel("t")
    plt.ylabel("x")
    plt.title(r"$E$")
    plt.colorbar()
    
    plt.tight_layout()
    plt.savefig("figures/" + output_filename[0:-8] + "-fields.pdf")
    plt.close()

    # plot max density and electric field
    
    if False:
        max_density = [max(charge_density[:,i]) for i in range(len(dumps))]
        max_field = [max(electric_field[:,i]) for i in range(len(dumps))]
        
        plt.figure()
        plt.plot(time_values, np.log(max_density))
        plt.xlabel(r"$t$")
        plt.ylabel(r"$\ln(\max(\rho))$")
        plt.savefig("figures/" + output_filename[0:-8] + "-max_rho.pdf")
        plt.close()
        
        logfield = np.log(max_field)
        plt.figure()
        plt.plot(time_values, logfield)
        plt.xlabel(r"$t$")
        plt.ylabel(r"$\ln(\max(\delta\rho))$")
        plt.savefig("figures/" + output_filename[0:-8] + "-max_drho.pdf")
        plt.close()

    # plot max electric field
    
    if False:
        wp = float(input_xml.findtext("particles/wp"))
        vth = float(input_xml.findtext("particles/thermal_velocity"))
        kp = float(input_xml.findtext("particles/kp"))
        la = vth / wp
        kla = kp * la
        gam = -0.5*np.sqrt(0.5 * np.pi) / np.power(kla, 3) * np.exp(- 0.5 / np.power(kla, 2))
        
        plt.figure()
        plt.plot(time_values, max_field, label="numeric")
        plt.plot(time_values, max_field[0] * np.exp(gam * time_values), label="analytic")
        plt.xlabel(r"$t$")
        plt.ylabel(r"$\max(E)$")
        plt.legend()
        plt.savefig("figures/" + output_filename[0:-8] + "-max_e.pdf")
        plt.close()
        
    # plot fourier electric field

    if False:
        num_modes_x = 10
        num_modes_t = 20
        fourier_electric = np.abs(np.fft.fft2(electric_field))
        fourier_electric = np.fft.fftshift(fourier_electric)
        fourier_electric = fourier_electric[number_of_cells/2:, number_of_times/2:]
        fourier_electric = fourier_electric[0:num_modes_x, 0:num_modes_t]
        
        x, y = np.meshgrid(np.arange(num_modes_t), np.arange(num_modes_x))
        plt.figure()
        plt.contourf(x, y, fourier_electric)
        plt.xlabel(r"$k_x$")
        plt.ylabel(r"$k_t$")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig("figures/" + output_filename[0:-8] + "-fourier.pdf")
        plt.close()
    
################
# Begin script #
################

if(len(sys.argv) != 2):
    sys.exit("usage: plot_results.py [output.xml]")

output_filename = sys.argv[1]

plot_fields(output_filename)
