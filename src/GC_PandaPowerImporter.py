import scipy.io
import pandapower as pp
import GridCalEngine as gce
import GC_utils
import os

# The initial proces had the following steps :
# pp.converter.to_mpc(gridPP,"network.mat")
# the network.mat file is opened with matlab/matpower
# am=load('network.mat')
# and saved as .m
# savecase('network.m',am.mpc)
# gridGCmp = gce.open_file("network.m")
# but this is optimized by converting the .mat file to .m directly in python

def __DeleteTmpFiles(filename):
    try:
        os.remove(filename)
        #print(f"File '{filename}' has been deleted.")
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
    except PermissionError:
        print(f"Permission denied: unable to delete '{filename}'.")
    except Exception as e:
        print(f"An error occurred while deleting '{filename}': {e}")

def PP2GC(gridPP):
    #by the moment, I still use two temporary files, which I need to delete at the end of the process
    tmp_file_in = "tmp.mat"
    tmp_file_out = "tmp.m"
    pp.converter.to_mpc(gridPP,tmp_file_in)
    # Load the .mat file
    mat_data = scipy.io.loadmat(tmp_file_in)   
    convert_mat_to_m(mat_data, tmp_file_out)
    gridGC = gce.open_file(tmp_file_out)
    gridGC.Sbase=gridPP.sn_mva
    __DeleteTmpFiles(tmp_file_in)
    __DeleteTmpFiles(tmp_file_out)
    return gridGC


def convert_mat_to_m(mat_data, output_m_file):
    """
    Convert a MATPOWER .mat file to a .m file.

    Parameters:
        input_mat_file (str): Path to the input .mat file containing the MATPOWER case.
        output_m_file (str): Path to the output .m file to be created.
    """
    
    # Extract the MATPOWER case structure (assuming it is stored under 'mpc')
    mpc = mat_data['mpc'][0, 0]

    # Open the output .m file for writing
    with open(output_m_file, 'w') as f:
        # Write the function definition and version
        f.write(f"function mpc = {output_m_file.split('.')[0]}\n")
        f.write("%% MATPOWER Case Format : Version 2\n")
        f.write("mpc.version = '2';\n\n")

        # Write the system MVA base
        f.write("%% system MVA base\n")
        f.write(f"mpc.baseMVA = {mpc['baseMVA'][0, 0]};\n\n")

        # Write the bus data
        f.write("%% bus data\n")
        f.write("%% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin\n")
        f.write("mpc.bus = [\n")
        for row in mpc['bus']:
            txtline=""
            for idx, element in enumerate(row):
                if idx<13:
                    txtline =  txtline + str(element) + " "
            #f.write("    " + " ".join(txtline) + ";\n")
            f.write("    " + txtline + ";\n")
        f.write("];\n\n")

        # Write the generator data
        f.write("%% generator data\n")
        f.write("%% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin ...\n")
        f.write("mpc.gen = [\n")
        for row in mpc['gen']:
            txtline=""
            for idx, element in enumerate(row):
                if idx<25:
                    txtline =  txtline + str(element) + " "
            #f.write("    " + " ".join(txtline) + ";\n")
            f.write("    " + txtline + ";\n")
        f.write("];\n\n")

        # Write the branch data
        f.write("%% branch data\n")
        f.write("%% fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax\n")
        f.write("mpc.branch = [\n")
        for row in mpc['branch']:
            f.write("    " + " ".join(map(str, row)) + ";\n")
        f.write("];\n\n")

        # Write the generator cost data (if it exists)
        if 'gencost' in mpc.dtype.names:
            f.write("%% generator cost data\n")
            f.write("%% 1 startup shutdown n x1 y1 ... xn yn\n")
            f.write("%% 2 startup shutdown n c(n-1) ... c0\n")
            f.write("mpc.gencost = [\n")
            for row in mpc['gencost']:
                f.write("    " + " ".join(map(str, row)) + ";\n")
            f.write("];\n\n")

        # Write the bus names (if they exist)
        if 'bus_name' in mpc.dtype.names:
            f.write("%% bus names\n")
            f.write("mpc.bus_name = {\n")
            for name in mpc['bus_name']:
                f.write(f"    '{name[0][0]}';\n")
            f.write("};\n\n")

    #print(f"MATPOWER case saved to {output_m_file}")


if __name__ == '__main__':

    gridPP=pp.networks.case33bw()
    pp.runpp(gridPP)
    print("powerflow of the original pandapower 33 buses network")
    print("solved in pandapower")
    print("   bus -5: ", gridPP.res_bus.tail(4))
    print("   power losses:", gridPP.res_line.pl_mw.sum(), gridPP.res_line.ql_mvar.sum())

    gridGC = PP2GC(gridPP)
    for line in gridGC.lines:
        line.active = True
    options = gce.PowerFlowOptions(gce.SolverType.NR, initialize_with_existing_solution=False,control_q=False, verbose=False)
    power_flowPP2GC = gce.PowerFlowDriver(gridGC, options)
    power_flowPP2GC.run()
    print("solved in gridcal")
    print("   ", power_flowPP2GC.results.get_bus_df().tail(4))
    print("   power losses:", power_flowPP2GC.results.losses.sum())


    gridPP=pp.networks.case118()
    pp.runpp(gridPP)
    print("powerflow of the original pandapower 118 buses network")
    print("solved in pandapower")
    print("   bus -5: ", gridPP.res_bus.tail(4))
    print("   power losses:", gridPP.res_line.pl_mw.sum(), gridPP.res_line.ql_mvar.sum())

    gridGC = PP2GC(gridPP)
    for line in gridGC.lines:
        line.active = True
    options = gce.PowerFlowOptions(gce.SolverType.NR, initialize_with_existing_solution=False,control_q=False, verbose=False)
    power_flowPP2GC = gce.PowerFlowDriver(gridGC, options)
    power_flowPP2GC.run()
    print("solved in gridcal")
    print("   ", power_flowPP2GC.results.get_bus_df().tail(4))
    print("   power losses:", power_flowPP2GC.results.losses.sum())

    import simbench as sb
    sb_code1 = "1-HVMV-urban-2.203-0-no_sw"
    gridPP = sb.get_simbench_net(sb_code1)
    gridPP.switch.drop([232,234,236,238,240, 242,244,246], inplace=True)
    gridPP.ext_grid.at[0,'name']="grid_ext"
    gridPP.line['in_service'] = True
    pp.runpp(gridPP)
    print("powerflow of the original simbench",sb_code1,"network")
    print("solved in pandapower")
    print("   bus -5: ", gridPP.res_bus.tail(4))
    print("   power losses:", gridPP.res_line.pl_mw.sum(), gridPP.res_line.ql_mvar.sum())

    gridGC = PP2GC(gridPP)
    for line in gridGC.lines:
        line.active = True
    options = gce.PowerFlowOptions(gce.SolverType.NR, initialize_with_existing_solution=False,control_q=False, verbose=False)
    power_flowPP2GC = gce.PowerFlowDriver(gridGC, options)
    power_flowPP2GC.run()
    print("solved in gridcal")
    print("   ", power_flowPP2GC.results.get_bus_df().tail(4))
    print("   power losses:", power_flowPP2GC.results.losses.sum())