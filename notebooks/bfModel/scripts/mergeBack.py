import re
import os
import glob
import shutil
import ruamel.yaml as yaml
from shutil import copyfile


def createOutputDir(output=None, makedir=False):
    """
    Create a directory to store outputs.
    """

    # Get output directory
    if output is not None:
        output = os.getcwd() + "/" + output
    else:
        output = os.getcwd() + "/output"

    if makedir:
        if os.path.exists(output):
            output += "_" + str(len(glob.glob(output + str("*"))) - 1)
    else:
        if os.path.exists(output):
            shutil.rmtree(output, ignore_errors=True)

    os.makedirs(output)
    os.makedirs(output + "/h5")
    os.makedirs(output + "/xmf")


def readOutputs(filename=None):

    # Check input file exists
    try:
        with open(filename) as finput:
            pass
    except IOError:
        print("Unable to open file: ", filename)
        raise IOError("The input file is not found...")

    # Open YAML file
    with open(filename, "r") as finput:
        input = yaml.load(finput, Loader=yaml.Loader)

    tStart, tEnd, tOut, outputDir = inputParser(input)

    nbCPUs = len(glob.glob1(outputDir + "/h5/", "topology.p*"))

    return tStart, tEnd, tOut, outputDir, nbCPUs


def inputParser(input):

    try:
        timeDict = input["time"]
    except KeyError:
        print("Key 'time' is required and is missing in the input file!")
        raise KeyError("Key time is required in the input file!")

    try:
        tStart = timeDict["start"]
    except KeyError:
        print("Key 'start' is required and is missing in the 'time' declaration!")
        raise KeyError("Simulation start time needs to be declared.")

    try:
        tEnd = timeDict["end"]
    except KeyError:
        print("Key 'end' is required and is missing in the 'time' declaration!")
        raise KeyError("Simulation end time needs to be declared.")

    try:
        tOut = timeDict["tout"]
    except KeyError:
        print("Key 'tout' is required and is missing in the 'time' declaration!")
        raise KeyError("Simulation output time needs to be declared.")

    try:
        outDict = input["output"]
        try:
            outputDir = outDict["dir"]
        except KeyError:
            outputDir = "output"
    except KeyError:
        outputDir = "output"

    return tStart, tEnd, tOut, outputDir


def save_DMPlex_XMF(
    output=None,
    cpus=None,
    step=None,
    start=None,
    end=None,
    dt=None,
    elems=None,
    nodes=None,
):
    """
    Saves mesh local information stored in the HDF5 to XMF file
    to visualise in Paraview.
    """

    time = start
    for k in range(step + 1):
        xmf_file = output + "/xmf/gospl" + str(k) + ".xmf"

        f = open(xmf_file, "w")

        # Header for xml file
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(" <Domain>\n")
        f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
        f.write('      <Time Type="Single" Value="%0.02f"/>\n' % time)

        for p in range(cpus):
            pfile = "h5/gospl." + str(k) + ".p" + str(p) + ".h5"
            tfile = "h5/topology.p" + str(p) + ".h5"
            f.write('      <Grid Name="Block.%s">\n' % (str(p)))
            f.write(
                '         <Topology Type="Triangle" NumberOfElements="%d" BaseOffset="1">\n'
                % elems[p]
            )
            f.write('          <DataItem Format="HDF" DataType="Int" ')
            f.write('Dimensions="%d 3">%s:/cells</DataItem>\n' % (elems[p], tfile))
            f.write("         </Topology>\n")

            f.write('         <Geometry Type="XYZ">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write('Dimensions="%d 3">%s:/coords</DataItem>\n' % (nodes[p], tfile))
            f.write("         </Geometry>\n")

            f.write('         <Attribute Type="Scalar" Center="Node" Name="Z">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write('Dimensions="%d 1">%s:/elev</DataItem>\n' % (nodes[p], pfile))
            f.write("         </Attribute>\n")
            f.write('         <Attribute Type="Vector" Center="Node" Name="hTec">\n')
            f.write(
                '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
            )
            f.write('Dimensions="%d 3">%s:/hdisp</DataItem>\n' % (nodes[p], pfile))
            f.write("         </Attribute>\n")
            f.write("      </Grid>\n")

        f.write("    </Grid>\n")
        f.write(" </Domain>\n")
        f.write("</Xdmf>\n")
        f.close()

        time += dt

    return


def save_XDMF(step=None, output=None):
    """
    This function writes the XDmF file which is calling the XmF file.
    """

    xdmf_file = output + "/gospl.xdmf"
    f = open(xdmf_file, "w")

    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
    f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write(" <Domain>\n")
    f.write('    <Grid GridType="Collection" CollectionType="Temporal">\n')

    for s in range(step + 1):
        xmf_file = "xmf/gospl" + str(s) + ".xmf"
        f.write(
            '      <xi:include href="%s" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n'
            % xmf_file
        )

    f.write("    </Grid>\n")
    f.write(" </Domain>\n")
    f.write("</Xdmf>\n")
    f.close()

    return


def getNodesElems(ndir=None, output=None, cpus=None):

    pathfile = ndir + "/xmf/gospl0.xmf"

    def findWholeWord(w):
        return re.compile(r"\b({0})\b".format(w), flags=re.IGNORECASE).search

    elems = []
    nodes = []
    with open(pathfile, "r") as file:
        for line in file:
            m = re.search('NumberOfElements="(\\d+)', line)
            if m:
                elems.append(int(m.group(1)))

            if '<Grid Name="Block.' in line:
                newBlock = True

            m = re.search('Precision="4" Dimensions="(\\d+)', line)
            if m:
                if len(nodes) == 0:
                    nodes.append(int(m.group(1)))
                    newBlock = False
                elif newBlock:
                    nodes.append(int(m.group(1)))
                    newBlock = False
                else:
                    if nodes[-1] != int(m.group(1)):
                        nodes.append(int(m.group(1)))

    # Copy topology files
    for nn in range(cpus):
        copyfile(
            ndir + "/h5/topology.p" + str(nn) + ".h5",
            output + "/h5/topology.p" + str(nn) + ".h5",
        )

    return elems, nodes


def mergeBackModels(inputs=None, output=None):

    odir = []
    instep = []
    for k in range(len(inputs)):
        if k == 0:
            tStart, tEnd, tOut, outputDir, nbCPUs = readOutputs(inputs[k])
            start = tStart
            end = tEnd
            dt = tOut
            odir.append(outputDir)
            cpus = nbCPUs
            instep.append(int((end - start) / dt))
        else:
            tStart, tEnd, tOut, outputDir, nbCPUs = readOutputs(inputs[k])
            if tOut != dt:
                raise ValueError(
                    "The simulations have to be ran with the same output time step..."
                )
            if nbCPUs != cpus:
                raise ValueError(
                    "The simulations have to be ran with the same number of CPUs..."
                )
            odir.append(outputDir)
            instep.append(int((tEnd - tStart) / tOut))
            start = min(tStart, start)
            end = max(tEnd, end)

    print("Read each input files done...")
    createOutputDir(output)
    step = int((end - start) / dt)
    save_XDMF(step, output)
    print("Create output directory and save XDMF file...")

    elems, nodes = getNodesElems(odir[0], output, cpus)

    save_DMPlex_XMF(output, cpus, step, start, end, dt, elems, nodes)
    print("Save XMF files...")

    k = -1
    for p in range(len(odir)):
        ndir = odir[p]
        nbstep = instep[p]
        if p == len(odir) - 1:
            nbstep += 1
        for s in range(nbstep):
            k += 1
            for n in range(cpus):
                copyfile(
                    ndir + "/h5/gospl." + str(s) + ".p" + str(n) + ".h5",
                    output + "/h5/gospl." + str(k) + ".p" + str(n) + ".h5",
                )

    print("Merging models done...")

    return
